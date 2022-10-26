#include "webgpu.h"
#include "../app.h"

#include <string.h>

WGPUDevice device;
WGPUQueue queue;
WGPUSwapChain swapchain;

WGPURenderPipeline pipeline;
WGPUBuffer vertBuf; // vertex buffer with triangle position and colours
WGPUBuffer indxBuf; // index buffer
WGPUBuffer uRotBuf; // uniform buffer (containing the rotation angle)
WGPUBindGroup bindGroup;

static char const blockrender_vert_wgsl[] = R"(
	struct VertexIn {
		@location(0) position : vec3<f32>;
		@location(1) normal : vec3<f32>;
	};
	struct VertexOut {
		builtin(position) position;
		@location(0) positionModel : vec3<f32>;
		@location(1) normalWorld : vec3<f32>;
		@location(2) albedo : vec3<f32>;
	};
	struct DrawConstantBuffer {
		@location(0) mWorld : mat4x4<f32>;
		@location(1) mViewProjection : mat4x4<f32>;
		@location(2) mSurfaceColor : mat4x4<f32>;
		@location(3) mDeepColor : mat4x4<f32>;
		@location(4) uint : mTextureIndex;
	};
	@group(0) @binding(0) var<uniform> uCB : DrawConstantBuffer;

	fn saturate(x: f32) -> f32 {
		return clamp(x, 0.0, 1.0);
	}

	fn linstep(min: f32, max: f32, s: f32) -> f32 {
		return saturate((s - min) / (max - min));
	}

	@stage(vertex)
	fn main(input : VertexIn) -> VertexOut {

		var output : VertexOut;

		var positionWorld : vec3<f32> = mul(mWorld, vec4<f32>(input.position, 1.0f)).xyz;
		output.position = mul(mViewProjection, vec4<f32>(positionWorld, 1.0f));

		output.positionModel = input.position;
		output.normalWorld = mul(mWorld, vec4(input.normal, 0.0f)).xyz; 

		var depth : f32 = linstep(0.5f, 0.7f, length(input.position.xyz));
		output.albedo = mix(mDeepColor.xyz, mSurfaceColor.xyz, depth);

		return output;
	}
)";


static char const blockrender_frag_wgsl[] = R"(
	fn saturate(x: f32) -> f32 {
		return clamp(x, 0.0, 1.0);
	}
	@stage(fragment)
	fn main( @location(0) positionModel : vec3<f32>, 
			 @location(1) normalWorld : vec3<f32>;
			 @location(2) albedo : vec3<f32>;) -> @location(0) vec4<f32> {

	// Tweaking
	var lightPos : vec3<f32>  = vec3<f32>(0.5, -0.25, -1);
	var applyNoise : bool = true;
	var applyLight : bool = true;
	var applyCoverage : bool = true;

	var normal : vec3 = normalize(input.normalWorld);

	// Triplanar projection
	var blendWeights : vec3<f32> = abs(normalize(input.positionModel));
	var uvw : vec3<f32> = input.positionModel * 0.5f + 0.5f;
	// Tighten up the blending zone
	blendWeights = saturate((blendWeights - 0.2f) * 7.0f);
	blendWeights /= (blendWeights.x + blendWeights.y + blendWeights.z).xxx;

	var coords1 : vec3<f32> = float3(uvw.yz, 0);
	var coords2 : vec3<f32> = float3(uvw.zx, 1);
	var coords3 : vec3<f32> = float3(uvw.xy, 2);

	// TODO: Should really branch out zero'd weight ones, but FXC is being a pain
	// and forward substituting the above and then refusing to compile "divergent"
	// coordinates...
	// Just disable the texture for now
	var detailTex : vec3<f32> = 1.0f;
	// float3 detailTex = 0.0f;
	// detailTex += blendWeights.x * Tex[mTextureIndex].Sample(Sampler, coords1).xyz;
	// detailTex += blendWeights.y * Tex[mTextureIndex].Sample(Sampler, coords2).xyz;
	// detailTex += blendWeights.z * Tex[mTextureIndex].Sample(Sampler, coords3).xyz;

	var wrap : f32 = 0.0f;
	var wrap_diffuse : f32 = saturate((dot(normal, normalize(lightPos)) + wrap) / (1.0f + wrap));
	var light : f32 = 3.0f * wrap_diffuse + 0.06f;

	// Approximate partial coverage on distant asteroids (by fading them out)
	var coverage : f32 = saturate(input.position.z * 4000.0f);

	var color : vec3<f32> = input.albedo;
	// flatten 3 if?
	if (applyNoise)    color = color * (2.0f * detailTex);
	if (applyLight)    color = color * light;
	if (applyCoverage) color = color * coverage;
	return vec4<f32>(color, 1.0f);
}
)";

/**
 * Helper to create a shader from WGSL source.
 *
 * \param[in] code WGSL shader source
 * \param[in] label optional shader name
 */
static WGPUShaderModule createShader(const char* const code, const char* label = nullptr) {
	WGPUShaderModuleWGSLDescriptor wgsl = {};
	wgsl.chain.sType = WGPUSType_ShaderModuleWGSLDescriptor;
	wgsl.source = code;
	WGPUShaderModuleDescriptor desc = {};
	desc.nextInChain = reinterpret_cast<WGPUChainedStruct*>(&wgsl);
	desc.label = label;
	return wgpuDeviceCreateShaderModule(device, &desc);
}

/**
 * Helper to create a buffer.
 *
 * \param[in] data pointer to the start of the raw data
 * \param[in] size number of bytes in \a data
 * \param[in] usage type of buffer
 */
static WGPUBuffer createBuffer(const void* data, size_t size, WGPUBufferUsage usage) {
	WGPUBufferDescriptor desc = {};
	desc.usage = WGPUBufferUsage_CopyDst | usage;
	desc.size  = size;
	WGPUBuffer buffer = wgpuDeviceCreateBuffer(device, &desc);
	wgpuQueueWriteBuffer(queue, buffer, 0, data, size);
	return buffer;
}

/**
 * Bare minimum pipeline to draw a triangle using the above shaders.
 */
static void createPipelineAndBuffers() {
	// compile shaders
	WGPUShaderModule vertMod = createShader(blockrender_vert_wgsl);
	WGPUShaderModule fragMod = createShader(blockrender_frag_wgsl);
	
	WGPUBufferBindingLayout buf = {};
	buf.type = WGPUBufferBindingType_Uniform;

	// bind group layout (used by both the pipeline layout and uniform bind group, released at the end of this function)
	WGPUBindGroupLayoutEntry bglEntry = {};
	bglEntry.binding = 0;
	bglEntry.visibility = WGPUShaderStage_Vertex;
	bglEntry.buffer = buf;

	WGPUBindGroupLayoutDescriptor bglDesc = {};
	bglDesc.entryCount = 1;
	bglDesc.entries = &bglEntry;
	WGPUBindGroupLayout bindGroupLayout = wgpuDeviceCreateBindGroupLayout(device, &bglDesc);

	// pipeline layout (used by the render pipeline, released after its creation)
	WGPUPipelineLayoutDescriptor layoutDesc = {};
	layoutDesc.bindGroupLayoutCount = 1;
	layoutDesc.bindGroupLayouts = &bindGroupLayout;
	WGPUPipelineLayout pipelineLayout = wgpuDeviceCreatePipelineLayout(device, &layoutDesc);

	// describe buffer layouts
	WGPUVertexAttribute vertAttrs[2] = {};
	vertAttrs[0].format = WGPUVertexFormat_Float32x2;
	vertAttrs[0].offset = 0;
	vertAttrs[0].shaderLocation = 0;
	vertAttrs[1].format = WGPUVertexFormat_Float32x3;
	vertAttrs[1].offset = 2 * sizeof(float);
	vertAttrs[1].shaderLocation = 1;
	WGPUVertexBufferLayout vertexBufferLayout = {};
	vertexBufferLayout.arrayStride = 5 * sizeof(float);
	vertexBufferLayout.attributeCount = 2;
	vertexBufferLayout.attributes = vertAttrs;

	// Fragment state
	WGPUBlendState blend = {};
	blend.color.operation = WGPUBlendOperation_Add;
	blend.color.srcFactor = WGPUBlendFactor_One;
	blend.color.dstFactor = WGPUBlendFactor_One;
	blend.alpha.operation = WGPUBlendOperation_Add;
	blend.alpha.srcFactor = WGPUBlendFactor_One;
	blend.alpha.dstFactor = WGPUBlendFactor_One;

	WGPUColorTargetState colorTarget = {};
	colorTarget.format = webgpu::getSwapChainFormat(device);
	colorTarget.blend = &blend;
	colorTarget.writeMask = WGPUColorWriteMask_All;

	WGPUFragmentState fragment = {};
	fragment.module = fragMod;
	fragment.entryPoint = "main";
	fragment.targetCount = 1;
	fragment.targets = &colorTarget;

	WGPURenderPipelineDescriptor desc = {};
	desc.fragment = &fragment;

	// Other state
	desc.layout = pipelineLayout;
	desc.depthStencil = nullptr;

	desc.vertex.module = vertMod;
	desc.vertex.entryPoint = "main";
	desc.vertex.bufferCount = 1;//0;
	desc.vertex.buffers = &vertexBufferLayout;

	desc.multisample.count = 1;
	desc.multisample.mask = 0xFFFFFFFF;
	desc.multisample.alphaToCoverageEnabled = false;

	desc.primitive.frontFace = WGPUFrontFace_CCW;
	desc.primitive.cullMode = WGPUCullMode_None;
	desc.primitive.topology = WGPUPrimitiveTopology_TriangleList;
	desc.primitive.stripIndexFormat = WGPUIndexFormat_Undefined;

	pipeline = wgpuDeviceCreateRenderPipeline(device, &desc);

	// partial clean-up (just move to the end, no?)
	wgpuPipelineLayoutRelease(pipelineLayout);

	wgpuShaderModuleRelease(fragMod);
	wgpuShaderModuleRelease(vertMod);

	// create the buffers (x, y, r, g, b)
	float const vertData[] = {
		-0.8f, -0.8f, 0.0f, 0.0f, 1.0f, // BL
		 0.8f, -0.8f, 0.0f, 1.0f, 0.0f, // BR
		-0.0f,  0.8f, 1.0f, 0.0f, 0.0f, // top
	};
	uint16_t const indxData[] = {
		0, 1, 2,
		0 // padding (better way of doing this?)
	};
	vertBuf = createBuffer(vertData, sizeof(vertData), WGPUBufferUsage_Vertex);
	indxBuf = createBuffer(indxData, sizeof(indxData), WGPUBufferUsage_Index);

	// create the uniform bind group (note 'rotDeg' is copied here, not bound in any way)
	uRotBuf = createBuffer(&rotDeg, sizeof(rotDeg), WGPUBufferUsage_Uniform);

	WGPUBindGroupEntry bgEntry = {};
	bgEntry.binding = 0;
	bgEntry.buffer = uRotBuf;
	bgEntry.offset = 0;
	bgEntry.size = sizeof(rotDeg);

	WGPUBindGroupDescriptor bgDesc = {};
	bgDesc.layout = bindGroupLayout;
	bgDesc.entryCount = 1;
	bgDesc.entries = &bgEntry;

	bindGroup = wgpuDeviceCreateBindGroup(device, &bgDesc);

	// last bit of clean-up
	wgpuBindGroupLayoutRelease(bindGroupLayout);
}

/**
 * Draws using the above pipeline and buffers.
 */
static bool redraw() {
	WGPUTextureView backBufView = wgpuSwapChainGetCurrentTextureView(swapchain);			// create textureView

	WGPURenderPassColorAttachment colorDesc = {};
	colorDesc.view    = backBufView;
	colorDesc.loadOp  = WGPULoadOp_Clear;
	colorDesc.storeOp = WGPUStoreOp_Store;
#ifdef __EMSCRIPTEN__
	// Dawn has both clearValue/clearColor but only Color works; Emscripten only has Value
	colorDesc.clearValue.r = 0.3f;
	colorDesc.clearValue.g = 0.3f;
	colorDesc.clearValue.b = 0.3f;
	colorDesc.clearValue.a = 1.0f;
#else
	colorDesc.clearColor.r = 0.3f;
	colorDesc.clearColor.g = 0.3f;
	colorDesc.clearColor.b = 0.3f;
	colorDesc.clearColor.a = 1.0f;
#endif

	WGPURenderPassDescriptor renderPass = {};
	renderPass.colorAttachmentCount = 1;
	renderPass.colorAttachments = &colorDesc;

	WGPUCommandEncoder encoder = wgpuDeviceCreateCommandEncoder(device, nullptr);			// create encoder
	WGPURenderPassEncoder pass = wgpuCommandEncoderBeginRenderPass(encoder, &renderPass);	// create pass

	// update the rotation
	rotDeg += 0.1f;
	wgpuQueueWriteBuffer(queue, uRotBuf, 0, &rotDeg, sizeof(rotDeg));

	// draw the triangle (comment these five lines to simply clear the screen)
	wgpuRenderPassEncoderSetPipeline(pass, pipeline);
	wgpuRenderPassEncoderSetBindGroup(pass, 0, bindGroup, 0, 0);
	wgpuRenderPassEncoderSetVertexBuffer(pass, 0, vertBuf, 0, WGPU_WHOLE_SIZE);
	wgpuRenderPassEncoderSetIndexBuffer(pass, indxBuf, WGPUIndexFormat_Uint16, 0, WGPU_WHOLE_SIZE);
	wgpuRenderPassEncoderDrawIndexed(pass, 3, 1, 0, 0, 0);

	//wgpuRenderPassEncoderDrawIndexedIndirect(pass, indirectBuffer, indirectOffset);

	wgpuRenderPassEncoderEnd(pass);
	wgpuRenderPassEncoderRelease(pass);														// release pass
	WGPUCommandBuffer commands = wgpuCommandEncoderFinish(encoder, nullptr);				// create commands
	wgpuCommandEncoderRelease(encoder);														// release encoder

	wgpuQueueSubmit(queue, 1, &commands);
	wgpuCommandBufferRelease(commands);														// release commands
#ifndef __EMSCRIPTEN__
	/*
	 * TODO: wgpuSwapChainPresent is unsupported in Emscripten, so what do we do?
	 */
	wgpuSwapChainPresent(swapchain);
#endif
	wgpuTextureViewRelease(backBufView);													// release textureView

	return true;
}

extern "C" int __main__(int /*argc*/, char* /*argv*/[]) {
	if (window::Handle wHnd = window::create()) {
		if ((device = webgpu::create(wHnd))) {
			queue = wgpuDeviceGetQueue(device);
			swapchain = webgpu::createSwapChain(device);
			createPipelineAndBuffers();

			window::show(wHnd);
			window::loop(wHnd, redraw);

		#ifndef __EMSCRIPTEN__
			wgpuBindGroupRelease(bindGroup);
			wgpuBufferRelease(uRotBuf);
			wgpuBufferRelease(indxBuf);
			wgpuBufferRelease(vertBuf);
			wgpuRenderPipelineRelease(pipeline);
			wgpuSwapChainRelease(swapchain);
			wgpuQueueRelease(queue);
			wgpuDeviceRelease(device);
		#endif
		}
	#ifndef __EMSCRIPTEN__
		window::destroy(wHnd);
	#endif
	}
	return 0;
}
