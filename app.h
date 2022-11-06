#pragma once

#include <deque>
#include <random>
#include <map>

#include "settings.h"
#include "simulation.h"
#include "camera.h"

#include "webgpu.h"

struct DrawConstantBuffer {
    CMatrix4f mWorld;
    CMatrix4f mViewProjection;
    CVector3f mSurfaceColor;
    float unused0;
    CVector3f mDeepColor;
    float unused1;
    unsigned int mTextureIndex;
};

struct ExecuteIndirectArgs {
 //   D3D12_GPU_VIRTUAL_ADDRESS mConstantBuffer;
 //   D3D12_DRAW_INDEXED_ARGUMENTS mDrawIndexed;
};

struct DynamicUploadHeap {
    DrawConstantBuffer mDrawConstantBuffers[NUM_ASTEROIDS];
    ExecuteIndirectArgs mIndirectArgs[NUM_ASTEROIDS];
};

class Subset
{
public:
    Subset()
    {
    }

    ~Subset()
    {
    }

    void Begin()
    {
    }

    void End()
    {
    }
};

class Asteroids {
public:
    Asteroids(AsteroidsSimulation* asteroids, unsigned int minCmdLsts);
    ~Asteroids();

    void WaitForReadyToRender();
    void Render(float frameTime, const OrbitCamera& camera, const Settings& settings);

    void ReleaseSwapChain();
    void ResizeSwapChain(HWND outputWindow, unsigned int width, unsigned int height);

	static WGPUDevice device;
	static WGPUQueue queue;
	static WGPUSwapChain swapchain;

	static WGPURenderPipeline pipeline;
	static WGPUBuffer vertBuf; // vertex buffer with triangle position and colours
	static WGPUBuffer indxBuf; // index buffer
	static WGPUBuffer uRotBuf; // uniform buffer (containing the rotation angle)
	static WGPUBindGroup bindGroup;

	static char const blockrender_vert_wgsl[];
    static char const blockrender_frag_wgsl[];

    static WGPUShaderModule createShader(const char* const code, const char* label = nullptr);
    static WGPUBuffer createBuffer(const void* data, size_t size, WGPUBufferUsage usage);
    static void createPipelineAndBuffers();
private:
    void WaitForAll();

    void RenderSubset(
        size_t frameIndex, float frameTime,
        Subset* subset, UINT subsetIdx,
        CVector4f cameraEye, CMatrix4f viewProjection,
        const Settings& settings);

    void CreatePSOs();

    void CreateSubsets(UINT numHeapsPerFrame);
    void ReleaseSubsets();

    void CreateMeshes();

    struct Frame {
        std::vector<Subset*>   mSubsets;

        UINT64                      mFrameCompleteFence = 0;
    } mFrame[NUM_FRAMES_TO_BUFFER];

    // Swap chain resources
    struct SwapChainBuffer {
    } mSwapChainBuffer[NUM_SWAP_CHAIN_BUFFERS];

    // Swap chain stuff


    // Fences and synchronization
    size_t                      mCurrentFrameIndex = 0;
    HANDLE                      mFenceEventHandle = NULL;
    UINT64                      mCurrentFence = 0;

    AsteroidsSimulation* mAsteroids = nullptr;

    // Mesh
    UINT                        mIndexOffsets[MESH_MAX_SUBDIV_LEVELS + 2]; // inclusive
    UINT                        mNumVerticesPerMesh = 0;

    UINT                        mSubsetCount = 0;
    UINT                        mDrawsPerSubset = 0;
};
