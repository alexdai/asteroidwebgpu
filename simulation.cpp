
#include "simulation.h"
#include "settings.h"
//#include "texture.h"
//#include "util.h"

#include <random>
#include <limits>
#include <algorithm>
#include <iostream>
// Non windows interface
//#include <ppl.h> 

static int const COLOR_SCHEMES[] = {
    156, 139, 113,  55,  49,  40,
    156, 139, 113,  58,  38,  14,
    156, 139, 113,  98, 101, 104,
    156, 139, 113, 205, 197, 178,
    153, 146, 136,  88,  88,  88,
    189, 181, 164, 148, 108, 102,
};

static int const NUM_COLOR_SCHEMES = (int)(sizeof(COLOR_SCHEMES) / (6 * sizeof(int)));

static CVector4f RandomPointOnSphere(std::mt19937& rng)
{
    // TODO: convert to our math format
    
    std::normal_distribution<float> dist;

    CVector4f r;
    for (;;) {
        r = CVector4f(dist(rng), dist(rng), dist(rng), 0.0f);
        auto d2 = r.SqrLength();
        if (d2 > std::numeric_limits<float>::min()) {
            return CVector4f(r.Normalize(), r.Normalize(), r.Normalize(), r.Normalize());
        }
    }
    
    // Unreachable
}

// From http://guihaire.com/code/?p=1135
static inline float VeryApproxLog2f(float x)
{
    union { float f; uint32_t i; } ux;
    ux.f = x;
    return (float)ux.i * 1.1920928955078125e-7f - 126.94269504f;
}


AsteroidsSimulation::AsteroidsSimulation(unsigned int rngSeed, unsigned int asteroidCount,
    unsigned int meshInstanceCount, unsigned int subdivCount,
    unsigned int textureCount)
    : mAsteroidStatic(asteroidCount)
    , mAsteroidDynamic(asteroidCount)
    , mIndexOffsets(subdivCount + 2) // Mesh subdivs are inclusive on both ends and need forward differencing for count
    , mSubdivCount(subdivCount)
{
    std::mt19937 rng(rngSeed);

    // Create meshes
    std::cout
        << "Creating " << meshInstanceCount << " meshes, each with "
        << subdivCount << " subdivision levels..." << std::endl;

    CreateAsteroidsFromGeospheres(&mMeshes, mSubdivCount, meshInstanceCount,
        rng(), mIndexOffsets.data(), &mVertexCountPerMesh);

    //CreateTextures(textureCount, rng());

    // Constants
    std::normal_distribution<float> orbitRadiusDist(SIM_ORBIT_RADIUS, 0.6f * SIM_DISC_RADIUS);
    std::normal_distribution<float> heightDist(0.0f, 0.4f);
    std::uniform_real_distribution<float> angleDist(-M_PI, M_PI);
    std::uniform_real_distribution<float> radialVelocityDist(5.0f, 15.0f);
    std::uniform_real_distribution<float> spinVelocityDist(-2.0f, 2.0f);
    std::normal_distribution<float> scaleDist(1.3f, 0.7f);
    std::normal_distribution<float> colorSchemeDist(0, NUM_COLOR_SCHEMES - 1);
    std::uniform_int_distribution<unsigned int> textureIndexDist(0, textureCount - 1);

    auto instancesPerMesh = std::max(1U, asteroidCount / meshInstanceCount);

    // Approximate SRGB->Linear for colors
    float linearColorSchemes[NUM_COLOR_SCHEMES * 6];
    for (int i = 0; i < NUM_COLOR_SCHEMES * 6; ++i) {
        linearColorSchemes[i] = std::powf((float)COLOR_SCHEMES[i] / 255.0f, 2.2f);
    }

    // Create a torus of asteroids that spin around the ring
    for (unsigned int i = 0; i < asteroidCount; ++i) {
        auto scale = scaleDist(rng);
#if SIM_USE_GAMMA_DIST_SCALE
        scale = scale * 0.3f;
#endif
        scale = std::max(scale, SIM_MIN_SCALE);
        auto scaleMatrix = CMatrix4f::Scale(CVector3f::Zero, scale);

        auto orbitRadius = orbitRadiusDist(rng);
        auto discPosY = float(SIM_DISC_RADIUS) * heightDist(rng);

        auto disc = CMatrix4f::Translate(orbitRadius, discPosY, 0.0f);

        auto positionAngle = angleDist(rng);
        auto orbit = CMatrix4f::Rotate(CVector3f::Zero, CVector3f(0.f, 1.f, 0.f), positionAngle);

        auto meshInstance = (unsigned int)(i / instancesPerMesh); // Vcache friendly ordering

        // Static data
        mAsteroidStatic[i].spinVelocity = spinVelocityDist(rng) / scale; // Smaller asteroids spin faster
        mAsteroidStatic[i].orbitVelocity = radialVelocityDist(rng) / (scale * orbitRadius); // Smaller asteroids go faster, and use arc length
        mAsteroidStatic[i].vertexStart = mVertexCountPerMesh * meshInstance;
        mAsteroidStatic[i].spinAxis = RandomPointOnSphere(rng);
        mAsteroidStatic[i].spinAxis.Normalize();
        mAsteroidStatic[i].scale = scale;
        mAsteroidStatic[i].textureIndex = textureIndexDist(rng);

        auto colorScheme = ((int)abs(colorSchemeDist(rng))) % NUM_COLOR_SCHEMES;
        auto c = linearColorSchemes + 6 * colorScheme;
        mAsteroidStatic[i].surfaceColor = CVector3f(c[0], c[1], c[2]);
        mAsteroidStatic[i].deepColor = CVector3f(c[3], c[4], c[5]);

        // Initialize dynamic data
        mAsteroidDynamic[i].world = scaleMatrix * disc * orbit;

        assert(mAsteroidStatic[i].scale > 0.0f);
        assert(mAsteroidStatic[i].orbitVelocity > 0.0f);
    }
}


void AsteroidsSimulation::Update(float frameTime, CVector4f cameraEye, const Settings& settings,
    size_t startIndex, size_t count)
{
    bool animate = settings.animate;

    // TODO: This constant should really depend on resolution and/or be configurable...
    static const float minSubdivSizeLog2 = std::log2f(0.0019f);

    size_t last = count ? startIndex + count : mAsteroidDynamic.size();
    for (size_t i = startIndex; i < last; ++i) {
        const AsteroidStatic& staticData = mAsteroidStatic[i];
        AsteroidDynamic& dynamicData = mAsteroidDynamic[i];

        if (animate) {
            auto orbit = CMatrix4f::Rotate(CVector3f::Zero, CVector3f(0.f, 1.f, 0.f), staticData.orbitVelocity * frameTime);
            auto spin = CMatrix4f::Rotate(CVector3f::Zero, staticData.spinAxis.Vec3(), staticData.spinVelocity * frameTime);
            dynamicData.world = spin * dynamicData.world * orbit;
        }

        // Pick LOD based on approx screen area - can be very approximate
        auto position = dynamicData.world[3];
        auto distanceToEyeRcp = 1.f / (cameraEye - position).Length();
        // Add one subdiv for each factor of 2 past min
        auto relativeScreenSizeLog2 = VeryApproxLog2f(staticData.scale * distanceToEyeRcp);
        float subdivFloat = std::max(0.0f, relativeScreenSizeLog2 - minSubdivSizeLog2);
        auto subdiv = std::min(mSubdivCount, (unsigned int)subdivFloat);

        // TODO: Ignore/cull/force lowest subdiv if offscreen?

        dynamicData.indexStart = mIndexOffsets[subdiv];
        dynamicData.indexCount = mIndexOffsets[subdiv + 1] - dynamicData.indexStart;
    }
}