#pragma once
#include <vector>
#include <algorithm>
#include <random>

#include "mesh.h"

#include "GMatrix4.h"

//#include Settings
struct Settings
{

};

// We may want to ISPC-ify this down the road and just let it own the data structure in AoSoA format or similar
// For now we'll just do the dumb thing and see if it's fast enough
struct AsteroidDynamic
{
    CMatrix4f world;
    // These depend on chosen subdiv level, hence are not constant
    unsigned int indexStart;
    unsigned int indexCount;
};

struct AsteroidStatic
{
    CVector3f surfaceColor;
    CVector3f deepColor;
    CVector4f spinAxis;
    float scale;
    float spinVelocity;
    float orbitVelocity;
    unsigned int vertexStart;
    unsigned int textureIndex;
};

class AsteroidsSimulation
{
private:
    // NOTE: Memory could be optimized further for efficient cache traversal, etc.
    std::vector<AsteroidStatic> mAsteroidStatic;
    std::vector<AsteroidDynamic> mAsteroidDynamic;

    Mesh mMeshes;
    std::vector<unsigned int> mIndexOffsets;
    unsigned int mSubdivCount;
    unsigned int mVertexCountPerMesh;

    unsigned int mTextureDim;
    unsigned int mTextureCount;
    unsigned int mTextureArraySize;
    unsigned int mTextureMipLevels;
    //std::vector<BYTE> mTextureDataBuffer;
    //std::vector<D3D11_SUBRESOURCE_DATA> mTextureSubresources;

    //unsigned int SubresourceIndex(unsigned int texture, unsigned int arrayElement = 0, unsigned int mip = 0)
    //{
    //    return mip + mTextureMipLevels * (arrayElement + mTextureArraySize * texture);
    //}

    //void CreateTextures(unsigned int textureCount, unsigned int rngSeed);

public:
    AsteroidsSimulation(unsigned int rngSeed, unsigned int asteroidCount,
        unsigned int meshInstanceCount, unsigned int subdivCount,
        unsigned int textureCount);

    const Mesh* Meshes() { return &mMeshes; }
    //const D3D11_SUBRESOURCE_DATA* TextureData(unsigned int textureIndex)
    //{
    //    return mTextureSubresources.data() + SubresourceIndex(textureIndex);
    //}

    const AsteroidStatic* StaticData() const { return mAsteroidStatic.data(); }
    const AsteroidDynamic* DynamicData() const { return mAsteroidDynamic.data(); }

    // Can optionall provide a range of asteroids to update; count = 0 => to the end
    // This is useful for multithreading
    void Update(float frameTime, CVector4f cameraEye, const Settings& settings,
        size_t startIndex = 0, size_t count = 0);
};
