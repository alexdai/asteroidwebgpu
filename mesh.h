#pragma once

#include <vector>
#include "GVec3.h"

typedef unsigned short IndexType;

struct Vertex
{
    float x;
    float y;
    float z;
    float nx;
    float ny;
    float nz;
};

struct Mesh
{
    void clear()
    {
        vertices.clear();
        indices.clear();
    }

    std::vector<Vertex> vertices;
    std::vector<IndexType> indices;
};

void CreateIcosahedron(Mesh* outMesh);

// 1 face -> 4 faces
void SubdivideInPlace(Mesh* outMesh);

void SpherifyInPlace(Mesh* outMesh, float radius = 1.0f);

void ComputeAvgNormalsInPlace(Mesh* outMesh);

// subdivIndexOffset array should be [subdivLevels+2] in size
void CreateGeospheres(Mesh* outMesh, unsigned int subdivLevelCount, unsigned int* outSubdivIndexOffsets);

// Returns a combined "mesh" that includes:
// - A set of indices for each subdiv level (outSubdivIndexOffsets for offsets/counts)
// - A set of vertices for each mesh instance (base vertices per mesh computed from vertexCountPerMesh)
// - Indices already have the vertex offsets for the correct subdiv level "baked-in", so only need the mesh offset
void CreateAsteroidsFromGeospheres(Mesh* outMesh,
    unsigned int subdivLevelCount, unsigned int meshInstanceCount,
    unsigned int rngSeed,
    unsigned int* outSubdivIndexOffsets, unsigned int* vertexCountPerMesh);