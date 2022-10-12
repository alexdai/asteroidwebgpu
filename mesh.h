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