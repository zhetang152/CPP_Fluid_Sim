#pragma once

#include <vector>

#include "Grid_Construction/vecmath.h"

struct TriangleMesh {
    std::vector<Vector3f> vertices;
    std::vector<Vector3f> normals;
    std::vector<int> faces;
};
