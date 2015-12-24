// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>



// TODO: reference additional headers your program requires here
#include <string>
#include <vector>
#include <array>
#include <deque>
#include <tuple>
#include <algorithm>
#include <limits>
#include <random>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#define _USE_MATH_DEFINES
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <PCG/pcg_random.hpp>
#include <tbb/blocked_range.h>

typedef OpenMesh::TriMesh_ArrayKernelT<> TriMesh;
static const float finf = std::numeric_limits<float>::infinity();
static const std::uniform_real_distribution<float> unif;