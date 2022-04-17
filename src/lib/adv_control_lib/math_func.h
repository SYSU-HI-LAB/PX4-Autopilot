#ifndef __MATH_FUNC_H__
#define __MATH_FUNC_H__

#include "matrix/math.hpp"
namespace MathFunc
{
    float floatLimit(float u, float limit);
    float floatConstrain(float u, float min, float max);
    uint8_t isZero(float u);
}

#endif
