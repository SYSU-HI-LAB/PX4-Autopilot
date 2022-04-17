#include "math_func.h"

#define EPS (2.2204e-8f)

float MathFunc::floatLimit(float u, float limit)
{
    if (u < -limit)
    {
        return -limit;
    }
    else if (u > limit)
    {
        return limit;
    }
    else
    {
        return u;
    }
}

float MathFunc::floatConstrain(float u, float min, float max)
{
    if (u < min)
    {
        return min;
    }
    else if (u > max)
    {
        return max;
    }
    else
    {
        return u;
    }
}

uint8_t MathFunc::isZero(float u)
{
    if (fabsf(u) < EPS)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
