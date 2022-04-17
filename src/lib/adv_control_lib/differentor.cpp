#include "differentiator.h"

void Differentiator1st::initDifferentiator(float fc, float fs)
{
    setBasicCoef(_coef_p, 1, fc, fs);
    uint8_t i = 0;
    for (i = 0; i < 2; i++)
    {
        _coef_z[i] = 0.0f;
    }
    _coef_z[i - 2] = 1.0f;

    initTFcontinuous(_coef_p, _coef_z, fs);
}

void Differentiator2nd::initDifferentiator(float fc, float fs)
{
    setBasicCoef(_coef_p, 2, fc, fs);
    uint8_t i = 0;
    for (i = 0; i < 3; i++)
    {
        _coef_z[i] = 0.0f;
    }
    _coef_z[i - 2] = 1.0f;

    initTFcontinuous(_coef_p, _coef_z, fs);
}

void Differentiator3rd::initDifferentiator(float fc, float fs)
{
    setBasicCoef(_coef_p, 3, fc, fs);
    uint8_t i = 0;
    for (i = 0; i < 4; i++)
    {
        _coef_z[i] = 0.0f;
    }
    _coef_z[i - 2] = 1.0f;

    initTFcontinuous(_coef_p, _coef_z, fs);
}

void Differentiator4th::initDifferentiator(float fc, float fs)
{
    setBasicCoef(_coef_p, 4, fc, fs);
    uint8_t i = 0;
    for (i = 0; i < 5; i++)
    {
        _coef_z[i] = 0.0f;
    }
    _coef_z[i - 2] = 1.0f;

    initTFcontinuous(_coef_p, _coef_z, fs);
}

void Differentiator5th::initDifferentiator(float fc, float fs)
{
    setBasicCoef(_coef_p, 5, fc, fs);
    uint8_t i = 0;
    for (i = 0; i < 6; i++)
    {
        _coef_z[i] = 0.0f;
    }
    _coef_z[i - 2] = 1.0f;

    initTFcontinuous(_coef_p, _coef_z, fs);
}

void Differentiator6th::initDifferentiator(float fc, float fs)
{
    setBasicCoef(_coef_p, 6, fc, fs);
    uint8_t i = 0;
    for (i = 0; i < 7; i++)
    {
        _coef_z[i] = 0.0f;
    }
    _coef_z[i - 2] = 1.0f;

    initTFcontinuous(_coef_p, _coef_z, fs);
}
