#include "butterworth_filter.h"

#define BUTTER_COEF_1ST     {1.0f,     1.0f}
#define BUTTER_COEF_2ND     {1.0f,   1.414213562373095f,   1.0f}
#define BUTTER_COEF_3RD     {1.0f,   2.0f,   2.0f,   1.0f}
#define BUTTER_COEF_4TH     {1.0f,   2.613125929752753f,   3.414213562373095f,   2.613125929752753f,   1.0f}
#define BUTTER_COEF_5TH     {1.0f,   3.236067977499789f,   5.236067977499789f,   5.236067977499789f,   3.236067977499789f,   1.0f}
#define BUTTER_COEF_6TH     {1.0f,   3.863703305156273f,   7.464101615137753f,   9.141620172685640f,   7.464101615137753f,   3.863703305156273f,   1.0f}
#define BUTTER_COEF_7TH     {1.0f,   4.493959207434933f,  10.097834679044610f,  14.591793886479543f,  14.591793886479543f,  10.097834679044610f,   4.493959207434933f,   1.0f}
#define BUTTER_COEF_8TH     {1.0f,   5.125830895483012f,  13.137071184544089f,  21.846150969207624f,  25.688355931461274f,  21.846150969207628f,  13.137071184544089f,   5.125830895483013f,   1.0f}

#ifndef C_PI_F
#define C_PI_F              (3.141592653589793f)
#endif

void ButterworthFilter::setBasicCoef(float *coef, uint8_t Np, float fc, float fs)
{
    /* predistoration to ensure accurate phase in fc */
    float fd = predistoration(fc, fs);

    float coef0[MAX_BUTTER_ORDER + 1] = {};

    switch (Np)
    {
    case 1:
        set1stCoef(coef0);
        break;
    case 2:
        set2ndCoef(coef0);
        break;
    case 3:
        set3rdCoef(coef0);
        break;
    case 4:
        set4thCoef(coef0);
        break;
    case 5:
        set5thCoef(coef0);
        break;
    case 6:
        set6thCoef(coef0);
        break;
    case 7:
        set7thCoef(coef0);
        break;
    case 8:
        set8thCoef(coef0);
        break;
    default:
        coef0[0] = 1.0f;
    }

    memcpy(coef, coef0, (Np + 1) * sizeof(float));

    for (uint8_t i = 0; i < Np; i++)
    {
        *(coef + i) = *(coef + i) / powf((2.0f * C_PI_F * fd), (float)((Np - i)));
    }
}

void ButterworthFilter::setTFCoef(uint8_t Np, uint8_t Nz, float fp, float fz, float fs)
{
    _Np = Np;
    _Nz = Nz;
    _fp = fp;
    _fz = fz;
    _fs = fs;
    setBasicCoef(_coef_p, Np, fp, fs);
    setBasicCoef(_coef_z, Nz, fz, fs);
}

void ButterworthFilter::getCoefPolesZeros(float *coef_p, float *coef_z)
{
    memcpy(coef_p, _coef_p, (_Np + 1) * sizeof(float));
    memcpy(coef_z, _coef_z, (_Nz + 1) * sizeof(float));
}

float ButterworthFilter::predistoration(float fc, float fs)
{
   return (fs / C_PI_F) * tanf(C_PI_F / fs * fc);
}

void ButterworthFilter::set1stCoef(float *coef)
{
    float coef_norm[] = BUTTER_COEF_1ST;
    memcpy(coef, coef_norm, 2 * sizeof(float));
}

void ButterworthFilter::set2ndCoef(float *coef)
{
    float coef_norm[] = BUTTER_COEF_2ND;
    memcpy(coef, coef_norm, 3 * sizeof(float));
}

void ButterworthFilter::set3rdCoef(float *coef)
{
    float coef_norm[] = BUTTER_COEF_3RD;
    memcpy(coef, coef_norm, 4 * sizeof(float));
}

void ButterworthFilter::set4thCoef(float *coef)
{
    float coef_norm[] = BUTTER_COEF_4TH;
    memcpy(coef, coef_norm, 5 * sizeof(float));
}

void ButterworthFilter::set5thCoef(float *coef)
{
    float coef_norm[] = BUTTER_COEF_5TH;
    memcpy(coef, coef_norm, 6 * sizeof(float));
}

void ButterworthFilter::set6thCoef(float *coef)
{
    float coef_norm[] = BUTTER_COEF_6TH;
    memcpy(coef, coef_norm, 7 * sizeof(float));
}

void ButterworthFilter::set7thCoef(float *coef)
{
    float coef_norm[] = BUTTER_COEF_7TH;
    memcpy(coef, coef_norm, 8 * sizeof(float));
}

void ButterworthFilter::set8thCoef(float *coef)
{
    float coef_norm[] = BUTTER_COEF_8TH;
    memcpy(coef, coef_norm, 9 * sizeof(float));
}

void ButterworthFilter1st::initButterSysLowpass(uint8_t Nz, float fp, float fz, float fs)
{
    if (Nz > 1)
    {
        Nz = 1;
    }
    setTFCoef(1, Nz, fp, fz, fs);

    if (_Nz < _Np)
    {
        for (uint8_t i = 0; i < _Nz + 1; i++)
        {
            _coef_z[_Np - i] = _coef_z[_Nz - i];
        }
        for (uint8_t i = 0; i < _Np - _Nz; i++)
        {
            _coef_z[i] = 0.0f;
        }
    }

    initTFcontinuous(_coef_p, _coef_z, fs);
}

void ButterworthFilter2nd::initButterSysLowpass(uint8_t Nz, float fp, float fz, float fs)
{
    if (Nz > 2)
    {
        Nz = 2;
    }

    setTFCoef(2, Nz, fp, fz, fs);

    if (_Nz < _Np)
    {
        for (uint8_t i = 0; i < _Nz + 1; i++)
        {
            _coef_z[_Np - i] = _coef_z[_Nz - i];
        }
        for (uint8_t i = 0; i < _Np - _Nz; i++)
        {
            _coef_z[i] = 0.0f;
        }
    }

    initTFcontinuous(_coef_p, _coef_z, fs);
}

void ButterworthFilter3rd::initButterSysLowpass(uint8_t Nz, float fp, float fz, float fs)
{
    if (Nz > 3)
    {
        Nz = 3;
    }

    setTFCoef(3, Nz, fp, fz, fs);

    if (_Nz < _Np)
    {
        for (uint8_t i = 0; i < _Nz + 1; i++)
        {
            _coef_z[_Np - i] = _coef_z[_Nz - i];
        }
        for (uint8_t i = 0; i < _Np - _Nz; i++)
        {
            _coef_z[i] = 0.0f;
        }
    }

    initTFcontinuous(_coef_p, _coef_z, fs);
}

void ButterworthFilter4th::initButterSysLowpass(uint8_t Nz, float fp, float fz, float fs)
{
    if (Nz > 4)
    {
        Nz = 4;
    }

    setTFCoef(4, Nz, fp, fz, fs);

    if (_Nz < _Np)
    {
        for (uint8_t i = 0; i < _Nz + 1; i++)
        {
            _coef_z[_Np - i] = _coef_z[_Nz - i];
        }
        for (uint8_t i = 0; i < _Np - _Nz; i++)
        {
            _coef_z[i] = 0.0f;
        }
    }

    initTFcontinuous(_coef_p, _coef_z, fs);
}

void ButterworthFilter5th::initButterSysLowpass(uint8_t Nz, float fp, float fz, float fs)
{
    if (Nz > 5)
    {
        Nz = 5;
    }

    setTFCoef(5, Nz, fp, fz, fs);

    if (_Nz < _Np)
    {
        for (uint8_t i = 0; i < _Nz + 1; i++)
        {
            _coef_z[_Np - i] = _coef_z[_Nz - i];
        }
        for (uint8_t i = 0; i < _Np - _Nz; i++)
        {
            _coef_z[i] = 0.0f;
        }
    }

    initTFcontinuous(_coef_p, _coef_z, fs);
}

void ButterworthFilter6th::initButterSysLowpass(uint8_t Nz, float fp, float fz, float fs)
{
    if (Nz > 6)
    {
        Nz = 6;
    }

    setTFCoef(6, Nz, fp, fz, fs);

    if (_Nz < _Np)
    {
        for (uint8_t i = 0; i < _Nz + 1; i++)
        {
            _coef_z[_Np - i] = _coef_z[_Nz - i];
        }
        for (uint8_t i = 0; i < _Np - _Nz; i++)
        {
            _coef_z[i] = 0.0f;
        }
    }

    initTFcontinuous(_coef_p, _coef_z, fs);
}
