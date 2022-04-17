#ifndef __BUTTERWORTH_FILTER_H__
#define __BUTTERWORTH_FILTER_H__

#include "statespace_siso_sys.h"

#define MAX_BUTTER_ORDER (8)

class ButterworthFilter
{
public:
    void setBasicCoef(float *coef, uint8_t Np, float fc, float fs);
    void setTFCoef(uint8_t Np, uint8_t Nz, float fp, float fz, float fs);
    void getCoefPolesZeros(float *coef_p, float *coef_z);

#ifndef DEBUG
protected:
#endif
    float           _fp {};
    float           _fz {};
    float           _fs {};
    float           _coef_p[MAX_BUTTER_ORDER] {};
    float           _coef_z[MAX_BUTTER_ORDER] {};
    uint8_t         _Np {};
    uint8_t         _Nz {};

    float predistoration(float fc, float fs);
    void setCoef(float *coef, uint8_t Np, float fd);
    void set1stCoef(float *coef);
    void set2ndCoef(float *coef);
    void set3rdCoef(float *coef);
    void set4thCoef(float *coef);
    void set5thCoef(float *coef);
    void set6thCoef(float *coef);
    void set7thCoef(float *coef);
    void set8thCoef(float *coef);
};

class ButterworthFilter1st: public StatespaceSISOSys1st, public ButterworthFilter
{
public:
    void initButterSysLowpass(uint8_t Nz, float fp, float fz, float fs); /* Np is 1 for first-order filter */

protected:

};

class ButterworthFilter2nd: public StatespaceSISOSys2nd, public ButterworthFilter
{
public:
    void initButterSysLowpass(uint8_t Nz, float fp, float fz, float fs); /* Np is 1 for first-order filter */

protected:

};

class ButterworthFilter3rd: public StatespaceSISOSys3rd, public ButterworthFilter
{
public:
    void initButterSysLowpass(uint8_t Nz, float fp, float fz, float fs); /* Np is 1 for first-order filter */

protected:

};
class ButterworthFilter4th: public StatespaceSISOSys4th, public ButterworthFilter
{
public:
    void initButterSysLowpass(uint8_t Nz, float fp, float fz, float fs); /* Np is 1 for first-order filter */

protected:

};

class ButterworthFilter5th: public StatespaceSISOSys5th, public ButterworthFilter
{
public:
    void initButterSysLowpass(uint8_t Nz, float fp, float fz, float fs); /* Np is 1 for first-order filter */

protected:

};

class ButterworthFilter6th: public StatespaceSISOSys6th, public ButterworthFilter
{
public:
    void initButterSysLowpass(uint8_t Nz, float fp, float fz, float fs); /* Np is 1 for first-order filter */

protected:

};

#endif
