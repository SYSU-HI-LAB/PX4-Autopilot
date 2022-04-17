#ifndef __DIFFERENTIATOR_H__
#define __DIFFERENTIATOR_H__

#include "butterworth_filter.h"

class Differentiator1st: public StatespaceSISOSys1st, public ButterworthFilter
{
public:
    void initDifferentiator(float fc, float fs);

#ifndef DEBUG
protected:
#endif
    float _fc;
};

class Differentiator2nd: public StatespaceSISOSys2nd, public ButterworthFilter
{
public:
    void initDifferentiator(float fc, float fs);

#ifndef DEBUG
protected:
#endif
    float _fc;
};

class Differentiator3rd: public StatespaceSISOSys3rd, public ButterworthFilter
{
public:
    void initDifferentiator(float fc, float fs);

#ifndef DEBUG
protected:
#endif
    float _fc;
};

class Differentiator4th: public StatespaceSISOSys4th, public ButterworthFilter
{
public:
    void initDifferentiator(float fc, float fs);

#ifndef DEBUG
protected:
#endif
    float _fc;
};

class Differentiator5th: public StatespaceSISOSys5th, public ButterworthFilter
{
public:
    void initDifferentiator(float fc, float fs);

#ifndef DEBUG
protected:
#endif
    float _fc;
};

class Differentiator6th: public StatespaceSISOSys6th, public ButterworthFilter
{
public:
    void initDifferentiator(float fc, float fs);

#ifndef DEBUG
protected:
#endif
    float _fc;
};

#endif
