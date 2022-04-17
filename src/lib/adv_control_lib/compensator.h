#ifndef __COMPENSATOR_H__
#define __COMPENSATOR_H__

#include "statespace_siso_sys.h"

/* basic ctrl sys: lead */
class LeadSys: public StatespaceSISOSys1st
{
public:
    void initLeadSys(float kp, float lead_fc, float lead_gain, float fs);
    void setOutBound(float out_min, float out_max);
    void updatewithSaturation(float u);

protected:
    float _kp;
    float _lead_fc;
    float _lead_gain;
    float _out_min;
    float _out_max;
};

/* basic ctrl sys: leadboost */
class LeadboostSys: public StatespaceSISOSys2nd
{
public:
    void initLeadboostSys(float kp, float lead_fc, float lead_gain, float boost_fc, float boost_gain, float fs);
    void setOutBound(float out_min, float out_max);
    void updatewithSaturation(float u);

protected:
    float _kp;
    float _lead_fc;
    float _lead_gain;
    float _boost_fc;
    float _boost_gain;
    float _out_min;
    float _out_max;
};

/* basic ctrl sys: integrade */
class IntSubsys: public StatespaceSISOSys1st
{
public:
    void initIntSys(float fi, float i_out_min, float i_out_max, float reset_gain, float reset_dead_zone, float fs);
    void updatewithSaturation(float u);
    void setDynamicIgain(float i_gain);
    void setInteOutbound(float i_out_min, float i_out_max);
    void setResetGain(float reset_gain);
    float getIoutMin();
    float getIoutMax();

#ifndef DEBUG
protected:
#endif

    float _fi;
    float _out_max;
    float _out_min;
    float _i_out_max;
    float _i_out_min;
    float _reset_gain;
    float _reset_dead_zone;
    float _i_gain;
    float _dynamic_i_gain;
};

/* derived ctrl sys: PI */
class PISys
{
public:
    void reset();
    void update(float u);
    void switchBuf(float u, float y);
    float getOutput();
    void initPISys(float kp, float fi, float out_min, float out_max, float i_out_min, float i_out_max, float reset_gain, float reset_dead_zone, float fs);
    float getkpTerm();
    float getkiTerm();

#ifndef DEBUG
protected:
#endif
    IntSubsys _intesys;
    float _kp;
    float _out_min;
    float _out_max;
    float _u;
    float _y;
    float _kp_term;
    float _ki_term;
};

/* derived ctrl sys: PI-lead */
class PILeadSys
{
public:
    void reset();
    void update(float u);
    void switchBuf(float u, float y);
    float getOutput();
    void initPILeadSys(float kp, float fi, float lead_fc, float lead_gain, float out_min, float out_max, float i_out_min, float i_out_max,
        float reset_gain, float reset_dead_zone, float fs);

#ifndef DEBUG
protected:
#endif
    PISys _pi_sys;
    LeadSys _lead_sys;
    float _kp;
    float _out_min;
    float _out_max;
    float _u;
    float _y;
    float _kcomp_term;
    float _kp_term;
    float _ki_term;
};

/* derived ctrl sys: PI-leadboost */
class PILeadboostSys
{
public:
    void reset();
    void update(float u);
    void switchBuf(float u, float y);
    float getOutput();
    void initPILeadboostSys(float kp, float fi, float lead_fc, float lead_gain, float boost_fc, float boost_gain, float out_min, float out_max,
        float i_out_min, float i_out_max, float reset_gain, float reset_dead_zone, float fs);

#ifndef DEBUG
protected:
#endif
    PISys _pi_sys;
    LeadboostSys _leadboost_sys;
    float _kp;
    float _out_min;
    float _out_max;
    float _u;
    float _y;
    float _kcomp_term;
    float _kp_term;
    float _ki_term;
};

class NotchSys: public StatespaceSISOSys2nd
{
public: 
    void initNotchSys(float fn, float bw, float gain, float fs);

#ifndef DEBUG
protected:
#endif
    float _fn;
    float _gain;
    float _bw;
};

#endif
