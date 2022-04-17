#include "compensator.h"
#include "math_func.h"

using namespace MathFunc;

#ifndef C_PI_F
#define C_PI_F              (3.141592653589793f)
#endif

#define SQUARE(X)           ((X) * (X))

void LeadSys::initLeadSys(float kp, float lead_fc, float lead_gain, float fs)
{
    _lead_fc = lead_fc;
    _lead_gain = lead_gain;
    _kp = kp;

    float coef_z[2] = {};
    float coef_p[2] = {};

    coef_z[0] = 1.0f / (2.0f * C_PI_F * lead_fc) * kp;
    coef_z[1] = 1.0f * kp;
    coef_p[0] = 1.0f / (2.0f * C_PI_F * lead_fc * lead_gain);
    coef_p[1] = 1.0f;

    initTFcontinuous(coef_p, coef_z, fs);
}

void LeadSys::setOutBound(float out_min, float out_max)
{
    _out_min = out_min;
    _out_max = out_max;
}

void LeadSys:: updatewithSaturation(float u)
{
    update(u);
    _y = floatConstrain(_y, _out_min, _out_max);
}

void LeadboostSys::initLeadboostSys(float kp, float lead_fc, float lead_gain, float boost_fc, float boost_gain, float fs)
{
    _kp = kp;
    _lead_fc = lead_fc;
    _lead_gain = lead_gain;
    _boost_fc = boost_fc;
    _boost_gain = boost_gain;

    float coef_z[3] = {};
    float coef_p[3] = {};

    coef_z[0] = kp * boost_gain * lead_gain;
    coef_z[1] = kp * 2.0f * C_PI_F * boost_fc * boost_gain * lead_gain + 2.0f * C_PI_F * boost_gain * lead_fc * lead_gain;
    coef_z[2] = kp * 4.0f * SQUARE(C_PI_F) * boost_fc * boost_gain * lead_fc * lead_gain;

    coef_p[0] = boost_gain;
    coef_p[1] = 2.0f * C_PI_F * boost_fc + 2.0f * C_PI_F * boost_gain * lead_fc * lead_gain;
    coef_p[2] = 4.0f * SQUARE(C_PI_F) * boost_fc * lead_fc * lead_gain;

    initTFcontinuous(coef_p, coef_z, fs);
}

void LeadboostSys::setOutBound(float out_min, float out_max)
{
    _out_min = out_min;
    _out_max = out_max;
}

void LeadboostSys::updatewithSaturation(float u)
{
    update(u);
    _y(0, 0) = floatConstrain(_y(0, 0), _out_min, _out_max);
}

void IntSubsys::initIntSys(float fi, float i_out_min, float i_out_max, float reset_gain, float reset_dead_zone, float fs)
{
    _i_out_min = i_out_min;
    _i_out_max = i_out_max;
    _reset_gain = reset_gain;
    _reset_dead_zone = reset_dead_zone;
    _i_gain = 2.0f * C_PI_F * fi;
    _dynamic_i_gain = 1.0f;
    float coef_z[2] = {0.0f, _i_gain};
    float coef_p[2] = {1.0f, 0.0f};
    initTFcontinuous(coef_p, coef_z, fs);
}

void IntSubsys::setDynamicIgain(float dynamic_i_gain)
{
    if (fabsf(_dynamic_i_gain - dynamic_i_gain) > FLT_EPSILON)
    {
        /* record original input and output */
        float u = _u;
        float y = _y;

        /* set sys by dynamic_i_gain */
        _dynamic_i_gain = dynamic_i_gain;
        float coef_z[2] = {0.0f, _dynamic_i_gain * _i_gain};
        float coef_p[2] = {1.0f, 0.0f};

        initTFcontinuous(coef_p, coef_z, _fs);

        /* reset original input and output */
        switchBuf(u, y);
    }
}

void IntSubsys::updatewithSaturation(float u)
{
    if ((u * (getOutput() + _reset_dead_zone) < 0.0f)
     && (u * (getOutput() - _reset_dead_zone) < 0.0f))
    {
        setDynamicIgain(_reset_gain);
    }
    else
    {
        setDynamicIgain(1.0f);
    }

    update(u);

    if (getOutput() > _i_out_max)
    {
        switchBuf(u, _i_out_max);
    }
    else if (getOutput() < _i_out_min)
    {
        switchBuf(u, _i_out_min);
    }
}

float IntSubsys::getIoutMin()
{
    return _i_out_min;
}
float IntSubsys::getIoutMax()
{
    return _i_out_max;
}

void PISys::initPISys(float kp, float fi, float out_min, float out_max, float i_out_min, float i_out_max, float reset_gain, float reset_dead_zone, float fs)
{
    _kp = kp;
    _out_min = out_min;
    _out_max = out_max;
    _intesys.initIntSys(fi, i_out_min, i_out_max, reset_gain, reset_dead_zone, fs);
    _u = 0.0f;
    _y = 0.0f;
    _kp_term = 0.0f;
    _ki_term = 0.0f;
}

void PISys::update(float u)
{
    _kp_term = _kp * u;
    _intesys.updatewithSaturation(u);
    _ki_term = _intesys.getOutput();
    _y = _kp_term + _ki_term;
    _y = floatConstrain(_y, _out_min, _out_max);
}

void PISys::reset()
{
    _u = 0.0f;
    _y = 0.0f;
    _kp_term = 0.0f;
    _ki_term = 0.0f;
    _intesys.reset();
}

void PISys::switchBuf(float u, float y)
{
    _u = u;
    _kp_term = _kp * u;
    float ki_term = y - _kp_term;
    _ki_term = floatConstrain(ki_term, _intesys.getIoutMin(), _intesys.getIoutMax());
    _intesys.switchBuf(u, _ki_term);
    _y = floatConstrain(_kp_term + _ki_term, _out_min, _out_max);
}

float PISys::getOutput()
{
    return _y;
}

float PISys::getkpTerm()
{
    return _kp_term;
}

float PISys::getkiTerm()
{
    return _ki_term;
}

void PILeadSys::initPILeadSys(float kp, float fi, float lead_fc, float lead_gain, float out_min, float out_max, float i_out_min, float i_out_max,
    float reset_gain, float reset_dead_zone, float fs)
{
    _kp = kp;
    _out_min = out_min;
    _out_max = out_max;
    _pi_sys.initPISys(kp, fi, out_min, out_max, i_out_min, i_out_max, reset_gain, reset_dead_zone, fs);
    _lead_sys.initLeadSys(1.0f, lead_fc, lead_gain, fs); /* lead kp = 1.0 */
    _u = 0.0f;
    _y = 0.0f;
    _kp_term = 0.0f;
    _ki_term = 0.0f;
    _kcomp_term = 0.0f;
}

void PILeadSys::reset()
{
    _u = 0.0f;
    _y = 0.0f;
    _kp_term = 0.0f;
    _ki_term = 0.0f;
    _kcomp_term = 0.0f;
    _pi_sys.reset();
    _lead_sys.reset();
}

void PILeadSys::update(float u)
{
    /* update lead sys first without saturation */
    _u = u;
    _lead_sys.update(u);
    _kcomp_term = _lead_sys.getOutput();

    /* update pi sys then */
    _pi_sys.update(_kcomp_term);
    _kp_term = _pi_sys.getkpTerm();
    _ki_term = _pi_sys.getkiTerm();

    _y = _pi_sys.getOutput();
}

void PILeadSys::switchBuf(float u, float y)
{
    /* when switch buf, kcomp_term gain can be 1.0 */
    _u = u;
    _lead_sys.reset();
    _lead_sys.switchBuf(u, u);
    _pi_sys.switchBuf(u, y);
    _kcomp_term = u;
    _kp_term = _pi_sys.getkpTerm();
    _ki_term = _pi_sys.getkiTerm();
    _y = _pi_sys.getOutput();
}

float PILeadSys::getOutput()
{
    return _y;
}

void PILeadboostSys::initPILeadboostSys(float kp, float fi, float lead_fc, float lead_gain, float boost_fc, float boost_gain, float out_min, float out_max,
    float i_out_min, float i_out_max, float reset_gain, float reset_dead_zone, float fs)
{
    _kp = kp;
    _out_min = out_min;
    _out_max = out_max;
    _pi_sys.initPISys(kp, fi, out_min, out_max, i_out_min, i_out_max, reset_gain, reset_dead_zone, fs);
    _leadboost_sys.initLeadboostSys(1.0f, lead_fc, lead_gain, boost_fc, boost_gain, fs); /* lead kp = 1.0 */
    _u = 0.0f;
    _y = 0.0f;
    _kp_term = 0.0f;
    _ki_term = 0.0f;
    _kcomp_term = 0.0f;
}

void PILeadboostSys::reset()
{
    _u = 0.0f;
    _y = 0.0f;
    _kp_term = 0.0f;
    _ki_term = 0.0f;
    _kcomp_term = 0.0f;
    _pi_sys.reset();
    _leadboost_sys.reset();
}

void PILeadboostSys::update(float u)
{
    /* update lead sys first without saturation */
    _u = u;
    _leadboost_sys.update(u);
    _kcomp_term = _leadboost_sys.getOutput();

    /* update pi sys then */
    _pi_sys.update(_kcomp_term);
    _kp_term = _pi_sys.getkpTerm();
    _ki_term = _pi_sys.getkiTerm();

    _y = _pi_sys.getOutput();
}

void PILeadboostSys::switchBuf(float u, float y)
{
    /* when switch buf, kcomp_term gain can be 1.0 */
    _u = u;
    _leadboost_sys.reset();
    _leadboost_sys.switchBuf(u, u);
    _pi_sys.switchBuf(u, y);
    _kcomp_term = u;
    _kp_term = _pi_sys.getkpTerm();
    _ki_term = _pi_sys.getkiTerm();
    _y = _pi_sys.getOutput();
}

float PILeadboostSys::getOutput()
{
    return _y;
}

void NotchSys::initNotchSys(float fn, float bw, float gain, float fs)
{
    _fn = fn;
    _gain = gain;
    _bw = bw;

    float coef_z[3] = {};
    float coef_p[3] = {};

    float ksi1 = bw;
    float ksi2 = bw / gain;
    float wn = 2.0f * C_PI_F * fn;

    coef_z[0] = 1.0f;
    coef_z[1] = 2.0f * ksi1 * wn;
    coef_z[2] = SQUARE(wn);

    coef_p[0] = 1.0f;
    coef_p[1] = 2.0f * ksi2 * wn;
    coef_p[2] = SQUARE(wn);

    initTFcontinuous(coef_p, coef_z, fs);
}
