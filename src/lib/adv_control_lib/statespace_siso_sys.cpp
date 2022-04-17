#include "statespace_siso_sys.h"

/* -------------------- for first-order system -------------------- */
void StatespaceSISOSys1st::reset()
{
    _u = 0.0f;
    _y = 0.0f;
    _x = 0.0f;
}

void StatespaceSISOSys1st::initSSdiscrete(float a, float b, float c, float d)
{
    if (!_inited_flag)
    {
        reset();
        _inited_flag = 1;
    }
    _Ad = a;
    _Bd = b;
    _Cd = c;
    _Dd = d;
    _Cd_pinv = 1 / c;
}

void StatespaceSISOSys1st::initSScontinuous(float a, float b, float c, float d, float fs)
{
    _Ac = a;
    _Bc = b;
    _Cc = c;
    _Dc = d;
    _fs = fs;

    float alpha = 2.0f * _fs;

    _Ad = (alpha + _Ac) / (alpha - _Ac);
    _Bd = _Bc / (alpha - _Ac);
    _Cd = _Cc * (_Ad + 1.0f);
    _Dd = _Cc * _Bd + _Dc;

    initSSdiscrete(_Ad, _Bd, _Cd, _Dd);
}


void StatespaceSISOSys1st::initTFcontinuous(float *coef_p, float *coef_z, float fs)
{
    for (uint8_t i = 0; i < 2; i++)
    {
        _num[i] = coef_z[i];
        _den[i] = coef_p[i];
    }

    float a0 = coef_p[1] / coef_p[0];
    float b0 = coef_z[1] / coef_p[0];
    float b1 = coef_z[0] / coef_p[0];
    float h0 = b1;
    float h1 = b0 - a0 * h0;

    _Ac = -a0;
    _Bc = h1;
    _Cc = 1.0f;
    _Dc = h0;

    initSScontinuous(_Ac, _Bc, _Cc, _Dc, fs);
}

void StatespaceSISOSys1st::update(float u)
{
    _u = u;
    _y = _Cd * _x + _Dd * _u;
    _x = _Ad * _x + _Bd * _u;
}

void StatespaceSISOSys1st::switchBuf(float u, float y)
{
    _u = u;
    _y = y;
    _x = _Cd_pinv * (_y - _Dd * _u);
    _x = _Ad * _x + _Bd * _u;
}

float StatespaceSISOSys1st::getOutput()
{
    return _y;
}



/* -------------------- for second-order system -------------------- */
void StatespaceSISOSys2nd::reset()
{
    _u.setZero();
    _y.setZero();
    _x.setZero();
}

void StatespaceSISOSys2nd::initSSdiscrete(matrix::SquareMatrix<float, 2> a, matrix::Matrix<float, 2, 1> b, matrix::Matrix<float, 1, 2> c, matrix::Matrix<float, 1, 1> d)
{
    if (!_inited_flag)
    {
        reset();
        _inited_flag = 1;
    }
    _Ad = a;
    _Bd = b;
    _Cd = c;
    _Dd = d;
    float tmp = _Cd(0, 0) * _Cd(0, 0) + _Cd(0, 1) * _Cd(0, 1);
    for (uint8_t i = 0; i < 2; i ++)
    {
        _Cd_pinv(i, 0) = c(0, i) / tmp;
    }
}

void StatespaceSISOSys2nd::initSScontinuous(matrix::SquareMatrix<float, 2> a, matrix::Matrix<float, 2, 1> b, matrix::Matrix<float, 1, 2> c, matrix::Matrix<float, 1, 1> d, float fs)
{
    _Ac = a;
    _Bc = b;
    _Cc = c;
    _Dc = d;
    _fs = fs;

    matrix::SquareMatrix<float, 2> Eye_alpha;
    matrix::SquareMatrix<float, 2> Eye;
    Eye.setIdentity();
    Eye_alpha.setIdentity();

    float alpha = 2.0f * _fs;

    for (uint8_t i = 0; i < 2; i++)
    {
        Eye_alpha(i,i) = alpha;
    }

    matrix::SquareMatrix<float, 2> tmp;
    tmp = Eye_alpha - _Ac;

    _Ad = inv(tmp) * (Eye_alpha + _Ac);
    _Bd = inv(tmp) * _Bc;
    _Cd = _Cc * (_Ad + Eye);
    _Dd = _Cc * _Bd + _Dc;

    initSSdiscrete(_Ad, _Bd, _Cd, _Dd);
}


void StatespaceSISOSys2nd::initTFcontinuous(float *coef_p, float *coef_z, float fs)
{
    for (uint8_t i = 0; i < 3; i++)
    {
        _num[i] = coef_z[i];
        _den[i] = coef_p[i];
    }

    float a[2] = {};
    uint8_t i = 0;
    for (i = 0; i < 2; i++)
    {
        a[i] = coef_p[2 - i] / coef_p[0];
    }

    float b[3] = {};
    for (i = 0; i < 3; i++)
    {
        b[i] = coef_z[2 - i] / coef_p[0];
    }

    float h[3] = {};
    h[0] = b[2];
    h[1] = b[1] - a[1] * h[0];
    h[2] = b[0] - a[1] * h[1] - a[0] * h[0];

    const float Ac[] ={0.0f,  1.0f,
                      -a[0], -a[1]};

    matrix::SquareMatrix<float, 2> Ac_value(Ac);
    _Ac = Ac_value;

    const float Bc[] = {h[1],
                        h[2]};
    matrix::Matrix<float, 2, 1> Bc_value(Bc);
    _Bc = Bc_value;

    _Cc.setZero();
    _Cc(0, 0) = 1.0f;

    _Dc(0, 0) = h[0];

    initSScontinuous(_Ac, _Bc, _Cc, _Dc, fs);
}

void StatespaceSISOSys2nd::update(float u)
{
    _u(0, 0) = u;
    _y = _Cd * _x + _Dd * _u;
    _x = _Ad * _x + _Bd * _u;
}

void StatespaceSISOSys2nd::switchBuf(float u, float y)
{
    _u(0, 0) = u;
    _y(0, 0) = y;
    _x = _Cd_pinv * (_y - _Dd * _u);
    _x = _Ad * _x + _Bd * _u;
}

float StatespaceSISOSys2nd::getOutput()
{
    return _y(0, 0);
}

/* -------------------- for third-order system -------------------- */
void StatespaceSISOSys3rd::reset()
{
    _u.setZero();
    _y.setZero();
    _x.setZero();
}

void StatespaceSISOSys3rd::initSSdiscrete(matrix::SquareMatrix<float, 3> a, matrix::Matrix<float, 3, 1> b, matrix::Matrix<float, 1, 3> c, matrix::Matrix<float, 1, 1> d)
{
    if (!_inited_flag)
    {
        reset();
        _inited_flag = 1;
    }
    _Ad = a;
    _Bd = b;
    _Cd = c;
    _Dd = d;
    float tmp = _Cd(0, 0) * _Cd(0, 0) + _Cd(0, 1) * _Cd(0, 1) + _Cd(0, 2) * _Cd(0, 2);
    for (uint8_t i = 0; i < 3; i ++)
    {
        _Cd_pinv(i, 0) = c(0, i) / tmp;
    }
}

void StatespaceSISOSys3rd::initSScontinuous(matrix::SquareMatrix<float, 3> a, matrix::Matrix<float, 3, 1> b, matrix::Matrix<float, 1, 3> c, matrix::Matrix<float, 1, 1> d, float fs)
{
    _Ac = a;
    _Bc = b;
    _Cc = c;
    _Dc = d;
    _fs = fs;

    matrix::SquareMatrix<float, 3> Eye_alpha;
    matrix::SquareMatrix<float, 3> Eye;
    Eye.setIdentity();
    Eye_alpha.setIdentity();

    float alpha = 2.0f * _fs;

    for (uint8_t i = 0; i < 3; i++)
    {
        Eye_alpha(i,i) = alpha;
    }

    matrix::SquareMatrix<float, 3> tmp;
    tmp = Eye_alpha - _Ac;

    _Ad = inv(tmp) * (Eye_alpha + _Ac);
    _Bd = inv(tmp) * _Bc;
    _Cd = _Cc * (_Ad + Eye);
    _Dd = _Cc * _Bd + _Dc;

    initSSdiscrete(_Ad, _Bd, _Cd, _Dd);
}

void StatespaceSISOSys3rd::initTFcontinuous(float *coef_p, float *coef_z, float fs)
{
    for (uint8_t i = 0; i < 4; i++)
    {
        _num[i] = coef_z[i];
        _den[i] = coef_p[i];
    }

    float a[3] = {};
    uint8_t i = 0;
    for (i = 0; i < 3; i++)
    {
        a[i] = coef_p[3 - i] / coef_p[0];
    }

    float b[4] = {};
    for (i = 0; i < 4; i++)
    {
        b[i] = coef_z[3 - i] / coef_p[0];
    }

    float h[4] = {};
    h[0] = b[3];
    h[1] = b[2] - a[2] * h[0];
    h[2] = b[1] - a[2] * h[1] - a[1] * h[0];
    h[3] = b[0] - a[2] * h[2] - a[1] * h[1] - a[0] * h[0];

    const float Ac[] ={0.0f,  1.0f,  0.0f,
                       0.0f,  0.0f,  1.0f,
                      -a[0], -a[1], -a[2]};

    matrix::SquareMatrix<float, 3> Ac_value(Ac);
    _Ac = Ac_value;

    const float Bc[] = {h[1],
                        h[2],
                        h[3]};
    matrix::Matrix<float, 3, 1> Bc_value(Bc);
    _Bc = Bc_value;

    _Cc.setZero();
    _Cc(0, 0) = 1.0f;

    _Dc(0, 0) = h[0];

    initSScontinuous(_Ac, _Bc, _Cc, _Dc, fs);
}

void StatespaceSISOSys3rd::update(float u)
{
    _u(0, 0) = u;
    _y = _Cd * _x + _Dd * _u;
    _x = _Ad * _x + _Bd * _u;
}

void StatespaceSISOSys3rd::switchBuf(float u, float y)
{
    _u(0, 0) = u;
    _y(0, 0) = y;
    _x = _Cd_pinv * (_y - _Dd * _u);
    _x = _Ad * _x + _Bd * _u;
}

float StatespaceSISOSys3rd::getOutput()
{
    return _y(0, 0);
}

/* -------------------- for fourth-order system -------------------- */
void StatespaceSISOSys4th::reset()
{
    _u.setZero();
    _y.setZero();
    _x.setZero();
}

void StatespaceSISOSys4th::initSSdiscrete(matrix::SquareMatrix<float, 4> a, matrix::Matrix<float, 4, 1> b, matrix::Matrix<float, 1, 4> c, matrix::Matrix<float, 1, 1> d)
{
    if (!_inited_flag)
    {
        reset();
        _inited_flag = 1;
    }
    _Ad = a;
    _Bd = b;
    _Cd = c;
    _Dd = d;
    float tmp = _Cd(0, 0) * _Cd(0, 0) + _Cd(0, 1) * _Cd(0, 1) + _Cd(0, 2) * _Cd(0, 2) + _Cd(0, 3) * _Cd(0, 3);
    for (uint8_t i = 0; i < 4; i ++)
    {
        _Cd_pinv(i, 0) = c(0, i) / tmp;
    }
}

void StatespaceSISOSys4th::initSScontinuous(matrix::SquareMatrix<float, 4> a, matrix::Matrix<float, 4, 1> b, matrix::Matrix<float, 1, 4> c, matrix::Matrix<float, 1, 1> d, float fs)
{
    _Ac = a;
    _Bc = b;
    _Cc = c;
    _Dc = d;
    _fs = fs;

    matrix::SquareMatrix<float, 4> Eye_alpha;
    matrix::SquareMatrix<float, 4> Eye;
    Eye.setIdentity();
    Eye_alpha.setIdentity();

    float alpha = 2.0f * _fs;

    for (uint8_t i = 0; i < 4; i++)
    {
        Eye_alpha(i,i) = alpha;
    }

    matrix::SquareMatrix<float, 4> tmp;
    tmp = Eye_alpha - _Ac;

    _Ad = inv(tmp) * (Eye_alpha + _Ac);
    _Bd = inv(tmp) * _Bc;
    _Cd = _Cc * (_Ad + Eye);
    _Dd = _Cc * _Bd + _Dc;

    initSSdiscrete(_Ad, _Bd, _Cd, _Dd);
}

void StatespaceSISOSys4th::initTFcontinuous(float *coef_p, float *coef_z, float fs)
{
    for (uint8_t i = 0; i < 5; i++)
    {
        _num[i] = coef_z[i];
        _den[i] = coef_p[i];
    }

    float a[4] = {};
    uint8_t i = 0;
    for (i = 0; i < 4; i++)
    {
        a[i] = coef_p[4 - i] / coef_p[0];
    }

    float b[5] = {};
    for (i = 0; i < 5; i++)
    {
        b[i] = coef_z[4 - i] / coef_p[0];
    }

    float h[5] = {};
    h[0] = b[4];
    h[1] = b[3] - a[3] * h[0];
    h[2] = b[2] - a[3] * h[1] - a[2] * h[0];
    h[3] = b[1] - a[3] * h[2] - a[2] * h[1] - a[1] * h[0];
    h[4] = b[0] - a[3] * h[3] - a[2] * h[2] - a[1] * h[1] - a[0] * h[0];

    const float Ac[] ={0.0f,  1.0f,  0.0f,  0.0f,
                       0.0f,  0.0f,  1.0f,  0.0f,
                       0.0f,  0.0f,  0.0f,  1.0f,
                      -a[0], -a[1], -a[2], -a[3]};

    matrix::SquareMatrix<float, 4> Ac_value(Ac);
    _Ac = Ac_value;

    const float Bc[] = {h[1],
                        h[2],
                        h[3],
                        h[4]};
    matrix::Matrix<float, 4, 1> Bc_value(Bc);
    _Bc = Bc_value;

    _Cc.setZero();
    _Cc(0, 0) = 1.0f;

    _Dc(0, 0) = h[0];

    initSScontinuous(_Ac, _Bc, _Cc, _Dc, fs);
}

void StatespaceSISOSys4th::update(float u)
{
    _u(0, 0) = u;
    _y = _Cd * _x + _Dd * _u;
    _x = _Ad * _x + _Bd * _u;
}

void StatespaceSISOSys4th::switchBuf(float u, float y)
{
    _u(0, 0) = u;
    _y(0, 0) = y;
    _x = _Cd_pinv * (_y - _Dd * _u);
    _x = _Ad * _x + _Bd * _u;
}

float StatespaceSISOSys4th::getOutput()
{
    return _y(0, 0);
}

/* -------------------- for fifth-order system -------------------- */
void StatespaceSISOSys5th::reset()
{
    _u.setZero();
    _y.setZero();
    _x.setZero();
}

void StatespaceSISOSys5th::initSSdiscrete(matrix::SquareMatrix<float, 5> a, matrix::Matrix<float, 5, 1> b, matrix::Matrix<float, 1, 5> c, matrix::Matrix<float, 1, 1> d)
{
    if (!_inited_flag)
    {
        reset();
        _inited_flag = 1;
    }
    _Ad = a;
    _Bd = b;
    _Cd = c;
    _Dd = d;
    float tmp = _Cd(0, 0) * _Cd(0, 0) + _Cd(0, 1) * _Cd(0, 1) + _Cd(0, 2) * _Cd(0, 2) + _Cd(0, 3) * _Cd(0, 3) + _Cd(0, 4) * _Cd(0, 4);
    for (uint8_t i = 0; i < 5; i ++)
    {
        _Cd_pinv(i, 0) = c(0, i) / tmp;
    }
}

void StatespaceSISOSys5th::initSScontinuous(matrix::SquareMatrix<float, 5> a, matrix::Matrix<float, 5, 1> b, matrix::Matrix<float, 1, 5> c, matrix::Matrix<float, 1, 1> d, float fs)
{
    _Ac = a;
    _Bc = b;
    _Cc = c;
    _Dc = d;
    _fs = fs;

    matrix::SquareMatrix<float, 5> Eye_alpha;
    matrix::SquareMatrix<float, 5> Eye;
    Eye.setIdentity();
    Eye_alpha.setIdentity();

    float alpha = 2.0f * _fs;

    for (uint8_t i = 0; i < 5; i++)
    {
        Eye_alpha(i,i) = alpha;
    }

    matrix::SquareMatrix<float, 5> tmp;
    tmp = Eye_alpha - _Ac;

    _Ad = inv(tmp) * (Eye_alpha + _Ac);
    _Bd = inv(tmp) * _Bc;
    _Cd = _Cc * (_Ad + Eye);
    _Dd = _Cc * _Bd + _Dc;

    initSSdiscrete(_Ad, _Bd, _Cd, _Dd);
}

void StatespaceSISOSys5th::initTFcontinuous(float *coef_p, float *coef_z, float fs)
{
    for (uint8_t i = 0; i < 6; i++)
    {
        _num[i] = coef_z[i];
        _den[i] = coef_p[i];
    }

    float a[5] = {};
    uint8_t i = 0;
    for (i = 0; i < 5; i++)
    {
        a[i] = coef_p[5 - i] / coef_p[0];
    }

    float b[6] = {};
    for (i = 0; i < 6; i++)
    {
        b[i] = coef_z[5 - i] / coef_p[0];
    }

    float h[6] = {};
    h[0] = b[5];
    h[1] = b[4] - a[4] * h[0];
    h[2] = b[3] - a[4] * h[1] - a[3] * h[0];
    h[3] = b[2] - a[4] * h[2] - a[3] * h[1] - a[2] * h[0];
    h[4] = b[1] - a[4] * h[3] - a[3] * h[2] - a[2] * h[1] - a[1] * h[0];
    h[5] = b[0] - a[4] * h[4] - a[3] * h[3] - a[2] * h[2] - a[1] * h[1] - a[0] * h[0];

    const float Ac[] ={0.0f,  1.0f,  0.0f,  0.0f,  0.0f,
                       0.0f,  0.0f,  1.0f,  0.0f,  0.0f,
                       0.0f,  0.0f,  0.0f,  1.0f,  0.0f,
                       0.0f,  0.0f,  0.0f,  0.0f,  1.0f,
                      -a[0], -a[1], -a[2], -a[3], -a[4]};

    matrix::SquareMatrix<float, 5> Ac_value(Ac);
    _Ac = Ac_value;

    const float Bc[] = {h[1],
                        h[2],
                        h[3],
                        h[4],
                        h[5]};
    matrix::Matrix<float, 5, 1> Bc_value(Bc);
    _Bc = Bc_value;

    _Cc.setZero();
    _Cc(0, 0) = 1.0f;

    _Dc(0, 0) = h[0];

    initSScontinuous(_Ac, _Bc, _Cc, _Dc, fs);
}

void StatespaceSISOSys5th::update(float u)
{
    _u(0, 0) = u;
    _y = _Cd * _x + _Dd * _u;
    _x = _Ad * _x + _Bd * _u;
}

void StatespaceSISOSys5th::switchBuf(float u, float y)
{
    _u(0, 0) = u;
    _y(0, 0) = y;
    _x = _Cd_pinv * (_y - _Dd * _u);
    _x = _Ad * _x + _Bd * _u;
}

float StatespaceSISOSys5th::getOutput()
{
    return _y(0, 0);
}

/* -------------------- for sixth-order system -------------------- */
void StatespaceSISOSys6th::reset()
{
    _u.setZero();
    _y.setZero();
    _x.setZero();
}

void StatespaceSISOSys6th::initSSdiscrete(matrix::SquareMatrix<float, 6> a, matrix::Matrix<float, 6, 1> b, matrix::Matrix<float, 1, 6> c, matrix::Matrix<float, 1, 1> d)
{
    if (!_inited_flag)
    {
        reset();
        _inited_flag = 1;
    }
    _Ad = a;
    _Bd = b;
    _Cd = c;
    _Dd = d;
    float tmp = _Cd(0, 0) * _Cd(0, 0) + _Cd(0, 1) * _Cd(0, 1) + _Cd(0, 2) * _Cd(0, 2) + _Cd(0, 3) * _Cd(0, 3) + _Cd(0, 4) * _Cd(0, 4) + _Cd(0, 5) * _Cd(0, 5);
    for (uint8_t i = 0; i < 6; i ++)
    {
        _Cd_pinv(i, 0) = c(0, i) / tmp;
    }
}

void StatespaceSISOSys6th::initSScontinuous(matrix::SquareMatrix<float, 6> a, matrix::Matrix<float, 6, 1> b, matrix::Matrix<float, 1, 6> c, matrix::Matrix<float, 1, 1> d, float fs)
{
    _Ac = a;
    _Bc = b;
    _Cc = c;
    _Dc = d;
    _fs = fs;

    matrix::SquareMatrix<float, 6> Eye_alpha;
    matrix::SquareMatrix<float, 6> Eye;
    Eye.setIdentity();
    Eye_alpha.setIdentity();

    float alpha = 2.0f * _fs;

    for (uint8_t i = 0; i < 6; i++)
    {
        Eye_alpha(i,i) = alpha;
    }

    matrix::SquareMatrix<float, 6> tmp;
    tmp = Eye_alpha - _Ac;

    _Ad = inv(tmp) * (Eye_alpha + _Ac);
    _Bd = inv(tmp) * _Bc;
    _Cd = _Cc * (_Ad + Eye);
    _Dd = _Cc * _Bd + _Dc;

    initSSdiscrete(_Ad, _Bd, _Cd, _Dd);
}

void StatespaceSISOSys6th::initTFcontinuous(float *coef_p, float *coef_z, float fs)
{
    for (uint8_t i = 0; i < 7; i++)
    {
        _num[i] = coef_z[i];
        _den[i] = coef_p[i];
    }

    float a[6] = {};
    uint8_t i = 0;
    for (i = 0; i < 6; i++)
    {
        a[i] = coef_p[6 - i] / coef_p[0];
    }

    float b[7] = {};
    for (i = 0; i < 7; i++)
    {
        b[i] = coef_z[6 - i] / coef_p[0];
    }

    float h[7] = {};
    h[0] = b[6];
    h[1] = b[5] - a[5] * h[0];
    h[2] = b[4] - a[5] * h[1] - a[4] * h[0];
    h[3] = b[3] - a[5] * h[2] - a[4] * h[1] - a[3] * h[0];
    h[4] = b[2] - a[5] * h[3] - a[4] * h[2] - a[3] * h[1] - a[2] * h[0];
    h[5] = b[1] - a[5] * h[4] - a[4] * h[3] - a[3] * h[2] - a[2] * h[1] - a[1] * h[0];
    h[6] = b[0] - a[5] * h[5] - a[4] * h[4] - a[3] * h[3] - a[2] * h[2] - a[1] * h[1] - a[0] * h[0];

    const float Ac[] ={0.0f,  1.0f,  0.0f,  0.0f,  0.0f,  0.0f,
                       0.0f,  0.0f,  1.0f,  0.0f,  0.0f,  0.0f,
                       0.0f,  0.0f,  0.0f,  1.0f,  0.0f,  0.0f,
                       0.0f,  0.0f,  0.0f,  0.0f,  1.0f,  0.0f,
                       0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  1.0f,
                      -a[0], -a[1], -a[2], -a[3], -a[4], -a[5]};

    matrix::SquareMatrix<float, 6> Ac_value(Ac);
    _Ac = Ac_value;

    const float Bc[] = {h[1],
                        h[2],
                        h[3],
                        h[4],
                        h[5],
                        h[6]};
    matrix::Matrix<float, 6, 1> Bc_value(Bc);
    _Bc = Bc_value;

    _Cc.setZero();
    _Cc(0, 0) = 1.0f;

    _Dc(0, 0) = h[0];

    initSScontinuous(_Ac, _Bc, _Cc, _Dc, fs);
}

void StatespaceSISOSys6th::update(float u)
{
    _u(0, 0) = u;
    _y = _Cd * _x + _Dd * _u;
    _x = _Ad * _x + _Bd * _u;
}

void StatespaceSISOSys6th::switchBuf(float u, float y)
{
    _u(0, 0) = u;
    _y(0, 0) = y;
    _x = _Cd_pinv * (_y - _Dd * _u);
    _x = _Ad * _x + _Bd * _u;
}

float StatespaceSISOSys6th::getOutput()
{
    return _y(0, 0);
}
