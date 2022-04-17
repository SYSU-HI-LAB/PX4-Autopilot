#ifndef __STATESPACE_SISO_SYS_H__
#define __STATESPACE_SISO_SYS_H__

#include <matrix/matrix/math.hpp>

typedef enum
{
    SISO_SYS_ORDER_ZERO,
    SISO_SYS_ORDER_1ST,
    SISO_SYS_ORDER_2ND,
    SISO_SYS_ORDER_3RD,
    SISO_SYS_ORDER_4TH,
    SISO_SYS_ORDER_5TH,
    SISO_SYS_ORDER_6TH,
    SISO_SYS_ORDER_7TH,
    SISO_SYS_ORDER_8TH,
    SISO_SYS_ORDER_MAX,
} sys_order_t;
class StatespaceSISOSys1st
{
public:
    void reset();
    void update(float u);
    void switchBuf(float u, float y);
    float getOutput();

    void initSSdiscrete(float a, float b, float c, float d);
    void initSScontinuous(float a, float b, float c, float d, float fs);
    void initTFcontinuous(float *coef_p, float *coef_z, float fs);

#ifndef DEBUG
protected:
#endif
    /* contitous system state matrices and transfer function */
    float _Ac;
    float _Bc;
    float _Cc;
    float _Dc;
    float _num[SISO_SYS_ORDER_1ST + 1];
    float _den[SISO_SYS_ORDER_1ST + 1];

    /* discrete state matrices */
    float _Ad;
    float _Bd;
    float _Cd;
    float _Dd;
    float _Cd_pinv;
    float _fs;

    /* discrete states */
    uint8_t _inited_flag;
    float _x;
    float _y;
    float _u;
};

class StatespaceSISOSys2nd
{
public:
    void reset();
    void update(float u);
    void switchBuf(float u, float y);
    float getOutput();

    void initSSdiscrete(matrix::SquareMatrix<float, SISO_SYS_ORDER_2ND> a, matrix::Matrix<float, SISO_SYS_ORDER_2ND, 1> b, matrix::Matrix<float, 1, SISO_SYS_ORDER_2ND> c, matrix::Matrix<float, 1, 1> d);
    void initSScontinuous(matrix::SquareMatrix<float, SISO_SYS_ORDER_2ND> a, matrix::Matrix<float, SISO_SYS_ORDER_2ND, 1> b, matrix::Matrix<float, 1, SISO_SYS_ORDER_2ND> c, matrix::Matrix<float, 1, 1> d, float fs);
    void initTFcontinuous(float *coef_p, float *coef_z, float fs);

#ifndef DEBUG
protected:
#endif
    /* contitous system state matrices and transfer function */
    matrix::SquareMatrix<float, SISO_SYS_ORDER_2ND> _Ac;
    matrix::Matrix<float, SISO_SYS_ORDER_2ND, 1> _Bc;
    matrix::Matrix<float, 1, SISO_SYS_ORDER_2ND> _Cc;
    matrix::Matrix<float, 1, 1> _Dc;
    float _num[SISO_SYS_ORDER_2ND + 1];
    float _den[SISO_SYS_ORDER_2ND + 1];

    /* discrete state matrices */
    matrix::SquareMatrix<float, SISO_SYS_ORDER_2ND> _Ad;
    matrix::Matrix<float, SISO_SYS_ORDER_2ND, 1> _Bd;
    matrix::Matrix<float, 1, SISO_SYS_ORDER_2ND> _Cd;
    matrix::Matrix<float, 1, 1> _Dd;
    matrix::Matrix<float, SISO_SYS_ORDER_2ND, 1> _Cd_pinv;
    float _fs;

    /* discrete states */
    uint8_t _inited_flag;
    matrix::Matrix<float, SISO_SYS_ORDER_2ND, 1> _x;
    matrix::Matrix<float, 1, 1> _y;
    matrix::Matrix<float, 1, 1> _u;
};

class StatespaceSISOSys3rd
{
public:
    void reset();
    void update(float u);
    void switchBuf(float u, float y);
    float getOutput();

    void initSSdiscrete(matrix::SquareMatrix<float, SISO_SYS_ORDER_3RD> a, matrix::Matrix<float, SISO_SYS_ORDER_3RD, 1> b, matrix::Matrix<float, 1, SISO_SYS_ORDER_3RD> c, matrix::Matrix<float, 1, 1> d);
    void initSScontinuous(matrix::SquareMatrix<float, SISO_SYS_ORDER_3RD> a, matrix::Matrix<float, SISO_SYS_ORDER_3RD, 1> b, matrix::Matrix<float, 1, SISO_SYS_ORDER_3RD> c, matrix::Matrix<float, 1, 1> d, float fs);
    void initTFcontinuous(float *coef_p, float *coef_z, float fs);

#ifndef DEBUG
protected:
#endif
    /* contitous system state matrices and transfer function */
    matrix::SquareMatrix<float, SISO_SYS_ORDER_3RD> _Ac;
    matrix::Matrix<float, SISO_SYS_ORDER_3RD, 1> _Bc;
    matrix::Matrix<float, 1, SISO_SYS_ORDER_3RD> _Cc;
    matrix::Matrix<float, 1, 1> _Dc;
    float _num[SISO_SYS_ORDER_3RD + 1];
    float _den[SISO_SYS_ORDER_3RD + 1];

    /* discrete state matrices */
    matrix::SquareMatrix<float, SISO_SYS_ORDER_3RD> _Ad;
    matrix::Matrix<float, SISO_SYS_ORDER_3RD, 1> _Bd;
    matrix::Matrix<float, 1, SISO_SYS_ORDER_3RD> _Cd;
    matrix::Matrix<float, 1, 1> _Dd;
    matrix::Matrix<float, SISO_SYS_ORDER_3RD, 1> _Cd_pinv;
    float _fs;

    /* discrete states */
    uint8_t _inited_flag;
    matrix::Matrix<float, SISO_SYS_ORDER_3RD, 1> _x;
    matrix::Matrix<float, 1, 1> _y;
    matrix::Matrix<float, 1, 1> _u;
};

class StatespaceSISOSys4th
{
public:
    void reset();
    void update(float u);
    void switchBuf(float u, float y);
    float getOutput();

    void initSSdiscrete(matrix::SquareMatrix<float, SISO_SYS_ORDER_4TH> a, matrix::Matrix<float, SISO_SYS_ORDER_4TH, 1> b, matrix::Matrix<float, 1, SISO_SYS_ORDER_4TH> c, matrix::Matrix<float, 1, 1> d);
    void initSScontinuous(matrix::SquareMatrix<float, SISO_SYS_ORDER_4TH> a, matrix::Matrix<float, SISO_SYS_ORDER_4TH, 1> b, matrix::Matrix<float, 1, SISO_SYS_ORDER_4TH> c, matrix::Matrix<float, 1, 1> d, float fs);
    void initTFcontinuous(float *coef_p, float *coef_z, float fs);

#ifndef DEBUG
protected:
#endif
    /* contitous system state matrices and transfer function */
    matrix::SquareMatrix<float, SISO_SYS_ORDER_4TH> _Ac;
    matrix::Matrix<float, SISO_SYS_ORDER_4TH, 1> _Bc;
    matrix::Matrix<float, 1, SISO_SYS_ORDER_4TH> _Cc;
    matrix::Matrix<float, 1, 1> _Dc;
    float _num[SISO_SYS_ORDER_4TH + 1];
    float _den[SISO_SYS_ORDER_4TH + 1];

    /* discrete state matrices */
    matrix::SquareMatrix<float, SISO_SYS_ORDER_4TH> _Ad;
    matrix::Matrix<float, SISO_SYS_ORDER_4TH, 1> _Bd;
    matrix::Matrix<float, 1, SISO_SYS_ORDER_4TH> _Cd;
    matrix::Matrix<float, 1, 1> _Dd;
    matrix::Matrix<float, SISO_SYS_ORDER_4TH, 1> _Cd_pinv;
    float _fs;

    /* discrete states */
    uint8_t _inited_flag;
    matrix::Matrix<float, SISO_SYS_ORDER_4TH, 1> _x;
    matrix::Matrix<float, 1, 1> _y;
    matrix::Matrix<float, 1, 1> _u;

};

class StatespaceSISOSys5th
{
public:
    void reset();
    void update(float u);
    void switchBuf(float u, float y);
    float getOutput();

    void initSSdiscrete(matrix::SquareMatrix<float, SISO_SYS_ORDER_5TH> a, matrix::Matrix<float, SISO_SYS_ORDER_5TH, 1> b, matrix::Matrix<float, 1, SISO_SYS_ORDER_5TH> c, matrix::Matrix<float, 1, 1> d);
    void initSScontinuous(matrix::SquareMatrix<float, SISO_SYS_ORDER_5TH> a, matrix::Matrix<float, SISO_SYS_ORDER_5TH, 1> b, matrix::Matrix<float, 1, SISO_SYS_ORDER_5TH> c, matrix::Matrix<float, 1, 1> d, float fs);
    void initTFcontinuous(float *coef_p, float *coef_z, float fs);

#ifndef DEBUG
protected:
#endif
    /* contitous system state matrices and transfer function */
    matrix::SquareMatrix<float, SISO_SYS_ORDER_5TH> _Ac;
    matrix::Matrix<float, SISO_SYS_ORDER_5TH, 1> _Bc;
    matrix::Matrix<float, 1, SISO_SYS_ORDER_5TH> _Cc;
    matrix::Matrix<float, 1, 1> _Dc;
    float _num[SISO_SYS_ORDER_5TH + 1];
    float _den[SISO_SYS_ORDER_5TH + 1];

    /* discrete state matrices */
    matrix::SquareMatrix<float, SISO_SYS_ORDER_5TH> _Ad;
    matrix::Matrix<float, SISO_SYS_ORDER_5TH, 1> _Bd;
    matrix::Matrix<float, 1, SISO_SYS_ORDER_5TH> _Cd;
    matrix::Matrix<float, 1, 1> _Dd;
    matrix::Matrix<float, SISO_SYS_ORDER_5TH, 1> _Cd_pinv;
    float _fs;

    /* discrete states */
    uint8_t _inited_flag;
    matrix::Matrix<float, SISO_SYS_ORDER_5TH, 1> _x;
    matrix::Matrix<float, 1, 1> _y;
    matrix::Matrix<float, 1, 1> _u;
};

class StatespaceSISOSys6th
{
public:
    void reset();
    void update(float u);
    void switchBuf(float u, float y);
    float getOutput();

    void initSSdiscrete(matrix::SquareMatrix<float, SISO_SYS_ORDER_6TH> a, matrix::Matrix<float, SISO_SYS_ORDER_6TH, 1> b, matrix::Matrix<float, 1, SISO_SYS_ORDER_6TH> c, matrix::Matrix<float, 1, 1> d);
    void initSScontinuous(matrix::SquareMatrix<float, SISO_SYS_ORDER_6TH> a, matrix::Matrix<float, SISO_SYS_ORDER_6TH, 1> b, matrix::Matrix<float, 1, SISO_SYS_ORDER_6TH> c, matrix::Matrix<float, 1, 1> d, float fs);
    void initTFcontinuous(float *coef_p, float *coef_z, float fs);

#ifndef DEBUG
protected:
#endif
    /* contitous system state matrices and transfer function */
    matrix::SquareMatrix<float, SISO_SYS_ORDER_6TH> _Ac;
    matrix::Matrix<float, SISO_SYS_ORDER_6TH, 1> _Bc;
    matrix::Matrix<float, 1, SISO_SYS_ORDER_6TH> _Cc;
    matrix::Matrix<float, 1, 1> _Dc;
    float _num[SISO_SYS_ORDER_6TH + 1];
    float _den[SISO_SYS_ORDER_6TH + 1];

    /* discrete state matrices */
    matrix::SquareMatrix<float, SISO_SYS_ORDER_6TH> _Ad;
    matrix::Matrix<float, SISO_SYS_ORDER_6TH, 1> _Bd;
    matrix::Matrix<float, 1, SISO_SYS_ORDER_6TH> _Cd;
    matrix::Matrix<float, 1, 1> _Dd;
    matrix::Matrix<float, SISO_SYS_ORDER_6TH, 1> _Cd_pinv;
    float _fs;

    /* discrete states */
    uint8_t _inited_flag;
    matrix::Matrix<float, SISO_SYS_ORDER_6TH, 1> _x;
    matrix::Matrix<float, 1, 1> _y;
    matrix::Matrix<float, 1, 1> _u;
};

class StatespaceSISOSys7th
{
public:
    void reset();
    void update(float u);
    void switchBuf(float u, float y);
    float getOutput();

    void initSSdiscrete(matrix::SquareMatrix<float, SISO_SYS_ORDER_7TH> a, matrix::Matrix<float, SISO_SYS_ORDER_7TH, 1> b, matrix::Matrix<float, 1, SISO_SYS_ORDER_7TH> c, matrix::Matrix<float, 1, 1> d);
    void initSScontinuous(matrix::SquareMatrix<float, SISO_SYS_ORDER_7TH> a, matrix::Matrix<float, SISO_SYS_ORDER_7TH, 1> b, matrix::Matrix<float, 1, SISO_SYS_ORDER_7TH> c, matrix::Matrix<float, 1, 1> d, float fs);
    void initTFcontinuous(float *coef_p, float *coef_z, float fs);

#ifndef DEBUG
protected:
#endif
    /* contitous system state matrices and transfer function */
    matrix::SquareMatrix<float, SISO_SYS_ORDER_7TH> _Ac;
    matrix::Matrix<float, SISO_SYS_ORDER_7TH, 1> _Bc;
    matrix::Matrix<float, 1, SISO_SYS_ORDER_7TH> _Cc;
    matrix::Matrix<float, 1, 1> _Dc;
    float _num[SISO_SYS_ORDER_7TH + 1];
    float _den[SISO_SYS_ORDER_7TH + 1];

    /* discrete state matrices */
    matrix::SquareMatrix<float, SISO_SYS_ORDER_7TH> _Ad;
    matrix::Matrix<float, SISO_SYS_ORDER_7TH, 1> _Bd;
    matrix::Matrix<float, 1, SISO_SYS_ORDER_7TH> _Cd;
    matrix::Matrix<float, 1, 1> _Dd;
    matrix::Matrix<float, SISO_SYS_ORDER_7TH, 1> _Cd_pinv;
    float _fs;

    /* discrete states */
    uint8_t _inited_flag;
    matrix::Matrix<float, SISO_SYS_ORDER_7TH, 1> _x;
    matrix::Matrix<float, 1, 1> _y;
    matrix::Matrix<float, 1, 1> _u;
};

class StatespaceSISOSys8th
{
public:
    void reset();
    void update(float u);
    void switchBuf(float u, float y);
    float getOutput();

    void initSSdiscrete(matrix::SquareMatrix<float, SISO_SYS_ORDER_8TH> a, matrix::Matrix<float, SISO_SYS_ORDER_8TH, 1> b, matrix::Matrix<float, 1, SISO_SYS_ORDER_8TH> c, matrix::Matrix<float, 1, 1> d);
    void initSScontinuous(matrix::SquareMatrix<float, SISO_SYS_ORDER_8TH> a, matrix::Matrix<float, SISO_SYS_ORDER_8TH, 1> b, matrix::Matrix<float, 1, SISO_SYS_ORDER_8TH> c, matrix::Matrix<float, 1, 1> d, float fs);
    void initTFcontinuous(float *coef_p, float *coef_z, float fs);

#ifndef DEBUG
protected:
#endif
    /* contitous system state matrices and transfer function */
    matrix::SquareMatrix<float, SISO_SYS_ORDER_8TH> _Ac;
    matrix::Matrix<float, SISO_SYS_ORDER_8TH, 1> _Bc;
    matrix::Matrix<float, 1, SISO_SYS_ORDER_8TH> _Cc;
    matrix::Matrix<float, 1, 1> _Dc;
    float _num[SISO_SYS_ORDER_8TH + 1];
    float _den[SISO_SYS_ORDER_8TH + 1];

    /* discrete state matrices */
    matrix::SquareMatrix<float, SISO_SYS_ORDER_8TH> _Ad;
    matrix::Matrix<float, SISO_SYS_ORDER_8TH, 1> _Bd;
    matrix::Matrix<float, 1, SISO_SYS_ORDER_8TH> _Cd;
    matrix::Matrix<float, 1, 1> _Dd;
    matrix::Matrix<float, SISO_SYS_ORDER_8TH, 1> _Cd_pinv;
    float _fs;

    /* discrete states */
    uint8_t _inited_flag;
    matrix::Matrix<float, SISO_SYS_ORDER_8TH, 1> _x;
    matrix::Matrix<float, 1, 1> _y;
    matrix::Matrix<float, 1, 1> _u;
};

#endif
