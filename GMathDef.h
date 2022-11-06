#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib> 

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>


#ifndef _MATH_DEFINES_DEFINED
#define _MATH_DEFINES_DEFINED

#ifndef M_E
#define M_E        2.71828182845904523536
#endif

#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340736
#endif

#ifndef M_LOG10E
#define M_LOG10E   0.434294481903251827651
#endif

#ifndef M_LN2
#define M_LN2      0.693147180559945309417
#endif

#ifndef M_LN10
#define M_LN10     2.30258509299404568402
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923
#endif

#ifndef M_PI_4
#define M_PI_4     0.785398163397448309616
#endif

#ifndef M_1_PI
#define M_1_PI     0.318309886183790671538
#endif

#ifndef M_2_PI
#define M_2_PI     0.636619772367581343076
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257390
#endif

#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2  0.707106781186547524401
#endif

#endif //_MATH_DEFINES_DEFINED

#define M_3PI_2    4.71238898038468985769
#define M_2PI      6.28318530717958647692
#define M_PI_6     0.523598775598298873076
#define M_PI_12    0.261799387799149436538
#define M_PI_180   0.017453292519943295
#define M_180_PI   57.29577951308232
#define M_1_2PI    0.159154943091895335769
#define M_SQRT3    1.7320508075688772
#define M_GOLDEN   0.6180339887498949
#define SIN_5_DEGREE 0.08715574274765817355806427083747

#define MaxLoopCount        1000 
#define MaxIterativeCount   32   
#define MaxNewtonIterativeCount 10 

#define MinParamV            0.1   
#define MaxParamV            10.0  

union MSVC_EVIL_FLOAT_HACK
{
    unsigned char Bytes[8];
    double Value;
};

static union MSVC_EVIL_FLOAT_HACK INFINITY_HACK = {{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xF0, 0x7F}};

static union MSVC_EVIL_FLOAT_HACK NAN_HACK = {{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xF8, 0x7F}};


const double g_DoubleResolution = 1E-12;

const double g_RelaxedDoubleResolution = 1E-8;

const double g_SingleResolution = 1E-6;

const double g_MinDistEpsilon = 1E-4;

const double g_MaxDistEpsilon = 1.0;


const int g_MinNum = -2147483647;

const int g_MaxNum = 2147483647;


const double g_Infinity = 1E100;

const double g_NegInfinity = -1E100;

const double g_Nan = NAN_HACK.Value;

template <typename Type>
inline Type Min(const Type &a, const Type &b)
{
    return a < b ? a : b;
}

template <typename Type>
inline Type Max(const Type &a, const Type &b)
{
    //return std::max(a,b);
    return a > b ? a : b;
}

template <typename Type>
inline Type Min3(const Type &a, const Type &b, const Type &c)
{
    return Min(Min(a, b), c);
}
    
template <typename Type>
inline Type Min3(const Type &a, const Type &b, const Type &c, int &min_index)
{
    min_index = 0;

    auto min_v = a;

    if(b < min_v)
    {
        min_index = 1;
        min_v = b;
    }

    if(c < min_v)
    {
        min_index = 2;
        min_v = c;
    }

    return min_v;
}

template <typename Type>
inline Type Max3(const Type &a, const Type &b, const Type &c)
{
    return Max(Max(a, b), c);
}

#define maxValue(d1, d2) ((d1) > (d2) ? (d1) : (d2))
#define minValue(d1, d2) ((d1) < (d2) ? (d1) : (d2))

inline int getSign(double dVal)
{
    return dVal < 0 ? -1 : 1;
}

inline int getSign(double dVal, double dEpsilon)
{
    return dVal < -dEpsilon ? -1 : (dVal > dEpsilon ? 1 : 0);
}

inline bool sameValue(double dVal1, double dVal2, double dEpsilon)
{
    if (dEpsilon == 0.0)
    {
        dEpsilon = dVal1 * g_DoubleResolution;
        if (dEpsilon < 0.0)
        {
            dEpsilon = -dEpsilon;
        }
        if (dEpsilon < g_DoubleResolution)
        {
            dEpsilon = g_DoubleResolution;
        }
    }

    return dVal1 > dVal2
        ? dVal1 - dVal2 <= dEpsilon
        : dVal2 - dVal1 <= dEpsilon;
}

enum ValueRelationship 
{
    vrGreaterThan =  1,
    vrEqual       =  0,
    vrLessThan    = -1
};

inline ValueRelationship compareValue(double dVal1, double dVal2, double dEpsilon)
{
    if (dEpsilon == 0.0)
    {
        dEpsilon = dVal1 * g_DoubleResolution;
        if (dEpsilon < 0.0)
        {
            dEpsilon = -dEpsilon;
        }
        if (dEpsilon < g_DoubleResolution)
        {
            dEpsilon = g_DoubleResolution;
        }
    }

    return dVal1 > dVal2
        ? (dVal1 - dVal2 <= dEpsilon ? vrEqual : vrGreaterThan)
        : (dVal2 - dVal1 <= dEpsilon ? vrEqual : vrLessThan);
}

inline bool inRange(double dVal, double dMin, double dMax, double dEpsilon, bool bClosed)
{
    return bClosed
        ? dVal - dMin >= -dEpsilon && dVal - dMax <= dEpsilon
        : dVal - dMin >= dEpsilon && dVal - dMax <= -dEpsilon;
}

inline bool isZero(double dVal, double dEpsilon)
{
    return dVal <= dEpsilon && dVal >= -dEpsilon;
}

inline ValueRelationship compareZero(double dVal, double dEpsilon)
{
    return dVal > 0.0
        ? (dVal <= dEpsilon ? vrEqual : vrGreaterThan)
        : (dVal >= -dEpsilon ? vrEqual : vrLessThan);
}

inline double normalAngle(double dAngle)
{
    int n = static_cast<int>(dAngle * M_1_2PI);
    dAngle -= n * M_2PI;
    if (dAngle < 0.0)
    {
        dAngle = dAngle + M_2PI;
    }
    return dAngle;
}

inline double asin_safe(double dSin)
{
    return dSin > 1.0 ? M_PI_2 : (dSin < -1.0 ? -M_PI_2 : asin(dSin));
}

inline double pow_safe(double x, double y)
{
    return x < 0.0 ? -pow(-x, y) : pow(x, y);
}

inline double sqrt_safe(double x)
{
    return x < 0.0 ? 0.0 : sqrt(x);
}

inline double acos_safe(double dCos)
{
    assert(fabs(dCos) < 1.0 + 100.0 * HUGE_VAL);
    return dCos > 1.0 ? 0.0 : (dCos < -1.0 ? M_PI : acos(dCos));
}

enum ValuePosition 
{
    vpUnknown,
    vpIn,
    vpOn,
    vpOut
};

template <typename Type>
inline void Clamp(const Type &min, Type &value, const Type &max)
{
    value = (value < min) ? min : (value > max) ? max : value;
}

template <typename Type>
Type Round(Type value)
{
    return value > 0.4
        ? (Type)::floor(value + 0.5)
        : (value < -0.4
        ? -(Type)::floor(-value + 0.5)
        : (Type) 0.0);
}

template <typename Type>
Type RoundTo(Type x, int exponent)
{
    Type power = pow((Type) 10.0, -exponent);

    return x > 0.0
        ? floor(x * power + 0.5) / power
        : (x < 0.0
        ? ceil(x * power - 0.5) / power
        : x);
}

template <typename Type>
inline Type Rand(const Type min, const Type max)
{
    return (std::rand() * (max - min) / (Type) RAND_MAX) + min;
}

template <typename Type>
inline int Sign(Type a)
{
    return (Type(0) <= a) - (a < Type(0));
}

template <> inline int Sign<int>(int a)
{
    return a >= 0 ? 1 : -1;
}

template <typename Type>
inline int Sign(Type a, Type epsilon)
{
    return (epsilon < a) - (a < -epsilon);
}

template <typename Type>
inline bool Equals(const Type &a, const Type &b, Type eps = std::numeric_limits<Type>::epsilon())
{
    return std::abs(a - b) <= eps;
}

/*!
*@brief     根据角度计算弧度
*@param[in] a 角度
*@return    弧度
*/
template <typename Type>
inline Type DegToRad(Type a)
{ 
    return (Type) (a * M_PI_180); 
}

/*!
*@brief     根据弧度计算角度
*@param[in] a 弧度
*@return    角度
*/
template <typename Type>
inline Type RadToDeg(Type a)
{ 
    return (Type) (a * M_180_PI); 
}

/*!
* @brief     判断value是否是2的幂
* @param[in] value 给定的幂值
*@return    是否是2的幂
* - true 是
* - false 不是
*/
inline bool IsPowerOfTwo(int value)
{
    return (value & (value - 1)) == 0;
}

/*!
* @brief     判断数值是否在给定误差接近零值
* @param[in] typeData 数值
* @param[in] eps      误差
*@return    是否接近零值
* - true 是
* - false 不是
*/
template <typename Type>
inline bool IsNearZero(const Type &typeData, Type eps = std::numeric_limits<Type>::epsilon())
{
    return std::abs(typeData) < eps;
}

/*!
* @brief 两个浮点数的大于条件比较
* @param[in] value1    数值1
* @param[in] value2    数值2
* @param[in] tolerance 误差
*@return    大小关系
* - true 数值1大于数值2
* - false 数值2大于数值1
*/
inline bool IsGreaterThan(double value1, double value2, double tolerance)
{
    return value1 - value2 > tolerance;
}

/*!
* @brief     两个浮点数的大于等于条件比较
* @param[in] value1    数值1
* @param[in] value2    数值2
* @param[in] tolerance 误差
* @return    大小关系
* - true 数值1大于等于数值2
* - false 数值2大于等于数值1
*/
inline bool IsGreaterEqualThan(double value1, double value2, double tolerance)
{
    return tolerance >= value2 - value1;
}

/*!
* @brief     两个浮点数的小于条件比较
* @param[in] value1    数值1
* @param[in] value2    数值2
* @param[in] tolerance 误差
* @return    大小关系
* - true 数值1小于于数值2
* - false 数值2小于于数值1
*/
inline bool IsLessThan(double value1, double value2, double tolerance)
{
    //      return value1 < value2 - tolerance;    //“12371.862901728357 - 1e-12 = 12371.862901728355”
    return tolerance < value2 - value1;    //“12371.862901728355 < 12371.862901728357 - 1e-12”不成立，
                                            //“1e-12 < 12371.862901728357 - 12371.862901728355”成立。极大数减极小数误差
}

/*!
* @brief     两个浮点数的小于等于条件比较
* @param[in] value1    数值1
* @param[in] value2    数值2
* @param[in] tolerance 误差
*@return    大小关系
* - true 数值1小于等于数值2
* - false 数值2小于等于数值1
*/
inline bool IsLessEqualThan(double value1, double value2, double tolerance)
{
    return value1 - value2 <= tolerance;
}

/*!
* @brief     判断数值 dVal 是否在区间（dMin,dMax）的开区间或闭区间内
* @param[in] dVal    数值
* @param[in] dMin    区间最小值
* @param[in] dMax    区间最大值
* @param[in] eps     误差
* @param[in] bClosed true:  闭区间内判断， false: 开区间内判断
* @return    是否在区间内
* - true 是
* - false 不是
*/
template <typename Type>
inline bool IsInRange(Type dVal, Type dMin, Type dMax,Type eps = std::numeric_limits<Type>::epsilon(), bool bClosed = false) 
{
    return bClosed
        ? dVal - dMin >= -eps && dVal - dMax <= eps
        : dVal - dMin >= eps && dVal - dMax <= -eps;
}

/*!
* @brief     将弧度值限制在 [0,M_2PI) 之间
* @param[in] dAngle 弧度值
* @return    [0,M_2PI) 之间的弧度值
*/
template <typename Type>
inline Type NormalAngle(Type dAngle) 
{
    double dAngInput = double(dAngle);
    int n = int(dAngle / M_2PI);
    dAngInput -= n * M_2PI;
    if(dAngInput < 0.0) 
    {
        dAngInput += M_2PI;
    }
    return Type(dAngInput);
}

template <typename Type>
Type CircleAngleDistance(Type a1, Type a2)
{
    Type a = a1 - a2;
    return a - (Type) M_2PI * floor((Type) M_1_2PI * (a + (Type) M_PI));
}

template <typename Type>
Type HalfCycleInterval(Type v, Type cycle)
{
    return v - cycle * floor(v / cycle + 0.5);
}

template <typename Type>
Type PiInterval(Type angle)
{
    return angle - M_2PI * floor(angle * M_1_2PI + 0.5);
}

template <typename Type>
Type CycleInterval(Type v, Type cycle)
{
    return v - cycle * floor(v / cycle);
}

template <typename Type>
inline Type AsinSafe(Type dSin) 
{
    Clamp((Type) -1.0, dSin, (Type) 1.0);
    return asin(dSin);
}
    
template <typename Type>
inline Type AcosSafe(Type dCos)
{
    Clamp((Type) -1.0, dCos, (Type) 1.0);
    return acos(dCos);
}

template <typename Type>
inline Type PowSafe(Type x, Type y) 
{
    return x < (Type) 0 ? -pow(-x, y) : pow(x, y);
}

template <typename Type>
inline Type SqrtSafe(Type x) 
{
    return x < (Type) 0 ? (Type) 0 : sqrt(x);
}

#ifdef _WIN32
inline void sincos(double x, double *sinVal, double *cosVal)
{
#if defined(_M_IX86) && defined(_MSC_VER)
    __asm
    {
        fld qword ptr[x]
        fsincos
        mov edx,[cosVal]
        mov eax,[sinVal]
        fstp qword ptr[edx]
        fstp qword ptr[eax]
    }
#else
    *sinVal = sin(x);
    *cosVal = cos(x);
#endif
}
#endif

inline double cosapprox(double x)
{
    return 1.0-0.5 * x * x;
}

inline float cosapprox(float x)
{
    return 1.0f - 0.5f * x * x;
}

inline double rsqrt(double x)
{
#if defined(_M_IX86) && defined(_MSC_VER)
    return 1.0/sqrt(x);
#else
    return 1.0 / sqrt(x);
#endif
}

inline double fastlog2(double x)
{
    int e = ( ((int32_t*)&x)[1]>>20 ) - 0x3ff;
    ((int32_t*)&x)[1] = 0x3ff<<20 | (0x0fffff&((int*)&x)[1]); //x指数变为0，成为区间[1,2)的一个数

    x -= 1.5; //x变为[-0.5, 0.5) 绝对值更小，展开式收敛快

    return 0.5849625007211562 + 0.9617966939259756*x*(1.0-0.2222222222222222*x) + e; //log2(1.5+x) 展开式 + e
}
    
inline int log2(int x)
{
    //int v; // 32-bit integer to find the log base 2 of
    int r; // result of log_2(v) goes here
    union { unsigned int u[2]; double d; } t; // temp

    t.u[1] = 0x43300000;
    t.u[0] = x;
    t.d -= 4503599627370496.0;
    r = (t.u[1] >> 20) - 0x3ff;

    return r;
}


template <int iter_count>
double fast_rsqrt(double x);

template <>
inline double fast_rsqrt<1>(double x)
{
    int64_t i;
    double y;

    y = x;
    i = *(int64_t*) &y;
    i = 0x5fe6eb50c7b537a9 - (i >> 1);
    y = *(double*) &i;
    y = y * (1.5 - (0.5 * x * y * y));
    return y;
}

template <>
inline double fast_rsqrt<2>(double x)
{
    int64_t i;
    double y;

    y = x;
    i = *(int64_t*) &y;
    i = 0x5fe6eb50c7b537a9 - (i >> 1);
    y = *(double*) &i;
    y = y * (1.5 - (0.5 * x * y * y));
    y = y * (1.5 - (0.5 * x * y * y));
    return y;
}

template <typename T>
bool is_infinite(const T &value)
{    
#undef max
    T max_value = std::numeric_limits<T>::max();
    T min_value = -max_value;     
    return !(min_value <= value && value <= max_value);
} 

template <typename T>
bool is_nan(const T &value)
{
    return value != value; // True if NAN
} 

template <typename T>
bool is_valid(const T &value)
{    
    return !is_infinite(value) && !is_nan(value);
}

template <typename T>
T relative_error(T v1, T v2)
{
    auto v = v1;
    if (v1 == 0.0)
    {
        v = 1.0;
    }
    return std::abs((v1 - v2)/v);
}

template <typename T>
int round_int(T x)
{
    return (int) ((x > 0) ? (x + (T) 0.5) : (x - (T) 0.5)); 
}

template <typename T>
int64_t round_int64(T x)
{
    return (int64_t)((x > 0) ? (x + (T) 0.5) : (x - (T) 0.5));
}

template<typename T>
T round(T x)
{
    return std::floor(x + (T)0.5);
}

enum UnivariateFunctionPointAttribute
{
        
    PointMinimum    = 1,
        
    PointInflection = 2,
        
    PointMaximum    = 4,
        
    PointInterval   = 8,
        
    PointBreak      = 16,
};
    
struct ExtremumPointInfo
{
    double t;

    double f;

    UnivariateFunctionPointAttribute point_attributes;
};
