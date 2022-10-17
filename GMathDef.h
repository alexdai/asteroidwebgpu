/*!
* @file      GMathDef.h
* @brief     基本数学常用函数和常量定义
* @details
* @author
* @date      2021/09/10
* @copyright Copyright 2012-2022 GLODON
**********************************************************************************
* @par 修改日志:
* |  Date  | Description | Author     |
* | :----: | :----       | :----      |
* | 2021/09/10 | 袁冠杰<yuangj-a> | GGP-55193 完成算法剩余部分的注释 |
*
**********************************************************************************
*/
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

    /* <math.h>中定义的有用的数学常量
    * M_E        - e
    * M_LOG2E    - log2(e)
    * M_LOG10E   - log10(e)
    * M_LN2      - ln(2)
    * M_LN10     - ln(10)
    * M_PI       - pi
    * M_PI_2     - pi/2
    * M_PI_4     - pi/4
    * M_1_PI     - 1/pi
    * M_2_PI     - 2/pi
    * M_2_SQRTPI - 2/sqrt(pi)
    * M_SQRT2    - sqrt(2)
    * M_SQRT1_2  - 1/sqrt(2) */

    //数学常数(自然对数的底数)
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

    /* 补充常量
    * M_3PI_2    - 3pi/2
    * M_PI_6     - pi/6
    * M_PI_12    - pi/12
    * M_1_2PI    - 1/2pi
    * M_2_SQRTPI - 2/sqrt(pi)
    * M_SQRT3    - sqrt(3)
    * EPS        - smallest such that 1.0+EPS != 1.0 */

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
//#define EPS        std::numeric_limits<Type>::epsilon()

#define MaxLoopCount        1000    //最大循环次数
#define MaxIterativeCount   32      //最大迭代次数
#define MaxNewtonIterativeCount 10  //最大牛顿迭代次数

#define MinParamV            0.1     //最小合理参数速度
#define MaxParamV            10.0    //最大合理参数速度

    /*
    * @union MSVC_EVIL_FLOAT_HACK
    * @brief 用于定义无穷大或无穷小值和无效值
    */
    union MSVC_EVIL_FLOAT_HACK
    {
        unsigned char Bytes[8];
        double Value;
    };
    //! 无穷大
    static union MSVC_EVIL_FLOAT_HACK INFINITY_HACK = {{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xF0, 0x7F}};
    //! 无效值
    static union MSVC_EVIL_FLOAT_HACK NAN_HACK = {{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xF8, 0x7F}};

    //! 全局的双精度浮点数精度
    const double g_DoubleResolution = 1E-12;
    //! 全局的宽松浮点数精度
    const double g_RelaxedDoubleResolution = 1E-8;
    //! 全局的单精度浮点数精度
    const double g_SingleResolution = 1E-6;
    //! 最小距离误差
    const double g_MinDistEpsilon = 1E-4;
    //! 最大距离误差
    const double g_MaxDistEpsilon = 1.0;

    //! 全局的最小整数值
    const int g_MinNum = -2147483647;
    //! 全局的最大整数值
    const int g_MaxNum = 2147483647;

    /// 全局的正无穷大值
    const double g_Infinity = 1E100;
    /// 全局的负无穷大值
    const double g_NegInfinity = -1E100;

    /// 全局的无效值
    const double g_Nan = NAN_HACK.Value;

    /*!
    * @brief     求两者数值最小值
    * @param[in] a 数值1
    * @param[in] b 数值2
    * @return    最小值
    */
    template <typename Type>
    inline Type Min(const Type &a, const Type &b)
    {
        //return std::min(a,b);
        return a < b ? a : b;
    }

    /*!
    * @brief     求两者数值最大值
    * @param[in] a 数值1
    * @param[in] b 数值2
    * @return    最大值
    */
    template <typename Type>
    inline Type Max(const Type &a, const Type &b)
    {
        //return std::max(a,b);
        return a > b ? a : b;
    }

    /*!
    * @brief     求三者数值最小值
    * @param[in] a 数值1
    * @param[in] b 数值2
    * @param[in] c 数值3
    * @return    最小值
    */
    template <typename Type>
    inline Type Min3(const Type &a, const Type &b, const Type &c)
    {
        //return std::min(std::Min(a, b), c);
        return Min(Min(a, b), c);
    }
    
    /*!
    * @brief      求三者数值最小值，并获得最小值的序号
    * @param[in]  a         数值1
    * @param[in]  b         数值2
    * @param[in]  c         数值3
    * @param[out] min_index 最小值的序号
    * @return     最小值
    */
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

    /*!
    * @brief     求三者数值最大值
    * @param[in] a 数值1
    * @param[in] b 数值2
    * @param[in] c 数值3
    * @return    最大值
    */
    template <typename Type>
    inline Type Max3(const Type &a, const Type &b, const Type &c)
    {
        //return std::max(std::Max(a, b), c);
        return Max(Max(a, b), c);
    }

    //! 求 d1, d2 两个数的最大值
#define maxValue(d1, d2) ((d1) > (d2) ? (d1) : (d2))
    //! 求d1,d2两个数的最小值
#define minValue(d1, d2) ((d1) < (d2) ? (d1) : (d2))

    /*!
    * @brief     得到浮点数 dVal 的正负，只返回 -1 或 1
    * @param[in] dVal 浮点数
    * @return    -1：负，1：正
    */
    inline int getSign(double dVal)
    {
        return dVal < 0 ? -1 : 1;
    }

    /*!
    * @brief     根据传入误差判断，判断浮点数dVal为正，负或 0
    * @param[in] dVal     浮点数
    * @param[in] dEpsilon 误差
    * @return    -1：负，1：正，0: 0
    */
    inline int getSign(double dVal, double dEpsilon)
    {
        return dVal < -dEpsilon ? -1 : (dVal > dEpsilon ? 1 : 0);
    }

    /*!
    *@brief     两个数值是否相等
    *@param[in] dVal1    数值1
    *@param[in] dVal2    数值2
    *@param[in] dEpsilon 误差
    *@return    是否相等
    * - true 相等
    * - false 不相等
    */ 
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

    /*!    
    * @enum  ValueRelationship
    * @brief 数值间的关系
    */
    enum ValueRelationship 
    {
        //! 大于
        vrGreaterThan =  1,
        //! 等于
        vrEqual       =  0,
        //! 小于
        vrLessThan    = -1
    };

    /*!
    * @brief     比较两个数值
    * @param[in] dVal1    数值 1
    * @param[in] dVal2    数值 2
    * @param[in] dEpsilon 误差
    * @return    数值 1 和数值 2 的大小关系
    */
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

    /*!
    * @brief     判断数值dVal 是否在区间（dMin,dMax）的开区间或闭区间内
    * @param[in] dVal     数值
    * @param[in] dMin     区间最小值
    * @param[in] dMax     期间最大值
    * @param[in] dEpsilon 误差
    * @param[in] bClosed  true:  闭区间内判断， false: 开区间内判断
    *@return    是否在区间内
    * - true 在
    * - false 不在
    */
    inline bool inRange(double dVal, double dMin, double dMax, double dEpsilon, bool bClosed)
    {
        return bClosed
            ? dVal - dMin >= -dEpsilon && dVal - dMax <= dEpsilon
            : dVal - dMin >= dEpsilon && dVal - dMax <= -dEpsilon;
    }

    /*!
    * @brief     判断数值dVal在给定误差内，是否为0
    * @param[in] dVal     数值
    * @param[in] dEpsilon 误差
    *@return    在给定误差内是否为0
    * - true 是
    * - false 不是
    */
    inline bool isZero(double dVal, double dEpsilon)
    {
        return dVal <= dEpsilon && dVal >= -dEpsilon;
    }

    /*!
    * @brief     判断数值dVal和0的大小关系
    * @param[in] dVal     数值
    * @param[in] dEpsilon 误差
    * @return    dVal和0的大小关系
    */
    inline ValueRelationship compareZero(double dVal, double dEpsilon)
    {
        return dVal > 0.0
            ? (dVal <= dEpsilon ? vrEqual : vrGreaterThan)
            : (dVal >= -dEpsilon ? vrEqual : vrLessThan);
    }

    /*!
    * @brief     把弧度 dAngle，转为 0 到 2PI 之间
    * @param[in] dAngle 传入弧度
    * @return    0 到 2PI 之间弧度
    */
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

    /*!
    * @brief     求dSin得反正弦，此函数内部会对传入值是否合法进行判断，调用者不用判断
    * @param[in] dSin 传入值
    * @return    弧度
    */
    inline double asin_safe(double dSin)
    {
        return dSin > 1.0 ? M_PI_2 : (dSin < -1.0 ? -M_PI_2 : asin(dSin));
    }

    /*!
    * @brief     计算 x 的 y次幂，此函数内部会对传入值是否合法进行判断，调用者不用判断
    * @param[in] x 底数
    * @param[in] y 指数
    * @return    结果
    */
    inline double pow_safe(double x, double y)
    {
        return x < 0.0 ? -pow(-x, y) : pow(x, y);
    }

    /*!
    * @brief     计算x的开方，此函数内部会对传入值是否合法进行判断，调用者不用判断
    * @param[in] x 需要进行开方的值
    * @return    开方结果
    */
    inline double sqrt_safe(double x)
    {
        return x < 0.0 ? 0.0 : sqrt(x);
    }

    /*!
    * @brief     求dCos的反余弦，此函数内部会对传入值是否合法进行判断，调用者不用判断
    * @param[in] dCos 传入值
    * @return    弧度
    */
    inline double acos_safe(double dCos)
    {
        assert(fabs(dCos) < 1.0 + 100.0 * HUGE_VAL);
        return dCos > 1.0 ? 0.0 : (dCos < -1.0 ? M_PI : acos(dCos));
    }

    /*!
    * @enum  ValuePosition
    * @brief 值和区间的位置关系
    */
    enum ValuePosition 
    {
        //! 未知
        vpUnknown,
        //! 在区间内
        vpIn,
        //! 在区间边界上
        vpOn,
        //! 在区间外
        vpOut
    };

    /*!
    * @brief     限制数值在min和max之间
    * @param[in] min   区间最小值
    * @param[in] value 数值
    * @param[in] max   区间最大值
    */
    template <typename Type>
    inline void Clamp(const Type &min, Type &value, const Type &max)
    {
        value = (value < min) ? min : (value > max) ? max : value;
    }

    /*!
    * @brief     数值取整
    * @param[in] value 数值
    * @return    整数值
    */
    template <typename Type>
    Type Round(Type value)
    {
        return value > 0.4
            ? (Type)::floor(value + 0.5)
            : (value < -0.4
            ? -(Type)::floor(-value + 0.5)
            : (Type) 0.0);
    }

    /*!
    * @brief     四舍五入到 10^exponent，例如 exponent=2 则近似到 1e2，exponent=-3 近似到 1e-3
    * @param[in] x        传入值
    * @param[in] exponent 10的指数
    * @return    四舍五入后的结果
    */
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

    /*!
    * @brief     获取随机数值，要求处于 min 和 max 之间
    * @param[in] min 区间最小值
    * @param[in] max 区间最大值
    * @return    区间中随机数
    */
    template <typename Type>
    inline Type Rand(const Type min, const Type max)
    {
        return (std::rand() * (max - min) / (Type) RAND_MAX) + min;
    }

    /*!
    * @brief     获取数值的符号
    * @param[in] a 数值
    * @return    1:  正号，-1： 负号
    */
    template <typename Type>
    inline int Sign(Type a)
    {
        return (Type(0) <= a) - (a < Type(0));
    }

    /*!
    * @brief     获取数值的符号
    * @param[in] a 数值
    * @return    1:  正号，-1： 负号
    */
    template <> inline int Sign<int>(int a)
    {
        return a >= 0 ? 1 : -1;
    }

    /*!
    * @brief     在指定误差范围内获取数值的符号
    * @param[in] a 数值
    * @param[in] epsilon 给定误差
    * @return    1:正号，-1:负号，0:0
    */
    template <typename Type>
    inline int Sign(Type a, Type epsilon)
    {
        return (epsilon < a) - (a < -epsilon);
    }

    /*!
    *@brief     判断数值是否在给定误差下相等
    *@param[in] a 数值1
    *@param[in] b 数值2
    *@param[in] eps 容差
    *@return    是否相等
    * - true 相等
    * - false 不相等
    */
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

    /*!
    * @brief     圆周上的两个角的距离（绝对值不大于Pi)
    * @param[in] a1 角1
    * @param[in] a2 角2
    * @return    从角2绕向角1的最小的角距离，符号表示转向，逆时针为正，顺时针为负 
    */
    template <typename Type>
    Type CircleAngleDistance(Type a1, Type a2)
    {
        Type a = a1 - a2;
        return a - (Type) M_2PI * floor((Type) M_1_2PI * (a + (Type) M_PI));
    }

    /*!
    * @brief     把一个值转换到 [-cycle/2, cycle/2] 区间内
    * @param[in] v     给定值
    * @param[in] cycle 周期
    * @return    转换后的值
    */
    template <typename Type>
    Type HalfCycleInterval(Type v, Type cycle)
    {
        return v - cycle * floor(v / cycle + 0.5);
    }

    /*!
    * @brief     把一个弧度值转换到 [-Pi, Pi) 区间内
    * @param[in] angle 给定弧度值
    * @return    转换后的弧度值
    */
    template <typename Type>
    Type PiInterval(Type angle)
    {
        return angle - M_2PI * floor(angle * M_1_2PI + 0.5);
    }

    /*!
    * @brief     把一个值转换到 [0, cycle] 区间内
    * @param[in] v     给定值
    * @param[in] cycle 周期
    * @return    转换后的值
    */
    template <typename Type>
    Type CycleInterval(Type v, Type cycle)
    {
        return v - cycle * floor(v / cycle);
    }

    /*!
    * @brief     安全的反正弦函数
    * @param[in] dSin 正弦值
    * @return    结果
    */
    template <typename Type>
    inline Type AsinSafe(Type dSin) 
    {
        Clamp((Type) -1.0, dSin, (Type) 1.0);
        return asin(dSin);
    }
    
    /*!
    * @brief     安全的反余弦函数
    * @param[in] dCos 余弦值
    * @return    结果
    */
    template <typename Type>
    inline Type AcosSafe(Type dCos)
    {
        Clamp((Type) -1.0, dCos, (Type) 1.0);
        return acos(dCos);
    }

    /*!
    * @brief     安全的幂函数
    * @param[in] x 底数
    * @param[in] y 指数
    * @return    结果
    */
    template <typename Type>
    inline Type PowSafe(Type x, Type y) 
    {
        return x < (Type) 0 ? -pow(-x, y) : pow(x, y);
    }

    /*!
    *@brief      安全二次开根号函数
    * @param[in] x 被开方数
    *@return     结果
    */
    template <typename Type>
    inline Type SqrtSafe(Type x) 
    {
        return x < (Type) 0 ? (Type) 0 : sqrt(x);
    }

#ifdef _WIN32
    /*!
    * @brief      同时计算正弦余弦
    * @param[in]  x      传入角度
    * @param[out] sinVal 正弦值
    * @param[out] cosVal 余弦值
    */
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

    /*!
    * @brief     计算接近于 0 的小角度的近似余弦值
    * @param[in] x 角度
    * @return    余弦值
    */
    inline double cosapprox(double x)
    {
        return 1.0-0.5 * x * x;
    }

    /*!
    * @brief     计算接近于0的小角度的近似余弦值
    * @param[in] x 角度
    * @return    余弦值
    */
    inline float cosapprox(float x)
    {
        return 1.0f - 0.5f * x * x;
    }

    /*!
    * @brief     x 的二次方根的倒数
    * @param[in] x 传入值
    * @return    结果
    */
    inline double rsqrt(double x)
    {
#if defined(_M_IX86) && defined(_MSC_VER)
        return 1.0/sqrt(x);
#else
        return 1.0 / sqrt(x);
#endif
    }

    /*!
    * @brief     快速但不精确的以 2 为底的对数
    * @param[in] x 真数
    * @return    对数
    */
    inline double fastlog2(double x)
    {
        int e = ( ((int32_t*)&x)[1]>>20 ) - 0x3ff;
        ((int32_t*)&x)[1] = 0x3ff<<20 | (0x0fffff&((int*)&x)[1]); //x指数变为0，成为区间[1,2)的一个数

        x -= 1.5; //x变为[-0.5, 0.5) 绝对值更小，展开式收敛快

        return 0.5849625007211562 + 0.9617966939259756*x*(1.0-0.2222222222222222*x) + e; //log2(1.5+x) 展开式 + e
    }
    
    /*!
    * @brief     快速但不精确的以 2 为底的对数
    * @param[in] x 真数
    * @return    对数
    */
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


    /*!
    * @brief     快速但不精确的平方根倒数
    *            http://en.wikipedia.org/wiki/Fast_inverse_square_root
    * @param[in] x 数
    * @return    平方根的倒数
    */
    template <int iter_count>
    double fast_rsqrt(double x);

    /*!
    * @brief     快速但不精确的平方根倒数
    *            http://en.wikipedia.org/wiki/Fast_inverse_square_root
    * @param[in] x 数
    * @return    平方根的倒数
    */
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

    /*!
    * @brief     快速但不精确的平方根倒数
    *            http://en.wikipedia.org/wiki/Fast_inverse_square_root
    * @param[in] x 数
    * @return    平方根的倒数
    */
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

    /*!
    *@brief 判断给定数是否是无穷值
    *@param[in] value 给定数
    *@return 是否是无穷值
    * - true 是
    * - false 不是
    */
    template <typename T>
    bool is_infinite(const T &value)
    {    
#undef max
        T max_value = std::numeric_limits<T>::max();
        T min_value = -max_value;     
        return !(min_value <= value && value <= max_value);
    } 

    /*!
    *@brief 判断给定数是否是无效值
    *@param[in] value 给定数
    *@return 是否是无效值
    * - true 是
    * - false 不是
    */
    template <typename T>
    bool is_nan(const T &value)
    {
        return value != value; // True if NAN
    } 

    /*!
    *@brief 判断给定数是否有效
    *@param[in] value 给定数
    *@return 是否有效
    * - true 是
    * - false 不是
    */
    template <typename T>
    bool is_valid(const T &value)
    {    
        return !is_infinite(value) && !is_nan(value);
    }

    /*!
    *@brief 两个数的相对误差
    *@param[in] v1 第一个数
    *@param[in] v2 第二个数
    *@return 相对误差
    */
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

    /*!
    *@brief 对给定值四舍五入
    *@param[in] x 值
    *@return 四舍五入结果
    */
    template <typename T>
    int round_int(T x)
    {
        return (int) ((x > 0) ? (x + (T) 0.5) : (x - (T) 0.5)); 
    }

    /*!
    *@brief 对给定值四舍五入
    *@param[in]  x 值
    *@return 四舍五入结果
    */
    template <typename T>
    int64_t round_int64(T x)
    {
        return (int64_t)((x > 0) ? (x + (T) 0.5) : (x - (T) 0.5));
    }

    /*!
    *@brief 对给定值四舍五入
    *@param[in] x 值
    *@return 四舍五入结果
    */
    template<typename T>
    T round(T x)
    {
        return std::floor(x + (T)0.5);
    }

    /*!
    * @enum UnivariateFunctionPointAttribute
    * @brief 一元函数点的类型
    */
    enum UnivariateFunctionPointAttribute
    {
        
        PointMinimum    = 1,///< 局部极小值点
        
        PointInflection = 2,///< 拐点
        
        PointMaximum    = 4,///< 局部极大值点
        
        PointInterval   = 8,///< 是一段区间
        
        PointBreak      = 16,///< 是断点
    };
    
    /*
    * @struct ExtremumPointInfo
    * @brief 极值点信息
    */
    struct ExtremumPointInfo
    {
        //! 参数值
        double t;

        //! 函数值
        double f;

        //! 当 point_type 为 PointInterval 时，表示极值整段区间都是，且函数值都
        //! 相同，此时 t 与 f 表示区间的左右端点
        UnivariateFunctionPointAttribute point_attributes;
    };
