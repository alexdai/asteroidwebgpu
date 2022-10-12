/*!
* @file      GVec2.h
* @brief     二维向量的基本定义  
* @details
* @author
* @date      2021/09/14
* @copyright Copyright 2012-2022 GLODON
**********************************************************************************
* @par 修改日志:
* |  Date  | Description | Author     |
* | :----: | :----       | :----      |
* | 2021/09/14 | 袁冠杰<yuangj-a> | GGP-55193 完成算法剩余部分的注释 |
*
**********************************************************************************
*/
#ifndef G_VEC2_H
#define G_VEC2_H

#include "GMathDef.h"
#include <complex>

#pragma warning( disable:4251)

namespace ggp
{
    /*!\addtogroup GMath GMath
    * @{
    */

    template<typename Type>
    struct CVec3;
    /*!
    * @struct CVec2
    * @brief  二维向量定义
            该类给出二维向量的定义，提供相关的计算方法，被其它很多类使用。该模板类对整型数据类型不支持
    */
    template<typename Type>
    struct CVec2
    {
    public:
        /*!
        *@brief  默认构造函数
        */
        CVec2():X(0), Y(0){}

        /*!
        *@brief  通过长度为2的数组构造二维向量
        *@param[in] v 长度为2的数组
        */
        explicit CVec2(const Type v[2]): X(v[0]),Y(v[1]) {}

        /*!
        *@brief  通过常量x和常量y构造二维向量
        *@param[in] x 常量x
        *@param[in] y 常量y
        */
        CVec2(const Type &x, const Type &y):X(x), Y(y) {}

        /*!
        *@brief  拷贝构造函数
        *@param[in] vecSrc 另一个CVec2
        */
        CVec2(const CVec2<Type>& vecSrc): X(vecSrc.X),Y(vecSrc.Y) {}

        /*!
        *@brief  整型向量
        *@return 结果向量
        */
        CVec2<int> Vec2i() const {return CVec2<int>((int)X, (int)Y);}

        /*!
        *@brief  单精度向量
        *@return 结果向量
        */
        CVec2<float>  Vec2f() const {return CVec2<float>((float)X, (float)Y); }

        /*!
        *@brief  双精度向量
        *@return 结果向量
        */
        CVec2<double> Vec2d() const {return CVec2<double>((double)X, (double)Y); }

        /*!
        *@brief  当成复数使用
        *@return 复数
        */
        complex<Type> AsComplex() { return complex<Type>(X, Y); }

        /*!
        *@brief  交换分量得到的向量
        *@return 结果向量
        */
        CVec2<Type> YX() const { return CVec2<Type>(Y, X); }

        /*!
        *@brief  转成三维向量
        *@param[in] zValue 第三维的值
        *@return 结果向量
        */
        CVec3<Type> Vec3(Type zValue = 0) const { return CVec3<Type>((Type)X, (Type)Y, zValue); }

        /*!
        *@brief  共轭向量
        *@return 结果向量
        */
        CVec2<Type> Conjugate() const {return CVec2<Type>(X, -Y); }

        /*!
        *@brief 计算该点从全局坐标转换到在新局部坐标系下的点,并设置为自身
        *@param[in] localOrg 新坐标系原点
        *@param[in] localVx 新坐标系X轴单位化正方向
        *@param[in] localVy 新坐标系Y轴单位化正方向
        *@return 新坐标系下的点
        */
        CVec2<Type> & ToLocalPoint(const CVec2<Type> & localOrg, const CVec2<Type> & localVx, const CVec2<Type> & localVy) 
        {
            CVec2<Type> tempVec(m_xy[0] - localOrg[0],m_xy[1] - localOrg[1]);

            m_xy[0] = tempVec * localVx;
            m_xy[1] = tempVec * localVy;

            return *this;
        }

        /*!
        *@brief 计算该点从局部坐标系转换到全局局部坐标系下的点，并设置为自身
        *@param[in] localOrg 局部坐标系原点
        *@param[in] localVx 局部坐标系X轴单位化正方向
        *@param[in] localVy 局部坐标系Y轴单位化正方向
        *@return 全局坐标系下的点
        */
        CVec2<Type> & ToWorldPoint(const CVec2<Type> & localOrg, const CVec2<Type> & localVx, const CVec2<Type> & localVy) 
        {
            Type x = m_xy[0];
            Type y = m_xy[1];

            m_xy[0] = localOrg.X + x * localVx.X + y * localVy.X;
            m_xy[1] = localOrg.Y + x * localVx.Y + y * localVy.Y;

            return *this;
        }

        /*!
        *@brief 计算该点从全局坐标转换到在新局部坐标系下的点
        *@param[in] localOrg 局部坐标系原点
        *@param[in] localVx 局部坐标系X轴单位化正方向
        *@param[in] localVy 局部坐标系Y轴单位化正方向
        *@return 局部坐标系下的点
        */
        CVec2<Type>  GetLocalPt(const CVec2<Type> & localOrg, const CVec2<Type> & localVx, const CVec2<Type> & localVy) const
        {
            CVec2<Type> tempVec(m_xy[0] - localOrg[0],m_xy[1] - localOrg[1]);
            return CVec2(tempVec.Dot(localVx), tempVec.Dot(localVy));
        }

        /*!
        *@brief 计算该点从局部坐标转换到世界坐标系下的点
        *@param[in] localOrg 局部坐标系原点
        *@param[in] localVx 局部坐标系X轴单位化正方向
        *@param[in] localVy 局部坐标系Y轴单位化正方向
        *@return 世界坐标系下的点
        */
        CVec2<Type>  GetWorldPt(const CVec2<Type> & localOrg, const CVec2<Type> & localVx, const CVec2<Type> & localVy) const
        {
            return localOrg + m_xy[0] * localVx + m_xy[1] * localVy;
        }

        /*!
        *@brief 将向量旋转给定角度
        *@param[in]  angle 旋转角度的弧度值，正值表示逆时针方向
        *@return 旋转后的向量
        */
        CVec2<Type> & Rotate(Type angle)
        {
            double cosVal;
            double sinVal;

            sincos(angle, &sinVal, &cosVal); //x86下会使用fsincos指令，比单独调sin,cos更快

            Type x = m_xy[0];
            Type y = m_xy[1];

            m_xy[0] = Type(x * cosVal - y * sinVal);
            m_xy[1] = Type(x * sinVal + y * cosVal);

            return *this;
        }

        /*!
        *@brief 将向量旋转由给定向量确定的角度
        *@param[in] v 给定向量
        *@return 旋转后的向量
        */
        CVec2<Type> Rotate(const CVec2<Type>& v)
        {
            return CVec2<Type>(X*v.X - Y*v.Y, X*v.Y + Y*v.X);//X分量对应cos,Y对应sin
        }

        /*!
        *@brief 将向量按照给定方向旋转90度
        *@param[in] bClockwise 旋转方向，true为顺时针方向，false为逆时针方向
        *@return 旋转后的向量
        */
        CVec2<Type> & RotateHalfPI(bool bClockwise)
        {
            Type dTempX = m_xy[0];
            if (bClockwise)
            {
                m_xy[0] = m_xy[1];
                m_xy[1] = -dTempX;
            }
            else 
            {
                m_xy[0] = -m_xy[1];
                m_xy[1] = dTempX;
            }

            return *this;
        }

        /*!
        *@brief 将点绕着指定中心点旋转给定角度
        *@param[in] vecCenter 中心点
        *@param[in] angle 旋转角度的弧度值，正值表示逆时针方向
        *@return 旋转后的点
        */
        CVec2<Type> & RotateAroundPt(const CVec2<Type> & vecCenter, Type angle)
        {
            *this -= vecCenter;


            double cosVal;
            double sinVal;
            sincos(angle, &sinVal, &cosVal);

            Type x = m_xy[0];
            Type y = m_xy[1];

            m_xy[0] = Type(vecCenter[0] + x * cosVal - y * sinVal);
            m_xy[1] = Type(vecCenter[1] + x * sinVal + y * cosVal);

            return *this;
        }

        /*!
        *@brief 向量的赋旋转操作，将向量旋转给定角度,
        *@param[in] angle 旋转角度的弧度值，正值表示逆时针方向
        *@return 旋转后的向量
        */
        CVec2<Type> GetRotateVector(Type angle) const
        {
            double cosVal;
            double sinVal;
            sincos(angle, &sinVal, &cosVal);

            return CVec2(Type(m_xy[0] * cosVal - m_xy[1] * sinVal), Type(m_xy[0] * sinVal + m_xy[1] * cosVal));
        }

        /*!
        *@brief 将向量按照给定方向旋转Pi/2
        *@param[in] bClockwise 旋转方向，true为顺时针方向，false为逆时针方向
        *@return 旋转后的向量
        */
        CVec2<Type> GetRotateHalfPIVector(bool bClockwise) const
        {
            CVec2<Type> resultVec(*this);
            return resultVec.RotateHalfPI(bClockwise);
        }

        /*!
        *@brief 将点绕着指定中心点旋转给定角度
        *@param[in] vecCenter 中心点
        *@param[in] typeAngle 旋转角度
        *@return 旋转后的点
        */
        CVec2<Type> GetRotateAroundPoint(const CVec2<Type> & vecCenter, Type typeAngle) const
        {
            CVec2<Type> tempVec((*this) - vecCenter);
            return CVec2(vecCenter + tempVec.Rotate(typeAngle));
        }

        /*!
        *@brief     计算由X轴旋转到该向量的角度
        *@return    旋转角度，弧度表示
        */
        Type AngleFromXAxis() const
        {
            double dAngle = atan2(m_xy[1], m_xy[0]);
            return Type((dAngle < 0) ? dAngle + M_2PI : dAngle);
        }

        /*!
        *@brief     各个分量绝对值构成的向量
        *@return    结果向量
        */
        CVec2<Type> Abs() const
        {
            return CVec2<Type>(std::abs(X), std::abs(Y));
        }

        /*!
        * @brief      两向量夹角的余弦
        * @param [in] vec  另一个向量
        * @return     余弦
        */
        Type Cos(const CVec2<Type> & vec) const
        {
            Type value = ((*this)*vec)/sqrt(double((this->SqrLength()*vec.SqrLength())));
            Clamp((Type) -1.0, value, (Type) 1.0);
            return value;
        }


        /*!
        * @brief      两向量夹角的正弦
        * @param [in] vec
        * @return     正弦值
        */
        Type Sin(const CVec2<Type> & vec) const
        {
            Type value = abs((*this)^vec)/sqrt(double(SqrLength()*vec.SqrLength()));
            Clamp((Type) 0, value, (Type) 1.0);
            return value;
        }


        /*!
        * @brief      两单位向量夹角的正弦
        * @param [in] vec
        * @return     正弦值
        */
        Type SinUnit(const CVec2<Type> & vec) const
        {
            Type value = abs(double((*this)^vec));
            Clamp((Type) 0, value, (Type) 1.0);
            return value;
        }


        /*!
        *@brief  通过给定数据设置向量
        *@param[in] v 给定数组
        */
        void Set(const Type v[2])
        {
            m_xy[0] = v[0];
            m_xy[1] = v[1];
        }

        /*!
        *@brief  设置该向量的X和Y分量
        *@param[in] iX X分量
        *@param[in] iY Y分量
        */
        void Set(const Type &iX, const Type &iY)
        {
            m_xy[0] = iX;
            m_xy[1] = iY;
        }

        /*!
        *@brief  返回向量数据的指针，数据长度为2
        *@return 向量数据的指针
        */
        const Type * Value() const
        {
            return &m_xy[0];
        }

        /*!
        *@brief  获取该向量的X和Y分量
        *@param[out] oX X分量
        *@param[out] oY Y分量
        */
        void Value(Type & oX, Type & oY) const
        {
            oX = m_xy[0];
            oY = m_xy[1];
        }

        /*!
        *@brief 按照索引获取向量的分量
        *@param[in] index 索引，由调用者确保分量索引的有效性
        *@return 向量分量值
        */
        Type Value(const int& index) const
        {
            return m_xy[index];
        }

        /*!
        *@brief 按照索引获取向量的分量的引用
        *@param[in] index 索引，由调用者确保分量索引的有效性
        *@return 向量分量值
        */
        Type& Value(const int& index) 
        {
            return m_xy[index];
        }


        /*!
        *@brief  计算向量与给定向量的点积
        *@param[in] vecSrc 给定向量
        *@return 点积
        */
        Type Dot(const CVec2<Type> & vecSrc) const
        {
            return (m_xy[0]*vecSrc[0] + m_xy[1]*vecSrc[1]);
        }

        /*!
        *@brief  计算向量与给定向量的点积
        *@param[in] x 给定向量X分量值
        *@param[in] y 给定向量Y分量值
        *@return 点积
        */
        Type Dot(Type x, Type y) const
        {
            return (X*x + Y*y);
        }

        /*!
        *@brief  计算向量的模
        *@return 模
        */
        Type Length() const
        {
            return Type(std::sqrt( double(SqrLength()) ));
        }

        /*!
        *@brief  计算向量模的平方
        *@return 模的平方
        */
        Type SqrLength() const
        {
            return (m_xy[0]*m_xy[0])+(m_xy[1]*m_xy[1]);
        }

        /*!
        *@brief  计算向量的平方
        *@return 结果向量
        */
        CVec2<Type> Square() const
        {
            return CVec2<Type>(X*X, Y*Y);
        }


        /*!
        *@brief  向量的单位化
        *@param[in] eps 容差
        *@return 单位化前向量的模
        */
        Type Normalize(Type eps = std::numeric_limits<Type>::epsilon())
        {
            Type magnitude = Length();

            if(!IsNearZero(magnitude, eps))
            {
                (*this) *= Type(1.0 / magnitude);
            }
            else
            {
                magnitude = 0;
                X = 0; Y = 0;
            }

            return magnitude;
        }

        /*!
        *@brief  向量的单位化，不改变当前向量
        *@param[in] eps 容差
        *@return 单位化后的向量
        */
        CVec2<Type> Unit(Type eps = std::numeric_limits<Type>::epsilon())
        {
            Type magnitude = Length();

            if(!IsNearZero(magnitude, eps))
                return CVec2<Type>(X / magnitude, Y / magnitude);
            else
                return CVec2<Type>(0.0, 0.0);
        }


        /*!
        *@brief  计算向量与给定向量的叉积
        *@param[in] vecSrc 给定向量
        *@return 叉积
        */
        Type Cross(const CVec2<Type> & vecSrc) const
        {
            return (m_xy[0] * vecSrc[1] - m_xy[1] * vecSrc[0]);
        }

        /*!
        *@brief  对向量取反
        */
        void Negate()
        {
            Set(-m_xy[0], -m_xy[1]);
        }

        /*!
        *@brief  返回数组的常指针
        *@return 数组的常指针
        */
        inline const Type* Ptr() const
        {
            return &m_xy[0];
        }

        /*!
        *@brief  返回数组的指针
        *@return 数组的指针
        */
        inline Type* Ptr()
        {
            return &m_xy[0];
        }

        /*!
        *@brief 按照索引获取向量的分量
        *@param[in] nIndex 索引，由调用者确保分量索引的有效性
        *@return 向量分量值
        */
        Type &  operator[](int nIndex)
        { 
            return m_xy[nIndex]; 
        }

        /*!
        *@brief 按照索引获取向量的分量
        *@param[in] nIndex 索引，由调用者确保分量索引的有效性
        *@return 向量分量值
        */
        const Type & operator[](int nIndex) const
        { 
            return m_xy[nIndex]; 
        }

        /*!
        *@brief 赋值给该向量
        *@param[in] vecSrc 用于赋值的向量
        *@return 结果向量
        */
        CVec2<Type> & operator =(const CVec2<Type> & vecSrc)
        {
            m_xy[0] = vecSrc.m_xy[0];
            m_xy[1] = vecSrc.m_xy[1];

            return *this;
        }

        /*!
        *@brief 该向量乘给定数值
        *@param[in] d 乘数
        *@return 结果向量
        */
        CVec2<Type> & operator *=(const Type &d)
        {
            m_xy[0] *= d;
            m_xy[1] *= d;

            return *this;
        }

        /*!
        *@brief 该向量除以给定数值
        *@param[in] d 除数，调用者确保输入参数的有效性，不能为0
        *@return 结果向量
        */
        CVec2<Type> & operator /=(const Type &d)
        {
            Type inv = Type(1.0/d);

            m_xy[0] *= inv;
            m_xy[1] *= inv;

            return *this;
        }

        /*!
        *@brief 向量与给定向量求和
        *@param[in] vecSrc 给定向量
        *@return 和向量
        */
        CVec2<Type> & operator +=(const CVec2<Type> & vecSrc)
        {
            m_xy[0] += vecSrc.m_xy[0];
            m_xy[1] += vecSrc.m_xy[1];

            return *this;
        }

        /*!
        *@brief 向量与给定向量求差
        *@param[in] vecSrc 给定向量
        *@return 差向量
        */
        CVec2<Type> & operator -=(const CVec2<Type> & vecSrc)
        {
            m_xy[0] -= vecSrc.m_xy[0];
            m_xy[1] -= vecSrc.m_xy[1];

            return *this;
        }

        /*!
        *@brief 向量反向
        *@return 反向量
        */
        CVec2<Type> operator-() const
        {
            return CVec2<Type>(-m_xy[0], -m_xy[1]);
        }

        /*!
        *@brief 向量绝对值最大分量的索引
        *@return 索引值
        */
        int MaxDimension() const
        {
            return (std::abs(X) >= std::abs(Y)) ? 0 : 1;
        }

        //向量绝对值最小分量的索引，即使三个分量的绝对值相同，
        //能保证MinDimension和MaxDimension不同，这样可以用3 - min - max 来获得中间索引值
        /*!
        *@brief 向量绝对值最小分量的索引
        *@return 索引值
        */
        int MinDimension() const
        {
            return (std::abs(X) < std::abs(Y)) ? 0 : 1;
        }



        /*!
        *@brief 两向量各分量乘法
        *@param[in] v 给定向量
        *@return 结果向量
        */
        CVec2<Type> Mul (const CVec2<Type> & v)
        {
            return CVec2<Type>(X*v.X, Y*v.Y);
        }

        /*!
        *@brief 计算向量乘给定数值
        *@param[in] vecSrc 给定向量
        *@param[in] d 乘数
        *@return 结果向量
        */
        friend CVec2<Type> operator *(const CVec2<Type> & vecSrc, const Type &d)
        { 
            return CVec2<Type>(vecSrc.m_xy[0] * d, vecSrc.m_xy[1] * d);
        }

        /*!
        *@brief 计算数值乘给定向量
        *@param[in] d 乘数
        *@param[in] vecSrc 给定向量
        *@return 结果向量
        */
        friend CVec2<Type> operator *(const Type &d, const CVec2<Type> & vecSrc)
        { 
            return vecSrc * d; 
        }


        /*!
        *@brief  两个向量的点乘
        *@param[in] v1   向量1
        *@param[in] v2   向量2
        *@return    点乘结果
        */
        friend Type operator *(const CVec2<Type> & v1, const CVec2<Type> & v2)
        {
            return (v1.X*v2.X + v1.Y*v2.Y);
        }

        /*!
        *@brief  两个向量的叉乘
        *@param[in] v1   向量1
        *@param[in] v2   向量2
        *@return    叉乘结果
        */
        friend Type operator^(const CVec2<Type> & v1, const CVec2<Type> & v2)
        {
            return v1.Cross(v2); 
        }


        /*!
        *@brief     向量除以一个数，
        *@param[in] vecSrc  被除向量  
        *@param[in] d                除数
        *@return    结果向量
        */
        friend CVec2<Type> operator /(const CVec2<Type> & vecSrc, const Type &d)
        { 
            Type d1 = ((Type)(1))/d;
            return CVec2<Type>(vecSrc.m_xy[0] * d1, vecSrc.m_xy[1]*d1);
        }

        /*!
        *@brief     到另一个向量的转角
        *@param[in] v       另一向量
        *@return    转角[-π; π]
        */
        Type AngleTo(const CVec2<Type>& v) const
        {
            return std::atan2(Cross(v), Dot(v));
        }

        /*!
        *@brief 两向量各分量除法
        *@param[in] v 给定向量
        *@return 结果向量
        */
        CVec2<Type> Div(const CVec2<Type> & v)
        {
            return CVec2<Type>(X/v.X, Y/v.Y);
        }

        /*!
        *@brief  两个向量的加法，即向量的每个分量相加
        *@param[in] v1   向量1
        *@param[in] v2   向量2
        *@return    向量的和
        */
        friend CVec2<Type> operator +(const CVec2<Type> & v1, const CVec2<Type> & v2)
        {
            return CVec2<Type>(v1.m_xy[0] + v2.m_xy[0],v1.m_xy[1] + v2.m_xy[1]);
        }

        /*!
        *@brief  两个向量的减法，即向量1的每个分量减去向量2的每个分量
        *@param[in] v1   向量1
        *@param[in] v2   向量2
        *@return    向量的差
        */
        friend CVec2<Type> operator -(const CVec2<Type> & v1, const CVec2<Type> & v2)
        {
            return CVec2<Type>(v1.m_xy[0] - v2.m_xy[0], v1.m_xy[1] - v2.m_xy[1]);
        }

        /*!
        *@brief  向量每个分量减去某值
        *@param[in] v   向量
        *@param[in] d                  某值
        *@return    结果
        */
        friend CVec2<Type> operator -(const CVec2<Type> & v, Type d)
        {
            return CVec2<Type>(v.X - d, v.Y - d);
        }

        /*!
        *@brief 判断向量是否相同
        *@param[in] v1   向量1
        *@param[in] v2   向量2
        *@return    是否相同
        * - true 是
        * - false 不是
        */
        friend bool operator ==(const CVec2<Type> & v1, const CVec2<Type> & v2)
        { 
            return v1.m_xy[0]==v2.m_xy[0] && v1.m_xy[1]==v2.m_xy[1]; 
        }

        /*!
        *@brief 判断向量是否相同
        *@param[in] v1   向量1
        *@param[in] v2   向量2
        *@return    是否不相同
        * - true 是
        * - false 不是
        */
        friend bool operator !=(const CVec2<Type> & v1, const CVec2<Type> & v2)
        {
            return !(v1 == v2); 
        }
        
        /*!
        *@brief 计算传入向量和当前向量是否平行，同向或反向
        *@param[in] dir 传入向量
        *@param[in] typeTol 角度容差
        *@return    是否平行
        * - true 是
        * - false 不是
        */
        bool IsParallel(const CVec2<Type>& dir, const Type& typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            Type dis = Cross(dir) ;
            //return fabs(dis) < typeTol;
            return dis * dis <= SqrLength() * dir.SqrLength() * typeTol * typeTol;
        }

        /*!
        *@brief 计算传入向量和当前向量是否垂直
        *@param[in] dir 传入向量
        *@param[in] typeTol 角度容差
        *@return    是否垂直
        * - true 是
        * - false 不是
        */
        bool IsPerpendicular(const CVec2<Type>& dir, const Type& typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            Type dis = Dot(dir) ;
            //return fabs(dis) < typeTol;
            return dis * dis <= SqrLength() * dir.SqrLength() * typeTol * typeTol;
        }

        /*!
        *@brief 计算传入向量和当前向量在给定误差下是否相等
        *@param[in] vecSrc 传入向量
        *@param[in] typeTol 容差
        *@return    是否相等
        * - true 是
        * - false 不是
        */
        bool IsEqual(const CVec2<Type> & vecSrc, const Type typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            return ( (*this - vecSrc).SqrLength() <= typeTol*typeTol );
        }
                
        /*!
        *@brief 判断向量在给定误差下是否零向量
        *@param[in] typeTol 向量长度误差
        *@return    是否为零向量
        * - true 是
        * - false 不是
        */
        bool IsZero(const Type typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            return ( SqrLength() <= typeTol*typeTol );
        }

        /*!
        *@brief 判断向量在给定误差下是否单位向量
        *@param[in] typeTol 向量长度误差
        *@return    是否为单位向量
        * - true 是
        * - false 不是
        */
        bool IsUnit(const Type typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            return fabs(SqrLength() - 1) <= 2*typeTol;
        }

        /*!
        *@brief  坐标是否是有效的坐标值
        *@return    是否为有效坐标值
        * - true 是
        * - false 不是
        */

        bool IsValid() const { return (is_valid(X) && is_valid(Y)); }

        /*!
        *@brief 计算传入点和当前点的距离
        *@param[in] src 给定点
        *@return   距离
        */
        Type DistanceTo(const CVec2<Type> & src) const
        {
            double dx = (X - src.X);
            double dy = (Y - src.Y);
            return static_cast<Type>(sqrt(dx * dx + dy * dy));
        }

        /*!
        *@brief 计算传入点和当前点的距离平方
        *@param[in] src 给定点
        *@return   距离平方
        */
        Type SqrDistanceTo(const CVec2<Type> & src) const
        {
            return ((*this) - src).SqrLength();
        }


        /*!
        *@brief 将向量流化到输出流
        *@param[out] oStream 流对象
        *@param[in] vecData 给定向量
        */
        friend std::ostream & operator<<(std::ostream & oStream, const CVec2<Type> & vecData)
        {
            return oStream << vecData.X << " " << vecData.Y; 
        }

        friend std::istream & operator>>(std::istream & iStream, CVec2<Type> & vecData)
        {
            return iStream >> vecData.X >> vecData.Y; 
        }

        union
        {
            Type m_xy[2];
            struct
            {
                Type X;
                Type Y;
            };
        };

    public:
        static const CVec2<Type> UnitX;
        static const CVec2<Type> UnitY;
        static const CVec2<Type> Zero;
        static const CVec2<Type> NaN;

    };

    template<class Type> const CVec2<Type> CVec2<Type>::UnitX(1, 0);
    template<class Type> const CVec2<Type> CVec2<Type>::UnitY(0, 1);
    template<class Type> const CVec2<Type> CVec2<Type>::Zero(0, 0);
    template<class Type> const CVec2<Type> CVec2<Type>::NaN(g_Nan, g_Nan);


    typedef CVec2<int>     CVector2i; 
    typedef CVec2<float>   CVector2f; 
    typedef CVec2<double>  CVector2d;
    typedef CVec2<int64_t> CVector2l;

    /*!
    *@brief 计算给定向量逆时针转到该向量的转角,要求两个向量为单位向量
    *@param[in] vecUnitFrom 起始向量
    *@param[in] vecUnitTo 终止向量
    *@return 弧度值 0到2Pi
    */
    template<typename Type>
    inline Type GetRotateAngle(const CVec2<Type> & vecUnitFrom, const CVec2<Type> & vecUnitTo)
    {
        double dAngle = AcosSafe(vecUnitTo.Dot(vecUnitFrom));
        if ((vecUnitFrom ^ vecUnitTo) < -std::numeric_limits<Type>::epsilon())
        {
            dAngle = M_2PI - dAngle;
        }

        return Type(dAngle);
    }

    /*! @} */
} // namespace

#endif  //G_VEC2_H
