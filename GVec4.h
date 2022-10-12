/*!
* @file      GVec4.h
* @brief     四维向量的基本定义
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
#ifndef G_VEC4_H
#define G_VEC4_H

#include "GMathDef.h"

    /*!\addtogroup GMath GMath
    * @{
    */

    template<typename Type> 
    struct CVec3;
    /*!
    * @struct CVec4
    * @brief  四维向量定义
              该类给出四维向量的定义，提供相关的计算方法，被其它很多类使用。该模板类对整型数据类型不支持
    */
    template<typename Type>
    struct CVec4
    {
    public:

        /*!
        *@brief  默认构造函数
        */
        CVec4(): X(0), Y(0), Z(0), W(0)
        {
        }

        /*!
        *@brief  通过长度为4的数组构造四维向量
        *@param[in] v 长度为4的数组
        */
        CVec4(const Type v[4])
        {
            Set(v);
        }

        /*!
        *@brief  通过三维向量与常量w构造四维向量
        *@param[in] v 三维向量
        *@param[in] w 常量
        */
        CVec4(const CVec3<Type> & v, Type w)
        {
            Set(v.X, v.Y, v.Z, w);
        }

        /*!
        *@brief  通过向量的X,Y分量和Z分量构造四维向量
        *@param[in] iX X分量
        *@param[in] iY Y分量
        *@param[in] iZ Z分量
        *@param[in] iW W分量
        */
        CVec4(Type iX, Type iY, Type iZ, Type iW)
        {
            Set(iX, iY, iZ, iW);
        }

        /*!
        *@brief  四维向量的拷贝构造函数
        *@param[in] vecSrc 拷贝数据源向量
        */
        CVec4(const CVec4<Type> & vecSrc)
        {
            Set(vecSrc[0], vecSrc[1], vecSrc[2], vecSrc[3]);
        }

        /*!
        *@brief  向量的赋值操作
        *@param[in] vecSrc 数据源向量
        *@return 修改后的向量
        */
        CVec4<Type>& operator=(const CVec4<Type>& vecSrc)
        {
            Set(vecSrc[0], vecSrc[1], vecSrc[2], vecSrc[3]);

            return (*this);
        }

        /*!
        *@brief  通过给定数据设置向量
        *@param[in] v 给定数组
        *@return 设置结果
        */
        CVec4<Type> &  Set(const Type v[4])
        {
            m_xyzw[0] = v[0];
            m_xyzw[1] = v[1];
            m_xyzw[2] = v[2];
            m_xyzw[3] = v[3];

            return *this;
        }

        /*!
        *@brief  设置该向量的X和Y分量
        *@param[in] iX X分量
        *@param[in] iY Y分量
        *@param[in] iZ Z分量
        *@param[in] iW W分量
        *@return 设置后的当前向量
        */
        CVec4<Type> &  Set(Type iX, Type iY, Type iZ, Type iW)
        {
            m_xyzw[0] = iX;
            m_xyzw[1] = iY;
            m_xyzw[2] = iZ;
            m_xyzw[3] = iW;

            return *this;
        }

        /*!
        *@brief  返回向量数据的指针，数据长度为4
        *@return 向量数据的指针
        */
        const Type * Value() const
        {
            return &m_xyzw[0];
        }

        /*!
        *@brief  获取该向量的X和Y分量
        *@param[out] oX X分量
        *@param[out] oY Y分量
        *@param[out] oZ Z分量
        *@param[out] oW W分量
        */
        void Value(Type & oX, Type & oY, Type & oZ, Type & oW) const
        {
            oX = m_xyzw[0];
            oY = m_xyzw[1];
            oZ = m_xyzw[2];
            oW = m_xyzw[3];
        }

        /*!
        *@brief 按照索引获取向量的分量
        *@param[in] index 索引，由调用者确保分量索引的有效性
        *@return 向量分量值
        */
        Type Value(const int& index) const
        {
            if (index == 0)
            {
                return m_xyzw[0];
            }
            else if (index == 1)
            {
                return m_xyzw[1];
            }
            else if (index == 2)
            {
                return m_xyzw[2];
            }
            else
            {
                return m_xyzw[3];
            }

        }

        /*!
        *@brief 按照索引获取向量的分量的引用
        *@param[in] index 索引，由调用者确保分量索引的有效性
        *@return 向量分量值
        */
        Type& Value(const int& index) 
        {
            if (index == 0)
            {
                return m_xyzw[0];
            }
            else if (index == 1)
            {
                return m_xyzw[1];
            }
            else if (index == 2)
            {
                return m_xyzw[2];
            }
            else
            {
                return m_xyzw[3];
            }

        }

        /*!
        *@brief  整型向量
        *@return 结果向量
        */
        CVec4<int> Vec4i() const {return CVec4<int>((int)X, (int)Y, (int)Z, (int)W);}

        /*!
        *@brief  单精度向量
        *@return 结果向量
        */
        CVec4<float> Vec4f() const {return CVec4<float>((float)X, (float)Y, (float)Z, (float)W); }

        /*!
        *@brief  双精度向量
        *@return 结果向量
        */
        CVec4<double> Vec4d() const {return CVec4<double>((double)X, (double)Y, (double)Z, (double)W); }

        /*!
        *@brief  计算向量与给定向量的点积
        *@param[in] vecSrc 给定向量
        *@return 点积
        */
        Type Dot(const CVec4<Type> & vecSrc) const
        {
            return (m_xyzw[0]*vecSrc[0] + m_xyzw[1]*vecSrc[1] + m_xyzw[2]*vecSrc[2] + m_xyzw[3]*vecSrc[3]);
        }

        /*!
        *@brief  计算向量的模
        *@return 模
        */
        Type Length() const
        {
            return (Type)std::sqrt( SqrLength() );
        }

        /*!
        *@brief  计算向量模的平方
        *@return 模的平方
        */
        Type SqrLength() const
        {
            return (m_xyzw[0]*m_xyzw[0])+(m_xyzw[1]*m_xyzw[1])+(m_xyzw[2]*m_xyzw[2])+(m_xyzw[3]*m_xyzw[3]);
        }

        /*!
        *@brief  对向量的每个分量取平方
        *@return 结果向量
        */
        CVec4<Type> Square() const
        {
            return CVec4<Type>(X*X, Y*Y, Z*Z, W*W);
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
                (*this) *= (Type)(1.0 / magnitude);
            }
            else
            {
                magnitude = 0;
                Set(0.0, 0.0, 0.0, 0.0);
            }

            return magnitude;
        }

        /*!
        *@brief  向量的单位化，不改变当前向量
        *@param[in] eps 容差
        *@return 单位化后的向量
        */
        CVec4<Type> Unit(Type eps = std::numeric_limits<Type>::epsilon())
        {
            Type magnitude = Length();

            if(!IsNearZero(magnitude, eps))
                return CVec4<Type>(X / magnitude, Y / magnitude, Z / magnitude, W / magnitude);
            else
                return CVec4<Type>(0.0, 0.0, 0.0, 0.0);
        }

        /*!
        *@brief  对向量取反
        */
        void Negate()
        {
            Set(-m_xyzw[0], -m_xyzw[1], -m_xyzw[2], -m_xyzw[3]);
        }

        /*!
        *@brief  返回数组的常指针
        *@return 数组的常指针
        */
        inline const Type* Ptr() const
        {
            return &m_xyzw[0];
        }

        /*!
        *@brief  返回数组的指针
        *@return 数组的指针
        */
        inline Type* Ptr()
        {
            return &m_xyzw[0];
        }

        /*!
        *@brief 按照索引获取向量的分量
        *@param[in] nIndex 索引，由调用者确保分量索引的有效性
        *@return 向量分量值
        */
        Type & operator[](int nIndex) 
        { 
            return m_xyzw[nIndex]; 
        }

        /*!
        *@brief 按照索引获取向量的分量
        *@param[in] nIndex 索引，由调用者确保分量索引的有效性
        *@return 向量分量值
        */
        const Type & operator[](int nIndex) const 
        { 
            return m_xyzw[nIndex]; 
        }

        /*!
        *@brief 该向量乘给定数值
        *@param[in] d 乘数
        *@return 结果向量
        */
        CVec4<Type> & operator *=(const Type d)
        {
            m_xyzw[0] *= d;
            m_xyzw[1] *= d;
            m_xyzw[2] *= d;
            m_xyzw[3] *= d;

            return *this;
        }

        /*!
        *@brief 该向量除给定数值
        *@param[in] d 除数，调用者确保输入参数的有效性，不能为0
        *@return 结果向量
        */
        CVec4<Type> & operator /=(const Type d)
        {
            Type inv = 1.0f/d;

            m_xyzw[0] *= inv;
            m_xyzw[1] *= inv;
            m_xyzw[2] *= inv;
            m_xyzw[3] *= inv;

            return *this;
        }

        /*!
        *@brief 向量与给定向量求和
        *@param[in] vecSrc 给定向量
        *@return 和向量引用
        */
        CVec4<Type> & operator +=(const CVec4<Type> & vecSrc)
        {
            m_xyzw[0] += vecSrc.m_xyzw[0];
            m_xyzw[1] += vecSrc.m_xyzw[1];
            m_xyzw[2] += vecSrc.m_xyzw[2];
            m_xyzw[3] += vecSrc.m_xyzw[3];

            return *this;
        }

        /*!
        *@brief 向量与给定向量求差
        *@param[in] vecSrc 给定向量
        *@return 差向量引用
        */
        CVec4<Type> & operator -=(const CVec4<Type> & vecSrc)
        {
            m_xyzw[0] -= vecSrc.m_xyzw[0];
            m_xyzw[1] -= vecSrc.m_xyzw[1];
            m_xyzw[2] -= vecSrc.m_xyzw[2];
            m_xyzw[3] -= vecSrc.m_xyzw[3];

            return *this;
        }

        /*!
        *@brief 向量反向
        *@return 反向量
        */
        CVec4<Type> operator-() const
        {
            return CVec4<Type>(-m_xyzw[0], -m_xyzw[1], -m_xyzw[2], -m_xyzw[3]);
        }

        /*!
        *@brief 向量乘给定数值
        *@param[in] vecSrc 给定向量
        *@param[in] d 乘数
        *@return 结果向量
        */
        friend CVec4<Type> operator *(const CVec4<Type> & vecSrc, const Type d)
        { 
            return CVec4<Type>(vecSrc.m_xyzw[0] * d, 
                              vecSrc.m_xyzw[1] * d, 
                              vecSrc.m_xyzw[2] * d, 
                              vecSrc.m_xyzw[3] * d);
        }

        /*!
        *@brief 数值乘给定向量
        *@param[in] d 乘数
        *@param[in] vecSrc 给定向量
        *@return 结果向量
        */
        friend CVec4<Type> operator *(const Type d, const CVec4<Type> & vecSrc)
        { 
            return vecSrc * d; 
        }

        /*!
        *@brief  两个向量的点乘
        *@param[in] v1   向量1
        *@param[in] v2   向量2
        *@return    点乘结果
        */
        friend Type operator *(const CVec4<Type> & v1, const CVec4<Type> & v2)
        {
            return v1.Dot(v2);
        }


        /*!
        *@brief 两向量各分量乘法
        *@param[in] v 给定向量
        *@return 结果向量
        */
        CVec4<Type> Mul (const CVec4<Type> & v)
        {
            return CVec4<Type>(X*v.X, Y*v.Y, Z*v.Z, W*v.W);
        }

        /*!
        *@brief 两向量各分量除法
        *@param[in] v 给定向量
        *@return 结果向量
        */
        CVec4<Type> Div (const CVec4<Type> & v)
        {
            return CVec4<Type>(X/v.X, Y/v.Y, Z/v.Z, W/v.W);
        }

        CVec3<Type> & Vec3() const
        {
            return *(CVec3<Type>*)this;
        }


        /*!
        *@brief     向量除以一个数，
        *@param[in] vecSrc  被除向量  
        *@param[in] d                除数
        *@return                    结果向量
        */
        friend CVec4<Type> operator /(const CVec4<Type> & vecSrc, const Type d)
        { 
            return CVec4<Type>(vecSrc.m_xyzw[0] / d, 
                              vecSrc.m_xyzw[1] / d, 
                              vecSrc.m_xyzw[2] / d, 
                              vecSrc.m_xyzw[3] / d);
        }

        /*!
        *@brief  两个向量的加法，即向量的每个分量相加
        *@param[in] v1   向量1
        *@param[in] v2   向量2
        *@return    向量的和
        */
        friend CVec4<Type> operator +(const CVec4<Type> & v1, const CVec4<Type> & v2)
        {
            return CVec4<Type>(v1.m_xyzw[0] + v2.m_xyzw[0],
                              v1.m_xyzw[1] + v2.m_xyzw[1],
                              v1.m_xyzw[2] + v2.m_xyzw[2],
                              v1.m_xyzw[3] + v2.m_xyzw[3]);
        }

        /*!
        *@brief  两个向量的减法，即向量1的每个分量减去向量2的每个分量
        *@param[in] v1   向量1
        *@param[in] v2   向量2
        *@return    向量的差
        */
        friend CVec4<Type> operator -(const CVec4<Type> & v1, const CVec4<Type> & v2)
        {
            return CVec4<Type>(v1.m_xyzw[0] - v2.m_xyzw[0],
                              v1.m_xyzw[1] - v2.m_xyzw[1],
                              v1.m_xyzw[2] - v2.m_xyzw[2],
                              v1.m_xyzw[3] - v2.m_xyzw[3]);
        }

        /*!
        *@brief 判断向量是否相同
        *@param[in] v1   向量1
        *@param[in] v2   向量2
        *@return    是否相同
        * - true 是
        * - false 不是
        */ 
        friend bool operator ==(const CVec4<Type> & v1, const CVec4<Type> & v2)
        { 
            return (v1.m_xyzw[0]==v2.m_xyzw[0] && 
                    v1.m_xyzw[1]==v2.m_xyzw[1] && 
                    v1.m_xyzw[2]==v2.m_xyzw[2] && 
                    v1.m_xyzw[3]==v2.m_xyzw[3]); 
        }

        /*!
        *@brief 判断向量是否相同
        *@param[in] v1   向量1
        *@param[in] v2   向量2
        *@return    是否不相同
        * - true 是
        * - false 不是
        */
        friend bool operator !=(const CVec4<Type> & v1, const CVec4<Type> & v2)
        { 
            return !(v1 == v2); 
        }

        /*!
        *@brief 判断向量是否在给定误差下相等
        *@param[in] vecSrc 给定向量
        *@param[in] typeTol 误差
        *@return    是否相等
        * - true 是
        * - false 不是
        */
        bool IsEqual(const CVec4<Type> & vecSrc, const Type typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            return ( (m_xyzw - vecSrc).SqrLength() <= typeTol*typeTol );
        }

        /*!
        *@brief 判断向量在给定误差下是否零向量
        *@param[in] typeTol 误差
        *@return    是否为零向量
        * - true 是
        * - false 不是
        */
        bool IsZero(const Type typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            return ( SqrLength() <= typeTol*typeTol );
        }

        /*!
        *@brief     判断向量在给定误差下是否单位向量
        *@param[in] typeTol 误差
        *@return    是否为单位向量
        * - true 是
        * - false 不是
        */
        bool IsUnit(const Type typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            return fabs(SqrLength() - 1) <= 2*typeTol;
        }

        /*!
        *@brief  把Vec4按照规定格式字符串输出
        *@return     字符串
        */
        char* AsString()
        {
            return AsStringHelper::NewString("{X=%.20le Y=%.20le Z=%.20le W=%.20le}", X, Y, Z, W);
        }

        /*!
        *@brief 将向量流化到输出流
        *@param[out] oStream 流对象
        *@param[in] vecData 给定向量
        */
        friend std::ostream & operator<<(std::ostream & oStream, const CVec4<Type> & vecData)
        { 
            return oStream << vecData.X << " " << vecData.Y << " " << vecData.Z << " " << vecData.W; 
        }

        friend std::istream & operator>>(std::istream & iStream, CVec4<Type> & vecData)
        {
            return iStream >> vecData.X >> vecData.Y >> vecData.Z >> vecData.W; 
        }

        union
        {
            Type m_xyzw[4];
            struct
            {
                Type X; 
                Type Y;
                Type Z;
                Type W;
            };
        };

    public:
        static const CVec4<Type> UnitX;
        static const CVec4<Type> UnitY;
        static const CVec4<Type> UnitZ;
        static const CVec4<Type> UnitW;
        static const CVec4<Type> Zero;
    };

    template<class Type> const CVec4<Type> CVec4<Type>::UnitX(1, 0, 0, 0);
    template<class Type> const CVec4<Type> CVec4<Type>::UnitY(0, 1, 0, 0);
    template<class Type> const CVec4<Type> CVec4<Type>::UnitZ(0, 0, 1, 0);
    template<class Type> const CVec4<Type> CVec4<Type>::UnitW(0, 0, 0, 1);
    template<class Type> const CVec4<Type> CVec4<Type>::Zero(0, 0, 0, 0);

    typedef CVec4<int>    CVector4i;
    typedef CVec4<float>  CVector4f; 
    typedef CVec4<double> CVector4d;

    /*! @} */

#endif
