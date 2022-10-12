/*!
* @file      GVec1.h
* @brief     一维向量的基本定义
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
#pragma once

#include "GMathDef.h"

#pragma warning(disable:4251)

    /*!\addtogroup GMath GMath
    * @{
    */

    template<typename Type>
    struct CVec3;
    /*!
    * @struct CVec1
    * @brief  一维向量定义
              该类给出一维向量的定义，提供相关的计算方法，被其它很多类使用。该模板类对整型数据类型不支持
    */
    template<typename Type>
    struct CVec1
    {
    public:
        /*!
        *@brief  默认构造函数
        */
        CVec1():X(0){}

        /*!
        *@brief  通过常量构造二维向量
        *@param[in] x 常量
        */
        CVec1(const Type &x):X(x){}

        /*!
        *@brief  拷贝构造函数
        *@param[in] vecSrc 传入的CVec1
        */
        CVec1(const CVec1<Type>& vecSrc): X(vecSrc.X){}


        /*!
        *@brief  通过给定数据设置向量
        *@param[in] v 给定数
        */
        void Set(const Type v)
        {
            X = v;
        }


        /*!
        *@brief 按照索引获取向量的分量
        *@param[in] nIndex 索引，由调用者确保分量索引的有效性
        *@return 向量分量值
        */
        Type &  operator[](int nIndex)
        {
            return m_x[nIndex];
        }

        /*!
        *@brief 按照索引获取向量的分量
        *@param[in] nIndex 索引，由调用者确保分量索引的有效性
        *@return 向量分量值
        */
        const Type & operator[](int nIndex) const
        {
            return m_x[nIndex];
        }

        /*!
        *@brief 赋值给该向量
        *@param[in] vecSrc 用于赋值的向量
        *@return 结果向量
        */
        CVec1<Type> & operator =(const CVec1<Type> & vecSrc)
        {
            X = vecSrc.X;
            return *this;
        }

        /*!
        *@brief 该向量乘给定数值
        *@param[in] d 乘数
        *@return 结果向量
        */
        CVec1<Type> & operator *=(const Type &d)
        {
            X *= d;

            return *this;
        }

        /*!
        *@brief 该向量除以给定数值
        *@param[in] d 除数，调用者确保输入参数的有效性，不能为0
        *@return 结果向量
        */
        CVec1<Type> & operator /=(const Type &d)
        {
            X/= d;

            return *this;
        }

        /*!
        *@brief 向量与给定向量求差
        *@param[in] vecSrc 给定向量
        *@return 差向量
        */
        CVec1<Type> & operator -=(const CVec1<Type> & vecSrc)
        {
            X -= vecSrc.X;

            return *this;
        }

        /*!
        *@brief 向量反向
        *@return 反向量
        */
        CVec1<Type> operator-() const
        {
            return CVec1<Type>(-X);
        }


        /*!
        *@brief     向量除以一个数，
        *@param[in] vecSrc  被除向量
        *@param[in] d                除数
        *@return    结果向量
        */
        friend CVec1<Type> operator /(const CVec1<Type> & vecSrc, const Type &d)
        {
            return CVec1<Type>(vecSrc.X/d);
        }

        /*!
        *@brief 该向量乘给定向量
        *@param[in] v 给定向量
        *@return 结果向量
        */
        CVec1<Type> operator *(const CVec1 & v) const
        {
            return CVec1(X*v.X);
        }

        /*!
        *@brief 该向量乘给定数值
        *@param[in] d 乘数
        *@return 结果向量
        */
        CVec1<Type> operator *(Type d) const
        {
            return CVec1(X*d);
        }

        /*!
        *@brief 将输入向量与输入值相乘
        *@param[in] d 输入值
        *@param[in] v 输入向量
        *@return 结果向量
        */
        friend CVec1<Type> operator*(Type d, const CVec1<Type> & v)
        {
            return CVec1(d*v.X);
        }

        /*!
        *@brief 向量与给定向量求和
        *@param[in]  v 给定向量
        *@return 和向量
        */
        CVec1<Type> & operator+=(const CVec1<Type> &v)
        {
            X += v.X;
            return *this;
        }

        /*!
        *@brief  两个向量的加法，即向量的每个分量相加
        *@param[in] v1   向量1
        *@param[in] v2   向量2
        *@return    向量的和
        */
        friend CVec1<Type> operator +(const CVec1<Type> & v1, const CVec1<Type> & v2)
        {
            return CVec1<Type>(v1.X + v2.X);
        }

        /*!
        *@brief  两个向量的减法，即向量1的每个分量减去向量2的每个分量
        *@param[in] v1   向量1
        *@param[in] v2   向量2
        *@return    向量的差
        */
        friend CVec1<Type> operator -(const CVec1<Type> & v1, const CVec1<Type> & v2)
        {
            return CVec1<Type>(v1.X - v2.X);
        }


        /*!
        *@brief 判断向量是否相同
        *@param[in] v1 向量1
        *@param[in] v2 向量2
        *@return    是否相同
        * - true 是
        * - false 不是
        */
        friend bool operator ==(const CVec1<Type> & v1, const CVec1<Type> & v2)
        {
            return v1.X==v2.X;
        }

        /*!
        *@brief 判断向量是否不同
        *@param[in] v1 向量1
        *@param[in] v2 向量2
        *@return    是否不相同
        * - true 是
        * - false 不是
        */
        friend bool operator !=(const CVec1<Type> & v1, const CVec1<Type> & v2)
        {
            return !(v1 == v2);
        }


        union
        {
            Type m_x[1];
            struct
            {
                Type X;
            };
        };

    public:
    };


    typedef CVec1<int>    CVector1i;
    typedef CVec1<float>  CVector1f;
    typedef CVec1<double> CVector1d;

    /*! @} */
