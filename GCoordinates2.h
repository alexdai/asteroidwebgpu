/*!
* @file      GCoordinates2.h
* @brief     二维坐标系的基本定义
* @details
* @author
* @date      2021/09/10
* @copyright Copyright 2012-2022 GLODON
**********************************************************************************
* @par 修改日志:
* |  Date  | Description | Author     |
* | :----: | :----       | :----      |
* | 2021/09/09 | 袁冠杰<yuangj-a> | GGP-55193 完成算法剩余部分的注释 |
*
**********************************************************************************
*/
#ifndef G_COORDINATES2_H
#define G_COORDINATES2_H

#include "GMath/GVec2.h"
#include "GMath/Interval.h"
#include "GMath/GMatrix3.h"

namespace ggp
{
    /*!\addtogroup GMath GMath
    * @{
    */

    template<typename Type>
    class CMatrix3;

    /*!
    \class CCoord2
    \brief 坐标系定义（用于传递参数）

    该结构定义了三维坐标轴，通过给定的XYZ正方向的矢量
    \sa CVec2
    */
    template<typename Type>
    /*!
    * @class CCoord2
    * @brief 二维坐标系基本定义
    */
    class CCoord2
    {
    public:
        /*!   
        *@brief 坐标系的原点
        */
        CVec2<Type> Origin; 
        /*!   
        *@brief 坐标系的X轴
        */
        CVec2<Type> X;
        /*!   
        *@brief 坐标系的Y轴
        */
        CVec2<Type> Y;

    public :
        /*!   
        *@brief 默认的构造函数
        */
        CCoord2() : Origin(CVec2<Type>::Zero), X(CVec2<Type>::UnitX), Y(CVec2<Type>::UnitY) {}

         /*!
        *@brief 构造函数，由X轴，和原点生成坐标系
        *@param[in] pos 原点
        *@param[in] dirX X轴
        */
        CCoord2(const CVec2<Type>& pos, const CVec2<Type>& dirX)
        {
            Set(pos, dirX) ;
        }
        
        /*!
        *@brief 由X轴，Y轴和原点构建坐标系
        *@param[in] pos 原点
        *@param[in] dirX X轴
        */
        void Set(const CVec2<Type>& pos, const CVec2<Type>& dirX)
        {
            Origin = pos ;
            X = dirX ;
            Y = dirX.GetRotateHalfPIVector(false);
        }

        /*!
        *@brief 二点生成坐标系
        *@param[in] p0 第一点，原点
        *@param[in] p1 第二点，X轴上的点
        */
        void CreateFromPt(const CVec2<Type>& p0, const CVec2<Type>& p1)
        {
            Origin = p0 ;
            X = p1 - p0 ; 
            X.Normalize();
            Y = X.GetRotateHalfPIVector(false);
        }

        /*!
        *@brief 由局部坐标得到世界坐标
        *@param[in] u 局部坐标X
        *@param[in] v 局部坐标Y
        *@return 世界坐标点
        */
        CVec2<Type> PointAt(double u, double v) const
        {
            return Origin + X * u + Y * v;
        }

        /*!
        *@brief     由局部坐标得到世界坐标
        *@param[in] point 局部坐标点
        *@return    世界坐标点
        */
        CVec2<Type> PointAt(const CVec2<Type>& point) const
        {
            return Origin + X * point.X + Y * point.Y;
        }

        /*!
        *@brief 由局部向量得到世界向量
        *@param[in] u 向量局部坐标X
        *@param[in] v 向量局部坐标Y
        *@return 世界向量
        */
        CVec2<Type> VectorAt(double u, double v) const
        {
            return X * u + Y * v;
        }

        /*!
        *@brief 由世界坐标得到局部坐标
        *@param[in] point 世界坐标点
        *@param[out] u 局部坐标X
        *@param[out] v 局部坐标Y
        */
        void LocalPoint(const CVec2<Type>& point, double* u, double* v) const
        {
            CVec2<Type> vec = point - Origin ;
            if (u)
                *u = vec.Dot(X);
            if (v)
                *v = vec.Dot(Y);
        }

         /*!
        *@brief 由世界坐标得到局部坐标
        *@param[in] point 世界坐标点
        *@return 局部坐标
        */
        CVec2<Type> LocalPoint(const CVec2<Type>& point) const
        {
            CVec2<Type> vec = point - Origin ;

            return CVec2<Type>(vec.Dot(X), vec.Dot(Y));
        }
       
        /*!
        *@brief 由世界向量得到局部向量
        *@param[in] vector 世界向量
        *@param[out] u 向量局部坐标X
        *@param[out] v 向量局部坐标Y
        *@param[out] w 无效参数，不知用户是否使用，暂不删除
        */
        void LocalVector(const CVec2<Type>& vector, double* u, double* v, double* w) const
        {
            if (u)
                *u = vector.Dot(X) ;
            if (v)
                *v = vector.Dot(Y);
        }

        /*!
        *@brief 判断坐标系是否相同
        *@param[in] coord1  坐标系1
        *@param[in] coord2  坐标系2
        *@return 坐标系是否相同
        * - true  相同
        * - false 不同
        */ 
        friend bool operator ==(const CCoord2<Type>& coord1, const CCoord2<Type>& coord2)
        { 
            return coord1.Origin == coord2.Origin && coord1.X == coord2.X && coord1.Y == coord2.Y;
        }

    public:
        static const CCoord2<Type> GlobalCCoord2D;
    };

    template<class Type> const CCoord2<Type> CCoord2<Type>::GlobalCCoord2D(CVec2<Type>::Zero, CVec2<Type>::UnitX);

    typedef CCoord2<int>    CCoordinates2i; 
    typedef CCoord2<float>  CCoordinates2f; 
    typedef CCoord2<double> CCoordinates2d;

    /*! @} */
} // namespace

#endif
