/*!
* @file      GCoordinates3.h
* @brief     三维坐标系的基本定义
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
#ifndef G_COORDINATES3_H
#define G_COORDINATES3_H

#include "GMath/GVec2.h"
#include "GMath/GVec3.h"
#include "GMath/Interval.h"
#include "GMath/GMatrix4.h"

namespace ggp
{
    /*!\addtogroup GMath GMath
    * @{
    */

    template<typename Type>
    class CMatrix4;

    /*!
    \class CCoord3
    \brief 坐标系定义（用于传递参数）

    该结构定义了三维坐标轴，通过给定的XYZ正方向的矢量
    \sa CVec3  
    */
    template<typename Type>
    /*!
    * @class CCoord3
    * @brief 坐标系定义（用于传递参数）
             该结构定义了三维坐标轴，通过给定的XYZ正方向的矢量
    */
    class CCoord3
    {
    public:
        /*!   
        *@brief 坐标系的原点
        */

        CVec3<Type> Origin; 

        union
        {
            struct
            {
                /*!
                *@brief  坐标系基向量
                */
                CVec3<Type> Base[3];
            };

            struct
            {
                /*!
                *@brief 坐标系的X轴
                */
                CVec3<Type> X;
                /*!
                *@brief 坐标系的Y轴
                */
                CVec3<Type> Y;
                /*!
                *@brief 坐标系的Z轴
                */
                CVec3<Type> Z;
            };
        };

    public :
        /*!   
        *@brief 默认的构造函数
        */
        CCoord3() : Origin(CVec3<Type>::Zero),
            X(CVec3<Type>::UnitX),
            Y(CVec3<Type>::UnitY),
            Z(CVec3<Type>::UnitZ)
        {
        }

        /*!
        *@brief 拷贝构造函数
        *@param[in] rhs 传入坐标
        */
        CCoord3(const CCoord3<Type> & rhs)
        {
            *this = rhs;
        }

        /*!
        *@brief 构造函数，由原点，X轴，Y轴，Z轴生成右手仿射坐标系
        *@param[in] pos 原点
        *@param[in] dirX X轴，单位向量
        *@param[in] dirY Y轴，单位向量
        *@param[in] dirZ Z轴，单位向量
        */
        CCoord3(const CVec3<Type>& pos, const CVec3<Type>& dirX, const CVec3<Type>& dirY, const CVec3<Type>& dirZ)
            : Origin(pos), X(dirX), Y(dirY), Z(dirZ)
        {
            //g_AngleEpsilon is not defined in GMath
            //assert(X.IsUnit(g_AngleEpsilon));
            //assert(Y.IsUnit(g_AngleEpsilon));
            //assert(Z.IsUnit(g_AngleEpsilon));
            assert(IsRightHanded());
        }

         /*!
        *@brief 构造函数，由原点，X轴，Y轴生成右手正交坐标系
        *@param[in] pos 原点
        *@param[in] dirX X轴
        *@param[in] dirY Y轴
        */
        CCoord3(const CVec3<Type>& pos, const CVec3<Type>& dirX, const CVec3<Type>& dirY)
        {
            Set(pos, dirX, dirY) ;
        }

         /*!
        *@brief 构造函数，由原点和法向生成右手正交坐标系
        *@param[in] pos 原点
        *@param[in] normal 法向
        */
        CCoord3(const CVec3<Type>& pos, const CVec3<Type>& normal)
        {
            Set(pos, normal) ;
        }

         /*!
        *@brief 由原点，X轴，Y轴构建坐标系
        *@param[in] pos 原点
        *@param[in] dirX X轴
        *@param[in] dirY Y轴
        *@return 是否构建成功
        * - true 构建成功
        * - false 构建不成功
        */
        bool Set(const CVec3<Type>& pos, const CVec3<Type>& dirX, const CVec3<Type>& dirY)
        {
            Origin = pos ;
            X = dirX ;
            Y = dirY ;
            Z = dirX.Cross(dirY) ;
            bool rc = X.Normalize() != 0;
            rc = rc && Z.Normalize() != 0;
            Y = Z.Cross(X) ;
            return  rc;

        }

         /*!
        *@brief 构建水平的坐标系
        *@param[in] pos 原点
        *@param[in] dirX X轴
        *@param[in] dirY Y轴
        *@return 是否构建成功
        * - true 构建成功
        * - false 构建不成功
        */
        bool Set(const CVec2<Type>& pos, const CVec2<Type>& dirX, const CVec2<Type>& dirY)
        {
            Origin = CVec3<Type>(pos, 0.0) ;
            X = CVec3<Type>(dirX, 0.0) ;
            Y = CVec3<Type>(dirY, 0.0) ;
            Z = X.Cross(Y) ;
            bool rc = X.Normalize() != 0;
            rc = rc && Z.Normalize() != 0;
            Y = Z.Cross(X) ;
            return  rc;

        }

         /*!
        *@brief 由原点和法向构建坐标系
        *@param[in] pos 原点
        *@param[in] normal 法向
        *@return 是否构建成功
        * - true 构建成功
        * - false 构建不成功
        */
        bool Set(const CVec3<Type>& pos, const CVec3<Type>& normal)
        {
            ZAxisCoordinates(normal) ;

            Origin = pos ;

            return true ;
        }


        /*!
        *@brief     根据原点和三个角度构造坐标系
        *@param[in] pos 局部坐标系的原点
        *@param[in] axDeg 绕世界坐标系的X轴的旋转角度
        *@param[in] ayDeg 绕经第一步旋转后生成的y'轴的旋转角度
        *@param[in] azDeg 绕经第一步及第二步修改后生成的z"轴的旋转角度
        *@return 是否构建成功
        * - true 构建成功
        * - false 构建不成功
        */
        bool CreateFromAngle(const CVec3<Type>& pos, const Type axDeg, const Type ayDeg, const Type azDeg)
        {
            Origin = pos;
            
            double ax = DegToRad(axDeg);
            double ay = DegToRad(ayDeg);
            double az = DegToRad(azDeg);
            double sinx, cosx, siny, cosy, sinz, cosz;
            sincos(ax, &sinx, &cosx);
            sincos(ay, &siny, &cosy);
            sincos(az, &sinz, &cosz);

            //X' = CVec3<Type>::UnitX;
            //Y' = CVec3<Type>(0, cosx, sinx);
            //Z' = CVec3<Type>(0, -sinx, cosx);

            //Co = Mz * My * Mx * Ci
            X.Set(cosy * cosz, cosz * sinx * siny + cosx * sinz, -cosx * cosz * siny + sinx * sinz);
            Y.Set(-cosy * sinz, -sinx * siny * sinz + cosx * cosz, cosx * siny * sinz + cosz * sinx);
            Z.Set(siny, -sinx * cosy, cosx * cosy);

            return true;
        }

        /*!
        *@brief 三点生成坐标系
        *@param[in] p0 第一点，原点
        *@param[in] p1 第二点，X轴上的点
        *@param[in] p2 第三点，XY平面上的点
        *@return 是否构建成功
        * - true 构建成功
        * - false 构建不成功
        */
        bool CreateFromPt(const CVec3<Type>& p0, const CVec3<Type>& p1, const CVec3<Type>& p2)
        {
            Origin = p0 ;
            X = p1 - p0 ; 
            Y = p2 - p0 ; 
            Z = X.Cross(Y) ;

            bool rc = X.Normalize() != 0 ;
            rc = rc && Z.Normalize() != 0;

            Y = Z.Cross(X);
            return rc ;
        }

        /*!
        *@brief 根据向量z，生成坐标轴
        *@param[in] z 坐标轴Z
        */
        void ZAxisCoordinates(const CVec3<Type>& z)
        {
            Z = z;
            if(Z.Normalize() == 0)
                Z.Set(0, 0, 1);

            Y = Z ^ CVec3<Type>(0,0,1);
            if(Y.Normalize() == 0)
            {
                if(Z.Z > 0.0)
                {
                    Y.Set(0, 1, 0);
                }
                else
                {
                    Y.Set(0, -1, 0);
                }
            }
            X = Y ^ Z;
            X.Normalize();
        }


        //保证正交，但X,Y模长不够精确
        void FastOrthogonal(const CVec3<Type>& z)
        {
            auto imax = z.MaxDimension();

            const int next[] = {1,2,0,1};
            int iy = next[imax];
            int iz = next[imax+1];
            this->X[imax] = -(z[iy] + z[iz]);
            this->Z = z;
            this->X[iy] = z[imax];
            this->X[iz] = z[imax];
            auto rx = ggp::fast_rsqrt<1>(this->X.SqrLength());
            this->X *= rx;

            this->Y = this->Z.Cross(this->X);
        }

        /*!
        *@brief 由局部坐标得到世界坐标
        *@param[in] u 局部坐标X
        *@param[in] v 局部坐标Y
        *@param[in] w 局部坐标Z
        *@return 世界坐标点
        */
        CVec3<Type> PointAt(double u, double v, double w) const
        {
            return Origin + X * u + Y * v + Z * w;
        }
 
        /*!
        *@brief     由局部坐标得到世界坐标
        *@param[in] point 局部坐标点
        *@return    世界坐标点
        */
        CVec3<Type> PointAt(const CVec3<Type>& point) const
        {
            return Origin + X * point.X + Y * point.Y + Z* point.Z;
        }


        /*!
        *@brief     由局部坐标得到世界坐标
        *@param[in] point 局部坐标点
        *@return    世界坐标点
        */
        CVec3<Type> PointAt(const CVec2<Type>& point) const
        {
            return Origin + X * point.X + Y * point.Y;
        }

        /*!
        *@brief 由局部向量得到世界向量
        *@param[in] u 向量局部坐标X
        *@param[in] v 向量局部坐标Y
        *@param[in] w 向量局部坐标Z
        *@return 世界向量
        */
        CVec3<Type> VectorAt(double u, double v, double w) const
        {
            return X * u + Y * v + Z * w;
        }

        /*!
        *@brief 由局部向量得到世界向量
        *@param[in] vec 局部向量
        *@return 世界向量
        */
        CVec3<Type> VectorAt(const CVec3<Type> & vec) const
        {
            return X * vec.X + Y * vec.Y + Z * vec.Z;
        }

        /*!
        *@brief 由世界坐标得到局部坐标
        *@param[in] point 世界坐标点
        *@param[out] u 局部坐标X
        *@param[out] v 局部坐标Y
        *@param[out] w 局部坐标Z
        */
        void LocalPoint(const CVec3<Type>& point, double* u, double* v, double* w) const
        {
            CVec3<Type> vec = point - Origin ;
            if (u)
                *u = vec.Dot(X);
            if (v)
                *v = vec.Dot(Y);
            if (w)
                *w = vec.Dot(Z);
        }

         /*!
        *@brief 由世界坐标得到局部坐标
        *@param[in] point 世界坐标点
        *@return 局部坐标
        */
        CVec3<Type> LocalPoint(const CVec3<Type>& point) const
        {
            auto v = point - Origin ;
            return CVec3<Type>(v.Dot(X), v.Dot(Y), v.Dot(Z));
        }

        /*!
        *@brief 由世界坐标得到局部XY坐标
        *@param[in]  point 世界坐标点
        *@return 局部XY坐标
        */
        CVec2<Type> LocalPointXY(const CVec3<Type>& point) const
        {
            CVec3<Type> v = point - Origin ;
            return CVec2<Type>(v.Dot(X), v.Dot(Y));
        }

        /*!
        *@brief 计算传入点和当前坐标系XY面的距离
        *@param[in] worldPt 传入点
        *@return   距离
        */
        Type Distance(const CVec3<Type>& worldPt) const
        {
            CVec3<Type> vec = worldPt - Origin ;
            return vec.Dot(Z) ;
        }

        /*!
        *@brief 判断是否右手坐标系
        *@return 是否是右手坐标系
        * - true  是
        * - false 不是
        */
        bool IsRightHanded() const
        {
            return ((X ^ Y) * Z > 0);
        }

        /*!
         *@brief 判断是否单位坐标系
         *@param[in] typeTol 误差
         *@return 是否是单位坐标系
         * - true  是
         * - false 不是
         */
        bool IsUnit(const Type typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            return (X.IsUnit(typeTol) && Y.IsUnit(typeTol) && Z.IsUnit(typeTol));
        }

        /*!
         *@brief 判断坐标系是否正交
         *@param[in] typeTol 误差
         *@return 坐标系是否正交
         * - true  是
         * - false 不是
         */
        bool IsPerpendicular(const Type& typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            return (compareZero(X * Y, typeTol) == vrEqual && 
                    compareZero(X * Z, typeTol) == vrEqual && 
                    compareZero(Y * Z, typeTol) == vrEqual);
        }

        /*!
         *@brief 判断坐标系是否是单位正交右手系
         *@param[in] typeTol 误差
         *@return 是否是单位正交右手坐标系
         * - true  是
         * - false 不是
         */
        bool IsValidFrame(const Type& typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            return IsUnit(typeTol) && IsPerpendicular(typeTol) && IsRightHanded();
        }

        /*!
        *@brief     判断当前坐标系与传入坐标系是否相等
        *@param[in] rCoord  传入坐标系
        *@param[in] typeTol 角度容差
        *@return    是否相等
        * - true 相等
        * - false  不相等
        */
        bool IsEqual(const CCoord3<Type> &rCoord, const Type& typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            return Origin.IsEqual(rCoord.Origin, typeTol) && 
                        X.IsEqual(rCoord.X, typeTol) &&
                        Y.IsEqual(rCoord.Y, typeTol) &&
                        Z.IsEqual(rCoord.Z, typeTol);
        }

        /*!
        *@brief 计算传入坐标系的XY面和当前坐标系XY面是否平行
        *@param[in] coord   传入坐标系
        *@param[in] typeTol 角度容差
        *@return    是否平行
        * - true 平行
        * - false  不平行
        */
        bool IsParallel(const CCoord3<Type>& coord, const Type& typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            return Z.IsParallel(coord.Z,typeTol) ;
        }

        /*!
        *@brief 计算传入向量和当前坐标系XY面是否平行
        *@param[in] dir 传入向量
        *@param[in] typeTol 角度容差
        *@return    是否平行
        * - true 平行
        * - false  不平行
        */
        bool IsParallel(const CVec3<Type>& dir,const Type& typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            return Z.IsPerpendicular(dir,typeTol) ;
        }

        /*!
        *@brief 计算传入向量和当前坐标系XY面是否垂直
        *@param[in] dir 传入向量
        *@param[in] typeTol 角度容差
        *@return    是否垂直
        * - true 垂直
        * - false  不垂直
        */
        bool IsPerpendicular(const CVec3<Type>& dir,const Type& typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            return Z.IsParallel(dir,typeTol) ;
        }

        /*!
        *@brief 计算传入点距离当前坐标系XY面的最近的点
        *@param[in] worldPt 传入点
        *@return   最近点
        */
        CVec3<Type> ClosedPoint(const CVec3<Type>& worldPt) const
        {
            CVec3<Type> uv = LocalPoint(worldPt) ;
            return PointAt(uv.X, uv.Y, 0.0) ;
        }

        /*!
        *@brief 重新计算坐标系三个坐标轴，令三个轴相互垂直
        */
        void ResetZaxis()
        {
            Y = Z.Cross(X) ;
            X = Y.Cross(Z) ;
            X.Normalize() ;
            Y.Normalize() ;
            Z.Normalize() ;
        }

       
        /*!
        *@brief 由世界向量得到局部向量
        *@param[in] vector 世界向量
        *@param[out] u 向量局部坐标X
        *@param[out] v 向量局部坐标Y
        *@param[out] w 向量局部坐标Z
        */
        void LocalVector(const CVec3<Type>& vector, double* u, double* v, double* w) const
        {
            if (u)
                *u = vector.Dot(X) ;
            if (v)
                *v = vector.Dot(Y);
            if (w)
                *w = vector.Dot(Z);
        }

        /*!
        *@brief 由世界向量得到局部向量
        *@param[in] vector 世界向量
        *@return 局部坐标
        */
        CVec3<Type> LocalVector(const CVec3<Type>& vector) const
        {
            return CVec3<Type>(vector.Dot(X), vector.Dot(Y), vector.Dot(Z));
        }

        /*!
        *@brief 获得由世界坐标到局部坐标的变换矩阵
        *@param[out] matrix 返回的变换矩阵
        *@return 是否获得成功
        * - true 成功
        * - false 不成功
        */
        bool GetLocalMatrix(CMatrix4<Type>& matrix) const
        {
            matrix[0][0] = X.X; matrix[1][0] = Y.X; matrix[2][0] = Z.X; matrix[3][0] = Origin.X;
            matrix[0][1] = X.Y; matrix[1][1] = Y.Y; matrix[2][1] = Z.Y; matrix[3][1] = Origin.Y;
            matrix[0][2] = X.Z; matrix[1][2] = Y.Z; matrix[2][2] = Z.Z; matrix[3][2] = Origin.Z;
            matrix[0][3] = 0.0;     matrix[1][3] = 0.0;     matrix[2][3] = 0.0;     matrix[3][3] = 1.0;

            matrix.Invert();

            return true;
        }
        /*!
        *@brief 坐标系矩阵变换
        *@param[in] matrix 输入的变换矩阵
        *@return 返回对象
        */
        CCoord3& TransformBy(const CMatrix4<Type>& matrix)
        {
            Origin=matrix.MultiPointLeft(Origin);
            X=matrix.MultiVecLeft(X);
            Y=matrix.MultiVecLeft(Y);
            Z=matrix.MultiVecLeft(Z);
            return *this;
        }
        /*!
        *@brief 获得由局部坐标到世界坐标的变换矩阵
        *@param[out] matrix 返回的变换矩阵
        *@return 是否获得成功
        * - true 成功
        * - false 不成功
        */
        bool GetWorldMatrix(CMatrix4<Type>& matrix) const
        {
            matrix[0][0] = X.X; matrix[1][0] = Y.X; matrix[2][0] = Z.X; matrix[3][0] = Origin.X;
            matrix[0][1] = X.Y; matrix[1][1] = Y.Y; matrix[2][1] = Z.Y; matrix[3][1] = Origin.Y;
            matrix[0][2] = X.Z; matrix[1][2] = Y.Z; matrix[2][2] = Z.Z; matrix[3][2] = Origin.Z;
            matrix[0][3] = 0.0;     matrix[1][3] = 0.0;     matrix[2][3] = 0.0;     matrix[3][3] = 1.0;

            return true ;
        }

        /*!
        *@brief 获得由当前坐标系到目标坐标系的变换矩阵
        *@param[in] target 目标坐标系
        *@param[out] matrix 返回的变换矩阵
        *@return 是否获得成功
        * - true 成功
        * - false 不成功
        */
        bool GetTransformMatrix(const CCoord3<Type> &target, CMatrix4<Type>& matrix) const
        {
            CMatrix4<Type> matrix1, matrix2;
            GetWorldMatrix(matrix1);
            target.GetLocalMatrix(matrix2);

            matrix = matrix1 * matrix2;

            return true ;
        }

        /*!
        *@brief 将当前坐标系从它所在的局部坐标系转换到世界坐标系定义下，并设置自身
        *@param[in] rLocalCoord 当前坐标系所在局部坐标系
        */
        void ToWorldCoord(const CCoord3<Type> &rLocalCoord)
        {
            Origin.ToWorldPoint(rLocalCoord);
            X.ToWorldVt(rLocalCoord);
            Y.ToWorldVt(rLocalCoord);
            Z.ToWorldVt(rLocalCoord);
        }

        /*!
        *@brief 计算当前坐标系从它所在的局部坐标系转换到世界坐标系定义下后的新坐标系
        *@param[in] rLocalCoord 当前坐标系所在局部坐标系
        *@return 返回转换到世界坐标系后的新坐标系
        */
        CCoord3<Type> WorldCoord(const CCoord3<Type> &rLocalCoord) const
        {
            CCoord3<Type> rCoord;
            rCoord.Origin = Origin.WorldPoint(rLocalCoord);
            rCoord.X = X.WorldVector(rLocalCoord);
            rCoord.Y = Y.WorldVector(rLocalCoord);
            rCoord.Z = Z.WorldVector(rLocalCoord);
            return rCoord;
        }

        /*!
        *@brief 将当前坐标系从世界坐标系转化到指定的局部坐标系下，并设置自身
        *@param[in] rLocalCoord 指定的局部坐标系
        */
        void ToLocalCoord(const CCoord3<Type> &rLocalCoord)
        {
            Origin.ToLocalPoint(rLocalCoord);
            X.ToLocalVt(rLocalCoord);
            Y.ToLocalVt(rLocalCoord);
            Z.ToLocalVt(rLocalCoord);
        }

        /*!
        *@brief 计算当前坐标系从世界坐标系转化到指定的局部坐标系下
        *@param[in] rLocalCoord 定的局部坐标系
        *@return 返回在局部坐标系下的新坐标系
        */
        CCoord3<Type> LocalCoord(const CCoord3<Type> &rLocalCoord) const
        {
            CCoord3<Type> rCoord;
            rCoord.Origin = Origin.LocalPoint(rLocalCoord);
            rCoord.X = X.LocalVector(rLocalCoord);
            rCoord.Y = Y.LocalVector(rLocalCoord);
            rCoord.Z = Z.LocalVector(rLocalCoord);
            return rCoord;
        }

        double GetDiscreteStartAngle()
        {
            const CVec3<Type> dirZ(X ^ Y);
            const CVec3<Type> &dirW(fabs(dirZ.X) > fabs(dirZ.Y) ? CVec3<Type>::UnitY : CVec3<Type>::UnitX);
            CVec2<Type> dirW2d(dirW.LocalVectorXY(X, Y));
            dirW2d.Normalize();
            return angle(dirW2d);
        }

        char* AsString() const
        {
        return AsStringHelper::NewString("Origin{X=%.20le Y=%.20le Z=%.20le}, X{X=%.20le Y=%.20le Z=%.20le},Y{X=%.20le Y=%.20le Z=%.20le},Z{X=%.20le Y=%.20le Z=%.20le}", 
            Origin.X, Origin.Y, Origin.Z,X.X, X.Y, X.Z,Y.X, Y.Y, Y.Z,Z.X, Z.Y, Z.Z);
        }

        /*!
        *@brief     根据字符串生成坐标系
        *@param[in] str   传入的字符串
        *@return   生成的坐标系
        */
         static CCoord3<Type>* LoadFromString(char* str)
         {
             CCoord3<Type>* pCoord = new CCoord3<Type>;
             g_sscanf(str, "Origin{X=%le Y=%le Z=%le}, X{X=%le Y=%le Z=%le},Y{X=%le Y=%le Z=%le},Z{X=%le Y=%le Z=%le}", 
                 &pCoord->Origin.X, &pCoord->Origin.Y, &pCoord->Origin.Z, &pCoord->X.X, &pCoord->X.Y, &pCoord->X.Z, &pCoord->Y.X, &pCoord->Y.Y, &pCoord->Y.Z,
                 &pCoord->Z.X, &pCoord->Z.Y, &pCoord->Z.Z);
             return pCoord;
         }

    public:
        static const CCoord3<Type> GlobalCCoord3D;
    };

    template<class Type> const CCoord3<Type> CCoord3<Type>::GlobalCCoord3D(CVec3<Type>::Zero, CVec3<Type>::UnitX, CVec3<Type>::UnitY, CVec3<Type>::UnitZ);

    typedef CCoord3<int>    CCoordinates3i; 
    typedef CCoord3<float>  CCoordinates3f; 
    typedef CCoord3<double> CCoordinates3d;

    /*! @} */
} // namespace

#endif
