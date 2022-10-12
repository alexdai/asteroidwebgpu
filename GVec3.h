/*!
* @file      GVec3.h
* @brief     三维向量的基本定义
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
#ifndef G_VEC3_H
#define G_VEC3_H

#include "GMathDef.h"
#include "GVec2.h"
#include "GCoordinates3.h"
#include "Matrix.h"

#pragma warning( disable:4251)

    /*!\addtogroup GMath GMath
    * @{
    */

    template<typename Type>
    class CCoord3;

    template<typename Type>
    struct CVec4;

    /*!
    * @struct CVec3
    * @brief  三维向量定义
              该类给出三维向量的定义，提供相关的计算方法，被其它很多类使用。该模板类对整型数据类型不支持
    */
    template<typename Type>
    struct CVec3
    {
    public:

        /*!
        *@brief  默认构造函数
        */
        CVec3(): X(0), Y(0), Z(0){}

        /*!
        *@brief  通过长度为3的数组构造三维向量
        *@param[in] v 长度为3的数组
        */
        CVec3(const Type v[3]):X(v[0]), Y(v[1]), Z(v[2]) {}

        /*!
        *@brief  通过向量的X,Y分量和Z分量构造三维向量
        *@param[in] x X分量
        *@param[in] y Y分量
        *@param[in] z Z分量
        */
        CVec3(const Type &x, const Type &y, const Type &z): X(x),Y(y),Z(z) {}

        /*!
        *@brief  通过向量的X,Y分量和Z分量以及权重w构造三维向量
        *@param[in] x X分量
        *@param[in] y Y分量
        *@param[in] z Z分量
        *@param[in] w 权重值
        */
        CVec3(const Type &x, const Type &y, const Type &z, const Type &w)
        {
            Type w1 = (Type)1/w;
            X = x*w1; Y = y*w1; Z = z*w1;
        }

        /*!
        *@brief  根据二维点和给定Z值组合三维点
        *@param[in] vecSrc2d 二维坐标点
        *@param[in] z Z分量
        */
        CVec3(const CVec2<Type> &vecSrc2d, const Type &z):X(vecSrc2d.X),Y(vecSrc2d.Y), Z(z) {}

        /*!
        *@brief  CVec3拷贝构造
        *@param[in] vecSrc3d  三维坐标点
        */
        CVec3(const CVec3<Type>& vecSrc3d):X(vecSrc3d.X),Y(vecSrc3d.Y), Z(vecSrc3d.Z) {}

        /*!
        *@brief  前两个分量的单精度向量
        *@return 结果向量
        */
        CVec2<float>  Vec2f() const {return CVec2<float>((float)X, (float)Y); }

        /*!
        *@brief  前两个分量的双精度向量
        *@return 结果向量
        */
        CVec2<double> Vec2d() const {return CVec2<double>((double)X, (double)Y); }

        /*!
        *@brief  XY构成的向量
        *@return XY分量向量的引用，可直接修改这两个分量
        */
        CVec2<Type> & XY() const { return *(CVec2<Type> *)this; }

        /*!
        *@brief  YZ构成的向量
        *@return YZ分量向量的引用，可直接修改这两个分量
        */
        CVec2<Type> & YZ() { return *(CVec2<Type> *)&Y; }

        /*!
        *@brief  按向量分量的排列,产生新向量
        *@param[in] dim 旧向量
        *@return 新向量
        */
        const CVec3<Type> Swizzle(const int dim[3]) const
        {
            return CVec3<Type>(m_xyz[dim[0]], m_xyz[dim[1]], m_xyz[dim[2]]);
        }

        /*!
        *@brief  ZX构成的向量
        *@return ZX分量向量
        */
        CVec2<Type> ZX() { return CVec2<Type>(Z, X); }

        /*!
        *@brief  整型向量
        *@return 结果向量
        */
        CVec3<int> Vec3i() const {return CVec3<int>((int)X, (int)Y, (int)Z);}

        /*!
        *@brief  单精度向量
        *@return 结果向量
        */
        CVec3<float>  Vec3f() const {return CVec3<float>((float)X, (float)Y, (float)Z); }

        /*!
        *@brief  双精度向量
        *@return 结果向量
        */
        CVec3<double> Vec3d() const {return CVec3<double>((double)X, (double)Y, (double)Z); }

        /*!
        *@brief  获得Vec2
        *@return Vec2
        */
        CVec2<Type> Vec2() const {return CVec2<Type>((Type)X, (Type)Y);}

        /*!
        *@brief  转成四维向量，第4个分量默认为1
        *@param[in] w 第四个分量
        *@return 结果向量
        */
        CVec4<Type>  Vec4(Type w = 1) const  {return CVec4<Type>((Type)X, (Type)Y, (Type)Z, (Type)w); }

        /*!
        *@brief  坐标是否是有效的坐标值
        *@return    是否是有效的坐标值
        * - true 是
        * - false 不是
        */
        bool IsValid() const { return (is_valid(X) && is_valid(Y) && is_valid(Z)); }

        /*!
        *@brief  通过给定数据设置向量
        *@param[in] v 给定数组
        */
        void Set(const Type v[3])
        {
            X = v[0];
            Y = v[1];
            Z = v[2];
        }

        /*!
        *@brief  设置该向量的X、Y分量以及Z分量
        *@param[in] iX X分量
        *@param[in] iY Y分量
        *@param[in] iZ Z分量
        */
        void Set(const Type &iX, const Type &iY, const Type &iZ)
        {
            X = iX;
            Y = iY;
            Z = iZ;
        }

        /*!
        *@brief  根据二维点和给定Z值组合三维点
        *@param[in] vecSrc2d 二维坐标点
        *@param[in] iZ Z分量
        */
        void Set(const CVec2<Type> &vecSrc2d, const Type &iZ)
        {
            m_xyz[0] = vecSrc2d.X;
            m_xyz[1] = vecSrc2d.Y;
            m_xyz[2] = iZ;
        }

        /*!
        *@brief  返回向量数据的指针，数据长度为3
        *@return 向量数据的指针
        */
        const Type * Value() const
        {
            return &m_xyz[0];
        }

        /*!
        *@brief  获取该向量的X和Y分量
        *@param[out] oX X分量
        *@param[out] oY Y分量
        *@param[out] oZ Z分量
        */
        void Value(Type & oX, Type & oY, Type & oZ) const
        {
            oX = X;
            oY = Y;
            oZ = Z;
        }


        /*!
        *@brief 按照索引获取向量的分量
        *@param[in] index 索引，由调用者确保分量索引的有效性
        *@return 向量分量值
        */
        Type Value(const int& index) const
        {
            return m_xyz[index];
        }


        /*!
        *@brief 按照索引获取向量的分量引用
        *@param[in] index 索引，由调用者确保分量索引的有效性
        *@return 向量分量值
        */
        Type& Value(const int& index) 
        {
            return m_xyz[index];
        }

        /*!
        *@brief  计算向量与给定向量的点积
        *@param[in] vecSrc 给定向量
        *@return 点积
        */
        Type Dot(const CVec3<Type> & vecSrc) const
        {
            return (m_xyz[0]*vecSrc[0] + m_xyz[1]*vecSrc[1] + m_xyz[2]*vecSrc[2]);
        }
        /*!
        *@brief  计算向量与给定向量的点积
        *@param[in] vecSrc 给定向量
        *@return 点积
        */
        Type DotProduct(const CVec3<Type> & vecSrc) const
        {
            return (m_xyz[0]*vecSrc[0] + m_xyz[1]*vecSrc[1] + m_xyz[2]*vecSrc[2]);
        }

        /*!
        *@brief 向量与给定向量相乘得到矩阵
        *@param[in] v 给定向量
        *@return 结果矩阵
        */
        Matrix<Type, 3, 3> Cov(const CVec3<Type> & v) const
        {
            Matrix<Type, 3, 3> M;

            M(0, 0) = X*v.X;  M(0, 1) = X*v.Y; M(0, 2) = X*v.Z;
            M(1, 0) = Y*v.X;  M(1, 1) = Y*v.Y; M(1, 2) = Y*v.Z;
            M(2, 0) = Z*v.X;  M(2, 1) = Z*v.Y; M(2, 2) = Z*v.Z;

            return M;
        }

        /*!
        *@brief  计算向量与给定向量的点积
        *@param[in] v 给定向量
        *@return 点积
        */
        Type operator * (const CVec3<Type> &v) const
        {
            return (X*v.X + Y*v.Y + Z*v.Z); 
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
            return (m_xyz[0]*m_xyz[0])+(m_xyz[1]*m_xyz[1])+(m_xyz[2]*m_xyz[2]);
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
                X = 0.0; Y = 0.0; Z = 0.0;
            }

            return magnitude;
        }

        /*!
        *@brief  向量的单位化，不改变当前向量
        *@param[in] eps 容差
        *@return 单位化后的向量
        */
        CVec3<Type> Unit(Type eps = std::numeric_limits<Type>::epsilon()) const
        {
            Type magnitude = Length();

            if(!IsNearZero(magnitude, eps))
                return CVec3<Type>(X / magnitude, Y / magnitude, Z / magnitude);
            else
                return CVec3<Type>(0.0, 0.0, 0.0);
        }

        /*!
        *@brief  计算向量与给定向量的叉积
        *@param[in] vecSrc 给定向量
        *@return 叉积
        */
        CVec3<Type> Cross(const CVec3<Type> & vecSrc) const
        {
            return CVec3<Type>(m_xyz[1] * vecSrc[2] - vecSrc[1] * m_xyz[2],
                              m_xyz[2] * vecSrc[0] - vecSrc[2] * m_xyz[0],
                              m_xyz[0] * vecSrc[1] - vecSrc[0] * m_xyz[1]);
        }
        /*!
        *@brief  计算向量与给定向量的叉积
        *@param[in] vecSrc 给定向量
        *@return 叉积
        */
        CVec3<Type> CrossProduct(const CVec3<Type> & vecSrc) const
        {
            return CVec3<Type>(m_xyz[1] * vecSrc[2] - vecSrc[1] * m_xyz[2],
                              m_xyz[2] * vecSrc[0] - vecSrc[2] * m_xyz[0],
                              m_xyz[0] * vecSrc[1] - vecSrc[0] * m_xyz[1]);
        }

        /*!
        *@brief  对向量取反
        */
        void Negate()
        {
            X = -X;
            Y = -Y;
            Z = -Z;
        }

        /*!
        *@brief 计算两个向量的夹角
        *@param[in] vec 另一个向量
        *@return 弧度值0到Pi
        */
        Type Angle(const CVec3<Type> & vec) const
        {
            Type cosVal = this->Dot(vec) /sqrt((this->SqrLength()*vec.SqrLength()));
            return AcosSafe(cosVal);
        }

        /*!
        *@brief 计算两个向量的夹角,要求两个向量为单位向量
        *@param[in] vec 另一个向量
        *@return 弧度值0到Pi
        */
        Type AngleUnit(const CVec3<Type> & vec) const
        {
            Type cosVal = this->Dot(vec);
            return AcosSafe(cosVal);
        }

        /*!
        * @brief      两向量夹角的余弦
        * @param [in] vec  另一个向量
        * @return     余弦值
        */
        Type Cos(const CVec3<Type> & vec) const
        {
            Type value = this->Dot(vec)/sqrt((this->SqrLength()*vec.SqrLength()));
            Clamp((Type) -1.0, value, (Type) 1.0);
            return value;
        }


        /*!
        * @brief      两向量夹角的正弦
        * @param [in] vec 另一个向量
        * @return     正弦值
        */
        Type Sin(const CVec3<Type> & vec) const
        {
            Type dot = (*this)*vec;
            Type cos2 = dot*dot/(SqrLength()*vec.SqrLength());
            Clamp((Type) 0, cos2, (Type) 1.0);
            return sqrt(1-cos2);
        }

        /*!
        *@brief 对向量的每个分量取绝对值
        *@return 结果向量
        */
        CVec3<Type> Abs() const 
        {
            return CVec3<Type>(abs(X), abs(Y), abs(Z));
        }

        /*!
        *@brief 对向量的每个分量取平方
        *@return 结果向量
        */
        CVec3<Type> Square() const
        {
            return CVec3<Type>(X*X, Y*Y, Z*Z);
        }

        /*!
        * @brief      两单位向量的夹角的正弦
        * @param [in] vec
        * @return     正弦值
        */
        Type SinUnit(const CVec3<Type> & vec)
        {
            Type dot = this->Dot(vec);
            Clamp((Type) -1.0, dot, (Type) 1.0);
            return sqrt(1-dot*dot);
        }

        /*!
        *@brief  返回数组的常指针
        *@return 数组的常指针
        */
        inline const Type* Ptr() const
        {
            return &m_xyz[0];
        }

        /*!
        *@brief  返回数组的指针
        *@return 数组的指针
        */
        inline Type* Ptr()
        {
            return &m_xyz[0];
        }

        /*!
        *@brief 按照索引获取向量的分量
        *@param[in] nIndex 索引，由调用者确保分量索引的有效性
        *@return 向量分量值
        */
        inline Type& operator[](int nIndex) 
        {
            return m_xyz[nIndex]; 
        }

        /*!
        *@brief 按照索引获取向量的分量
        *@param[in] nIndex 索引，由调用者确保分量索引的有效性
        *@return 向量分量值
        */
        inline const Type & operator[](int nIndex) const 
        { 
            return m_xyz[nIndex]; 
        }

        /*!
        *@brief  两个向量的加法，即向量的每个分量相加
        *@param[in] v   另一个向量
        *@return    向量的和
        */
        inline CVec3<Type> operator +(const CVec3<Type> & v) const
        {    
            return CVec3<Type>(X + v.X, Y + v.Y, Z + v.Z);
        }

        /*!
        *@brief  两个向量的减法，即向量1的每个分量减去向量2的每个分量
        *@param[in] v   向量
        *@return    向量的差
        */
        inline CVec3<Type> operator -(const CVec3<Type> & v) const
        {    
            return CVec3<Type>(X - v.X, Y - v.Y, Z - v.Z);
        }


        /*!
        *@brief 该向量乘给定数值
        *@param[in] d 乘数
        *@return 结果向量
        */
        CVec3<Type> & operator *=(const Type &d)
        {
            m_xyz[0] *= d;
            m_xyz[1] *= d;
            m_xyz[2] *= d;

            return *this;
        }

        /*!
        *@brief 该向量除以给定数值
        *@param[in] d 除数，调用者确保输入参数的有效性，不能为0
        *@return 结果向量
        */
        CVec3<Type> & operator /=(const Type &d)
        {
            *this *= ((Type)(1))/d;

            return *this;
        }


        /*!
        *@brief 向量与给定向量求和
        *@param[in] vecSrc 给定向量
        *@return 和向量引用
        */
        CVec3<Type> & operator +=(const CVec3<Type> & vecSrc)
        {
            m_xyz[0] += vecSrc.m_xyz[0];
            m_xyz[1] += vecSrc.m_xyz[1];
            m_xyz[2] += vecSrc.m_xyz[2];

            return *this;
        }

        /*!
        *@brief 向量与给定向量求差
        *@param[in] vecSrc 给定向量
        *@return 差向量引用
        */
        CVec3<Type> & operator -=(const CVec3<Type> & vecSrc)
        {
            m_xyz[0] -= vecSrc.m_xyz[0];
            m_xyz[1] -= vecSrc.m_xyz[1];
            m_xyz[2] -= vecSrc.m_xyz[2];

            return *this;
        }

        /*!
        *@brief 向量反向
        *@return 反向量
        */
        CVec3<Type> operator-() const
        {
            return CVec3<Type>(-m_xyz[0], -m_xyz[1], -m_xyz[2]);
        }

        /*!
        *@brief 将this三维点从局部坐标系转换到全局坐标系下，局部坐标系由Pos、vec3DirX、vec3DirY、vec3DirZ定义，并设置自身
        *@param[in] Pos 原点
        *@param[in] vec3DirX X正方向
        *@param[in] vec3DirY Y正方向
        *@param[in] vec3DirZ Z正方向
        *@return 转换以后的三维点
        */
        CVec3<Type> & ToWorldPoint(const CVec3<Type> & Pos, const CVec3<Type> & vec3DirX, const CVec3<Type> & vec3DirY, const CVec3<Type> & vec3DirZ)
        {
            *this = Pos + vec3DirX * m_xyz[0] + vec3DirY * m_xyz[1] + vec3DirZ * m_xyz[2];
            return *this;
        }

        /*!
        *@brief 将this三维点由rCoord3定义的局部坐标系转换到全局坐标系，并设置自身
        *@param[in] rLocalCoord 三维坐标系
        *@return 转换以后的三维点
        */
        CVec3<Type> & ToWorldPoint(const CCoord3<Type> & rLocalCoord)
        {
            ToWorldPoint(rLocalCoord.Origin, rLocalCoord.X, rLocalCoord.Y, rLocalCoord.Z);
            return *this;
        }

        /*!
        *@brief 将this三维点从全局坐标系转换到由Pos，vec3DirX、vec3DirY、vec3DirZ定义的局部坐标系下，并设置自身
        *@param[in] Pos 原点
        *@param[in] vec3DirX X正方向
        *@param[in] vec3DirY Y正方向
        *@param[in] vec3DirZ Z正方向
        *@return 转换以后的三维点
        */
        CVec3<Type>& ToLocalPoint(const CVec3<Type> & Pos, const CVec3<Type> & vec3DirX, const CVec3<Type> & vec3DirY, const CVec3<Type> & vec3DirZ)
        {
            CVec3<Type> vecThis = *this;
            CVec3<Type> vecTmp = vecThis - Pos;
            m_xyz[0] = vecTmp.Dot(vec3DirX);
            m_xyz[1] = vecTmp.Dot(vec3DirY);
            m_xyz[2] = vecTmp.Dot(vec3DirZ);

            return *this;
        }

        /*!
        *@brief 将this三维点从全局坐标系转换到由rCoord3定义的局部坐标系下，并设置自身
        *@param[in] rCoord3 三维坐标系
        *@return 转换以后的三维点
        */
        CVec3<Type>& ToLocalPoint(const CCoord3<Type> & rCoord3)
        {
            ToLocalPoint(rCoord3.Origin, rCoord3.X, rCoord3.Y, rCoord3.Z);
            return *this;
        }

        /*!
        *@brief 将this向量从局部坐标系转换到全局坐标系下，局部坐标系由vec3DirX、vec3DirY、vec3DirZ定义，并设置自身
        *@param[in] vec3DirX X正方向
        *@param[in] vec3DirY Y正方向
        *@param[in] vec3DirZ Z正方向
        *@return 转换以后的三维向量
        */
        CVec3<Type> & ToWorldVt(const CVec3<Type> & vec3DirX, const CVec3<Type> & vec3DirY, const CVec3<Type> & vec3DirZ)
        {
            *this = vec3DirX * m_xyz[0] + vec3DirY * m_xyz[1] + vec3DirZ * m_xyz[2];
            return *this;
        }

        /*!
        *@brief 将this向量由rCoord3定义的局部坐标系转换到全局坐标系，并设置自身
        *@param[in] rLocalCoord 三维坐标系
        *@return 转换以后的三维向量
        */
        CVec3<Type> & ToWorldVt(const CCoord3<Type> & rLocalCoord)
        {
            ToWorldVt(rLocalCoord.X, rLocalCoord.Y, rLocalCoord.Z);
            return *this;
        }

        /*!
        *@brief 将this向量从全局坐标系转换到由vec3DirX、vec3DirY、vec3DirZ定义的局部坐标系下，并设置自身
        *@param[in] vec3DirX X正方向
        *@param[in] vec3DirY Y正方向
        *@param[in] vec3DirZ Z正方向
        *@return 转换以后的三维向量
        */
        CVec3<Type>& ToLocalVt(const CVec3<Type> & vec3DirX, const CVec3<Type> & vec3DirY, const CVec3<Type> & vec3DirZ)
        {
            CVec3<Type> vec = (*this);
            m_xyz[0] = vec.Dot(vec3DirX);
            m_xyz[1] = vec.Dot(vec3DirY);
            m_xyz[2] = vec.Dot(vec3DirZ);
            return *this;
        }

        /*!
        *@brief 将this向量从全局坐标系转换到由rCoord3定义的局部坐标系下，并设置自身
        *@param[in] rLocalCoord 三维坐标系
        *@return 转换以后的三维向量
        */
        CVec3<Type>& ToLocalVt(const CCoord3<Type> & rLocalCoord)
        {
            ToLocalVt(rLocalCoord.X, rLocalCoord.Y, rLocalCoord.Z);
            return *this;
        }


        /*!
        *@brief 计算由Pos,vec3DirX、vec3DirY、vec3DirZ定义的局部坐标系下的三维点转换到全局坐标系下
        *@param[in] Pos 原点
        *@param[in] vec3DirX X正方向
        *@param[in] vec3DirY Y正方向
        *@param[in] vec3DirZ Z正方向
        *@return 转换以后的三维点
        */
        CVec3<Type> WorldPoint(const CVec3<Type> & Pos,const CVec3<Type> & vec3DirX, const CVec3<Type> & vec3DirY, const CVec3<Type> & vec3DirZ) const
        {
            return Pos + vec3DirX * m_xyz[0] + vec3DirY * m_xyz[1] + vec3DirZ * m_xyz[2];
        }

        /*!
        *@brief 计算给定局部坐标系下的点在世界坐标下的坐标
        *@param[in] coord 局部坐标系
        *@return 世界坐标系下的坐标
        */
        CVec3<Type> WorldPoint(const CCoord3<Type> & coord) const
        {
            return coord.Origin + coord.X*m_xyz[0] + coord.Y*m_xyz[1] + coord.Z*m_xyz[2];
        }

        /*!
        *@brief 计算将三维点从全局坐标系转换到由Pos，vec3DirX、vec3DirY、vec3DirZ定义的局部坐标系下
        *@param[in] Pos 原点
        *@param[in] vec3DirX X正方向
        *@param[in] vec3DirY Y正方向
        *@param[in] vec3DirZ Z正方向
        *@return 转换以后的三维点
        */
        CVec3<Type> LocalPoint(const CVec3<Type> & Pos,const CVec3<Type> & vec3DirX, const CVec3<Type> & vec3DirY, const CVec3<Type> & vec3DirZ) const
        {
            CVec3<Type> vec = *this - Pos;

            return CVec3<Type>(vec.Dot(vec3DirX), vec.Dot(vec3DirY), vec.Dot(vec3DirZ));
        }

        /*!
        *@brief 计算将向量由全局坐标投影到局部XY坐标下的结果
        *@param[in] Pos 局部坐标的基点
        *@param[in] vec3DirX 局部坐标的X正方向
        *@param[in] vec3DirY 局部坐标的Y正方向
        *@return 投影后的二维向量
        */
        CVec2<Type> LocalPointXY(const CVec3<Type> & Pos,const CVec3<Type> & vec3DirX, const CVec3<Type> & vec3DirY) const
        {
            CVec3<Type> v = *this - Pos;

            return CVec2<Type>(v.Dot(vec3DirX), v.Dot(vec3DirY));
        }

        /*!
        *@brief 计算将三维点从全局坐标系转换到局部坐标系下的结果
        *@param[in]  coord 局部坐标系
        *@return 转换以后的三维点
        */
        CVec3<Type> LocalPoint(const CCoord3<Type>& coord) const
        {
            return LocalPoint(coord.Origin, coord.X, coord.Y, coord.Z);
        }

        /*!
        *@brief 计算将向量由全局坐标投影到局部XY坐标下的结果
        *@param[in]  coord 局部坐标，取前两个分量得到局部XY坐标
        *@return 投影后的二维向量
        */
        CVec2<Type> LocalPointXY(const CCoord3<Type>& coord) const
        {
            return LocalPointXY(coord.Origin, coord.X, coord.Y);
        }


        /*!
        *@brief 计算由vec3DirX、vec3DirY、vec3DirZ定义的局部坐标系下的向量转换到全局坐标系下
        *@param[in] vec3DirX X正方向
        *@param[in] vec3DirY Y正方向
        *@param[in] vec3DirZ Z正方向
        *@return 转换以后的三维向量
        */
        CVec3<Type> WorldVector(const CVec3<Type> & vec3DirX, const CVec3<Type> & vec3DirY, const CVec3<Type> & vec3DirZ) const
        {
            return vec3DirX * m_xyz[0] + vec3DirY * m_xyz[1] + vec3DirZ * m_xyz[2];
        }

        /*!
        *@brief 计算由rCoord3定义的局部坐标系下的向量转换到全局坐标系下
        *@param[in]  rLocalCoord 三维坐标系
        *@return 转换以后的三维向量
        */
        CVec3<Type> WorldVector(const CCoord3<Type> & rLocalCoord) const
        {
            return WorldVector(rLocalCoord.X, rLocalCoord.Y, rLocalCoord.Z);
        }

        /*!
        *@brief 计算将向量从全局坐标系转换到由vec3DirX、vec3DirY、vec3DirZ定义的局部坐标系下
        *@param[in] vec3DirX X正方向
        *@param[in] vec3DirY Y正方向
        *@param[in] vec3DirZ Z正方向
        *@return 转换以后的三维向量
        */
        CVec3<Type> LocalVector(const CVec3<Type> & vec3DirX, const CVec3<Type> & vec3DirY, const CVec3<Type> & vec3DirZ) const
        {
            return CVec3<Type>(Dot(vec3DirX), Dot(vec3DirY), Dot(vec3DirZ));
        }

        /*!
        *@brief 计算将向量由全局坐标投影到局部XY坐标下的结果
        *@param[in] vec3DirX 局部坐标的X正方向
        *@param[in] vec3DirY 局部坐标的Y正方向
        *@return 投影后的二维向量
        */
        CVec2<Type> LocalVectorXY(const CVec3<Type> & vec3DirX, const CVec3<Type> & vec3DirY) const
        {
            return CVec2<Type>(Dot(vec3DirX), Dot(vec3DirY));
        }

        /*!
        *@brief 计算向量从全局坐标系转换到由rCoord3定义的局部坐标系下
        *@param[in] rLocalCoord 三维坐标系
        *@return 转换以后的三维向量
        */
        CVec3<Type> LocalVector(const CCoord3<Type> & rLocalCoord) const
        {
            return LocalVector(rLocalCoord.X, rLocalCoord.Y, rLocalCoord.Z);
        }


        /*!
        *@brief 计算将向量由全局坐标投影到局部XY坐标下的结果
        *@param[in] rLocalCoord 局部坐标，取前两个分量得到局部XY坐标
        *@return 投影后的二维向量
        */
        CVec2<Type> LocalVectorXY(const CCoord3<Type> & rLocalCoord) const
        {
            return LocalVectorXY(rLocalCoord.X, rLocalCoord.Y);
        }

        /*!
        *@brief 两向量各分量乘法
        *@param[in] v 给定向量
        *@return 结果向量
        */
        CVec3<Type> Mul (const CVec3<Type> & v) const
        {
            return CVec3<Type>(X*v.X, Y*v.Y, Z*v.Z);
        }

        /*!
        *@brief 向量绝对值最大分量的索引
        *@return 索引
        */
        int MaxDimension() const
        {
            int max_dim = 0;
            auto max_len = abs(X);

            if(abs(Y) > max_len)
            {
                max_dim = 1;
                max_len = abs(Y);
            }

            if(abs(Z) > max_len)
            {
                max_dim = 2;
            }

            return max_dim;
        }

        //向量绝对值最小分量的索引，即使三个分量的绝对值相同，
        //能保证MinDimension和MaxDimension不同，这样可以用3 - min - max 来获得中间索引值
        /*!
        *@brief 向量绝对值最小分量的索引
        *@return 索引
        */
        int MinDimension() const
        {
            int min_dim = 1;
            auto min_len = abs(Y);

            if(abs(X) < min_len)
            {
                min_dim = 0;
                min_len = abs(X);
            }

            if(abs(Z) < min_len)
            {
                min_dim = 2;
            }

            return min_dim;
        }

        /*!
        *@brief 计算向量乘给定数值
        *@param[in] vecSrc 给定向量
        *@param[in] d 乘数
        *@return 结果向量
        */
        friend CVec3<Type> operator *(const CVec3<Type> & vecSrc, const Type d)
        { 
            return CVec3<Type>(vecSrc.m_xyz[0] * d, vecSrc.m_xyz[1] * d, vecSrc.m_xyz[2] * d);
        }

        /*!
        *@brief 计算数值乘给定向量
        *@param[in] d 乘数
        *@param[in] vecSrc 给定向量
        *@return 结果向量
        */
        friend CVec3<Type> operator *(const Type &d, const CVec3<Type> & vecSrc)
        { 
            return vecSrc * d; 
        }

        /*!
        *@brief  两个向量的叉乘
        *@param[in] v1   向量1
        *@param[in] v2   向量2
        *@return    叉乘结果
        */
        friend CVec3<Type> operator^(const CVec3<Type> & v1, const CVec3<Type> & v2)
        {
            return v1.Cross(v2); 
        }

        /*!
        *@brief     向量除以一个数，
        *@param[in] vecSrc  被除向量  
        *@param[in] d                除数
        *@return      结果向量
        */
        friend CVec3<Type> operator /(const CVec3<Type> & vecSrc, Type d)
        { 
            Type d1 = ((Type)1)/d;
            return CVec3<Type>(vecSrc.X * d1, vecSrc.Y * d1, vecSrc.Z * d1);
        }

        /*!
        *@brief 两向量各分量除法
        *@param[in] v 给定向量
        *@return 结果向量
        */
        CVec3<Type> Div(const CVec3<Type> & v)
        {
            return CVec3<Type>(X/v.X, Y/v.Y, Z/v.Z);
        }

        /*!
        *@brief 判断向量是否相同
        *@param[in] v1 向量1
        *@param[in] v2 向量2
        *@return    是否相同
        * - true 是
        * - false 不是
        */ 
        friend bool operator ==(const CVec3<Type> & v1, const CVec3<Type> & v2)
        { 
            return (v1.m_xyz[0]==v2.m_xyz[0] && v1.m_xyz[1]==v2.m_xyz[1] && v1.m_xyz[2]==v2.m_xyz[2]); 
        }

        /*!
        *@brief 判断向量是否相同
        *@param[in] v1 向量1
        *@param[in] v2 向量2
        *@return    是否不相同
        * - true 是
        * - false 不是
        */
        friend bool operator !=(const CVec3<Type> & v1, const CVec3<Type> & v2)
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
        bool IsEqual(const CVec3<Type> & vecSrc, const Type typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            return ( (*this - vecSrc).SqrLength() <= typeTol*typeTol );
        }
        
        /*!
        *@brief 判断向量在给定误差下是否零向量
        *@param[in] typeTol 误差
        *@return    是否是零向量
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
        *@return    是否是单位向量
        * - true 是
        * - false 不是
        */
        bool IsUnit(const Type typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            //return fabs(SqrLength() - 1) <= 2*typeTol;
            return abs(Length()-1.0) <= typeTol;
        }

        /*!
        *@brief 根据绝对值最大的分量确定要退化的维度
        *@param[out] reversed 退化的维度对应分量是否为正值
        *@return 要退化的维度
        */
        int GetMainIndex(bool& reversed) const
        {
            int mainIndex = 0; // drop x
            reversed = m_xyz[0] < 0;

            Type abs = fabs(m_xyz[0]);
            Type max = abs;
            if ((abs = fabs(m_xyz[1])) > max)
            {
                max = abs;
                mainIndex = 1; // drop y
                reversed = m_xyz[1] < 0;
            }
            if ((abs = fabs(m_xyz[2])) > max)
            {
                mainIndex = 2;  // drop z
                reversed = m_xyz[2] < 0;
            }

            return mainIndex;
        }

        /*!
        *@brief 根据绝对值最大的分量确定要退化的维度
        *@return 要退化的维度
        */
        int GetMainIndex() const
        {
            bool reversed;
            return GetMainIndex(reversed);
        }

        /*!
        *@brief 获取退化后的向量
        *@param[in] dropIndex 退化的维度
        *@return 退化后得到的二维向量
        */
        CVec2<Type> GetDroppedVector(int dropIndex) const
        {
            switch (dropIndex)
            {
            case 0: return CVec2<Type>(m_xyz[1], m_xyz[2]);
            case 1: return CVec2<Type>(m_xyz[2], m_xyz[0]);
            case 2: return CVec2<Type>(m_xyz[0], m_xyz[1]);
            default:
                assert(false);
            }
            return CVec2<Type>(0, 0);
        }

        /*!
        *@brief 计算传入点和当前点的距离
        *@param[in] src 给定向量
        *@return   距离
        */
        Type DistanceTo(const CVec3<Type> & src) const
        {
            auto dx = (X - src.X);
            auto dy = (Y - src.Y);
            auto dz = (Z - src.Z);
            return sqrt(dx * dx + dy * dy + dz * dz);
        }

        /*!
        *@brief 计算传入点和当前点的距离平方
        *@param[in] src 给定向量
        *@return   距离平方
        */
        Type SqrDistanceTo(const CVec3<Type> & src) const
        {
            return ((*this) - src).SqrLength();
        }

        /*!
        *@brief 计算当前向量投影到给定向量后的向量及其垂直向量
        *@param[in] src 给定向量
        *@param[out] projVec 投影向量
        *@param[out] projVertVec 与投影向量相垂直向量
        *@return    是否投影有效
        * - true 是
        * - false 不是
        */
        bool ProjectTo(const CVec3<Type>& src, CVec3<Type>& projVec, CVec3<Type>& projVertVec) const
        {
            Type sqrLen = src.SqrLength();
            if (sqrLen > 0)
            {
                projVec = src * (this->Dot(src) / sqrLen);
                projVertVec = *this - projVec;
                return true;
            }
            else
            {
                return false;
            }
        }

        /*!
        *@brief 计算传入向量和当前向量是否平行，同向或反向
        *@param[in] dir 传入向量
        *@param[in] typeTol 角度容差
        *@return    是否平行
        * - true 是
        * - false 不是
        */
        bool IsParallel(const CVec3<Type>& dir, const Type& typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            CVec3<Type> vec = Cross(dir) ;
            return vec.SqrLength() <= SqrLength() * dir.SqrLength() * typeTol * typeTol;
        }

        /*!
        *@brief 计算传入向量和当前向量是否垂直
        *@param[in] dir 传入向量
        *@param[in] typeTol 角度容差
        *@return    是否垂直
        * - true 是
        * - false 不是
        */
        bool IsPerpendicular(const CVec3<Type>& dir, const Type& typeTol=std::numeric_limits<Type>::epsilon()) const
        {
            Type dis = Dot(dir) ;
            //return fabs(dis) < typeTol;
            return dis * dis <= SqrLength() * dir.SqrLength() * typeTol * typeTol;
        }

        bool AskPerpendicularVector(CVec3<Type>& perpendicularVec)const
        {
            if (IsZero(g_DoubleResolution))
            {
                perpendicularVec.Set(0.0,0.0,0.0);
                return false;
            }

            const double fabs_x = fabs(X);
            const double fabs_y = fabs(Y);
            const double fabs_z = fabs(Z);

            if (IsLessThan(fabs_x,fabs_y,g_DoubleResolution))
            {
                if (IsLessThan(fabs_x,fabs_z,g_DoubleResolution))
                {
                    perpendicularVec.X = 0.0;
                    if (IsLessThan(fabs_y,fabs_z,g_DoubleResolution))
                    {
                        perpendicularVec.Y = Z;
                        perpendicularVec.Z = -Y;
                    }
                    else
                    {
                        perpendicularVec.Y = -Z;
                        perpendicularVec.Z = Y;
                    }
                }
                else
                {
                    perpendicularVec.Z = 0.0;
                    if (IsLessThan(fabs_y,fabs_x,g_DoubleResolution))
                    {
                        perpendicularVec.X = -Y;
                        perpendicularVec.Y = X;
                    }
                    else
                    {
                        perpendicularVec.X = Y;
                        perpendicularVec.Y = -X;
                    }
                }
            }
            else
            {
                if (IsLessThan(fabs_y,fabs_z,g_DoubleResolution))
                {
                    perpendicularVec.Y = 0.0;
                    if (IsLessThan(fabs_x,fabs_z,g_DoubleResolution)) 
                    {
                        perpendicularVec.X = Z;
                        perpendicularVec.Z = -X;
                    }
                    else
                    {
                        perpendicularVec.X = -Z;
                        perpendicularVec.Z = X;
                    }
                }
                else
                {
                    perpendicularVec.Z = 0.0;
                    if (IsLessThan(fabs_y,fabs_x,g_DoubleResolution))
                    {
                        perpendicularVec.X = -Y;
                        perpendicularVec.Y = X;
                    }
                    else
                    {
                        perpendicularVec.X = Y;
                        perpendicularVec.Y = -X;
                    }
                }
            }

            return true;
        }

        /*!
        *@brief  把Vec3按照规定格式字符串输出
        *@return    字符串
        */
        char* AsString() const
        {
            return AsStringHelper::NewString("{X=%.20le Y=%.20le Z=%.20le}", X, Y, Z);
        }

        /*!
        *@brief 将向量流化到输出流
        *@param[out] oStream 流对象
        *@param[in] vecData 给定向量
        */
        friend std::ostream & operator<<(std::ostream & oStream, const CVec3<Type> & vecData)
        { 
            return oStream << vecData.X << " " << vecData.Y << " " << vecData.Z; 
        }

        friend std::istream & operator>>(std::istream & iStream, CVec3<Type> & vecData)
        { 
            return iStream >> vecData.X >> vecData.Y >> vecData.Z; 
        }

        friend std::wostream& operator<< (std::wostream& s, const CVec3<Type>& v)
        {
            return s << L"[X: " << v.X << L", Y: " << v.Y << L", Z: " << v.Z;
        }

        union
        {
            Type m_xyz[3];
            struct
            {
                Type X; 
                Type Y;
                Type Z;
            };
        };

    public:
        static const CVec3<Type> UnitX;
        static const CVec3<Type> UnitY;
        static const CVec3<Type> UnitZ;
        static const CVec3<Type> NegaUnitX;
        static const CVec3<Type> NegaUnitY;
        static const CVec3<Type> NegaUnitZ;
        static const CVec3<Type> Zero;
        static const CVec3<Type> NaN;
    };

    template<class Type> const CVec3<Type> CVec3<Type>::UnitX(1, 0, 0);
    template<class Type> const CVec3<Type> CVec3<Type>::UnitY(0, 1, 0);
    template<class Type> const CVec3<Type> CVec3<Type>::UnitZ(0, 0, 1);
    template<class Type> const CVec3<Type> CVec3<Type>::NegaUnitX(-1, 0, 0);
    template<class Type> const CVec3<Type> CVec3<Type>::NegaUnitY(0, -1, 0);
    template<class Type> const CVec3<Type> CVec3<Type>::NegaUnitZ(0, 0, -1);
    template<class Type> const CVec3<Type> CVec3<Type>::Zero(0, 0, 0);
    template<class Type> const CVec3<Type> CVec3<Type>::NaN(g_Nan, g_Nan, g_Nan);

    typedef CVec3<int>    CVector3i;
    typedef CVec3<float>  CVector3f; 
    typedef CVec3<double> CVector3d;

    /*! @} */

#endif
