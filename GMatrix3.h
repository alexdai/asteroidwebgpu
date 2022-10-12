/*!
* @file     GMatrix3.h
* @brief    3*3矩阵的基本定义
* @details  矩阵类是一个3x3的模板类，提供对矩阵的基本操作方法外，同时提供平移，缩放，旋转和镜像等几何变化的常用矩阵的设置，向量的变换和矩阵的操作方法
            预定了类型CMatrix3f和CMatrix3d可以直接使用。
* @author
* @date      2021/09/13
* @copyright Copyright 2012-2022 GLODON
**********************************************************************************
* @par 修改日志:
* |  Date  | Description | Author     |
* | :----: | :----       | :----      |
* | 2021/09/13 | 袁冠杰<yuangj-a> | GGP-55193 完成算法剩余部分的注释 |
*
**********************************************************************************
*/
#ifndef G_MATRIX3_H
#define G_MATRIX3_H

#include "GMathDef.h"
#include "GVec3.h"

#pragma warning( disable:4251)

namespace ggp
{
    /*!\addtogroup GMath GMath
    * @{
    */

    /*!
    * @class CMatrix3
    * @brief 3*3矩阵的基本定义
             该类给出3*3矩阵的定义，提供相关的计算方法，被其它很多类使用。该模板类对整型数据类型不支持
    */
    template<typename Type>
    class  CMatrix3
    {
    public:
        /*!
        * @enum  ConstructMethod
        * @brief 
        */
        enum ConstructMethod { CM_ZERO, CM_IDENTITY };

        CMatrix3(ConstructMethod method = CM_IDENTITY)
        {
            if (method == CM_IDENTITY)
            {
                MakeIdentity();
            }
            else if (method == CM_ZERO)
            {
                MakeZero();
            }
        }

        /*!
        *@brief  根据9个数值，构造矩阵
        *@param[in] a11 第一行
        *@param[in] a12
        *@param[in] a13
        *@param[in] a21 第二行
        *@param[in] a22
        *@param[in] a23
        *@param[in] a31 第三行
        *@param[in] a32
        *@param[in] a33  
        */
        CMatrix3(const Type& a11, const Type& a12, const Type& a13,
            const Type& a21, const Type& a22, const Type& a23,
            const Type& a31, const Type& a32, const Type& a33)
        {
            m_data[0][0] = a11; m_data[0][1] = a12; m_data[0][2] = a13;
            m_data[1][0] = a21; m_data[1][1] = a22; m_data[1][2] = a23;
            m_data[2][0] = a31; m_data[2][1] = a32; m_data[2][2] = a33;
        }

        /*!
        *@brief  用三个行向量，构造矩阵    
        *@param[in] v1 行向量1
        *@param[in] v2 行向量2
        *@param[in] v3 行向量3
        */
        CMatrix3(const CVec3<Type> & v1, const CVec3<Type> & v2, const CVec3<Type> & v3)
        {
            m_data[0][0] = v1[0]; m_data[0][1] = v1[1]; m_data[0][2] = v1[2];
            m_data[1][0] = v2[0]; m_data[1][1] = v2[1]; m_data[1][2] = v2[2];
            m_data[2][0] = v3[0]; m_data[2][1] = v3[1]; m_data[2][2] = v3[2];
        }
        
        /*!
        *@brief  矩阵的拷贝构造
        *@param[in] matrix  传入的矩阵
        */
        CMatrix3(const CMatrix3<Type> & matrix)
        {
            *this = matrix;
        }

        /*!
        *@brief  令矩阵本身为单位阵
        */
        void MakeIdentity()
        {
            memset(m_data, 0, sizeof(Type) * 9);
            m_data[0][0]=m_data[1][1]=m_data[2][2]= (Type)1.0;
        }

        /*!
        *@brief  令矩阵本身为零矩阵
        */
        void MakeZero()
        {
            memset(m_data, 0, sizeof(Type) * 9);
        }
        
        /*!
        *@brief  判断矩阵是否为单位阵
        *@return    是否为单位阵
        * - true 是
        * - false 不是        
        */
        bool IsIdentity() const
        {
            CMatrix3<Type> Mat;
            return IsEqual(Mat, g_DoubleResolution);
        }

        /*!
        *@brief    用代数余子式计算3阶矩阵的行列式
        *@return   行列式的值
        */
        Type Det3() const
        {
            Type a11 = m_data[0][0];
            Type a12 = m_data[0][1];
            Type a13 = m_data[0][2];
            Type a21 = m_data[1][0];
            Type a22 = m_data[1][1];
            Type a23 = m_data[1][2];
            Type a31 = m_data[2][0];
            Type a32 = m_data[2][1];
            Type a33 = m_data[2][2];

            Type M11 = a22 * a33 - a32 * a23;
            Type M21 = -(a12 * a33 - a32 * a13);
            Type M31 = a12 * a23 - a22 * a13;

            return (a11 * M11 + a21 * M21 + a31 * M31);
        }

        // Inverse the current matrix.    
        /*!
        *@brief 求矩阵的逆
        *@return 逆矩阵
        */
        CMatrix3<Type> & Invert()
        {
            Type tmp[3][3];
            // Invert using cofactors.
            tmp[0][0] = m_data[1][1]*m_data[2][2] - m_data[1][2]*m_data[2][1];
            tmp[0][1] = m_data[0][2]*m_data[2][1] - m_data[0][1]*m_data[2][2];
            tmp[0][2] = m_data[0][1]*m_data[1][2] - m_data[0][2]*m_data[1][1];
            tmp[1][0] = m_data[1][2]*m_data[2][0] - m_data[1][0]*m_data[2][2];
            tmp[1][1] = m_data[0][0]*m_data[2][2] - m_data[0][2]*m_data[2][0];
            tmp[1][2] = m_data[0][2]*m_data[1][0] - m_data[0][0]*m_data[1][2];
            tmp[2][0] = m_data[1][0]*m_data[2][1] - m_data[1][1]*m_data[2][0];
            tmp[2][1] = m_data[0][1]*m_data[2][0] - m_data[0][0]*m_data[2][1];
            tmp[2][2] = m_data[0][0]*m_data[1][1] - m_data[0][1]*m_data[1][0];

            Type mInvDet = (Type)1.0 / (m_data[0][0]*tmp[0][0]+m_data[0][1]*tmp[1][0]+m_data[0][2]*tmp[2][0]);

            m_data[0][0] = tmp[0][0] * mInvDet; m_data[0][1] = tmp[0][1] * mInvDet; m_data[0][2] = tmp[0][2] * mInvDet;
            m_data[1][0] = tmp[1][0] * mInvDet; m_data[1][1] = tmp[1][1] * mInvDet; m_data[1][2] = tmp[1][2] * mInvDet;
            m_data[2][0] = tmp[2][0] * mInvDet; m_data[2][1] = tmp[2][1] * mInvDet; m_data[2][2] = tmp[2][2] * mInvDet;

            return *this;
        }


        //Return a new matrix which is the inverse matrix of this.    
        /*!
        *@brief 求矩阵的逆
        *@return 逆矩阵
        */
        CMatrix3<Type> Inverse() const
        {
            CMatrix3<Type> mat(*this);

            return mat.Invert();
        }

        /*!
        *@brief  把获得其转置矩阵
        *@return 转置矩阵
        */
        CMatrix3<Type> GetTranspose() const
        {
            return CMatrix3<Type>(m_data[0][0], m_data[1][0], m_data[2][0], 
                m_data[0][1], m_data[1][1], m_data[2][1], 
                m_data[0][2], m_data[1][2], m_data[2][2]);
        }

        /*!
        *@brief  把自身矩阵进行转置
        *@return 转置矩阵
        */
        CMatrix3<Type> & Transpose()
        {
            Type temp;
#define Swap(i,j) temp = m_data[i][j]; m_data[i][j] = m_data[j][i]; m_data[j][i] = temp;
            Swap(0, 1); Swap(0, 2);
            Swap(1, 2); 

            return *this;
#undef Swap
        }

       
        /*!
        *@brief  判断自身矩阵和传入矩阵是否在给定误差下相等
        *@param[in] matrix 传入矩阵
        *@param[in] a_tolerance  给定误差
        *@return    是否相等
        * - true 是
        * - false 不是
        */
        bool IsEqual(const CMatrix3<Type> & matrix, double a_tolerance  = g_DoubleResolution) const
        {
            for(int i = 0; i < 3; i++)
            {
                for(int j = 0; j < 3; j++)
                {
                    if (std::abs(m_data[i][j] - matrix.m_data[i][j]) > a_tolerance)
                        return false;
                }
            }
            return true;
        }

         
        /*!
        *@brief  返回给定行和列的值
        *@param[in] row  给定行
        *@param[in] col  给定列
        *@return    值
        */
        inline const Type& Value(int row, int col) const 
        { 
            return m_data[row][col]; 
        }

        /*!
        *@brief  设置给定行和列的值
        *@param[in] row  给定行
        *@param[in] col  给定列
        *@param[in] Val 值 
        */ 
        inline void Set(int row, int col,const Type& Val)  
        { 
            m_data[row][col] = Val; 
        }


        /*!
        *@brief  返回给定行和列的值
        *@param[in] row  给定行
        *@param[in] col  给定列
        *@return    值
        */
        inline const Type& operator()(int row, int col) const 
        { 
            return m_data[row][col];
        }

        /*!
        *@brief  设置给定行和列的值
        *@param[in] row  给定行
        *@param[in] col  给定列
        *@return 值
        */ 
        inline Type& operator()(int row, int col)  
        { 
            return m_data[row][col];
        }

       
        /*!
        *@brief  返回二维数组的常指针
        *@return 二维数组的常指针
        */ 
        inline const Type * Ptr() const
        { 
            return &(m_data[0][0]); 
        }

        /*!
        *@brief  返回二维数组的指针
        *@return 二维数组的指针
        */
        inline Type* Ptr()
        {
            return &(m_data[0][0]); 
        }

       /*!
        *@brief  返回二维数组中某行的指针
        *@param[in] i  给定行
        *@return  行指针
        */ 
        Type * operator [](int i)
        { 
            return m_data[i]; 
        }

         /*!
        *@brief  返回二维数组中某行的指针
        *@param[in] i  给定行
        *@return  行指针
        */ 
        const Type * operator [](int i) const
        { 
            return m_data[i]; 
        }

        
       
        /*!
        *@brief  矩阵赋值操作
        *@param[in] m  传入矩阵
        *@return    自身矩阵
        */
        CMatrix3<Type> & operator =(const CMatrix3<Type> & m)
        {
            memmove(m_data, m.m_data, sizeof(Type) * 9);
            return *this;
        }

        /*!
        *@brief  自身矩阵乘以传入矩阵，结果设置为自身
        *@param[in] matrix   传入矩阵
        *@return    相乘后结果
        *@sa MultRight(const CMatrix3<Type> & m)
        */
        CMatrix3<Type> & operator *=(const CMatrix3<Type> & matrix)
        {
            return MultRight(matrix);
        }

        /*!
        *@brief  两矩阵相乘，返回结果矩阵
        *@param[in] m1
        *@param[in] m2
        *@return    结果矩阵
        */
        friend CMatrix3<Type> operator *(const CMatrix3<Type> & m1, const CMatrix3<Type> & m2)
        { 
            CMatrix3<Type> m = m1; 
            m *= m2; 
            return m; 
        }

        /*!
        *@brief  矩阵乘以列向量，返回列向量
        *@param[in] m   3阶矩阵 
        *@param[in] v      列向量
        *@return    结果向量
        */
        friend CVec3<Type> operator *(const CMatrix3<Type> & m, const CVec3<Type> & v)
        {
            CVec3<Type> vRet;

            vRet[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
            vRet[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
            vRet[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];

            return vRet;
        }

        /*!
        *@brief  行向量乘以矩阵，返回行向量
        *@param[in] v      行向量
        *@param[in] m   3阶矩阵 
        *@return    结果向量
        */
        friend CVec3<Type> operator *(const CVec3<Type> & v, const CMatrix3<Type> & m)
        {
            CVec3<Type> vRet;

            vRet[0] = v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0];
            vRet[1] = v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1];
            vRet[2] = v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2];

            return vRet;
        }

        /*!
        *@brief  判断两矩阵是否相等
        *@param[in] m1   矩阵1
        *@param[in] m2   矩阵2
        *@return    是否相等
        * - true 是
        * - false 不是
        */
        friend bool operator ==(const CMatrix3<Type> & m1, const CMatrix3<Type> & m2)
        {
            return (
                m1.m_data[0][0] == m2.m_data[0][0] &&
                m1.m_data[0][1] == m2.m_data[0][1] &&
                m1.m_data[0][2] == m2.m_data[0][2] &&

                m1.m_data[1][0] == m2.m_data[1][0] &&
                m1.m_data[1][1] == m2.m_data[1][1] &&
                m1.m_data[1][2] == m2.m_data[1][2] &&

                m1.m_data[2][0] == m2.m_data[2][0] &&
                m1.m_data[2][1] == m2.m_data[2][1] &&
                m1.m_data[2][2] == m2.m_data[2][2]
            );
        }

        /*!
        *@brief  判断两矩阵是否不相等
        *@param[in] m1   矩阵1
        *@param[in] m2   矩阵2
        *@return    是否不相等
        * - true 是
        * - false 不是
        */
        friend bool operator !=(const CMatrix3<Type> & m1, const CMatrix3<Type> & m2)
        { 
            return !(m1 == m2);
        }


        /*!
        *@brief  构造平移矩阵
        *@param[in] t  平移量
        *@return    平移矩阵
        */
        static inline CMatrix3<Type> Translate(const CVec2<Type> & t)
        {
            CMatrix3<Type> m;
            m.SetTranslate(t);

            return m;
        }

        /*!
        *@brief  根据平移向量的二个分量返回构造的平移矩阵
        *@param[in] tx  平移X分量
        *@param[in] ty  平移Y分量
        *@return    平移矩阵
        */
        static inline CMatrix3<Type> Translate(const Type& tx,const Type& ty)
        {
            CMatrix3<Type> m;
            m.SetTranslate(tx,ty);

            return m;
        }
        
        /*!
        *@brief  返回绕任意点的旋转矩阵
        *@param[in] Pos   中心点
        *@param[in] Angle        旋转角度
        *@return    旋转矩阵
        */
        static inline CMatrix3<Type> Rotate(const CVec2<Type>& Pos, const Type& Angle)
        {
            CMatrix3<Type> m;
            m.SetRotate(Pos, Angle);

            return m;
        }

        /*!
        *@brief  返回构造的等比缩放矩阵
        *@param[in] Pos    基点
        *@param[in] Scale         缩放比例
        *@return    比例矩阵
        */
        static inline CMatrix3<Type> Scale(const CVec2<Type>& Pos, const Type& Scale)
        {
            CMatrix3<Type> m;
            m.SetScale(Pos,Scale);

            return m;
        }

        /*!
        *@brief  返回构造的非等比缩放矩阵
        *@param[in] Pos    基点
        *@param[in] Scales 缩放比例
        *@return    比例矩阵
        */
        static inline CMatrix3<Type> Scale(const CVec2<Type>& Pos, const CVec2<Type>& Scales)
        {
            CMatrix3<Type> m;
            m.SetScale(Pos,Scales);

            return m;
        }

        /*!
        *@brief  返回构造任意轴的镜像矩阵
        *@param[in] Pos   镜像轴的基点
        *@param[in] Dir   镜像轴的方向
        *@return    镜像矩阵
        */
        static inline CMatrix3<Type> Mirror(const CVec2<Type>& Pos, const CVec2<Type>& Dir)
        {
            CMatrix3<Type> m;
            m.SetMirror(Pos, Dir);

            return m;
        }


        /*!
        *@brief  设置自身矩阵为平移矩阵
        *@param[in] t  平移量
        */
        inline void SetTranslate(const CVec2<Type> & t)
        {
            MakeIdentity();
            m_data[2][0] = t[0];
            m_data[2][1] = t[1];
        }

        /*!
        *@brief  根据平移向量的二个分量设置自身矩阵为平移矩阵
        *@param[in] tx  平移X分量
        *@param[in] ty  平移Y分量
        */
        inline void SetTranslate(const Type& tx,const Type& ty)
        {
            MakeIdentity();
            m_data[2][0] = tx;
            m_data[2][1] = ty;
        }

        /*!
        *@brief  设置自身矩阵为绕任意点的旋转矩阵
        *@param[in] Pos   中心点
        *@param[in] Angle        旋转角度的弧度值，正值表示逆时针方向
        */
        void SetRotate(const CVec2<Type>& Pos, const Type& Angle)
        {
            Type cosVal = cos(Angle);
            Type sinVal = sin(Angle);

            Type dx = (1-cosVal)*Pos.X + sinVal*Pos.Y;
            Type dy = (1-cosVal)*Pos.Y - sinVal*Pos.X;

            m_data[0][0] = cosVal;      m_data[0][1] = sinVal;      m_data[0][2] = 0.0;
            m_data[1][0] = -sinVal;     m_data[1][1] = cosVal;      m_data[1][2] = 0.0;
            m_data[2][0] = dx;          m_data[2][1] = dy;          m_data[2][2] = 1.0;    
        }
       
        /*!
        *@brief  设置自身矩阵等比缩放矩阵
        *@param[in] Pos    基点
        *@param[in] Scale         缩放比例
        */
        void SetScale(const CVec2<Type>& Pos, const Type& Scale)
        {
            assert(Scale > 0.0);

//             if(Scale<0)
//                 Scale *= -1;

            MakeIdentity();

            m_data[0][0] = Scale;
            m_data[1][1] = Scale;
            m_data[2][0] = Pos.X * ((Type)1.0 - Scale);
            m_data[2][1] = Pos.Y * ((Type)1.0 - Scale);
        }

        /*!
        *@brief  设置自身矩阵非等比缩放矩阵
        *@param[in] Pos   基点
        *@param[in] Scales 缩放比例
        */
        void SetScale(const CVec2<Type>& Pos, const CVec2<Type>& Scales)
        {
            assert(Scales.X != 0.0 && Scales.Y!=0.0);

            MakeIdentity();

            m_data[0][0] = Scales.X;
            m_data[1][1] = Scales.Y;
            m_data[2][0] = Pos.X * ((Type)1.0 - Scales.X);   
            m_data[2][1] = Pos.Y * ((Type)1.0 - Scales.Y);     
        }

        /*!
        *@brief  设置自身矩阵为任意轴的镜像矩阵
        *@param[in] Pos   镜像轴的基点
        *@param[in] Dir   镜像轴的方向
        */
        void SetMirror(const CVec2<Type>& Pos, const CVec2<Type>& Dir)
        {
            CVec2<Type> DirX = Dir;
            CVec2<Type> DirY;
            Type Length = DirX.Normalize();
            if (IsNearZero(Length))
            {
                DirX.Set((Type)1.0, (Type)0.0);
                DirY.Set((Type)0.0, (Type)1.0);
            }
            else
            {
                DirY.Set(-DirX.Y, DirX.X);
            }

            CMatrix3<Type> MatTmp(DirX.X, DirX.Y, (Type)0.0,
                                  DirY.X, DirY.Y, (Type)0.0,
                                  Pos.X,  Pos.Y,  (Type)1.0);

            CMatrix3<Type> MatTmpInverse = MatTmp.Inverse();

            CMatrix3<Type> MirrorXMatTmp;
            MirrorXMatTmp.Set(1,1,(Type)(-1.0));

            SetProduct(MatTmpInverse,MirrorXMatTmp);
            MultRight(MatTmp);
        }

        /*!
        *@brief  两个矩阵相乘，结果设置为自身矩阵
        *@param[in] m1  矩阵1
        *@param[in] m2  矩阵2
        *@return    矩阵相乘结果
        */
        CMatrix3<Type>&  SetProduct(const CMatrix3<Type> & m1,const CMatrix3<Type> & m2)
        {
            (*this) = m1;
            return MultRight(m2);
        }

       /*!
        *@brief  传入矩阵右乘当前矩阵，结果设置为自身矩阵
        *@param[in] m  传入矩阵
        *@return    结果
        */
        CMatrix3<Type> & MultRight(const CMatrix3<Type> & m)
        {
            // Trivial cases
            if(m.IsIdentity()) 
                return *this;
            else if(IsIdentity())
                return (*this = m);

            CMatrix3<Type> tmp;
            

            tmp[0][0] = m_data[0][0]*m.m_data[0][0] + m_data[0][1]*m.m_data[1][0] + m_data[0][2]*m.m_data[2][0];
            tmp[0][1] = m_data[0][0]*m.m_data[0][1] + m_data[0][1]*m.m_data[1][1] + m_data[0][2]*m.m_data[2][1];
            tmp[0][2] = m_data[0][0]*m.m_data[0][2] + m_data[0][1]*m.m_data[1][2] + m_data[0][2]*m.m_data[2][2];

            tmp[1][0] = m_data[1][0]*m.m_data[0][0] + m_data[1][1]*m.m_data[1][0] + m_data[1][2]*m.m_data[2][0];
            tmp[1][1] = m_data[1][0]*m.m_data[0][1] + m_data[1][1]*m.m_data[1][1] + m_data[1][2]*m.m_data[2][1];
            tmp[1][2] = m_data[1][0]*m.m_data[0][2] + m_data[1][1]*m.m_data[1][2] + m_data[1][2]*m.m_data[2][2];

            tmp[2][0] = m_data[2][0]*m.m_data[0][0] + m_data[2][1]*m.m_data[1][0] + m_data[2][2]*m.m_data[2][0];
            tmp[2][1] = m_data[2][0]*m.m_data[0][1] + m_data[2][1]*m.m_data[1][1] + m_data[2][2]*m.m_data[2][1];
            tmp[2][2] = m_data[2][0]*m.m_data[0][2] + m_data[2][1]*m.m_data[1][2] + m_data[2][2]*m.m_data[2][2];


            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    m_data[i][j] = tmp[i][j];
                }
            }
            
            return *this;
        }

        /*!
        *@brief  二维点src左乘矩阵
        *@param[in] src  传入的二维点
        *@return    返回的二维点
        */
        CVec2<Type> MultiPointLeft(const CVec2<Type> & src) const
        {
            Type x = src[0]*m_data[0][0] + src[1]*m_data[1][0] + m_data[2][0];
            Type y = src[0]*m_data[0][1] + src[1]*m_data[1][1] + m_data[2][1];
            Type w = src[0]*m_data[0][2] + src[1]*m_data[1][2] + m_data[2][2];

            CVec2<Type> dst;
            dst.Set(x/w, y/w);

            return dst;
        }

        /*!
        *@brief  二维行向量src左乘矩阵,无视偏移信息
        *@param[in] src  传入的行向量
        *@return  返回的向量
        */
        CVec2<Type> MultiVecLeft(const CVec2<Type> & src) const
        {
            Type x = src[0]*m_data[0][0] + src[1]*m_data[1][0];
            Type y = src[0]*m_data[0][1] + src[1]*m_data[1][1];

            CVec2<Type> dst;
            dst.Set(x, y);

            return dst;
        }

        /*!
        *@brief     三维点src右乘矩阵，注意若是二维点转三维点应该用齐次坐标！！！
        *@param[in] v 传入的三维点
        *@return    右乘结果
        */
        CVec3<Type> MultRight(const CVec3<Type> & v)
        {
            Type x = m_data[0][0]*v[0] + m_data[0][1]*v[1] + m_data[0][2]*v[2];
            Type y = m_data[1][0]*v[0] + m_data[1][1]*v[1] + m_data[1][2]*v[2];
            Type z = m_data[2][0]*v[0] + m_data[2][1]*v[1] + m_data[2][2]*v[2];
            return CVec3<Type>(x, y, z);
        }

        /*!
        *@brief  三维点src左乘矩阵，注意若是二维点转三维点应该用齐次坐标！！！
        *@param[in] src  传入的三维点
        *@return    返回的三维点
        */
        CVec3<Type> MultiLeft(const CVec3<Type> & src) const
        {
            Type x = src[0]*m_data[0][0] + src[1]*m_data[1][0] + src[2]*m_data[2][0];
            Type y = src[0]*m_data[0][1] + src[1]*m_data[1][1] + src[2]*m_data[2][1];
            Type z = src[0]*m_data[0][2] + src[1]*m_data[1][2] + src[2]*m_data[2][2];

            return CVec3<Type> (x, y, z);
        }

        /*!
        *@brief 将矩阵流化到输出流
        *@param[out] os 流对象
        *@param[in] mat 给定矩阵
        */
        friend std::ostream & operator<<(std::ostream & os, const CMatrix3<Type> & mat)
        { 
            for(unsigned int i=0; i<3; i++){
                os << mat[i][0] <<'\t'<< mat[i][1] <<'\t'<< mat[i][2] << std::endl;
            }
            return os;
        }

        friend std::istream & operator>>(std::istream & iStream, CMatrix3<Type> & mat)
        {
            for(unsigned int i=0; i<3; i++){
                iStream >> mat[i][0] >> mat[i][1] >> mat[i][2];
            }
            return iStream;
        }

        /*!
        *@brief  矩阵对象按照规定格式字符串输出
        *@return    字符串
        */
        char* AsString()
        {
            char * a = AsStringHelper::NewString("Matrix3=[(%.20le %.20le %.20le),(%.20le %.20le %.20le),(%.20le %.20le %.20le)]", 
                m_data[0][0], m_data[0][1], m_data[0][2], 
                m_data[1][0], m_data[1][1], m_data[1][2], 
                m_data[2][0], m_data[2][1], m_data[2][2]);
            return a;
        }

        /*!
        *@brief     根据字符串生成矩阵
        *@param[in] str   传入的字符串
        *@return    生成的矩阵
        */
        static CMatrix3<double>* LoadFromString(char* str)
        {
            CMatrix3<double>* pMatrix = new CMatrix3<double>;
            CMatrix3<double>& matrix = *pMatrix;
            g_sscanf(str, "Matrix3=[(%le %le %le),(%le %le %le),(%le %le %le)]", 
                &matrix(0, 0), &matrix(0, 1), &matrix(0, 2), 
                &matrix(1, 0), &matrix(1, 1), &matrix(1, 2), 
                &matrix(2, 0), &matrix(2, 1), &matrix(2, 2));

            return pMatrix;
        }

        /*!
        * @brief 返回矩阵中的旋转部分对应的坐标轴分量。
        * @param[out] xAxis x轴
        * @param[out] yAxis y轴
        */
        void GetRotationAxes(CVec2<Type>& xAxis, CVec2<Type>& yAxis) const
        {
            xAxis = CVec3<Type>(m_data[0][0], m_data[0][1]);
            yAxis = CVec3<Type>(m_data[1][0], m_data[1][1]);
        }

        /*!
        * @brief 基于向量的1范数
        * @return 1范数
        */
        Type Norm1() const
        {
            Type s1 = abs(m_data[0][0]) + abs(m_data[1][0]) + abs(m_data[2][0]);
            Type s2 = abs(m_data[0][1]) + abs(m_data[1][1]) + abs(m_data[2][1]);
            Type s3 = abs(m_data[0][2]) + abs(m_data[1][2]) + abs(m_data[2][2]);
            return (s1 + s2 + s3);
        }

        /*!
        * @brief 两个三维向量相乘获得3 x 3矩阵。
        * @param[in]  u 表示3 x 1向量
        * @param[in]  v 表示1 x 3向量
        * @return 相乘后结果
        */

        static CMatrix3<Type> DirectProduct(const CVec3<Type>& u, const CVec3<Type>& v)
        {
            return CMatrix3<Type>(v * u.X, v * u.Y, v * u.Z);
        }


    private:
        Type m_data[3][3];

    public:
        static const CMatrix3<Type> Zero;
        static const CMatrix3<Type> Identity;
    };

    template<class T> const CMatrix3<T> CMatrix3<T>::Zero(CMatrix3<T>::CM_ZERO);
    template<class T> const CMatrix3<T> CMatrix3<T>::Identity(CMatrix3<T>::CM_IDENTITY);

    typedef CMatrix3<float>  CMatrix3f;
    typedef CMatrix3<double> CMatrix3d;

    /*! @} */
} // namespace

#endif
