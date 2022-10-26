#ifndef G_MATRIX4_H
#define G_MATRIX4_H

#include "GMathDef.h"
#include "GVec3.h"
#include "GVec4.h"
#include "GMatrix3.h"
#include <cfloat> // for FLT_EPSILON

#pragma warning( disable:4251)


    template<typename Type>
    class  CMatrix4
    {
    public:
        enum ConstructMethod { CM_ZERO, CM_IDENTITY };

        CMatrix4(ConstructMethod method = CM_IDENTITY)
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

        CMatrix4(const Type& a11, const Type& a12, const Type& a13, const Type& a14,
                 const Type& a21, const Type& a22, const Type& a23, const Type& a24,
                 const Type& a31, const Type& a32, const Type& a33, const Type& a34,
                 const Type& a41, const Type& a42, const Type& a43, const Type& a44)
        {
            m_data[0][0] = a11; m_data[0][1] = a12; m_data[0][2] = a13; m_data[0][3] = a14;
            m_data[1][0] = a21; m_data[1][1] = a22; m_data[1][2] = a23; m_data[1][3] = a24;
            m_data[2][0] = a31; m_data[2][1] = a32; m_data[2][2] = a33; m_data[2][3] = a34;
            m_data[3][0] = a41; m_data[3][1] = a42; m_data[3][2] = a43; m_data[3][3] = a44;
        }


        CMatrix4(const CVec3<Type> & r1, const CVec3<Type> & r2, const CVec3<Type> & r3, const CVec3<Type> & t)
        {
            m_data[0][0] = r1[0]; m_data[0][1] = r1[1]; m_data[0][2] = r1[2]; m_data[0][3] = 0;
            m_data[1][0] = r2[0]; m_data[1][1] = r2[1]; m_data[1][2] = r2[2]; m_data[1][3] = 0;
            m_data[2][0] = r3[0]; m_data[2][1] = r3[1]; m_data[2][2] = r3[2]; m_data[2][3] = 0;
            m_data[3][0] = t[0];  m_data[3][1] = t[1];  m_data[3][2] = t[2];  m_data[3][3] = 1;//左乘
        }

        CMatrix4(const CVec4<Type> & r1, const CVec4<Type> & r2, const CVec4<Type> & r3, const CVec4<Type> & r4)
        {
            m_data[0][0] = r1[0]; m_data[0][1] = r1[1]; m_data[0][2] = r1[2]; m_data[0][3] = r1[3];
            m_data[1][0] = r2[0]; m_data[1][1] = r2[1]; m_data[1][2] = r2[2]; m_data[1][3] = r2[3];
            m_data[2][0] = r3[0]; m_data[2][1] = r3[1]; m_data[2][2] = r3[2]; m_data[2][3] = r3[3];
            m_data[3][0] = r4[0]; m_data[3][1] = r4[1]; m_data[3][2] = r4[2]; m_data[3][3] = r4[3];
        }

        CMatrix4(const CMatrix3<Type> & R, const CVec3<Type> & t)
        {
            m_data[0][0] = R[0][0]; m_data[0][1] = R[0][1]; m_data[0][2] = R[0][2]; m_data[0][3] = 0;
            m_data[1][0] = R[1][0]; m_data[1][1] = R[1][1]; m_data[1][2] = R[1][2]; m_data[1][3] = 0;
            m_data[2][0] = R[2][0]; m_data[2][1] = R[2][1]; m_data[2][2] = R[2][2]; m_data[2][3] = 0;
            m_data[3][0] = t[0];    m_data[3][1] = t[1];    m_data[3][2] = t[2];    m_data[3][3] = 1;
        }

        CMatrix4(const CMatrix3<Type> & m)
        {
            m_data[0][0] = m[0][0]; m_data[0][1] = m[0][1]; m_data[0][2] = 0; m_data[0][3] = m[0][2];
            m_data[1][0] = m[1][0]; m_data[1][1] = m[1][1]; m_data[1][2] = 0; m_data[1][3] = m[1][2];
            m_data[2][0] = 0;       m_data[2][1] = 0;       m_data[2][2] = 1; m_data[2][3] = 0;
            m_data[3][0] = m[2][0]; m_data[3][1] = m[2][1]; m_data[3][2] = 0; m_data[3][3] = m[2][2];
        }

        CMatrix4(const CMatrix4<Type> & matrix)
        {
            *this = matrix;
        }

        CMatrix3<Type> Matrix3() const
        {
            CMatrix3<Type> mat3;
            for(int i = 0 ; i < 3; ++i)
                for(int j = 0; j < 3; ++j)
                    mat3[i][j] = m_data[i][j];

            return mat3;
        }

        CMatrix4<double> Matrix4d() const
        {
            CMatrix4<double> mat;
            mat[0][0] = m_data[0][0];
            mat[0][1] = m_data[0][1];
            mat[0][2] = m_data[0][2];
            mat[0][3] = m_data[0][3];

            mat[1][0] = m_data[1][0];
            mat[1][1] = m_data[1][1];
            mat[1][2] = m_data[1][2];
            mat[1][3] = m_data[1][3];

            mat[2][0] = m_data[2][0];
            mat[2][1] = m_data[2][1];
            mat[2][2] = m_data[2][2];
            mat[2][3] = m_data[2][3];

            mat[3][0] = m_data[3][0];
            mat[3][1] = m_data[3][1];
            mat[3][2] = m_data[3][2];
            mat[3][3] = m_data[3][3];

            return mat;
        }

        CMatrix4<float> Matrix4f() const
        {
            CMatrix4<float> mat;
            mat[0][0] = (float)m_data[0][0];
            mat[0][1] = (float)m_data[0][1];
            mat[0][2] = (float)m_data[0][2];
            mat[0][3] = (float)m_data[0][3];

            mat[1][0] = (float)m_data[1][0];
            mat[1][1] = (float)m_data[1][1];
            mat[1][2] = (float)m_data[1][2];
            mat[1][3] = (float)m_data[1][3];

            mat[2][0] = (float)m_data[2][0];
            mat[2][1] = (float)m_data[2][1];
            mat[2][2] = (float)m_data[2][2];
            mat[2][3] = (float)m_data[2][3];

            mat[3][0] = (float)m_data[3][0];
            mat[3][1] = (float)m_data[3][1];
            mat[3][2] = (float)m_data[3][2];
            mat[3][3] = (float)m_data[3][3];

            return mat;
        }


        void MakeIdentity()
        {
            memset(m_data, 0, sizeof(Type) * 16);
            m_data[0][0]=m_data[1][1]=m_data[2][2]=m_data[3][3] = (Type)1.0;
        }

        void MakeZero()
        {
            memset(m_data, 0, sizeof(Type) * 16);
        }

        bool IsIdentity() const
        {
            auto eps = 2.0*std::numeric_limits<Type>::epsilon();
            if(std::abs(m_data[0][0] - 1.0) > eps)
            {
                return false;
            }

            if(std::abs(m_data[1][1] - 1.0) > eps)
            {
                return false;
            }

            if(std::abs(m_data[2][2] - 1.0) > eps)
            {
                return false;
            }

            if(std::abs(m_data[3][3] - 1.0) > eps)
            {
                return false;
            }

            auto sum = m_data[0][1]*m_data[0][1] + m_data[0][2]*m_data[0][2] + m_data[0][3]*m_data[0][3]
                     + m_data[1][0]*m_data[1][0] + m_data[1][2]*m_data[1][2] + m_data[1][3]*m_data[1][3]
                     + m_data[2][0]*m_data[2][0] + m_data[2][1]*m_data[2][1] + m_data[2][3]*m_data[2][3]
                     + m_data[3][0]*m_data[3][0] + m_data[3][1]*m_data[3][1] + m_data[3][2]*m_data[3][2];

            if(sum > eps*eps)
            {
                return false;
            }

            return true;
        }

        inline bool IsMirror()const
        {
            CVec3<Type> &vx = *(CVec3<Type>*)m_data[0];
            CVec3<Type> &vy = *(CVec3<Type>*)m_data[1];
            CVec3<Type> &vz = *(CVec3<Type>*)m_data[2];;

            CVec3<Type> vc = vx.Cross(vy);

            Type d = vc.Dot(vz);

            if (d < (Type)0.0)
            {
                return true;
            }

            return false;
        }

        inline bool IsXOYMirror() const
        {
            CVec2<Type> vx (m_data[0][0], m_data[0][1]);
            CVec2<Type> vy (m_data[1][0], m_data[1][1]);

            return vx.Cross(vy) < (Type)0.0;
        }

        CMatrix4<Type> & Invert(bool* pSuccess = NULL)
        {

            Type tmp[12];
            Type dst[16];

            Transpose();

            tmp[0]  = m_data[2][2] * m_data[3][3];
            tmp[1]  = m_data[3][2] * m_data[2][3];
            tmp[2]  = m_data[1][2] * m_data[3][3];
            tmp[3]  = m_data[3][2] * m_data[1][3];
            tmp[4]  = m_data[1][2] * m_data[2][3];
            tmp[5]  = m_data[2][2] * m_data[1][3];
            tmp[6]  = m_data[0][2] * m_data[3][3];
            tmp[7]  = m_data[3][2] * m_data[0][3];
            tmp[8]  = m_data[0][2] * m_data[2][3];
            tmp[9]  = m_data[2][2] * m_data[0][3];
            tmp[10] = m_data[0][2] * m_data[1][3];
            tmp[11] = m_data[1][2] * m_data[0][3];


            dst[0]  = tmp[0]*m_data[1][1] + tmp[3]*m_data[2][1] + tmp[4]*m_data[3][1];
            dst[0] -= tmp[1]*m_data[1][1] + tmp[2]*m_data[2][1] + tmp[5]*m_data[3][1];
            dst[1]  = tmp[1]*m_data[0][1] + tmp[6]*m_data[2][1] + tmp[9]*m_data[3][1];
            dst[1] -= tmp[0]*m_data[0][1] + tmp[7]*m_data[2][1] + tmp[8]*m_data[3][1];
            dst[2]  = tmp[2]*m_data[0][1] + tmp[7]*m_data[1][1] + tmp[10]*m_data[3][1];
            dst[2] -= tmp[3]*m_data[0][1] + tmp[6]*m_data[1][1] + tmp[11]*m_data[3][1];
            dst[3]  = tmp[5]*m_data[0][1] + tmp[8]*m_data[1][1] + tmp[11]*m_data[2][1];
            dst[3] -= tmp[4]*m_data[0][1] + tmp[9]*m_data[1][1] + tmp[10]*m_data[2][1];
            dst[4]  = tmp[1]*m_data[1][0] + tmp[2]*m_data[2][0] + tmp[5]*m_data[3][0];
            dst[4] -= tmp[0]*m_data[1][0] + tmp[3]*m_data[2][0] + tmp[4]*m_data[3][0];
            dst[5]  = tmp[0]*m_data[0][0] + tmp[7]*m_data[2][0] + tmp[8]*m_data[3][0];
            dst[5] -= tmp[1]*m_data[0][0] + tmp[6]*m_data[2][0] + tmp[9]*m_data[3][0];
            dst[6]  = tmp[3]*m_data[0][0] + tmp[6]*m_data[1][0] + tmp[11]*m_data[3][0];
            dst[6] -= tmp[2]*m_data[0][0] + tmp[7]*m_data[1][0] + tmp[10]*m_data[3][0];
            dst[7]  = tmp[4]*m_data[0][0] + tmp[9]*m_data[1][0] + tmp[10]*m_data[2][0];
            dst[7] -= tmp[5]*m_data[0][0] + tmp[8]*m_data[1][0] + tmp[11]*m_data[2][0];


            tmp[0]  = m_data[2][0]*m_data[3][1];
            tmp[1]  = m_data[3][0]*m_data[2][1];
            tmp[2]  = m_data[1][0]*m_data[3][1];
            tmp[3]  = m_data[3][0]*m_data[1][1];
            tmp[4]  = m_data[1][0]*m_data[2][1];
            tmp[5]  = m_data[2][0]*m_data[1][1];
            tmp[6]  = m_data[0][0]*m_data[3][1];
            tmp[7]  = m_data[3][0]*m_data[0][1];
            tmp[8]  = m_data[0][0]*m_data[2][1];
            tmp[9]  = m_data[2][0]*m_data[0][1];
            tmp[10] = m_data[0][0]*m_data[1][1];
            tmp[11] = m_data[1][0]*m_data[0][1];


            dst[8]   = tmp[0]*m_data[1][3] + tmp[3]*m_data[2][3] + tmp[4]*m_data[3][3];
            dst[8]  -= tmp[1]*m_data[1][3] + tmp[2]*m_data[2][3] + tmp[5]*m_data[3][3];
            dst[9]   = tmp[1]*m_data[0][3] + tmp[6]*m_data[2][3] + tmp[9]*m_data[3][3];
            dst[9]  -= tmp[0]*m_data[0][3] + tmp[7]*m_data[2][3] + tmp[8]*m_data[3][3];
            dst[10]  = tmp[2]*m_data[0][3] + tmp[7]*m_data[1][3] + tmp[10]*m_data[3][3];
            dst[10] -= tmp[3]*m_data[0][3] + tmp[6]*m_data[1][3] + tmp[11]*m_data[3][3];
            dst[11]  = tmp[5]*m_data[0][3] + tmp[8]*m_data[1][3] + tmp[11]*m_data[2][3];
            dst[11] -= tmp[4]*m_data[0][3] + tmp[9]*m_data[1][3] + tmp[10]*m_data[2][3];
            dst[12]  = tmp[2]*m_data[2][2] + tmp[5]*m_data[3][2] + tmp[1]*m_data[1][2];
            dst[12] -= tmp[4]*m_data[3][2] + tmp[0]*m_data[1][2] + tmp[3]*m_data[2][2];
            dst[13]  = tmp[8]*m_data[3][2] + tmp[0]*m_data[0][2] + tmp[7]*m_data[2][2];
            dst[13] -= tmp[6]*m_data[2][2] + tmp[9]*m_data[3][2] + tmp[1]*m_data[0][2];
            dst[14]  = tmp[6]*m_data[1][2] + tmp[11]*m_data[3][2] + tmp[3]*m_data[0][2];
            dst[14] -= tmp[10]*m_data[3][2] + tmp[2]*m_data[0][2] + tmp[7]*m_data[1][2];
            dst[15]  = tmp[10]*m_data[2][2] + tmp[4]*m_data[0][2] + tmp[9]*m_data[1][2];
            dst[15] -= tmp[8]*m_data[1][2] + tmp[11]*m_data[2][2] + tmp[5]*m_data[0][2];

            Type d = m_data[0][0]*dst[0] + m_data[1][0]*dst[1] + m_data[2][0]*dst[2] + m_data[3][0]*dst[3];
            //if (!IsNearZero(d))
            if ( d != 0)
            {
                if (pSuccess)
                {
                    *pSuccess = true;
                }

                Type det = (Type)1.0 / d;

                m_data[0][0] = dst[0] * det; m_data[0][1] = dst[4] * det; m_data[0][2] = dst[8] * det; m_data[0][3] = dst[12] * det;
                m_data[1][0] = dst[1] * det; m_data[1][1] = dst[5] * det; m_data[1][2] = dst[9] * det; m_data[1][3] = dst[13] * det;
                m_data[2][0] = dst[2] * det; m_data[2][1] = dst[6] * det; m_data[2][2] = dst[10] * det; m_data[2][3] = dst[14] * det;
                m_data[3][0] = dst[3] * det; m_data[3][1] = dst[7] * det; m_data[3][2] = dst[11] * det; m_data[3][3] = dst[15] * det;

                return *this;
            }

            if (pSuccess)
            {
                *pSuccess = false; 
            }

            Transpose();
            return *this;

        }




        CMatrix4<Type> Inverse(bool* pSuccess = NULL) const
        {
            CMatrix4<Type> mat(*this);

            return mat.Invert(pSuccess);
        }


        CMatrix4<Type> & Transpose()
        {
            Type temp;
#define Swap(i,j) temp = m_data[i][j]; m_data[i][j] = m_data[j][i]; m_data[j][i] = temp;
            Swap(0, 1); Swap(0, 2); Swap(0, 3);
            Swap(1, 2); Swap(1, 3); 
            Swap(2, 3);

            return *this;
#undef Swap
        }

        CMatrix4<Type> GetTranspose() const
        {
            return CMatrix4<Type>(m_data[0][0], m_data[1][0], m_data[2][0], m_data[3][0],
                m_data[0][1], m_data[1][1], m_data[2][1], m_data[3][1],
                m_data[0][2], m_data[1][2], m_data[2][2], m_data[3][2],
                m_data[0][3], m_data[1][3], m_data[2][3], m_data[3][3]);
        }


        Type Det3(int r1, int r2, int r3, int c1, int c2, int c3) const
        {

            Type a11 = m_data[r1][c1];
            Type a12 = m_data[r1][c2];
            Type a13 = m_data[r1][c3];
            Type a21 = m_data[r2][c1];
            Type a22 = m_data[r2][c2];
            Type a23 = m_data[r2][c3];
            Type a31 = m_data[r3][c1];
            Type a32 = m_data[r3][c2];
            Type a33 = m_data[r3][c3];

            Type M11 = a22 * a33 - a32 * a23;
            Type M21 = -(a12 * a33 - a32 * a13);
            Type M31 = a12 * a23 - a22 * a13;

            return (a11 * M11 + a21 * M21 + a31 * M31);
        }



        Type Det4() const
        {
            Type det = Type(0.0);
            det += m_data[0][0] * Det3(1, 2, 3, 1, 2, 3);
            det -= m_data[1][0] * Det3(0, 2, 3, 1, 2, 3);
            det += m_data[2][0] * Det3(0, 1, 3, 1, 2, 3);
            det -= m_data[3][0] * Det3(0, 1, 2, 1, 2, 3);
            return det;
        }

        
        Type Norm1() const
        {
            Type s1 = abs(m_data[0][0]) + abs(m_data[1][0]) + abs(m_data[2][0]) + abs(m_data[3][0]);
            Type s2 = abs(m_data[0][1]) + abs(m_data[1][1]) + abs(m_data[2][1]) + abs(m_data[3][1]);
            Type s3 = abs(m_data[0][2]) + abs(m_data[1][2]) + abs(m_data[2][2]) + abs(m_data[3][2]);
            Type s4 = abs(m_data[0][3]) + abs(m_data[1][3]) + abs(m_data[2][3]) + abs(m_data[3][3]);
            return (s1 + s2 + s3 + s4);
        }

        bool IsEqual(const CMatrix4<Type> & matrix, double a_tolerance  = g_DoubleResolution) const
        {
            for(int i = 0; i < 4; i++)
            {
                for(int j = 0; j < 4; j++)
                {
                    if (std::abs(m_data[i][j] - matrix.m_data[i][j]) > a_tolerance)
                        return false;
                }
            }
            return true;
        }


        inline const Type& Value(int row, int col) const 
        { 
            return m_data[row][col]; 
        }

        inline void Set(int row, int col,const Type& Val)  
        { 
            m_data[row][col] = Val; 
        }

        void Set(Type m[4][4])
        {
            m_data[0][0] = m[0][0]; m_data[0][1] = m[0][1]; m_data[0][2] = m[0][2]; m_data[0][3] = m[0][3];
            m_data[1][0] = m[1][0]; m_data[1][1] = m[1][1]; m_data[1][2] = m[1][2]; m_data[1][3] = m[1][3];
            m_data[2][0] = m[2][0]; m_data[2][1] = m[2][1]; m_data[2][2] = m[2][2]; m_data[2][3] = m[2][3];
            m_data[3][0] = m[3][0]; m_data[3][1] = m[3][1]; m_data[3][2] = m[3][2]; m_data[3][3] = m[3][3];
        }

        inline const Type& operator()(int row, int col) const 
        { 
            return m_data[row][col];
        }

        inline Type& operator()(int row, int col)  
        { 
            return m_data[row][col];
        }

        inline const Type * Ptr() const
        { 
            return &(m_data[0][0]); 
        }

        inline Type* Ptr()
        {
            return &(m_data[0][0]); 
        }

        Type * operator [](int i)
        { 
            return m_data[i]; 
        }

        inline const Type * operator [](int i) const
        { 
            return m_data[i]; 
        }

        inline CMatrix4<Type> & operator =(const CMatrix4<Type> & m)
        {
            memmove(m_data, m.m_data, sizeof(Type) * 16);
            return *this;
        }

        CMatrix4<Type> & operator *=(const CMatrix4<Type> & matrix)
        {
            return MultRight(matrix);
        }

        friend CMatrix4<Type> operator+(const CMatrix4<Type> & m1, const CMatrix4<Type> & m2)
        {
            CMatrix4<Type> m;
            m[0][0] = m1[0][0] + m2[0][0];
            m[1][0] = m1[1][0] + m2[1][0];
            m[2][0] = m1[2][0] + m2[2][0];
            m[3][0] = m1[3][0] + m2[3][0];
            
            m[0][1] = m1[0][1] + m2[0][1];
            m[1][1] = m1[1][1] + m2[1][1];
            m[2][1] = m1[2][1] + m2[2][1];
            m[3][1] = m1[3][1] + m2[3][1];

            m[0][2] = m1[0][2] + m2[0][2];
            m[1][2] = m1[1][2] + m2[1][2];
            m[2][2] = m1[2][2] + m2[2][2];
            m[3][2] = m1[3][2] + m2[3][2];

            m[0][3] = m1[0][3] + m2[0][3];
            m[1][3] = m1[1][3] + m2[1][3];
            m[2][3] = m1[2][3] + m2[2][3];
            m[3][3] = m1[3][3] + m2[3][3];

            return m;
        }

        friend CMatrix4<Type> operator-(const CMatrix4<Type> & m1, const CMatrix4<Type> & m2)
        {
            CMatrix4<Type> m;
            m[0][0] = m1[0][0] - m2[0][0];
            m[1][0] = m1[1][0] - m2[1][0];
            m[2][0] = m1[2][0] - m2[2][0];
            m[3][0] = m1[3][0] - m2[3][0];

            m[0][1] = m1[0][1] - m2[0][1];
            m[1][1] = m1[1][1] - m2[1][1];
            m[2][1] = m1[2][1] - m2[2][1];
            m[3][1] = m1[3][1] - m2[3][1];

            m[0][2] = m1[0][2] - m2[0][2];
            m[1][2] = m1[1][2] - m2[1][2];
            m[2][2] = m1[2][2] - m2[2][2];
            m[3][2] = m1[3][2] - m2[3][2];

            m[0][3] = m1[0][3] - m2[0][3];
            m[1][3] = m1[1][3] - m2[1][3];
            m[2][3] = m1[2][3] - m2[2][3];
            m[3][3] = m1[3][3] - m2[3][3];

            return m;
        }

        friend CMatrix4<Type> operator *(const CMatrix4<Type> & m1, const CMatrix4<Type> & m2)
        { 
            CMatrix4<Type> m = m1; 
            m *= m2; 
            return m; 
        }

        friend CVec4<Type> operator *(const CMatrix4<Type> & m, const CVec4<Type> & v)
        {
            CVec4<Type> vRet;

            vRet[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2] + m[0][3] * v[3];
            vRet[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2] + m[1][3] * v[3];
            vRet[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2] + m[2][3] * v[3];
            vRet[3] = m[3][0] * v[0] + m[3][1] * v[1] + m[3][2] * v[2] + m[3][3] * v[3];

            return vRet;
        }

        friend CVec4<Type> operator *(const CVec4<Type> & v, const CMatrix4<Type> & m)
        {
            CVec4<Type> vRet;

            vRet[0] = v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0] + v[3] * m[3][0];
            vRet[1] = v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1] + v[3] * m[3][1];
            vRet[2] = v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2] + v[3] * m[3][2];
            vRet[3] = v[0] * m[0][3] + v[1] * m[1][3] + v[2] * m[2][3] + v[3] * m[3][3];

            return vRet;
        }

        friend bool operator ==(const CMatrix4<Type> & m1, const CMatrix4<Type> & m2)
        {
            return (
                m1.m_data[0][0] == m2.m_data[0][0] &&
                m1.m_data[0][1] == m2.m_data[0][1] &&
                m1.m_data[0][2] == m2.m_data[0][2] &&
                m1.m_data[0][3] == m2.m_data[0][3] &&

                m1.m_data[1][0] == m2.m_data[1][0] &&
                m1.m_data[1][1] == m2.m_data[1][1] &&
                m1.m_data[1][2] == m2.m_data[1][2] &&
                m1.m_data[1][3] == m2.m_data[1][3] &&

                m1.m_data[2][0] == m2.m_data[2][0] &&
                m1.m_data[2][1] == m2.m_data[2][1] &&
                m1.m_data[2][2] == m2.m_data[2][2] &&
                m1.m_data[2][3] == m2.m_data[2][3] &&

                m1.m_data[3][0] == m2.m_data[3][0] &&
                m1.m_data[3][1] == m2.m_data[3][1] &&
                m1.m_data[3][2] == m2.m_data[3][2] &&
                m1.m_data[3][3] == m2.m_data[3][3]
            );
        }

        friend bool operator !=(const CMatrix4<Type> & m1, const CMatrix4<Type> & m2)
        { 
            return !(m1 == m2);
        }


        static inline CMatrix4<Type> Translate(const CVec3<Type> & t)
        {
            CMatrix4<Type> m;
            m.SetTranslate(t);

            return m;
        }

        static inline CMatrix4<Type> Translate(const Type& tx,const Type& ty,const Type& tz)
        {
            CMatrix4<Type> m;
            m.SetTranslate(tx,ty,tz);

            return m;
        }

        static inline CMatrix4<Type> Rotate(const CVec3<Type>& Pos,const CVec3<Type>& Dir, const Type& Angle)
        {
            CMatrix4<Type> m;
            m.SetRotate(Pos,Dir,Angle);

            return m;
        }

        static inline CMatrix4<Type> Scale(const CVec3<Type>& Pos, const Type& Scale)
        {
            CMatrix4<Type> m;
            m.SetScale(Pos,Scale);

            return m;
        }

        static inline CMatrix4<Type> Scale(const CVec3<Type>& Pos, const CVec3<Type>& Scales)
        {
            CMatrix4<Type> m;
            m.SetScale(Pos,Scales);

            return m;
        }

        static inline CMatrix4<Type> Mirror(const CVec3<Type>& Pos, const CVec3<Type>& Normal)
        {
            CMatrix4<Type> m;
            m.SetMirror(Pos,Normal);

            return m;
        }


        void SetRotateAndTranslate(const CMatrix3<Type> & R, const CVec3<Type> & t)
        {
            m_data[0][0] = R[0][0]; m_data[0][1] = R[0][1]; m_data[0][2] = R[0][2];
            m_data[1][0] = R[1][0]; m_data[1][1] = R[1][1]; m_data[1][2] = R[1][2];
            m_data[2][0] = R[2][0]; m_data[2][1] = R[2][1]; m_data[2][2] = R[2][2];

            m_data[3][0] = t.X; m_data[3][1] = t.Y; m_data[3][2] = t.Z;
        }

        void SetTranslate(const CVec3<Type> & t)
        {
            MakeIdentity();
            m_data[3][0] = t[0];
            m_data[3][1] = t[1];
            m_data[3][2] = t[2];
        }

        void AddTranslate(const CVec3<Type> & t)
        {
            m_data[3][0] += t[0];
            m_data[3][1] += t[1];
            m_data[3][2] += t[2];
        }

        CMatrix4<Type>& SetToTranslation(const CVec3<Type> & t)
        {
            MakeIdentity();
            m_data[3][0] = t[0];
            m_data[3][1] = t[1];
            m_data[3][2] = t[2];
            return *this;
        }

        inline void SetTranslate(const Type& tx,const Type& ty,const Type& tz)
        {
            MakeIdentity();
            m_data[3][0] = tx;
            m_data[3][1] = ty;
            m_data[3][2] = tz;
        }

        void SetRotate(const CVec3<Type>& Pos,const CVec3<Type>& Dir, const Type& Angle)
        {
            Type a,b,c;
            Dir.Value(a,b,c);

            Type Cos = cos(Angle);
            Type Sin = sin(Angle);

            if (IsNearZero( b * b + c * c))
            {
                if (a < 0)
                {
                    Cos = cos((Type)(M_2PI - Angle));
                    Sin = sin((Type)(M_2PI - Angle));
                }

                Type dy = Pos.Y * ((Type)1.0 - Cos) + Pos.Z * Sin;
                Type dz = Pos.Z * ((Type)1.0 - Cos) - Pos.Y * Sin; 

                m_data[0][0] = 1.0;        m_data[0][1] = 0.0;        m_data[0][2] = 0.0;       m_data[0][3] = 0.0;    
                m_data[1][0] = 0.0;        m_data[1][1] = Cos;       m_data[1][2] = Sin;      m_data[1][3] = 0.0;    
                m_data[2][0] = 0.0;        m_data[2][1] = -Sin;      m_data[2][2] = Cos;      m_data[2][3] = 0.0;    
                m_data[3][0] = 0.0;        m_data[3][1] = dy;         m_data[3][2] = dz;       m_data[3][3] = 1.0;    

                return;
            }

            if ( IsNearZero(a * a + c * c) )//绕Y轴旋转
            {
                if (b < 0)
                {
                    Cos = cos( ((Type)M_2PI) - Angle);
                    Sin = sin( ((Type)M_2PI) - Angle);
                }

                Type dx = Pos.X * ((Type)1.0-Cos) - Pos.Z * Sin;
                Type dz = Pos.Z * ((Type)1.0-Cos) + Pos.X * Sin; 

                m_data[0][0] = Cos;      m_data[0][1] = 0.0;     m_data[0][2] = -Sin;    m_data[0][3] = 0.0;    
                m_data[1][0] = 0.0;        m_data[1][1] = 1.0;     m_data[1][2] = 0.0;      m_data[1][3] = 0.0;    
                m_data[2][0] = Sin;       m_data[2][1] = 0.0;      m_data[2][2] = Cos;     m_data[2][3] = 0.0;    
                m_data[3][0] = dx;        m_data[3][1] = 0.0;       m_data[3][2] = dz;        m_data[3][3] = 1.0;    

                return;
            }

            if ( IsNearZero(a * a + b * b) )//绕Z轴旋转
            {
                if (c < 0)
                {
                    Cos = cos( ((Type)M_2PI) - Angle);
                    Sin = sin( ((Type)M_2PI) - Angle);
                }

                Type dx = Pos.X * ((Type)1.0 - Cos) + Pos.Y * Sin;
                Type dy = Pos.Y * ((Type)1.0 - Cos) - Pos.X * Sin; 

                m_data[0][0] = Cos;      m_data[0][1] = Sin;   m_data[0][2] = 0.0;   m_data[0][3] = 0.0; 
                m_data[1][0] = -Sin;     m_data[1][1] = Cos;   m_data[1][2] = 0.0;   m_data[1][3] = 0.0;    
                m_data[2][0] = 0.0;      m_data[2][1] = 0.0;    m_data[2][2] = 1.0;   m_data[2][3] = 0.0;    
                m_data[3][0] = dx;       m_data[3][1] = dy;     m_data[3][2] = 0.0;   m_data[3][3] = 1.0;    

                return;
            }

            Type v = sqrt( b * b + c * c);
            Type s = sqrt(a * a + b * b + c * c);



            CMatrix4<Type> MatTa;
            MatTa.Set(3,0,-Pos.X);
            MatTa.Set(3,1,-Pos.Y);
            MatTa.Set(3,2,-Pos.Z);
            MatTa.Set(3,3,(Type)1.0);

            CMatrix4<Type> InverseMatTa;
            InverseMatTa.Set(3,0,Pos.X);
            InverseMatTa.Set(3,1,Pos.Y);
            InverseMatTa.Set(3,2,Pos.Z);
            InverseMatTa.Set(3,3,(Type)1.0);

            Type Cos1 = c/v;
            Type Sin1 = b/v;

            CMatrix4<Type> MatRx;
            MatRx.Set(1,1,Cos1);
            MatRx.Set(1,2,Sin1);
            MatRx.Set(2,1,-Sin1);
            MatRx.Set(2,2,Cos1);

            CMatrix4<Type> InverseMatRx ;

            InverseMatRx.Set(1,1,Cos1);
            InverseMatRx.Set(1,2,-Sin1);
            InverseMatRx.Set(2,1,Sin1);
            InverseMatRx.Set(2,2,Cos1);

            Type Cos2 = v/s;
            Type Sin2 = a/s;

            CMatrix4<Type> MatRy;
            MatRy.Set(0,0,Cos2);
            MatRy.Set(0,2,Sin2);
            MatRy.Set(2,0,-Sin2);
            MatRy.Set(2,2,Cos2);

            CMatrix4<Type> InverseMatRy ;

            InverseMatRy.Set(0,0,Cos2);
            InverseMatRy.Set(0,2,-Sin2);
            InverseMatRy.Set(2,0,Sin2);
            InverseMatRy.Set(2,2,Cos2);



            CMatrix4<Type> MatRz;
            MatRz.Set(0,0,Cos);
            MatRz.Set(0,1,Sin);
            MatRz.Set(1,0,-Sin);
            MatRz.Set(1,1,Cos);

            SetProduct(MatTa,MatRx);
            MultRight(MatRy);
            MultRight(MatRz);
            MultRight(InverseMatRy);
            MultRight(InverseMatRx);
            MultRight(InverseMatTa);

        }

        void SetScale(const CVec3<Type>& Pos, const Type& Scale)
        {
            assert(Scale > 0.0);
//             if(Scale < 0.0)
//                 Scale *= -1;

            MakeIdentity();

            m_data[0][0] = Scale;                         
            m_data[1][1] = Scale;                         
            m_data[2][2] = Scale;                 
            m_data[3][0] = Pos.X * ((Type)1.0 - Scale);   
            m_data[3][1] = Pos.Y * ((Type)1.0 - Scale);  
            m_data[3][2] = Pos.Z * ((Type)1.0 - Scale);      

        }

        void SetScale(const CVec3<Type>& Pos, const CVec3<Type>& Scales)
        {
            assert(Scales.X != 0.0 && Scales.Y!=0.0 && Scales.Z!=0);

            MakeIdentity();

            m_data[0][0] = Scales.X;
            m_data[1][1] = Scales.Y;                         
            m_data[2][2] = Scales.Z;                 
            m_data[3][0] = Pos.X * ((Type)1.0 - Scales.X);   
            m_data[3][1] = Pos.Y * ((Type)1.0 - Scales.Y);  
            m_data[3][2] = Pos.Z * ((Type)1.0 - Scales.Z);      
        }

        void SetMirror(const CVec3<Type>& Pos, const CVec3<Type>& Normal)
        {
            CVec3<Type> v1;
            v1.Set((Type)0.0, (Type)0.0, (Type)1.0);
            CVec3<Type> DirU = Normal ^ v1;
            Type Length = DirU.Length();
            if ( IsNearZero(Length))
            {
                DirU.Set((Type)1.0, (Type)0.0, (Type)0.0);   
            }
            else
            {
                DirU /= Length;     
            }

            CVec3<Type> PlnNormal = Normal;
            PlnNormal.Normalize();

            CVec3<Type> CrossDir = PlnNormal ^ DirU;

            Type cx,cy,cz,nx,ny,nz,ux,uy,uz,px,py,pz;
            DirU.Value(ux,uy,uz);
            CrossDir.Value(cx,cy,cz);
            PlnNormal.Value(nx,ny,nz);
            Pos.Value(px,py,pz);

            CMatrix4<Type> MatTmp(ux,uy,uz,(Type)0.0,cx,cy,cz,(Type)0.0,nx,ny,nz,(Type)0.0,px,py,pz,(Type)1.0);

            CMatrix4<Type> MatTmpInverse = MatTmp.Inverse();

            CMatrix4<Type> MirrorXMatTmp;
            MirrorXMatTmp.Set(2,2,(Type)(-1.0));

            SetProduct(MatTmpInverse,MirrorXMatTmp);
            MultRight(MatTmp);

        }


        CMatrix4<Type>&  SetProduct(const CMatrix4<Type> & m1,const CMatrix4<Type> & m2)
        {
            (*this) = m1;
            return MultRight(m2);
        }

        CMatrix4<Type> MultRight(const CMatrix3<Type> & R, const CVec3<Type> & t) const
        {
            CMatrix4<Type> m;
            m[0][0] = m_data[0][0] * R[0][0] + m_data[0][1] * R[1][0] + m_data[0][2] * R[2][0];
            m[0][1] = m_data[0][0] * R[0][1] + m_data[0][1] * R[1][1] + m_data[0][2] * R[2][1];
            m[0][2] = m_data[0][0] * R[0][2] + m_data[0][1] * R[1][2] + m_data[0][2] * R[2][2];

            m[1][0] = m_data[1][0] * R[0][0] + m_data[1][1] * R[1][0] + m_data[1][2] * R[2][0];
            m[1][1] = m_data[1][0] * R[0][1] + m_data[1][1] * R[1][1] + m_data[1][2] * R[2][1];
            m[1][2] = m_data[1][0] * R[0][2] + m_data[1][1] * R[1][2] + m_data[1][2] * R[2][2];

            m[2][0] = m_data[2][0] * R[0][0] + m_data[2][1] * R[1][0] + m_data[2][2] * R[2][0];
            m[2][1] = m_data[2][0] * R[0][1] + m_data[2][1] * R[1][1] + m_data[2][2] * R[2][1];
            m[2][2] = m_data[2][0] * R[0][2] + m_data[2][1] * R[1][2] + m_data[2][2] * R[2][2];

            CVec3<Type> & v = *(CVec3<Type> *)m_data[3];

            CVec3<Type> & mv = *(CVec3<Type> *)m[3];
            mv = R.MultiLeft(v) + t;
            return m;
        }

        CMatrix4<Type> & MultRight(const CMatrix4<Type> & m)
        {
            // Trivial cases
            if(m.IsIdentity()) 
                return *this;
            else if(IsIdentity())
                return (*this = m);

            Type tmp[4][4];


            tmp[0][0] = m_data[0][0]*m.m_data[0][0] + m_data[0][1]*m.m_data[1][0] + m_data[0][2]*m.m_data[2][0] + m_data[0][3]*m.m_data[3][0];
            tmp[0][1] = m_data[0][0]*m.m_data[0][1] + m_data[0][1]*m.m_data[1][1] + m_data[0][2]*m.m_data[2][1] + m_data[0][3]*m.m_data[3][1];
            tmp[0][2] = m_data[0][0]*m.m_data[0][2] + m_data[0][1]*m.m_data[1][2] + m_data[0][2]*m.m_data[2][2] + m_data[0][3]*m.m_data[3][2];
            tmp[0][3] = m_data[0][0]*m.m_data[0][3] + m_data[0][1]*m.m_data[1][3] + m_data[0][2]*m.m_data[2][3] + m_data[0][3]*m.m_data[3][3];
            tmp[1][0] = m_data[1][0]*m.m_data[0][0] + m_data[1][1]*m.m_data[1][0] + m_data[1][2]*m.m_data[2][0] + m_data[1][3]*m.m_data[3][0];
            tmp[1][1] = m_data[1][0]*m.m_data[0][1] + m_data[1][1]*m.m_data[1][1] + m_data[1][2]*m.m_data[2][1] + m_data[1][3]*m.m_data[3][1];
            tmp[1][2] = m_data[1][0]*m.m_data[0][2] + m_data[1][1]*m.m_data[1][2] + m_data[1][2]*m.m_data[2][2] + m_data[1][3]*m.m_data[3][2];
            tmp[1][3] = m_data[1][0]*m.m_data[0][3] + m_data[1][1]*m.m_data[1][3] + m_data[1][2]*m.m_data[2][3] + m_data[1][3]*m.m_data[3][3];
            tmp[2][0] = m_data[2][0]*m.m_data[0][0] + m_data[2][1]*m.m_data[1][0] + m_data[2][2]*m.m_data[2][0] + m_data[2][3]*m.m_data[3][0];
            tmp[2][1] = m_data[2][0]*m.m_data[0][1] + m_data[2][1]*m.m_data[1][1] + m_data[2][2]*m.m_data[2][1] + m_data[2][3]*m.m_data[3][1];
            tmp[2][2] = m_data[2][0]*m.m_data[0][2] + m_data[2][1]*m.m_data[1][2] + m_data[2][2]*m.m_data[2][2] + m_data[2][3]*m.m_data[3][2];
            tmp[2][3] = m_data[2][0]*m.m_data[0][3] + m_data[2][1]*m.m_data[1][3] + m_data[2][2]*m.m_data[2][3] + m_data[2][3]*m.m_data[3][3];
            tmp[3][0] = m_data[3][0]*m.m_data[0][0] + m_data[3][1]*m.m_data[1][0] + m_data[3][2]*m.m_data[2][0] + m_data[3][3]*m.m_data[3][0];
            tmp[3][1] = m_data[3][0]*m.m_data[0][1] + m_data[3][1]*m.m_data[1][1] + m_data[3][2]*m.m_data[2][1] + m_data[3][3]*m.m_data[3][1];
            tmp[3][2] = m_data[3][0]*m.m_data[0][2] + m_data[3][1]*m.m_data[1][2] + m_data[3][2]*m.m_data[2][2] + m_data[3][3]*m.m_data[3][2];
            tmp[3][3] = m_data[3][0]*m.m_data[0][3] + m_data[3][1]*m.m_data[1][3] + m_data[3][2]*m.m_data[2][3] + m_data[3][3]*m.m_data[3][3];

            Set(tmp);

            return *this;
        }

        CVec3<Type> MultiPointLeft(const CVec3<Type> & src) const
        {
            Type x = src[0]*m_data[0][0] + src[1]*m_data[1][0] + src[2]*m_data[2][0] + m_data[3][0];
            Type y = src[0]*m_data[0][1] + src[1]*m_data[1][1] + src[2]*m_data[2][1] + m_data[3][1];
            Type z = src[0]*m_data[0][2] + src[1]*m_data[1][2] + src[2]*m_data[2][2] + m_data[3][2];
            Type w = src[0]*m_data[0][3] + src[1]*m_data[1][3] + src[2]*m_data[2][3] + m_data[3][3];

            CVec3<Type> dst;
            auto w1 = (Type(1))/w;
            return CVec3<Type>(x*w1, y*w1, z*w1);
        }

        CVec3<Type> MultiVecLeft(const CVec3<Type> & src) const
        {
            Type x = src[0]*m_data[0][0] + src[1]*m_data[1][0] + src[2]*m_data[2][0];
            Type y = src[0]*m_data[0][1] + src[1]*m_data[1][1] + src[2]*m_data[2][1];
            Type z = src[0]*m_data[0][2] + src[1]*m_data[1][2] + src[2]*m_data[2][2];

            CVec3<Type> dst;
            dst.Set(x, y, z);

            return dst;
        }

        CVec4<Type> MultiVecLeft(const CVec4<Type> & src) const
        {
            Type x = src[0]*m_data[0][0] + src[1]*m_data[1][0] + src[2]*m_data[2][0] + src[3]*m_data[3][0];
            Type y = src[0]*m_data[0][1] + src[1]*m_data[1][1] + src[2]*m_data[2][1] + src[3]*m_data[3][1];
            Type z = src[0]*m_data[0][2] + src[1]*m_data[1][2] + src[2]*m_data[2][2] + src[3]*m_data[3][2];
            Type w = src[0]*m_data[0][3] + src[1]*m_data[1][3] + src[2]*m_data[2][3] + src[3]*m_data[3][3];

            CVec4<Type> dst;
            dst.Set(x, y, z,w);

            return dst;
        }


        bool IsOrthoScale(double & dScale, double dTol = g_DoubleResolution) const
        {
            CVec3<Type> a3(m_data[0][3], m_data[1][3], m_data[2][3]);
            if (a3.IsZero(dTol))
            {
                CVec3<Type> r0;
                CVec3<Type> r1; 
                CVec3<Type> r2;
                GetRotationAxes(r0, r1, r2);
                double dTol2 = dTol * dTol;
                if (isZero(r0.Dot(r1), dTol2) &&
                    isZero(r1.Dot(r2), dTol2) &&
                    isZero(r2.Dot(r1), dTol2))
                {
                    double r0sqr = r0.Dot(r0);
                    double r1sqr = r1.Dot(r1);
                    double r2sqr = r2.Dot(r2);

                    if (sameValue(r0sqr, r1sqr, dTol2) &&
                        sameValue(r1sqr, r2sqr, dTol2) &&
                        sameValue(r2sqr, r0sqr, dTol2))
                    {
                        dScale = sqrt_safe(r0sqr);
                        return true;
                    }
                }
                
            }

            return false;
        }

        bool IsXOYOrthoScale(double & dScale, double dTol = g_DoubleResolution) const
        {
            CVec2<Type> a3(m_data[0][3], m_data[1][3]);
            if (a3.IsZero(dTol))
            {
                CVec2<Type> r0(m_data[0][0], m_data[0][1]);
                CVec2<Type> r1(m_data[1][0], m_data[1][1]); 
                if (isZero(r0.Dot(r1), dTol))
                {
                    //double r0sqr = r0.Dot(r0);
                    //double r1sqr = r1.Dot(r1);
                    double r0sqr = r0.Length();
                    double r1sqr = r1.Length();
                    if (sameValue(r0sqr, r1sqr, dTol))
                    {
                        dScale = r0sqr;
                        return true;
                    }
                }
                
            }

            return false;
        }

        bool IsXYZOrthoScale(double & dScale, double dTol = g_RelaxedDoubleResolution) const
        {
            CVec3<Type> r0;
            CVec3<Type> r1;
            CVec3<Type> r2;
            GetRotationAxes(r0, r1, r2);
            if (isZero(r0.Dot(r1), dTol) &&
                isZero(r1.Dot(r2), dTol) &&
                isZero(r0.Dot(r2), dTol))
            {
                dScale = r0.Length();
                if (sameValue(r1.Length(), dScale, dTol) &&
                    sameValue(r2.Length(), dScale, dTol) &&
                    !isZero(dScale, dTol))
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }

        friend std::ostream & operator<<(std::ostream & os, const CMatrix4<Type> & mat)
        { 
            for(unsigned int i=0; i<4; i++){
                os << mat[i][0] <<'\t'<< mat[i][1] <<'\t'<< mat[i][2] <<'\t'<< mat[i][3] <<std::endl;
            }
            return os;
        }

        friend std::istream & operator>>(std::istream & iStream, CMatrix4<Type> & mat)
        {
            for(unsigned int i=0; i<4; i++){
                iStream >> mat[i][0] >> mat[i][1] >> mat[i][2] >> mat[i][3];
            }
            return iStream;
        }

        char* AsString()
        {
            char * a = AsStringHelper::NewString("Matrix4=[(%.20le %.20le %.20le %.20le),(%.20le %.20le %.20le %.20le),(%.20le %.20le %.20le %.20le),(%.20le %.20le %.20le %.20le)]", 
                m_data[0][0], m_data[0][1], m_data[0][2], m_data[0][3], 
                m_data[1][0], m_data[1][1], m_data[1][2], m_data[1][3], 
                m_data[2][0], m_data[2][1], m_data[2][2], m_data[2][3],
                m_data[3][0], m_data[3][1], m_data[3][2], m_data[3][3]);
            return a;
        }

        static CMatrix4<double>* LoadFromString(char* str)
        {
            CMatrix4<double>* pMatrix = new CMatrix4<double>;
            CMatrix4<double>& matrix = *pMatrix;
            g_sscanf(str, "Matrix4=[(%le %le %le %le),(%le %le %le %le),(%le %le %le %le),(%le %le %le %le)]", 
                &matrix(0, 0), &matrix(0, 1), &matrix(0, 2), &matrix(0, 3), 
                &matrix(1, 0), &matrix(1, 1), &matrix(1, 2), &matrix(1, 3), 
                &matrix(2, 0), &matrix(2, 1), &matrix(2, 2), &matrix(2, 3), 
                &matrix(3, 0), &matrix(3, 1), &matrix(3, 2), &matrix(3, 3));

            return pMatrix;
        }

        void GetRotationAxes(CVec3<Type>& xAxis, CVec3<Type>& yAxis, CVec3<Type>& zAxis) const
        {
            xAxis = CVec3<Type>(m_data[0][0], m_data[0][1], m_data[0][2]);
            yAxis = CVec3<Type>(m_data[1][0], m_data[1][1], m_data[1][2]);
            zAxis = CVec3<Type>(m_data[2][0], m_data[2][1], m_data[2][2]);
        }

        CVec3<Type> GetTranslation() const
        {
            return CVec3<Type>(m_data[3][0], m_data[3][1], m_data[3][2]);    
        }

        void Decompose(CVec3<Type>& vTrans, CVec3<Type>& vScale, CMatrix3<Type>& mRot) const
        {
            vTrans.X = m_data[3][0];
            vTrans.Y = m_data[3][1];
            vTrans.Z = m_data[3][2];

            CVec3<Type> vCols[3] = {
                CVec3<Type>(m_data[0][0],m_data[0][1],m_data[0][2]),
                CVec3<Type>(m_data[1][0],m_data[1][1],m_data[1][2]),
                CVec3<Type>(m_data[2][0],m_data[2][1],m_data[2][2])
            };

            vScale.X = vCols[0].Length();
            vScale.Y = vCols[1].Length();
            vScale.Z = vCols[2].Length();

            if(vScale.X != 0)
            {
                vCols[0].X /= vScale.X;
                vCols[0].Y /= vScale.X;
                vCols[0].Z /= vScale.X;
            }
            if(vScale.Y != 0)
            {
                vCols[1].X /= vScale.Y;
                vCols[1].Y /= vScale.Y;
                vCols[1].Z /= vScale.Y;
            }
            if(vScale.Z != 0)
            {
                vCols[2].X /= vScale.Z;
                vCols[2].Y /= vScale.Z;
                vCols[2].Z /= vScale.Z;
            }

            mRot.Set(0,0, vCols[0].X);
            mRot.Set(1,0, vCols[0].Y);
            mRot.Set(2,0, vCols[0].Z);
            //mRot.Set(0,3,0);
            //mRot.Set(3,0,0);
            mRot.Set(0,1, vCols[1].X);
            mRot.Set(1,1, vCols[1].Y);
            mRot.Set(2,1, vCols[1].Z);
            //mRot.Set(1,3,0);
            //mRot.Set(3,1,0);
            mRot.Set(0,2, vCols[2].X);
            mRot.Set(1,2, vCols[2].Y);
            mRot.Set(2,2, vCols[2].Z);
            //mRot.Set(2,3,0);
            //mRot.Set(3,2,0);  
            //mRot.Set(3,3,1);
        }

        void Compose(const CVec3<Type>& vTrans, const CVec3<Type>& vScale, const CMatrix3<Type>& mRot)
        {
            m_data[3][0] = vTrans.X;
            m_data[3][1] = vTrans.Y;
            m_data[3][2] = vTrans.Z;

            m_data[0][0] = mRot.Value(0,0) * vScale.X;
            m_data[0][1] = mRot.Value(1,0) * vScale.X;
            m_data[0][2] = mRot.Value(2,0) * vScale.X;

            m_data[1][0] = mRot.Value(0,1) * vScale.Y;
            m_data[1][1] = mRot.Value(1,1) * vScale.Y;
            m_data[1][2] = mRot.Value(2,1) * vScale.Y;

            m_data[2][0] = mRot.Value(0,2) * vScale.Z;
            m_data[2][1] = mRot.Value(1,2) * vScale.Z;
            m_data[2][2] = mRot.Value(2,2) * vScale.Z;
        }

        void GetBlock3(CMatrix3<Type> & A, CVec3<Type> & b, CVec3<Type> & c, Type & d) const
        {
            A[0][0] = m_data[0][0];
            A[0][1] = m_data[0][1];
            A[0][2] = m_data[0][2];

            A[1][0] = m_data[1][0];
            A[1][1] = m_data[1][1];
            A[1][2] = m_data[1][2];

            A[2][0] = m_data[2][0];
            A[2][1] = m_data[2][1];
            A[2][2] = m_data[2][2];

            b.Set(m_data[0][3], m_data[1][3], m_data[2][3]);
            c.Set(m_data[3][0], m_data[3][1], m_data[3][2]);

            d = m_data[3][3];
        }

        bool IsValidTransform() const
        {
            CVec3<Type> translation;
            CVec3<Type> scale;
            CMatrix3<Type> rotation;

            Decompose(translation, scale, rotation);

            // check rotation axis.
            bool bValidRotation = true;
            for (int i=0; i<3; ++i)
            {
                CVec3<Type> axis(rotation[i][0], rotation[i][1], rotation[i][2]);
                // invalid vector or non uniform axis length.
                if (!axis.IsValid()/* || fabs(axis.Length() - 1.0f) > FLT_EPSILON*/)
                {
                    bValidRotation = false;
                    break;
                }
            }

            // non dim of scale is zero.
            bool bValidScale = scale.IsValid() && fabs(scale[scale.MinDimension()]) > FLT_EPSILON;

            return translation.IsValid() && bValidRotation && bValidScale; 
        }

        static CMatrix4<Type> MakeLookAtMatrix(const CVec3<Type>& eye, const CVec3<Type>& center,
            const CVec3<Type>& up)
        {
            CVec3<Type> f(center - eye);
            f.Normalize();
            CVec3<Type> s(f ^ up);
            s.Normalize();
            CVec3<Type> u(s ^ f);
            u.Normalize();

            CMatrix4<Type> mat(
                s[0], u[0], -f[0], 0.0,
                s[1], u[1], -f[1], 0.0,
                s[2], u[2], -f[2], 0.0,
                0.0,  0.0,  0.0,   1.0 );

            CVec3<Type> pos = mat.MultiPointLeft(-eye);
            mat.Set(3, 0, pos.X);
            mat.Set(3, 1, pos.Y);
            mat.Set(3, 2, pos.Z);

            return mat;
        }

        static CMatrix4<Type> MakePerspectiveMatrix(Type fovy, Type aspectRatio, Type zNear, Type zFar)
        {
            Type t = 1.0f / tan(DegToRad(fovy * 0.5f));
            Type A = zNear - zFar;
            Type B = zNear + zFar;
            Type C = 2.0f * zFar * zNear;
            Type D = t / aspectRatio;

            CMatrix4<Type> mat(
                D,   0.0, 0.0,    0.0,
                0.0, t,   0.0,    0.0,
                0.0, 0.0, B / A, -1.0,
                0.0, 0.0, C / A,  0.0
                );

            return mat;
        }

        static CMatrix4<Type> MakeFrustumMatrix(Type left, Type right, Type bottom, Type top,
            Type zNear, Type zFar)
        {
            Type A =  (right + left) / (right - left);
            Type B =  (top + bottom) / (top - bottom);
            Type C = -(zFar + zNear) / (zFar - zNear);
            Type D = -2.0f * zFar * zNear / (zFar - zNear);
            Type E =  2.0f * zNear / (right - left);
            Type F =  2.0f * zNear / (top - bottom);

            CMatrix4<Type> mat(
                E,    0.0f, 0.0f,  0.0f,
                0.0f, F,    0.0f,  0.0f,
                A,    B,    C,    -1.0f,
                0.0f, 0.0f, D,     0.0f );

            return mat;
        }

        static CMatrix4<Type> MakeOrthoMatrix(Type left, Type right, Type bottom, Type top, 
            Type zNear, Type zFar)
        {
            Type tx = -(right + left) / (right - left);
            Type ty = -(top + bottom) / (top - bottom);
            Type tz = -(zFar + zNear) / (zFar - zNear);
            Type A =  2.0f / (right - left);
            Type B =  2.0f / (top - bottom);
            Type C = -2.0f / (zFar - zNear);

            CMatrix4<Type> mat(
                A,   0.0, 0.0, 0.0,
                0.0, B,   0.0, 0.0,
                0.0, 0.0, C,   0.0,
                tx,  ty,  tz,  1.0 );

            return mat;
        }

        static CMatrix4<Type> MakeOrtho2DMatrix(Type left, Type right, Type bottom, Type top)
        {
            return MakeOrthoMatrix(left, right, bottom, top, (Type)-1.0, (Type)1.0);
        }

    private:
        Type m_data[4][4];

    public:
        static const CMatrix4<Type> Zero;
        static const CMatrix4<Type> Identity;
    };

    template<class T> const CMatrix4<T> CMatrix4<T>::Zero(CMatrix4<T>::CM_ZERO);
    template<class T> const CMatrix4<T> CMatrix4<T>::Identity(CMatrix4<T>::CM_IDENTITY);

    typedef CMatrix4<float>  CMatrix4f;
    typedef CMatrix4<double> CMatrix4d;

    /*! @} */

#endif
