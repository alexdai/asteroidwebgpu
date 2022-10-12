/*!
* @file      Matrix.h
* @brief     Matrix
* @details
* @author
* @date      2021.7.30
* @copyright Copyright 2012-2022 GLODON
**********************************************************************************
*/

#ifndef Matrix_h__
#define Matrix_h__

//#include <initializer_list>
#include "GMath/GVec2.h"
#include "GMath/GVec3.h"
#include "GMath/GVec4.h"

namespace ggp
{
    /*!
    * @class Matrix
    * @brief 
    */
    template<typename Type, int Rows, int Cols>
    class Matrix
    {
    public:
        Matrix() {};
        Matrix(const CVec2<Type> & v)
        {
            if (Rows == 1 && Cols == 2)
            {
                SetRow(0, v);
            }
            else if(Rows == 2 && Cols == 1)
            {
                SetColumn(0, v);
            }
            else
            {
                assert(false);
            }
        }

        Matrix(const CVec3<Type> & v)
        {
            if (Rows == 1 && Cols == 3)
            {
                SetRow(0, v);
            }
            else if (Rows == 3 && Cols == 1)
            {
                SetColumn(0, v);
            }
            else
            {
                assert(false);
            }
        }

        Matrix(const CVec4<Type> & v)
        {
            if (Rows == 1 && Cols == 4)
            {
                SetRow(0, v);
            }
            else if (Rows == 4 && Cols == 1)
            {
                SetColumn(0, v);
            }
            else
            {
                assert(false);
            }
        }

        //Matrix(std::initializer_list<std::initializer_list<Type> > list)
        //{
        //    int r = 0;
        //    for (auto irow = list.begin(); irow != list.end(); ++irow, ++r)
        //    {
        //        auto & row = *irow;
        //        int c = 0;
        //        for (auto icol = row.begin(); icol != row.end(); ++icol, ++c)
        //        {
        //            data_[r][c] = *icol;
        //        }

        //    }

        //}

        Matrix(Type M[Rows][Cols])
        {
            for (int i = 0; i < Rows; ++i)
            {
                for (int j = 0; j < Cols; ++j)
                {
                    data_[i][j] = M[i][j];
                }
            }
        }

        Matrix(const Matrix<Type, Rows, Cols> & M)
        {
            for (int i = 0; i < Rows; ++i)
            {
                for (int j = 0; j < Cols; ++j)
                {
                    data_[i][j] = M.data_[i][j];
                }
            }
        }

        Matrix<Type, Rows, Cols> & SetIdentity()
        {
            assert(Rows == Cols);

            for (int i = 0; i < Rows; ++i)
            {
                for (int j = 0; j < Cols; ++j)
                {
                    data_[i][j] = 0.0;
                }

                data_[i][i] = 1.0;
            }

            return *this;
        }

        Matrix<Type, Rows, Cols> & SetRow(int iRow, const CVec2<Type> & v)
        {
            assert(Cols >= (sizeof(v) / sizeof(v.X)));
            data_[iRow][0] = v.X;
            data_[iRow][1] = v.Y;

            return *this;
        }

        Matrix<Type, Rows, Cols> & SetRow(int iRow, const CVec3<Type> & v)
        {
            assert(Cols >= (sizeof(v) / sizeof(v.X)));
            data_[iRow][0] = v.X;
            data_[iRow][1] = v.Y;
            data_[iRow][2] = v.Z;

            return *this;
        }

        Matrix<Type, Rows, Cols> & SetRow(int iRow, const CVec4<Type> & v)
        {
            assert(Cols >= (sizeof(v) / sizeof(v.X)));
            data_[iRow][0] = v.X;
            data_[iRow][1] = v.Y;
            data_[iRow][2] = v.Z;
            data_[iRow][3] = v.W;

            return *this;
        }

        Matrix<Type, Rows, Cols> & SetColumn(int iCol, const CVec2<Type> & v)
        {
            assert(Rows >= (sizeof(v) / sizeof(v.X)));
            data_[0][iCol] = v.X;
            data_[1][iCol] = v.Y;

            return *this;
        }

        Matrix<Type, Rows, Cols> & SetColumn(int iCol, const CVec3<Type> & v)
        {
            assert(Rows >= (sizeof(v) / sizeof(v.X)));
            data_[0][iCol] = v.X;
            data_[1][iCol] = v.Y;
            data_[2][iCol] = v.Z;

            return *this;
        }

        Matrix<Type, Rows, Cols> & SetColumn(int iCol, const CVec4<Type> & v)
        {
            assert(Rows >= (sizeof(v) / sizeof(v.X)));
            data_[0][iCol] = v.X;
            data_[1][iCol] = v.Y;
            data_[2][iCol] = v.Z;
            data_[3][iCol] = v.W;

            return *this;
        }


        Type * operator [](int i) const
        {
            return (Type*)data_[i];
        }


        Type & operator()(int i, int j)
        {
            return data_[i][j];
        }



        template<int n>
        Matrix<Type, Rows, n> operator*(Matrix<Type, Cols, n> & N) const
        {
            Matrix<Type, Rows, n> L;

            for (int i = 0; i < Rows; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    L(i,j) = 0;

                    for (int k = 0; k < Cols; ++k)
                    {
                        L(i,j) += data_[i][k] * N(k,j);
                    }
                }
            }

            return L;
        }


        template<int n>
        Matrix<Type, Rows, n> MultiplyTranspose(Matrix<Type, n, Cols> & N) // this * N' //矩阵乘另一个矩阵的转置
        {
            Matrix<Type, Rows, n> L;

            for (int i = 0; i < Rows; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    L(i, j) = 0;

                    for (int k = 0; k < Cols; ++k)
                    {
                        L(i, j) += data_[i][k] * N(j, k);
                    }
                }
            }

            return L;
        }

        Matrix<Type, Rows, 1> operator*(const CVec2<Type> & v) const
        {
            assert(sizeof(v) / sizeof(v[0]) == Cols);

            return (*this)*(*(Matrix<Type, 2, 1>*)(&v));
        }


        Matrix<Type, Rows, 1> operator*(const CVec3<Type> & v) const
        {
            assert(sizeof(v) / sizeof(v[0]) == Cols);

            return (*this)*(*(Matrix<Type, 3, 1>*)(&v));
        }

        Matrix<Type, Rows, 1> operator*(const CVec4<Type> & v) const
        {
            assert(sizeof(v) / sizeof(v[0]) == Cols);

            return (*this)*(*(Matrix<Type, 4, 1>*)(&v));
        }


        Matrix<Type, Rows, Cols> & operator+=(const Matrix<Type, Rows, Cols> & M)
        {
            for (int i = 0; i < Rows; ++i)
            {
                for (int j = 0; j < Cols; ++j)
                {
                    data_[i][j] += M.data_[i][j];
                }
            }

            return *this;
        }
        Matrix<Type, Rows, Cols> & operator-=(const Matrix<Type, Rows, Cols> & M)
        {
            for (int i = 0; i < Rows; ++i)
            {
                for (int j = 0; j < Cols; ++j)
                {
                    data_[i][j] -= M.data_[i][j];
                }
            }

            return *this;
        }

        Matrix<Type, Rows, Cols> & operator*=(Type k)
        {
            for (int i = 0; i < Rows; ++i)
            {
                for (int j = 0; j < Cols; ++j)
                {
                    data_[i][j] *= k;
                }
            }

            return *this;
        }

        CVec2<Type> & Vec2() const
        {
            assert((Rows == 2) && (Cols == 1) || (Rows == 1) && (Cols == 2));
            return *(CVec2<Type>*)&data_[0][0];
        }

        CVec3<Type> & Vec3() const
        {
            assert((Rows == 3) && (Cols == 1) || (Rows == 1) && (Cols == 3));
            return *(CVec3<Type>*)&data_[0][0];
        }

        CVec4<Type> & Vec4() const
        {
            assert((Rows == 4) && (Cols == 1) || (Rows == 1) && (Cols == 4));
            return *(CVec4<Type>*)&data_[0][0];
        }



    private:
        Type data_[Rows][Cols];


    };

    typedef Matrix<double, 4, 4> Mat4d;
    typedef Matrix<double, 3, 3> Mat3d;

    typedef Matrix<double, 2, 3> Mat2x3d;
    typedef Matrix<double, 3, 2> Mat3x2d;

    typedef Matrix<double, 2, 2> Mat2d;

    typedef Matrix<double, 1, 3> RowVector3d;
    typedef Matrix<double, 3, 1> ColumnVector3d;

    typedef Matrix<double, 1, 2> RowVector2d;
    typedef Matrix<double, 2, 1> ColumnVector2d;


}


#endif // Matrix_h__
