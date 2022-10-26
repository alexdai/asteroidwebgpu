#pragma once

#include "GMathDef.h"

#pragma warning(disable:4251)


template<typename Type>
struct CVec3;

template<typename Type>
struct CVec1
{
public:
    CVec1():X(0){}

    CVec1(const Type &x):X(x){}

    CVec1(const CVec1<Type>& vecSrc): X(vecSrc.X){}


    void Set(const Type v)
    {
        X = v;
    }


    Type &  operator[](int nIndex)
    {
        return m_x[nIndex];
    }

    const Type & operator[](int nIndex) const
    {
        return m_x[nIndex];
    }

    CVec1<Type> & operator =(const CVec1<Type> & vecSrc)
    {
        X = vecSrc.X;
        return *this;
    }

    CVec1<Type> & operator *=(const Type &d)
    {
        X *= d;

        return *this;
    }

    CVec1<Type> & operator /=(const Type &d)
    {
        X/= d;

        return *this;
    }

    CVec1<Type> & operator -=(const CVec1<Type> & vecSrc)
    {
        X -= vecSrc.X;

        return *this;
    }

    CVec1<Type> operator-() const
    {
        return CVec1<Type>(-X);
    }


    friend CVec1<Type> operator /(const CVec1<Type> & vecSrc, const Type &d)
    {
        return CVec1<Type>(vecSrc.X/d);
    }

    CVec1<Type> operator *(const CVec1 & v) const
    {
        return CVec1(X*v.X);
    }

    CVec1<Type> operator *(Type d) const
    {
        return CVec1(X*d);
    }

    friend CVec1<Type> operator*(Type d, const CVec1<Type> & v)
    {
        return CVec1(d*v.X);
    }

    CVec1<Type> & operator+=(const CVec1<Type> &v)
    {
        X += v.X;
        return *this;
    }

    friend CVec1<Type> operator +(const CVec1<Type> & v1, const CVec1<Type> & v2)
    {
        return CVec1<Type>(v1.X + v2.X);
    }

    friend CVec1<Type> operator -(const CVec1<Type> & v1, const CVec1<Type> & v2)
    {
        return CVec1<Type>(v1.X - v2.X);
    }

    friend bool operator ==(const CVec1<Type> & v1, const CVec1<Type> & v2)
    {
        return v1.X==v2.X;
    }

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
