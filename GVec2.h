#ifndef G_VEC2_H
#define G_VEC2_H

#include "GMathDef.h"
#include <complex>

#pragma warning( disable:4251)

template<typename Type>
struct CVec3;

template<typename Type>
struct CVec2
{
public:
    CVec2():X(0), Y(0){}

    explicit CVec2(const Type v[2]): X(v[0]),Y(v[1]) {}

    CVec2(const Type &x, const Type &y):X(x), Y(y) {}

    CVec2(const CVec2<Type>& vecSrc): X(vecSrc.X),Y(vecSrc.Y) {}

    CVec2<int> Vec2i() const {return CVec2<int>((int)X, (int)Y);}

    CVec2<float>  Vec2f() const {return CVec2<float>((float)X, (float)Y); }

    CVec2<double> Vec2d() const {return CVec2<double>((double)X, (double)Y); }

    std::complex<Type> AsComplex() { return complex<Type>(X, Y); }

    CVec2<Type> YX() const { return CVec2<Type>(Y, X); }

    CVec3<Type> Vec3(Type zValue = 0) const { return CVec3<Type>((Type)X, (Type)Y, zValue); }

    CVec2<Type> Conjugate() const {return CVec2<Type>(X, -Y); }

    CVec2<Type> & ToLocalPoint(const CVec2<Type> & localOrg, const CVec2<Type> & localVx, const CVec2<Type> & localVy) 
    {
        CVec2<Type> tempVec(m_xy[0] - localOrg[0],m_xy[1] - localOrg[1]);

        m_xy[0] = tempVec * localVx;
        m_xy[1] = tempVec * localVy;

        return *this;
    }

    CVec2<Type> & ToWorldPoint(const CVec2<Type> & localOrg, const CVec2<Type> & localVx, const CVec2<Type> & localVy) 
    {
        Type x = m_xy[0];
        Type y = m_xy[1];

        m_xy[0] = localOrg.X + x * localVx.X + y * localVy.X;
        m_xy[1] = localOrg.Y + x * localVx.Y + y * localVy.Y;

        return *this;
    }

    CVec2<Type>  GetLocalPt(const CVec2<Type> & localOrg, const CVec2<Type> & localVx, const CVec2<Type> & localVy) const
    {
        CVec2<Type> tempVec(m_xy[0] - localOrg[0],m_xy[1] - localOrg[1]);
        return CVec2(tempVec.Dot(localVx), tempVec.Dot(localVy));
    }

    CVec2<Type>  GetWorldPt(const CVec2<Type> & localOrg, const CVec2<Type> & localVx, const CVec2<Type> & localVy) const
    {
        return localOrg + m_xy[0] * localVx + m_xy[1] * localVy;
    }

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

    CVec2<Type> Rotate(const CVec2<Type>& v)
    {
        return CVec2<Type>(X*v.X - Y*v.Y, X*v.Y + Y*v.X);//X分量对应cos,Y对应sin
    }

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

    CVec2<Type> GetRotateVector(Type angle) const
    {
        double cosVal;
        double sinVal;
        sincos(angle, &sinVal, &cosVal);

        return CVec2(Type(m_xy[0] * cosVal - m_xy[1] * sinVal), Type(m_xy[0] * sinVal + m_xy[1] * cosVal));
    }

    CVec2<Type> GetRotateHalfPIVector(bool bClockwise) const
    {
        CVec2<Type> resultVec(*this);
        return resultVec.RotateHalfPI(bClockwise);
    }

    CVec2<Type> GetRotateAroundPoint(const CVec2<Type> & vecCenter, Type typeAngle) const
    {
        CVec2<Type> tempVec((*this) - vecCenter);
        return CVec2(vecCenter + tempVec.Rotate(typeAngle));
    }

    Type AngleFromXAxis() const
    {
        double dAngle = atan2(m_xy[1], m_xy[0]);
        return Type((dAngle < 0) ? dAngle + M_2PI : dAngle);
    }

    CVec2<Type> Abs() const
    {
        return CVec2<Type>(std::abs(X), std::abs(Y));
    }

    Type Cos(const CVec2<Type> & vec) const
    {
        Type value = ((*this)*vec)/sqrt(double((this->SqrLength()*vec.SqrLength())));
        Clamp((Type) -1.0, value, (Type) 1.0);
        return value;
    }


    Type Sin(const CVec2<Type> & vec) const
    {
        Type value = abs((*this)^vec)/sqrt(double(SqrLength()*vec.SqrLength()));
        Clamp((Type) 0, value, (Type) 1.0);
        return value;
    }


    Type SinUnit(const CVec2<Type> & vec) const
    {
        Type value = abs(double((*this)^vec));
        Clamp((Type) 0, value, (Type) 1.0);
        return value;
    }


    void Set(const Type v[2])
    {
        m_xy[0] = v[0];
        m_xy[1] = v[1];
    }

    void Set(const Type &iX, const Type &iY)
    {
        m_xy[0] = iX;
        m_xy[1] = iY;
    }

    const Type * Value() const
    {
        return &m_xy[0];
    }

    void Value(Type & oX, Type & oY) const
    {
        oX = m_xy[0];
        oY = m_xy[1];
    }

    Type Value(const int& index) const
    {
        return m_xy[index];
    }

    Type& Value(const int& index) 
    {
        return m_xy[index];
    }


    Type Dot(const CVec2<Type> & vecSrc) const
    {
        return (m_xy[0]*vecSrc[0] + m_xy[1]*vecSrc[1]);
    }

    Type Dot(Type x, Type y) const
    {
        return (X*x + Y*y);
    }

    Type Length() const
    {
        return Type(std::sqrt( double(SqrLength()) ));
    }

    Type SqrLength() const
    {
        return (m_xy[0]*m_xy[0])+(m_xy[1]*m_xy[1]);
    }

    CVec2<Type> Square() const
    {
        return CVec2<Type>(X*X, Y*Y);
    }


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

    CVec2<Type> Unit(Type eps = std::numeric_limits<Type>::epsilon())
    {
        Type magnitude = Length();

        if(!IsNearZero(magnitude, eps))
            return CVec2<Type>(X / magnitude, Y / magnitude);
        else
            return CVec2<Type>(0.0, 0.0);
    }


    Type Cross(const CVec2<Type> & vecSrc) const
    {
        return (m_xy[0] * vecSrc[1] - m_xy[1] * vecSrc[0]);
    }

    void Negate()
    {
        Set(-m_xy[0], -m_xy[1]);
    }

    inline const Type* Ptr() const
    {
        return &m_xy[0];
    }

    inline Type* Ptr()
    {
        return &m_xy[0];
    }

    Type &  operator[](int nIndex)
    { 
        return m_xy[nIndex]; 
    }

    const Type & operator[](int nIndex) const
    { 
        return m_xy[nIndex]; 
    }

     CVec2<Type> & operator =(const CVec2<Type> & vecSrc)
    {
        m_xy[0] = vecSrc.m_xy[0];
        m_xy[1] = vecSrc.m_xy[1];

        return *this;
    }

    CVec2<Type> & operator *=(const Type &d)
    {
        m_xy[0] *= d;
        m_xy[1] *= d;

        return *this;
    }

    CVec2<Type> & operator /=(const Type &d)
    {
        Type inv = Type(1.0/d);

        m_xy[0] *= inv;
        m_xy[1] *= inv;

        return *this;
    }

    CVec2<Type> & operator +=(const CVec2<Type> & vecSrc)
    {
        m_xy[0] += vecSrc.m_xy[0];
        m_xy[1] += vecSrc.m_xy[1];

        return *this;
    }

    CVec2<Type> & operator -=(const CVec2<Type> & vecSrc)
    {
        m_xy[0] -= vecSrc.m_xy[0];
        m_xy[1] -= vecSrc.m_xy[1];

        return *this;
    }

    CVec2<Type> operator-() const
    {
        return CVec2<Type>(-m_xy[0], -m_xy[1]);
    }

    int MaxDimension() const
    {
        return (std::abs(X) >= std::abs(Y)) ? 0 : 1;
    }

    int MinDimension() const
    {
        return (std::abs(X) < std::abs(Y)) ? 0 : 1;
    }


    CVec2<Type> Mul (const CVec2<Type> & v)
    {
        return CVec2<Type>(X*v.X, Y*v.Y);
    }

    friend CVec2<Type> operator *(const CVec2<Type> & vecSrc, const Type &d)
    { 
        return CVec2<Type>(vecSrc.m_xy[0] * d, vecSrc.m_xy[1] * d);
    }

    friend CVec2<Type> operator *(const Type &d, const CVec2<Type> & vecSrc)
    { 
        return vecSrc * d; 
    }


    friend Type operator *(const CVec2<Type> & v1, const CVec2<Type> & v2)
    {
        return (v1.X*v2.X + v1.Y*v2.Y);
    }

    friend Type operator^(const CVec2<Type> & v1, const CVec2<Type> & v2)
    {
        return v1.Cross(v2); 
    }


    friend CVec2<Type> operator /(const CVec2<Type> & vecSrc, const Type &d)
    { 
        Type d1 = ((Type)(1))/d;
        return CVec2<Type>(vecSrc.m_xy[0] * d1, vecSrc.m_xy[1]*d1);
    }

    Type AngleTo(const CVec2<Type>& v) const
    {
        return std::atan2(Cross(v), Dot(v));
    }

    CVec2<Type> Div(const CVec2<Type> & v)
    {
        return CVec2<Type>(X/v.X, Y/v.Y);
    }

    friend CVec2<Type> operator +(const CVec2<Type> & v1, const CVec2<Type> & v2)
    {
        return CVec2<Type>(v1.m_xy[0] + v2.m_xy[0],v1.m_xy[1] + v2.m_xy[1]);
    }

    friend CVec2<Type> operator -(const CVec2<Type> & v1, const CVec2<Type> & v2)
    {
        return CVec2<Type>(v1.m_xy[0] - v2.m_xy[0], v1.m_xy[1] - v2.m_xy[1]);
    }

    friend CVec2<Type> operator -(const CVec2<Type> & v, Type d)
    {
        return CVec2<Type>(v.X - d, v.Y - d);
    }

    friend bool operator ==(const CVec2<Type> & v1, const CVec2<Type> & v2)
    { 
        return v1.m_xy[0]==v2.m_xy[0] && v1.m_xy[1]==v2.m_xy[1]; 
    }

    friend bool operator !=(const CVec2<Type> & v1, const CVec2<Type> & v2)
    {
        return !(v1 == v2); 
    }
        
    bool IsParallel(const CVec2<Type>& dir, const Type& typeTol=std::numeric_limits<Type>::epsilon()) const
    {
        Type dis = Cross(dir) ;
        //return fabs(dis) < typeTol;
        return dis * dis <= SqrLength() * dir.SqrLength() * typeTol * typeTol;
    }

    bool IsPerpendicular(const CVec2<Type>& dir, const Type& typeTol=std::numeric_limits<Type>::epsilon()) const
    {
        Type dis = Dot(dir) ;
        //return fabs(dis) < typeTol;
        return dis * dis <= SqrLength() * dir.SqrLength() * typeTol * typeTol;
    }

    bool IsEqual(const CVec2<Type> & vecSrc, const Type typeTol=std::numeric_limits<Type>::epsilon()) const
    {
        return ( (*this - vecSrc).SqrLength() <= typeTol*typeTol );
    }
                
    bool IsZero(const Type typeTol=std::numeric_limits<Type>::epsilon()) const
    {
        return ( SqrLength() <= typeTol*typeTol );
    }

    bool IsUnit(const Type typeTol=std::numeric_limits<Type>::epsilon()) const
    {
        return fabs(SqrLength() - 1) <= 2*typeTol;
    }


    bool IsValid() const { return (is_valid(X) && is_valid(Y)); }

    Type DistanceTo(const CVec2<Type> & src) const
    {
        double dx = (X - src.X);
        double dy = (Y - src.Y);
        return static_cast<Type>(sqrt(dx * dx + dy * dy));
    }

    Type SqrDistanceTo(const CVec2<Type> & src) const
    {
        return ((*this) - src).SqrLength();
    }

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


#endif  //G_VEC2_H
