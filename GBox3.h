/*!
* @file      GBox3.h
* @brief     三维包围盒的基本定义
* @details
* @author
* @date      2021/09/09
* @copyright Copyright 2012-2022 GLODON
**********************************************************************************
* @par 修改日志:
* |  Date  | Description | Author     |
* | :----: | :----       | :----      |
* | 2021/09/09 | 袁冠杰<yuangj-a> | GGP-55193 完成算法剩余部分的注释 |
*
**********************************************************************************
*/
#ifndef G_BOX3_H
#define G_BOX3_H

#include "GMath/GVec3.h"
#include "GMath/GBox2.h"
#include "GMath/GMatrix4.h"

#pragma warning( disable:4251)

namespace ggp
{
    /*!\addtogroup GMath GMath
    * @{
    */

    /*!
    * \class CBox3 
    * \brief 三维包围盒的基本定义
    * 该类给出三维包围盒的定义，提供相关的计算方法，被其它很多类使用
    */
    template<typename Type>
    /*!
    * @class CBox3
    * @brief 
    */
    class CBox3
    {
    public:
        
         /*!
        *@brief  构造缺省包围盒    
        */
        CBox3()
        {
            MakeEmpty();
        }

         /*!
        *@brief  根据包围盒的最小点和最大点构造包围盒
        *@param[in] min  最小点
        *@param[in] max  最大点    
        */
        CBox3(const CVec3<Type> & min, const CVec3<Type> & max)
        {
            Set(min, max);
        }

        CBox3(const CVec3<Type> & pt): m_min(pt), m_max(pt){}

        /*!
        *@brief  根据包围盒的最小点的坐标和最大点坐标构造包围盒
        *@param[in] MinX  最小点X坐标
        *@param[in] MinY  最小点Y坐标
        *@param[in] MinZ  最小点Z坐标
        *@param[in] MaxX  最大点X坐标
        *@param[in] MaxY  最大点Y坐标
        *@param[in] MaxZ  最大点Z坐标    
        */
        CBox3(const Type& MinX, const Type& MinY, const Type& MinZ, const Type& MaxX, const Type& MaxY, const Type& MaxZ)
        {
            Set(MinX, MinY, MinZ, MaxX, MaxY, MaxZ);
        }

         /*!
        *@brief  拷贝构造包围盒
        *@param[in] box  传入包围盒    
        */
        CBox3(const CBox3<Type> & box)
        {
            *this = box;
        }

        /*!
        *@brief  设置包围盒的最小点和最大点
        *@param[in] min  最小点
        *@param[in] max  最大点
        */
        void Set(const CVec3<Type> & min, const CVec3<Type> & max)
        {
            m_min = min;
            m_max = max;
        }

        /*!
        *@brief  设置包围盒的最小点和最大点
        *@param[in] pt  最小点与最大点
        */
        void Set(const CVec3<Type> & pt)
        {
            m_min = pt;
            m_max = pt;
        }

        /*!
        *@brief  设置包围盒的最小点坐标和最大点坐标
        *@param[in] MinX  最小点X坐标
        *@param[in] MinY  最小点Y坐标
        *@param[in] MinZ  最小点Z坐标
        *@param[in] MaxX  最大点X坐标
        *@param[in] MaxY  最大点Y坐标
        *@param[in] MaxZ  最大点Z坐标 
        */
         void Set(const Type& MinX, const Type& MinY, const Type& MinZ, const Type& MaxX, const Type& MaxY, const Type& MaxZ) 
         {
             m_min.Set(MinX,MinY,MinZ);
             m_max.Set(MaxX,MaxY,MaxZ);
         }
         /*!
         *@brief  根据数组数据设置包围盒
         *@param[in] p0   第一点X,Y,Z的值
         *@param[in] p1   第二点X,Y,Z的值
         */
         void Set(const Type* p0, const Type* p1) 
         {
             if (p0[0] < p1[0])   
                 m_min.X = p0[0], m_max.X = p1[0] ;
             else                               
                 m_min.X = p1[0], m_max.X = p0[0] ;

             if (p0[1] < p1[1])   
                 m_min.Y = p0[1], m_max.Y = p1[1] ;
             else                             
                 m_min.Y = p1[1], m_max.Y = p0[1] ;

             if (p0[2] < p1[2])   
                 m_min.Z = p0[2], m_max.Z = p1[2] ;
             else                               
                 m_min.Z = p1[2], m_max.Z = p0[2] ;
         }

         /*!
         *@brief  根据两点设置包围盒
         *@param[in] p0  第一点
         *@param[in] p1   第二点
         */
         void SafeSet(const CVec3<Type>& p0, const CVec3<Type>& p1) 
         {
             if (p0.X < p1.X) 
                 m_min.X = p0.X, m_max.X = p1.X ;
             else                           
                 m_min.X = p1.X, m_max.X = p0.X ;

             if (p0.Y < p1.Y)   
                 m_min.Y = p0.Y, m_max.Y = p1.Y ;
             else                           
                 m_min.Y = p1.Y, m_max.Y = p0.Y ;

             if (p0.Z < p1.Z)   
                 m_min.Z = p0.Z, m_max.Z = p1.Z ;
             else                           
                 m_min.Z = p1.Z, m_max.Z = p0.Z ;
         }
        
         /*!
         *@brief  根据点数组设置包围盒
         *@param[in] count        点数组的个数
         *@param[in] p  点数组
         */
         void Set(int count, CVec3<Type>* p)
         {
             if (count>0) 
             {
                 m_min = p[0] ; 
                 m_max = p[0] ;

                 for (int i=1; i<count; i++)
                 {
                     if (m_min.X>p[i].X)
                         m_min.X = p[i].X ;
                     else if (m_max.X<p[i].X)
                         m_max.X = p[i].X ;

                     if (m_min.Y>p[i].Y)
                         m_min.Y = p[i].Y ;
                     else if (m_max.Y<p[i].Y) 
                         m_max.Y = p[i].Y ;

                     if (m_min.Z>p[i].Z)
                         m_min.Z = p[i].Z ;
                     else if (m_max.Z<p[i].Z) 
                         m_max.Z = p[i].Z ;
                 }
             }
         }

         /*!
         *@brief  获得包围盒的8个角点
         *@param[out] corner  返回角点数组
         *@param[in] bSortByMinMax 是否按照XYZ的min到max的顺序，从000~111进行排列，若为否，则为逆时针顺序
         */
         void GetCorners(CVec3<Type> corner[8], bool bSortByMinMax = true) const 
         {
             if (bSortByMinMax)
             {
                 corner[0].Set(m_min.X, m_min.Y, m_min.Z);
                 corner[1].Set(m_max.X, m_min.Y, m_min.Z);
                 corner[2].Set(m_min.X, m_max.Y, m_min.Z);
                 corner[3].Set(m_max.X, m_max.Y, m_min.Z);
                 corner[4].Set(m_min.X, m_min.Y, m_max.Z);
                 corner[5].Set(m_max.X, m_min.Y, m_max.Z);
                 corner[6].Set(m_min.X, m_max.Y, m_max.Z);
                 corner[7].Set(m_max.X, m_max.Y, m_max.Z);
             }
             else
             {
                 corner[0].Set(m_min.X, m_min.Y, m_min.Z);
                 corner[1].Set(m_max.X, m_min.Y, m_min.Z);
                 corner[2].Set(m_max.X, m_max.Y, m_min.Z);
                 corner[3].Set(m_min.X, m_max.Y, m_min.Z);
                 corner[4].Set(m_min.X, m_min.Y, m_max.Z);
                 corner[5].Set(m_max.X, m_min.Y, m_max.Z);
                 corner[6].Set(m_max.X, m_max.Y, m_max.Z);
                 corner[7].Set(m_min.X, m_max.Y, m_max.Z);
             }
         }

        /*!
        *@brief  如果最小点 m_min 的分量大于最大点 m_max 的分量则交换
        */
         void Normalize()
         {
#define NORMALIZE(a, b) if(a > b){ auto temp = a; a = b; b = temp; }
             NORMALIZE(m_min.X, m_max.X);
             NORMALIZE(m_min.Y, m_max.Y);
             NORMALIZE(m_min.Z, m_max.Z);
#undef NORMALIZE
         }

       /*!
        *@brief  判断包围盒是否合法
        *@return 是否合法
        * - true 合法
        * - false 不合法
        */
        inline bool NotEmpty() const
        {
            return (m_min.X <= m_max.X + g_DoubleResolution 
                   && m_min.Y <= m_max.Y + g_DoubleResolution 
                   && m_min.Z <= m_max.Z + g_DoubleResolution);
        }

        /*!
        *@brief  初始化包围盒
        */
        void MakeEmpty()    //移除inline关键字（函数过长，没有inline必要）
        {
            if(std::numeric_limits<Type>::has_infinity)
            {
                const Type & c_infinity = std::numeric_limits<Type>::infinity();
                m_min.Set(c_infinity, c_infinity, c_infinity);
                m_max.Set(-c_infinity, -c_infinity, -c_infinity);
            }
            else
            {
#undef min
#undef max
                const Type & c_min = std::numeric_limits<Type>::min();
                const Type & c_max = std::numeric_limits<Type>::max();
                m_min.Set(c_max, c_max, c_max);
                m_max.Set(c_min, c_min, c_min);
            }
        }

        /*!
        *@brief  获得包围盒的最小点
        *@return 最小点
        */
        inline const CVec3<Type> & MinPt() const
        { 
            return m_min; 
        }

        /*!
        *@brief  获得包围盒的最大点
        *@return 最大点
        */
        inline const CVec3<Type> & MaxPt() const
        { 
            return m_max; 
        }

        /*!
        *@brief  返回包围盒的最小点的引用
        *@return 最小点
        */
        inline CVec3<Type> & MinPt() 
        { 
            return m_min; 
        }

        /*!
        *@brief  返回包围盒的最大点的引用
        *@return 最大点
        */
        inline CVec3<Type> & MaxPt() 
        { 
            return m_max; 
        }

        /*!
        *@brief  获得包围盒的长度，宽度和高度
        *@return 长度，宽度和高度,X坐标为长度，Y坐标为宽度，Z坐标为高度
        */
        CVec3<Type> GetSize() const
        {
            assert(NotEmpty());
            return m_max - m_min;
        }

        /*!
        *@brief  获得Box2
        *@return Box2
        */
        CBox2<Type> Box2() const
        {
            CBox2<Type> box2;
            box2.Set(m_min.Vec2(), m_max.Vec2());
            return box2;
        }

        /*!
        *@brief  获得Box3i
        *@return Box3i
        */
        CBox3<int> Box3i() const
        {
            return CBox3<int>((int)m_min.X, (int)m_min.Y, (int)m_min.Z, (int)m_max.X, (int)m_max.Y, (int)m_max.Z);
        }

        /*!
        *@brief  获得Box3f
        *@return Box3f
        */
        CBox3<float> Box3f() const
        {
            return CBox3<float>((float)m_min.X, (float)m_min.Y, (float)m_min.Z, (float)m_max.X, (float)m_max.Y, (float)m_max.Z);
        }

        /*!
        *@brief  获得Box3d
        *@return Box3d
        */
        CBox3<double> Box3d() const
        {
            return CBox3<double>((double)m_min.X, (double)m_min.Y, (double)m_min.Z, (double)m_max.X, (double)m_max.Y, (double)m_max.Z);
        }

        /*!
        *@brief  包围盒的中心点
        *@return 中心点
        */
        CVec3<Type> CenterPt() const
        {
            assert(NotEmpty());
            return CVec3<Type>(static_cast<Type>((m_max[0] + m_min[0]) * 0.5),
                static_cast<Type>((m_max[1] + m_min[1]) * 0.5),
                static_cast<Type>((m_max[2] + m_min[2]) * 0.5));
        }

        /*!
        *@brief  缩放外包盒，以盒的中心为中心进行缩放
        *@param[in] x  X方向缩放量
        *@param[in] y  Y方向缩放量
        *@param[in] z  Z方向缩放量
        *@param[in] asRatio  是否按照比例进行缩放，若否则按照绝对大小
        */
        void Expand(Type x, Type y, Type z, bool asRatio)
        {
            assert(NotEmpty());

            CVec3<Type> expandVec;
            if (asRatio)
            {
                CVec3<Type> boxDir = m_max - m_min;
                expandVec.Set(boxDir.X * x, boxDir.Y * y, boxDir.Z * z);
                
            }
            else
            {
                expandVec.Set(x, y, z);
            }

            expandVec.Set(expandVec.X / ((Type)(2)), expandVec.Y / ((Type)(2)), expandVec.Z / ((Type)(2)));
            m_max += expandVec;
            m_min -= expandVec;
        }

        /*!
        *@brief  当前包围盒和传入点合并，并设置当前包围盒为合并的结果
        *@param[in] point  传入点
        */
        void MergeBox(const CVec3<Type> & point)
        {
            if(point[0] < m_min[0]) 
                m_min[0] = point[0];
            if(point[1] < m_min[1])
                m_min[1] = point[1];
            if(point[2] < m_min[2]) 
                m_min[2] = point[2];

            if(point[0] > m_max[0]) 
                m_max[0] = point[0];
            if(point[1] > m_max[1]) 
                m_max[1] = point[1];
            if(point[2] > m_max[2])
                m_max[2] = point[2];
        }

        /*!
        *@brief  当前包围盒和传入包围盒合并，并设置当前包围盒为合并后的结果。操作类似Union。
        *@param[in] box 传入BOX
        */
        void MergeBox(const CBox3<Type> & box)
        {
            if(box.m_min[0] < m_min[0]) 
                m_min[0] = box.m_min[0];
            if(box.m_min[1] < m_min[1]) 
                m_min[1] = box.m_min[1];
            if(box.m_min[2] < m_min[2]) 
                m_min[2] = box.m_min[2];

            if (box.m_max[0] > m_max[0])
                m_max[0] = box.m_max[0];
            if (box.m_max[1] > m_max[1])
                m_max[1] = box.m_max[1];
            if (box.m_max[2] > m_max[2])
                m_max[2] = box.m_max[2];
        }

        /*!
        *@brief  获得当前包围盒和传入包围盒的交集
        *@param[in] Box       传入包围盒
        *@param[out] IntersectBox   交集
        *@param[in] tolerance                给定误差
        */
        void GetIntersectBox(const CBox3<Type> & Box, CBox3<Type> & IntersectBox, Type tolerance = std::numeric_limits<Type>::epsilon())     //移除inline关键字（函数过长，没有inline必要）
        {
            assert(NotEmpty());
            assert(Box.NotEmpty());
            IntersectBox = *this;

            if (MinPt().X < Box.MinPt().X)
            {
                IntersectBox.MinPt().X = Box.MinPt().X;
            }

            if (MinPt().Y < Box.MinPt().Y) 
            {
                IntersectBox.MinPt().Y = Box.MinPt().Y;
            }

            if (MinPt().Z < Box.MinPt().Z) 
            {
                IntersectBox.MinPt().Z = Box.MinPt().Z;
            }

            if (MaxPt().X > Box.MaxPt().X)
            {
                IntersectBox.MaxPt().X = Box.MaxPt().X;
            }

            if (MaxPt().Y > Box.MaxPt().Y) 
            {
                IntersectBox.MaxPt().Y = Box.MaxPt().Y;
            }

            if (MaxPt().Z > Box.MaxPt().Z) 
            {
                IntersectBox.MaxPt().Z = Box.MaxPt().Z;
            }
        }

        /*!
        *@brief  在给定误差内，判断当前包围盒是否包含传入点（在包围盒内或在包围盒上）
        *@param[in] point 传入点
        *@param[in] tolerance            给定误差
        *@return    是否包含
        * - true    包含
        * - false   不包含
        */
        bool Contains(const CVec3<Type> & point,Type tolerance = std::numeric_limits<Type>::epsilon()) const
        {
            assert(NotEmpty());
            
            return point.X >= m_min.X - tolerance 
                && point.X <= m_max.X + tolerance 
                && point.Y >= m_min.Y - tolerance 
                && point.Y <= m_max.Y + tolerance 
                && point.Z >= m_min.Z - tolerance 
                && point.Z <= m_max.Z + tolerance;
        }

        /*!
        *@brief  在给定误差内，判断传入包围盒是否和当前包围盒相等
        *@param[in] box  传入包围盒
        *@param[in] typeTol           给定误差
        *@return    是否相等
        * - true    相等 
        * - false   不相等
        */
        bool IsEqual(const CBox3<Type> & box, const Type typeTol = std::numeric_limits<Type>::epsilon()) const
        {
            return m_max.IsEqual(box.MaxPt(), typeTol) && m_min.IsEqual(box.MinPt(), typeTol);
        }

        /*!
        *@brief  在给定误差内，判断传入包围盒是否和当前包围盒相交
        *@param[in] box  传入包围盒
        *@param[in] tolerance           给定误差
        *@return    是否相交
        * - true    相交 
        * - false   不相交
        */
        bool IsIntersect(const CBox3<Type> & box,Type tolerance = std::numeric_limits<Type>::epsilon()) const
        {
            assert(NotEmpty());
            assert(box.NotEmpty());
            
            return (m_max.X >= box.MinPt().X - tolerance 
                && m_min.X <= box.MaxPt().X + tolerance 
                && m_max.Y >= box.MinPt().Y - tolerance 
                && m_min.Y <= box.MaxPt().Y + tolerance 
                && m_max.Z >= box.MinPt().Z  - tolerance 
                && m_min.Z <= box.MaxPt().Z + tolerance);
        }

        /*!
        *@brief  判断线段是否和包围盒相交
        *@param[in] rStartPt  线段起点
        *@param[in] rEndPt    线段终点
        *@param[in] tolerance                给定误差
        *@return    是否相交     
        * - true    相交 
        * - false   不相交
        */
        bool IsIntersect( const CVec3<Type> & rStartPt, const CVec3<Type> & rEndPt, Type tolerance = std::numeric_limits<Type>::epsilon()) const    //移除inline关键字（函数过长，没有inline必要）
        {
            assert(NotEmpty());

            CVec3<Type> c = ((Type)0.5)*(m_min + m_max);     //包围盒中心
            CVec3<Type> e = m_max - c;                       //包围盒半长向量
            CVec3<Type> m = ((Type)0.5)*(rStartPt + rEndPt); //线段中点
            CVec3<Type> d = rEndPt - m;                      //线段半长向量
            m = m - c;                                       //变换 box 和线段到原点

            //尝试世界坐标轴作为分离轴
            Type adx = abs(d.X);
            e.X += tolerance;
            if(abs(m.X) > e.X + adx) return false;

            Type ady = abs(d.Y);
            e.Y += tolerance;
            if(abs(m.Y) > e.Y + ady) return false;

            Type adz = abs(d.Z);
            e.Z += tolerance;
            if(abs(m.Z) > e.Z + adz) return false;

            //增加容差项。当线段几乎平行于坐标轴的时候误差会增大
            if(abs(m.Y*d.Z - m.Z*d.Y) > e.Y*adz + e.Z*ady) return false;

            if(abs(m.Z*d.X - m.X*d.Z) > e.Z*adx + e.X*adz) return false;

            if(abs(m.X*d.Y - m.Y*d.X) > e.X*ady + e.Y*adx) return false;

            return true;
        }

        /*!
        *@brief  判断线段是否和包围盒相交，并求出裁剪点
        *@param[in] rStartPt  线段起点
        *@param[in] rEndPt    线段终点
        *@param[in] isComputeClipPoints      是否计算出裁剪点
        *@param[out] clipedStartPt   被裁减的起点
        *@param[out] clipedEndPt     被裁减的终点
        *@param[in] tolerance                给定误差
        *@return    是否相交
        * - true    相交 
        * - false   不相交
        */
        bool ClipSegmentLine(const CVec3<Type>& rStartPt, const CVec3<Type>& rEndPt, bool isComputeClipPoints, 
            CVec3<Type>& clipedStartPt, CVec3<Type>& clipedEndPt,Type tolerance = std::numeric_limits<Type>::epsilon()) const        //移除inline关键字（函数过长，没有inline必要）
        {
            clipedStartPt = rStartPt;
            clipedEndPt = rEndPt;

            if (clipedStartPt.X - clipedEndPt.X <= tolerance)
            {
                if (clipedEndPt.X - m_min.X < tolerance) return false;
                if (clipedStartPt.X - m_max.X > -tolerance) return false;

                if (isComputeClipPoints)
                {
                    if (clipedStartPt.X - m_min.X < tolerance)
                    {
                        clipedStartPt = clipedStartPt + (clipedEndPt - clipedStartPt) * 
                            (m_min.X - clipedStartPt.X) / (clipedEndPt.X - clipedStartPt.X);
                    }

                    if (rEndPt.X - m_max.X > -tolerance)
                    {
                        clipedEndPt = clipedStartPt + (clipedEndPt - clipedStartPt) * 
                            (m_max.X - clipedStartPt.X) / (clipedEndPt.X - clipedStartPt.X);
                    }
                }
            }
            else
            {
                if (clipedStartPt.X - m_min.X < tolerance) return false;
                if (clipedEndPt.X - m_max.X > -tolerance) return false;

                if (isComputeClipPoints)
                {
                    if (clipedEndPt.X - m_min.X < tolerance)
                    {
                        clipedEndPt = clipedStartPt + (clipedEndPt - clipedStartPt) * 
                            (m_min.X - rStartPt.X) / (clipedEndPt.X - clipedStartPt.X);
                    }

                    if (clipedStartPt.X - m_max.X > -tolerance)
                    {
                        clipedStartPt = clipedStartPt + (clipedEndPt - clipedStartPt) * 
                            (m_max.X - rStartPt.X) / (clipedEndPt.X - clipedStartPt.X);
                    }
                }
            }

            if (clipedStartPt.Y - clipedEndPt.Y <= tolerance)
            {
                if (clipedEndPt.Y - m_min.Y < tolerance) return false;
                if (clipedStartPt.Y - m_max.Y > -tolerance) return false;

                if (isComputeClipPoints)
                {
                    if (clipedStartPt.Y - m_min.Y < tolerance)
                    {
                        clipedStartPt = clipedStartPt + (clipedEndPt - clipedStartPt) * 
                            (m_min.Y - clipedStartPt.Y) / (clipedEndPt.Y - clipedStartPt.Y);
                    }

                    if (rEndPt.Y - m_max.Y > -tolerance)
                    {
                        clipedEndPt = clipedStartPt + (clipedEndPt - clipedStartPt) * 
                            (m_max.Y - clipedStartPt.Y) / (clipedEndPt.Y - clipedStartPt.Y);
                    }
                }
            }
            else
            {
                if (clipedStartPt.Y - m_min.Y < tolerance) return false;
                if (clipedEndPt.Y - m_max.Y > -tolerance) return false;

                if (isComputeClipPoints)
                {
                    if (clipedEndPt.Y - m_min.Y < tolerance)
                    {
                        clipedEndPt = clipedStartPt + (clipedEndPt - clipedStartPt) * 
                            (m_min.Y - rStartPt.Y) / (clipedEndPt.Y - clipedStartPt.Y);
                    }

                    if (clipedStartPt.Y - m_max.Y > -tolerance)
                    {
                        clipedStartPt = clipedStartPt + (clipedEndPt - clipedStartPt) * 
                            (m_max.Y - rStartPt.Y) / (clipedEndPt.Y - clipedStartPt.Y);
                    }
                }
            }

            if (clipedStartPt.Z - clipedEndPt.Z <= tolerance)
            {
                if (clipedEndPt.Z - m_min.Z < tolerance) return false;
                if (clipedStartPt.Z - m_max.Z > -tolerance) return false;

                if (isComputeClipPoints)
                {
                    if (clipedStartPt.Z - m_min.Z < tolerance)
                    {
                        clipedStartPt = clipedStartPt + (clipedEndPt - clipedStartPt) * 
                            (m_min.Z - clipedStartPt.Z) / (clipedEndPt.Z - clipedStartPt.Z);
                    }

                    if (rEndPt.Z - m_max.Z > -tolerance)
                    {
                        clipedEndPt = clipedStartPt + (clipedEndPt - clipedStartPt) * 
                            (m_max.Z - clipedStartPt.Z) / (clipedEndPt.Z - clipedStartPt.Z);
                    }
                }
            }
            else
            {
                if (clipedStartPt.Z - m_min.Z < tolerance) return false;
                if (clipedEndPt.Z - m_max.Z > -tolerance) return false;

                if (isComputeClipPoints)
                {
                    if (clipedEndPt.Z - m_min.Z < tolerance)
                    {
                        clipedEndPt = clipedStartPt + (clipedEndPt - clipedStartPt) * 
                            (m_min.Z - rStartPt.Z) / (clipedEndPt.Z - clipedStartPt.Z);
                    }

                    if (clipedStartPt.Z - m_max.Z > -tolerance)
                    {
                        clipedStartPt = clipedStartPt + (clipedEndPt - clipedStartPt) * 
                            (m_max.Z - rStartPt.Z) / (clipedEndPt.Z - clipedStartPt.Z);
                    }
                }
            }

            return true;
        }

        /*!
        *@brief  在给定误差下，包围盒是否在传入包围内，包括内相切
        *@param[in] box   传入包围盒
        *@param[in] tolerance            给定误差
        *@return    是否在传入包围盒内    
        * - true    当前包围盒在传入包围盒内
        * - false   当前包围盒不在传入包围盒内
        */
        inline bool In(const CBox3<Type> & box, Type tolerance = std::numeric_limits<Type>::epsilon()) const
        {  
            assert(NotEmpty());
            assert(box.NotEmpty());

            return (box.MinPt().X <= m_min.X + tolerance
                && box.MaxPt().X >=  m_max.X - tolerance
                && box.MinPt().Y <=  m_min.Y + tolerance 
                && box.MaxPt().Y >=  m_max.Y - tolerance 
                && box.MinPt().Z <=  m_min.Z + tolerance 
                && box.MaxPt().Z >=  m_max.Z - tolerance);
        }

         /*!
         *@brief  平移当前包围盒
         *@param[in] Vt  平移量
         */
         inline void Translate( const CVec3<Type> & Vt) 
         {
              assert(NotEmpty());
              m_min += Vt;
              m_max += Vt;
         }
      
         /*!
         *@brief  将rBox3d投影到局部坐标系中, 获得局部坐标系下的包围盒
         *@param[in] localCoord  局部坐标系
         *@return  局部坐标系下的包围盒
         */
         CBox3<Type> LocalBox(const CCoord3<Type>& localCoord)
         {
             CVec3<Type> corners[8];
             GetCorners(corners);

             CBox3<Type> localBox;
             for (int i = 0; i < 8; ++i)
             {
                 localBox.MergeBox(corners[i].LocalPoint(localCoord));
             }

             return localBox;
         }

         /*!
         *@brief  将rBox3d投影到局部坐标系中, 获得局部坐标系下的包围盒
         *@param[in] localCoord  局部坐标系
         *@return  局部坐标系下的包围盒
         */
         CBox3<Type> WorldBox(const CCoord3<Type>& localCoord)
         {
             CVec3<Type> corners[8];
             GetCorners(corners);

             CBox3<Type> worldBox;
             for (int i = 0; i < 8; ++i)
             {
                 worldBox.MergeBox(corners[i].WorldPoint(localCoord));
             }

             return worldBox;
         }

         void Transform(const CMatrix4<Type> & mat)
         {
             auto m = CenterPt();
             auto h = m_max - m;

             m = mat.MultiPointLeft(m);
             auto v0 = mat.MultiVecLeft(CVec3<Type>(h[0], 0.0, 0.0));
             auto v1 = mat.MultiVecLeft(CVec3<Type>(0.0, h[1], 0.0));
             auto v2 = mat.MultiVecLeft(CVec3<Type>(0.0, 0.0, h[2]));

             h = v0.Abs() + v1.Abs() + v2.Abs();

             m_min = m - h;
             m_max = m + h;
         }

         /*!
        *@brief  判断当前包围盒外一点到当前包围盒的最短距离
        *@param[in] Pt   包围盒外一点
        *@return    最短距离
        */
         inline Type Distance( const CVec3<Type> & Pt) const
         {
             assert(NotEmpty());
             CVec3<Type> vMindist;
             MindistVector(Pt, vMindist);
             return vMindist.Length();
         }


        /*!
        *@brief  计算当前包围盒外一点到当前包围盒的最短距离的平方
        *@param[in] pt   包围盒外一点
        *@return    最短距离的平方
        */
         inline Type Distance2( const CVec3<Type> & pt) const
         {
             assert(NotEmpty());
             CVec3<Type> vMindist;
             MindistVector(pt, vMindist);
             return vMindist.SqrLength();
         }

         /*!
         *@brief  计算当前包围盒外一点到当前包围盒的最长距离的平方
         *@param[in] pt   包围盒外一点
         *@return    最长距离的平方
         */
         inline Type MaxDistance2(const CVec3<Type>& pt) const
         {
             CVec3<Type> v_max;
             MaxdistVector(pt, v_max);
             return v_max.SqrLength();
         }

         /*!
         *@brief  计算一个点到包围盒最小距离点距离等长的向量
         *@param[in]  pt  点
         *@param[in,out] v_maxdist 点到包围盒最小距离点距离等长的向量
         */
         void MaxdistVector(const CVec3<Type> & pt, CVec3<Type> & v_maxdist) const    //移除inline关键字（函数过长，没有inline必要）
         {
             CVec3<Type> center = CenterPt(); //box中心
             CVec3<Type> half = m_max - center;

             v_maxdist = pt - center;

             //将v向量映射到右上角
             v_maxdist.X = abs(v_maxdist.X);
             v_maxdist.Y = abs(v_maxdist.Y);
             v_maxdist.Z = abs(v_maxdist.Z);

             v_maxdist += half;
         }

         /*!
        *@brief  把包围盒对象按照规定格式字符串输出
        *@return 字符串
        */
         char* AsString()
         {
             char * a = AsStringHelper::NewString("Box3d=[(%.20le %.20le %.20le),(%.20le %.20le %.20le)]", 
                 m_min.X, m_min.Y, m_min.Z, m_max.X, m_max.Y, m_max.Z);
             return a;

         }

         /*!
        *@brief     根据字符串生成包围盒
        *@param[in] str   传入的字符串
        *@return    生成的包围盒
        */
         static CBox3<Type>* LoadFromString(char* str)
         {
             CBox3<Type>* pBox = new CBox3<Type>;
             g_sscanf(str, "Box3d=[(%le %le %le),(%le %le %le)]", 
                 &pBox->MinPt().X, &pBox->MinPt().Y, &pBox->MinPt().Z, &pBox->MaxPt().X, &pBox->MaxPt().Y, &pBox->MaxPt().Z);
             return pBox;

         }

        /*
        *@brief     赋值操作符
        *@param[in] box  传入的包围盒
        *@return    
        */
        /*CBox3<Type> & operator = (const CBox3<Type> & box)
        {
            Set(box.MinPt(), box.MaxPt());
            return *this;
        }*/

         /*!
        *@brief     逻辑等号操作符，判断两个BOX是否相等
        *@param[in] b1  传入的包围盒1
        *@param[in] b2  传入的包围盒2
        *@return    是否相等   
        * - true    相等
        * - false   不相等
        */
        friend bool operator ==(const CBox3<Type> & b1, const CBox3<Type> & b2)
        { 
            assert(b1.NotEmpty());
            assert(b2.NotEmpty());
            return b1.MinPt() == b2.MinPt() && b1.MaxPt() == b2.MaxPt(); 
        }

        /*!
        *@brief     逻辑不等号操作符，判断两个BOX是否不相等
        *@param[in] b1  传入的包围盒1
        *@param[in] b2  传入的包围盒2
        *@return    是否不相等   
        * - true    不相等
        * - false   相等
        */
        friend bool operator !=(const CBox3<Type> & b1, const CBox3<Type> & b2)
        { 
            assert(b1.NotEmpty());
            assert(b2.NotEmpty());
            return !(b1 == b2); 
        }

        /*!
        *@brief    Minkowski差
        *@param[in] b1  包围盒1
        *@param[in] b2  包围盒2
        *@return    Minkowski差结果
        */
        friend CBox3<Type> operator -(const CBox3<Type> & b1, const CBox3<Type> & b2)
        {
            auto center1 = b1.CenterPt();
            auto center2 = b2.CenterPt();

            auto center = center1 - center2;
            auto half1 = b1.MaxPt() - center1;
            auto half2 = b2.MaxPt() - center2;

            auto half = half1 + half2;

            return CBox3<Type>(center - half, center + half);
        }

        /*!
        *@brief     Minkowski和
        *@param[in] b1   包围盒1
        *@param[in] b2   包围盒2
        *@return    Minkowski和结果
        */
        friend CBox3<Type> operator +(const CBox3<Type> & b1, const CBox3<Type> & b2)
        {
            auto center1 = b1.CenterPt();
            auto center2 = b2.CenterPt();

            auto center = center1 + center2;
            auto half1 = b1.MaxPt() - center1;
            auto half2 = b2.MaxPt() - center2;

            auto half = half1 + half2;

            return CBox3<Type>(center - half, center + half);
        }


    private:
        //计算一个点到包围盒最小距离点距离等长的向量
        void MindistVector(const CVec3<Type> & pt, CVec3<Type> & vMinDist) const    //移除inline关键字（函数过长，没有inline必要）
        {
            CVec3<Type> center = CenterPt(); //box中心
            CVec3<Type> half = m_max - center;

            vMinDist = pt - center;

            //将v向量映射到右上角
            vMinDist.X = abs(vMinDist.X);
            vMinDist.Y = abs(vMinDist.Y);
            vMinDist.Z = abs(vMinDist.Z);

            vMinDist -= half;

            if(vMinDist.X < (Type)0.0) 
                vMinDist.X = (Type)0.0;

            if(vMinDist.Y < (Type)0.0)
                vMinDist.Y = (Type)0.0;

            if(vMinDist.Z < (Type)0.0)
                vMinDist.Z = (Type)0.0;
        }

    private:
        CVec3<Type> m_min;
        CVec3<Type> m_max;
    };

    typedef CBox3<int>    CBox3i;
    typedef CBox3<float>  CBox3f;
    typedef CBox3<double> CBox3d;

    /*! @} */
} 

#endif
