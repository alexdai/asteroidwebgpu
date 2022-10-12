/*!
* @file      GBox2.h
* @brief     二维包围盒的基本定义
* @details   该类GBox2.h是一个模板类，预定义CBox2i,CBox2f和CBox2d三个类型，提供平行于坐标轴的包围盒的描述和相关算法。
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
#ifndef G_BOX2_H
#define G_BOX2_H

#include "GVec2.h"

#pragma warning( disable:4251)

namespace ggp
{
    /*!\addtogroup GMath GMath
    * 基础类型定义， 包括包围盒、坐标系、矩阵、向量、区间、四元数、三角形等基本类型
    * @{
    */

    template<typename Type>
    class CBox3;

    /*!
    * \class CBox2 
    * \brief 二维包围盒的基本定义
    * 该类给出二维包围盒的定义，提供相关的计算方法，被其它很多类使用
    */
    template<typename Type>
    class CBox2
    {
    public:
        
        /*!
        *@brief  缺省构造包围盒    
        */
        CBox2()
        {
            MakeEmpty();
        }
       
        CBox2(const CVec2<Type>& pt)
        {
            Set(pt,pt);
        }

        /*!
        *@brief  通过BOX的最小点和最大点构造包围盒
        *@param[in] min  最小点
        *@param[in] max  最大点
        */
        CBox2(const CVec2<Type> & min, const CVec2<Type> & max)
        {
            Set(min,max);
        }

        /*!
        *@brief  通过包围盒的最小点的坐标和最大点坐标构造包围盒
        *@param[in] minx  最小点X坐标
        *@param[in] miny  最小点Y坐标
        *@param[in] maxx  最大点X坐标
        *@param[in] maxy  最大点Y坐标
        */
        CBox2(const Type& minx, const Type& miny,const Type& maxx, const Type& maxy)
        {
            Set(minx,miny,maxx,maxy);
        }
        
        /*!
        *@brief  通过点集构造包围盒
        *@param[in] pPts 点数组指针
        *@param[in] nCount 点的个数
        */
        CBox2(const CVec2<Type>* pPts, const int nCount)
        {
            assert(nCount > 0);
            MakeEmpty();
            for (int i = 0; i < nCount; i++)
            {
                MergeBox(pPts[i]);
            }
        }

        CBox2(const vector<vector<CVector2d>> &ptGroups)
        {
            MakeEmpty();
            for (int i = 0; i < ptGroups.size(); i++)
            {
                for (int j = 0; j < ptGroups[i].size(); j++)
                {
                    MergeBox(ptGroups[i][j]);
                }
            }
        }

        /*!
        *@brief  通过包围盒中心点的坐标和X与Y方向上的范围构造包围盒
        *@param[in] centerPos  最小点X坐标
        *@param[in] extentX  最小点Y坐标
        *@param[in] extentY  最大点X坐标
        */
        CBox2(const CVec2<Type>& centerPos, const Type& extentX, const Type& extentY)
        {
            m_min.X = centerPos.X - extentX / (Type)2;
            m_min.Y = centerPos.Y - extentY / (Type)2;
            m_max.X = centerPos.X + extentX / (Type)2;
            m_max.Y = centerPos.Y + extentY / (Type)2;
        }

        /*!
        *@brief  获得Box2i
        *@return Box2i
        */
        CBox2<int> Box2i() const
        {
            return CBox2<int>((int)m_min.X, (int)m_min.Y, (int)m_max.X, (int)m_max.Y);
        }

        /*!
        *@brief  获得Box2f
        *@return Box2f
        */
        CBox2<float> Box2f() const
        {
            return CBox2<float>((float)m_min.X, (float)m_min.Y, (float)m_max.X, (float)m_max.Y);
        }

        /*!
        *@brief  获得Box2d
        *@return Box2d
        */
        CBox2<double> Box2d() const
        {
            return CBox2<double>((double)m_min.X, (double)m_min.Y, (double)m_max.X, (double)m_max.Y);
        }

        /*!
        *@brief  通过输入二维向量获得Box3
        *@param[in] minMaxZ 输入二维向量
        *@return Box3
        */
        CBox3<Type> Box3(const CVec2<Type>& minMaxZ = CVec2<Type>::Zero) const
        {
            CBox3<Type> box3;
            box3.Set(m_min.Vec3(minMaxZ.X), m_max.Vec3(minMaxZ.Y));
            return box3;
        }

        /*!
        *@brief  盒子移动
        *@param[in] moveVector 移动向量
        */
        inline void Translate(const CVec2<Type>& moveVector)
        {
            assert(NotEmpty());
            m_min += moveVector;
            m_max += moveVector;
        }


        /*!
        *@brief  拷贝构造包围盒
        *@param[in] box  传入的包围盒
        */
        CBox2(const CBox2<Type> & box)
        {
            *this = box;
        }

       
        /*!
        *@brief  设置包围盒的最小点坐标和最大点坐标
        *@param[in] minx  最小点X坐标
        *@param[in] miny  最小点Y坐标
        *@param[in] maxx  最大点X坐标
        *@param[in] maxy  最大点Y坐标
        */
        inline void Set(const Type& minx, const Type& miny,const Type& maxx, const Type& maxy)
        {
            m_min.Set(minx,miny);
            m_max.Set(maxx,maxy);    
        }

        /*!
        *@brief  设置包围盒的最小点和最大点
        *@param[in] min  最小点
        *@param[in] max  最大点   
        */
        void Set(const CVec2<Type> & min, const CVec2<Type> & max)
        {
            m_min = min;
            m_max = max;
        }

        /*!
        *@brief  设置包围盒的最小点和最大点
        *@param[in] pt  最小点与最大点
        */
        void Set(const CVec2<Type> & pt)
        {
            m_min = pt;
            m_max = pt;
        }
       
        /*!
        *@brief  判断包围盒是否合法
        *@return 是否合法   
        * - true 合法
        * - false 不合法
        */
        inline bool NotEmpty() const
        {
            return (m_min.X <= m_max.X + g_DoubleResolution) && (m_min.Y <= m_max.Y + g_DoubleResolution);
        }

        /*!
        *@brief  初始化包围盒
        */
        void MakeEmpty()    //移除inline关键字（函数过长，没有inline必要）
         {
            if(std::numeric_limits<Type>::has_infinity)
            {
                const Type & c_infinity = std::numeric_limits<Type>::infinity();
                m_min.Set(c_infinity, c_infinity);
                m_max.Set(-c_infinity, -c_infinity);
            }
            else
            {
#undef min
#undef max
                const Type & c_min = std::numeric_limits<Type>::min();
                const Type & c_max = std::numeric_limits<Type>::max();
                m_min.Set(c_max, c_max);
                m_max.Set(c_min, c_min);
            }
        }

        /*!
        *@brief  获得包围盒的最小点
        *@return 最小点
        */
        inline const CVec2<Type> & MinPt() const
        { 
            return m_min; 
        }

        /*!
        *@brief  获得包围盒的最大点
        *@return 最大点
        */
        inline const CVec2<Type> & MaxPt() const
        { 
            return m_max; 
        }

         /*!
        *@brief  返回包围盒的最小点的引用
        *@return 最小点
        */
        inline CVec2<Type> & MinPt() 
        { 
            return m_min; 
        }

       /*!
        *@brief  返回包围盒的最大点的引用
        *@return 最大点
        */
        inline CVec2<Type> & MaxPt() 
        { 
            return m_max; 
        }

       
        /*!
        *@brief  获得包围盒的宽度和高度
        *@return 宽度和高度，X坐标为包围盒的宽度，Y坐标为包围盒的高度
        */
        inline CVec2<Type> GetSize() const
        {
            assert(NotEmpty());
            return m_max - m_min;
        }

        
        /*!
        *@brief     获得包围盒的中心点
        *@return    中心点
        */
        inline CVec2<Type> CenterPt() const
        {
            assert(NotEmpty());
            return CVec2<Type>(static_cast<Type>((m_max[0] + m_min[0]) * 0.5),
                static_cast<Type>((m_max[1] + m_min[1]) * 0.5));
        }

        /*!
        *@brief  缩放外包盒，以盒的中心为中心进行缩放
        *@param[in] x  X方向缩放量
        *@param[in] y  Y方向缩放量
        *@param[in] asRatio  是否按照比例进行缩放，若否则按照绝对大小
        *@return  自身的引用
        */
        CBox2<Type>& Expand(Type x, Type y, bool asRatio)
        {
            assert(NotEmpty());

            CVec2<Type> expandVec;
            if (asRatio)
            {
                CVec2<Type> boxDir = m_max - m_min;
                expandVec.Set(boxDir.X * x, boxDir.Y * y);
            }
            else
            {
                expandVec.Set(x, y);
            }

            expandVec.Set(expandVec.X / ((Type)(2)), expandVec.Y / ((Type)(2)));
            m_max += expandVec;
            m_min -= expandVec;

            return *this;
        }

        /*!
        *@brief  当前包围盒扩大使其包含传入点
        *@param[in] point  传入点
        */
        void  MergeBox(const CVec2<Type> & point)    //移除inline关键字（函数过长，没有inline必要）
        {
            if(point[0] < m_min[0]) 
                m_min[0] = point[0];
            if(point[1] < m_min[1]) 
                m_min[1] = point[1];

            if(point[0] > m_max[0])
                m_max[0] = point[0];
            if(point[1] > m_max[1])
                m_max[1] = point[1];
        }

        /*!
        *@brief  当前包围盒和传入包围盒合并，并设置当前包围盒为合并的结果
        *@param[in] box   传入包围盒
        */
        void MergeBox(const CBox2<Type> & box)    //移除inline关键字（函数过长，没有inline必要）
        {
            if(box.m_min[0] < m_min[0]) 
                m_min[0] = box.m_min[0];
            if(box.m_min[1] < m_min[1]) 
                m_min[1] = box.m_min[1];

            if (box.m_max[0] > m_max[0])
                m_max[0] = box.m_max[0];
            if (box.m_max[1] > m_max[1])
                m_max[1] = box.m_max[1];
        }

        /*!
        *@brief  在给定误差内，判断当前包围盒是否包含传入点（在包围盒内或在包围盒上）
        *@param[in] point 传入点
        *@param[in] tolerance            给定误差
        *@return    是否包含传入点
        * - true    包含
        * - false   不包含
        */
        inline bool Contains(const CVec2<Type> & point,Type tolerance = std::numeric_limits<Type>::epsilon()) const
        {
            assert(NotEmpty());
            return ( point.X >= MinPt().X - tolerance 
                  && point.X <= MaxPt().X + tolerance 
                  && point.Y >= MinPt().Y - tolerance 
                  && point.Y <= MaxPt().Y + tolerance );
        }

        /*!
        *@brief  判断点是否在包围盒内
        *@param[in] point 传入点
        *@param[in] tolerance            给定误差
        *@return    点是否在包围盒内
        * - true    在
        * - false   不在
        */
        inline bool IsIn(const CVec2<Type> & point, Type tolerance = std::numeric_limits<Type>::epsilon()) const
        {
            assert(NotEmpty());
            return (point.X >= MinPt().X + tolerance
                 && point.X <= MaxPt().X - tolerance
                 && point.Y >= MinPt().Y + tolerance
                 && point.Y <= MaxPt().Y - tolerance);
        }

        /*!
        *@brief  判断点是否在包围盒上
        *@param[in] point 传入点
        *@param[in] tolerance            给定误差
        *@return    点是否在包围盒上
        * - true    在
        * - false   不在
        */
        inline bool IsOn(const CVec2<Type> & point, Type tolerance = std::numeric_limits<Type>::epsilon()) const
        {
            assert(NotEmpty());
            return Contains(point, tolerance) && !IsIn(point, tolerance);
        }

        /*!
        * @brief                              若点在矩形域外, 则根据参考点将其拉回到矩形域的边界
        * @param[in] point      点
        * @param[in] lastPoint  矩形域内的参考点
        * @return  是否拉回
        * - true  点被拉回到矩形域的边界
        * - false  点在矩形域内, 未被拉回  
        */
        bool Clamp(CVec2<Type> & point, const CVec2<Type> & lastPoint)
        {
            bool clamped = false;
            if (!Contains(point))
            {
                clamped = true;

                CVec2<Type> pts[4];
                pts[0] = MinPt();
                pts[1] = CVec2<Type>(MaxPt().X, MinPt().Y);
                pts[2] = MaxPt();
                pts[3] = CVec2<Type>(MinPt().X, MaxPt().Y);

                for (int i = 0; i < 4; ++i)
                {
                    if ((pts[i] - lastPoint).Cross(point - lastPoint) >0 
                        && (point - lastPoint).Cross(pts[(i + 1) % 4] - lastPoint) > 0 )
                    { 
                        if (i % 2 == 0)
                        {
                            point.X = point.X + (pts[i].Y - point.Y) * (lastPoint.X - point.X) / (lastPoint.Y - point.Y);
                            point.Y = pts[i].Y;
                        }
                        else
                        {
                            point.Y = point.Y + (pts[i].X - point.X) * (lastPoint.Y - point.Y) / (lastPoint.X - point.X);
                            point.X = pts[i].X;
                        }
                    }       
                }
            }

            return clamped;
        }

//         CVec2<Type> Clamp(const CVec2<Type> & point, bool & clamped)
//         {
//             bool beyondIndex = false;
//             if (value < Min)
//             {
//                 value = Min;
//                 beyondIndex = true;
//             }
//             if (value > Max)
//             {
//                 value = Max;
//                 beyondIndex = true;
//             }
// 
//             return beyondIndex;
//         }

        /*!
        *@brief  在给定误差内，判断传入包围盒是否和当前包围盒相等
        *@param[in] box  传入包围盒
        *@param[in] typeTol          给定误差
        *@return    是否相等
        * - true  相等 
        * - false  不相等
        */
        bool IsEqual(const CBox2<Type> & box, const Type typeTol = std::numeric_limits<Type>::epsilon()) const
        {
            return m_max.IsEqual(box.MaxPt(), typeTol) && m_min.IsEqual(box.MinPt(), typeTol);
        }

        /*!
        *@brief  在给定误差内，判断传入包围盒是否和当前包围盒相交
        *@param[in] box  传入包围盒
        *@param[in] tolerance           给定误差
        *@return    是否相交
        * - true  相交 
        * - false  不相交
        */
        inline bool IsIntersect(const CBox2<Type> & box,Type tolerance = std::numeric_limits<Type>::epsilon()) const
        {
            //assert(NotEmpty());
            //assert(box.NotEmpty());
            return (m_max.X >= box.MinPt().X - tolerance
                 && m_min.X <= box.MaxPt().X + tolerance
                 && m_max.Y >= box.MinPt().Y - tolerance 
                 && m_min.Y <= box.MaxPt().Y + tolerance
                 && NotEmpty()    //有空盒，必然不是相交的
                 && box.NotEmpty());
        }

        /*!
        *@brief  获得当前包围盒和传入包围盒的交集
        *@param[in] Box       传入包围盒
        *@param[out] IntersectBox    交集
        *@param[in] tolerance                给定误差
        *@return    是否存在交集
        * - true 存在相交包围盒
        * - false  不相交
        */
        bool GetIntersectBox(const CBox2<Type> & Box, CBox2<Type> & IntersectBox, Type tolerance  = std::numeric_limits<Type>::epsilon())    //移除inline关键字（函数过长，没有inline必要）
         {
            IntersectBox.MinPt().X = Max(Box.MinPt().X,MinPt().X);
            IntersectBox.MaxPt().X = Min(Box.MaxPt().X,MaxPt().X);
            if (IntersectBox.MaxPt().X < IntersectBox.MinPt().X - tolerance) 
            {
                return false;
            }

            IntersectBox.MinPt().Y = Max(Box.MinPt().Y, MinPt().Y);
            IntersectBox.MaxPt().Y = Min(Box.MaxPt().Y, MaxPt().Y);
            if (IntersectBox.MaxPt().Y < IntersectBox.MinPt().Y - tolerance) 
            {
                return false;
            }

            return true;
        }

        
        /*!
        *@brief  在给定误差下，包围盒是否在传入包围内，包括内相切
        *@param[in] box   传入包围盒
        *@param[in] tolerance            给定误差
        *@return   包围盒是否在传入包围内
        * - true  当前包围盒在传入包围盒内
        * - false  当前包围盒不在传入包围盒内
        */
        inline bool In(const CBox2<Type> & box, Type tolerance = std::numeric_limits<Type>::epsilon()) const
        {   
            assert(NotEmpty());
            assert(box.NotEmpty());
            return MinPt().X >= box.MinPt().X - tolerance 
                && MaxPt().X <= box.MaxPt().X + tolerance 
                && MinPt().Y >= box.MinPt().Y - tolerance 
                && MaxPt().Y <= box.MaxPt().Y + tolerance;
        }


        /*!
        *@brief  计算点到包围盒所有内点的最短距离，如果在包围盒内，距离为0
        *@param[in] pt   给定点
        *@return     最短距离
        */
        inline Type Distance(const CVec2<Type>& pt) const
        {
            assert(NotEmpty());
            CVec2<Type> vMindist;
            MindistVector(pt, vMindist);
            return vMindist.Length();
        }


        /*!
        *@brief  计算点到包围盒所有内点的最短距离平方，如果在包围盒内为0
        *@param[in] pt   给定点
        *@return   最短距离的平方
        */
        inline Type Distance2( const CVec2<Type> & pt) const
        {
            assert(NotEmpty());
            CVec2<Type> vMindist;
            MindistVector(pt, vMindist);
            return vMindist.SqrLength();
        }

        /*!
        *@brief  计算点到包围盒所有内点的最长距离
        *@param[in] pt   给定点
        *@return     最长距离
        */
        Type MaxDistance(const CVec2<Type>& pt) const
        {
            CVec2<Type> v_max;
            MaxdistVector(pt, v_max);
            return v_max.Length();
        }

        /*!
        *@brief  计算点到包围盒所有内点的最长距离平方
        *@param[in] pt   给定点
        *@return   最长距离的平方
        */
        Type MaxDistance2(const CVec2<Type>& pt) const
        {
            CVec2<Type> v_max;
            MaxdistVector(pt, v_max);
            return v_max.SqrLength();
        }

        /*!
         *@brief  获得包围盒的4个角点
         *@param[out] corner  返回角点数组
         *@param[in] bSortByMinMax 是否按照XY的min到max的顺序，从00~11进行排列，若为否则逆时针排列
         */
         void GetCorners(CVec2<Type> corner[4], bool bSortByMinMax = false) const 
         {
             if (bSortByMinMax)
             {
                 corner[0].Set(m_min.X, m_min.Y);
                 corner[1].Set(m_max.X, m_min.Y);
                 corner[2].Set(m_min.X, m_max.Y);
                 corner[3].Set(m_max.X, m_max.Y);
             }
             else
             {
                 corner[0].Set(m_max.X, m_max.Y);
                 corner[1].Set(m_min.X, m_max.Y);
                 corner[2].Set(m_min.X, m_min.Y);
                 corner[3].Set(m_max.X, m_min.Y);
             }
         }

        /*!
        *@brief  包围盒对象按照规定格式字符串输出
        *@return   字符串
        */
        char* AsString()
        {
            char * a = AsStringHelper::NewString("Box2d=[(%.20le %.20le),(%.20le %.20le)]", 
                m_min.X, m_min.Y, m_max.X, m_max.Y);
            return a;

        }
        /*!
        *@brief     根据字符串生成包围盒
        *@param[in] str   传入的字符串
        *@return     生成的包围盒
        */
        static CBox2<Type>* LoadFromString(char* str)
        {
            CBox2<Type>* pBox = new CBox2<Type>;
            g_sscanf(str, "Box2d=[(%le %le),(%le %le)]", 
                &pBox->MinPt().X, &pBox->MinPt().Y, &pBox->MaxPt().X, &pBox->MaxPt().Y);

            return pBox;
        }


        /*!
        *@brief     赋值操作符
        *@param[in] box 传入的包围盒
        *@return 当前包围盒
        */
        CBox2<Type> & operator = (const CBox2<Type> & box)
        {
            Set(box.MinPt(), box.MaxPt());

            return *this;
        }

         /*!
        *@brief     逻辑等号操作符，判断两个包围盒的是否相等
        *@param[in] b1   传入的包围盒1
        *@param[in] b2   传入的包围盒2
        *@return    是否相等
        * - true  相等
        * - false  不相等
        */
        friend bool operator ==(const CBox2<Type> & b1, const CBox2<Type> & b2)
        { 
            assert(b1.NotEmpty());
            assert(b2.NotEmpty());
            return b1.MinPt() == b2.MinPt() && b1.MaxPt() == b2.MaxPt();
        }

        /*!
        *@brief     逻辑不等号操作符，判断两个包围盒是否不相等
        *@param[in] b1   传入的包围盒1
        *@param[in] b2   传入的包围盒2
        *@return    是否不相等
        * - true  不相等
        * - false  相等
        */
        friend bool operator !=(const CBox2<Type> & b1, const CBox2<Type> & b2)
        { 
            assert(b1.NotEmpty());
            assert(b2.NotEmpty());
            return !(b1 == b2); 
        }

        /*!
        *@brief     -= 操作符，用于包围盒平移
        *@param[in] v   移动的向量
        *@return    *this
        */
        CBox2<Type> & operator -= (const CVec2<Type> & v)
        {
            m_min -= v;
            m_max -= v;
            return *this;
        }


        /*!
        *@brief     += 操作符，用于包围盒平移
        *@param[in] v   移动的向量
        *@return    *this
        */
        CBox2<Type> & operator += (const CVec2<Type> & v)
        {
            m_min += v;
            m_max += v;
            return *this;
        }

        CBox2<Type> & operator *= (Type scale)
        {
            if(scale > 0)
            {
                m_min *= scale;
                m_max *= scale;
            }
            else
            {
                CVec2<Type> temp = m_min; 
                m_min = scale*m_max;
                m_max = scale*temp;
            }

            return *this;
        }


    private:
        //计算一个点到包围盒的最小距离点距离等长的向量
        void MindistVector(const CVec2<Type> & pt, CVec2<Type> & vMinDist) const    //移除inline关键字（函数过长，没有inline必要）
        {
            CVec2<Type> center = CenterPt(); //box中心
            CVec2<Type> half = m_max - center;

            vMinDist = pt - center;

            //将v向量映射到右上角
            vMinDist.X = std::abs(vMinDist.X);
            vMinDist.Y = std::abs(vMinDist.Y);

            vMinDist -= half;

            if(vMinDist.X < (Type)0.0) 
                vMinDist.X = (Type)0.0;

            if(vMinDist.Y < (Type)0.0)
                vMinDist.Y= (Type)0.0;
        }

        inline void MaxdistVector(const CVec2<Type> & pt, CVec2<Type> & v_maxdist) const
        {
            CVec2<Type> center = CenterPt(); //box中心
            CVec2<Type> half = m_max - center;

            v_maxdist = pt - center;

            //将v向量映射到右上角
            v_maxdist.X = abs(v_maxdist.X);
            v_maxdist.Y = abs(v_maxdist.Y);

            v_maxdist += half;
        }

    private:
        CVec2<Type> m_min;
        CVec2<Type> m_max;

    public:
        static const CBox2<Type> Zero;
    };

    template<class Type> const CBox2<Type> CBox2<Type>::Zero(0, 0, 0, 0);

    typedef CBox2<int>    CBox2i;
    typedef CBox2<float>  CBox2f;
    typedef CBox2<double> CBox2d;

    /*! @} */
} 

#endif
