/*!
* @file      interval.h
* @brief     区间的基本定义
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
#ifndef __interval_h__
#define __interval_h__

#include "GMathDef.h"
#include "GVec2.h"

#pragma warning( disable:4251)

    /*!\addtogroup GMath GMath
    * @{
    */

    /*!
    * @struct CInterval
    * @brief  区间定义
              该类给出区间的定义，提供相关的计算方法，被其它很多类使用
    */
    template<typename Type>
    class CInterval
    {
    public:
        /*!
        *@brief 默认构造函数
        */
        CInterval(): Min(0), Max(0){}

        /*!
        *@brief 通过给定区间端点值构造区间
        *@param[in] a 区间左端点
        *@param[in] b 区间右端点
        */
        CInterval(const Type &a, const Type &b): Min(a), Max(b) {}

        /*!
        *@brief 设置区间左端点
        *@param[in] a 左端点值
        */
        void SetMin(const Type &a)
        {
            Min = a;
        }

        /*!
        *@brief 设置区间右端点
        *@param[in] a 右端点值
        */
        void SetMax(const Type &a)
        {
            Max = a;
        }

        
        /*!
        * @brief      设置区间
        * @param[in] a 区间左端点
        * @param[in] b 区间右端点
        */
        void Set(const Type &a, const Type &b)
        {
            Min = a;
            Max = b;
        }

        
        /*!
        * @brief      区间分割
        * @param [in] a 分割点
        * @param [in] bLow  保留左半段还是右半段
        */
        void CutTo(Type a, bool bLow)
        {
            if(bLow)
            {
                Max = a;
            }
            else
            {
                Min = a;
            }
        }
        
        /*!
        * @brief      区间对原点的镜像
        */
        void Mirror()
        {
            double temp = Min;
            Min = -Max;
            Max = -temp;
        }
        
        /*!
        * @brief      初始化区间成无穷区间
        */
        void SetInfinite()
        {
            Min = g_NegInfinity;
            Max = g_Infinity;
        }

        /*!
        * @brief      初始化区间成空区间
        */

        void SetEmpty()
        {
            if(std::numeric_limits<Type>::has_infinity)
            {
                const Type & c_infinity = std::numeric_limits<Type>::infinity();
                Min = c_infinity;
                Max = -c_infinity;
            }
            else
            {
#undef min
#undef max
                const Type & c_min = std::numeric_limits<Type>::min();
                const Type & c_max = std::numeric_limits<Type>::max();
                Min = c_max;
                Max = c_min;
            }
        }

        /*!
        * @brief      判断区间是否为空
        *@return    是否为空
        * - true 是
        * - false 不是
        */

        bool IsEmpty() const
        {
            return (Min > Max);
        }

        /*!
        * @brief      区间正常化，如果下限大于上限，则交换它们
        */
        void Normalize()
        {
            if(Min > Max)
            {
                Swap();
            }
        }

        
        /*!
        * @brief      交换区间上下限
        */
        void Swap()
        {
            Type temp = Min;
            Min = Max;
            Max = temp;
        }

        
        /*!
        * @brief      区间膨胀
        * @param[in]  size 扩充大小
        */
        void Inflate(Type size)
        {
            Min -= size;
            Max += size;
        }

        
        /*!
        * @brief      合并区间，如果给定值a恰好在原来区间内，区间维持不变，如果不在，将区间扩充至恰好包含给定值a
        * @param[in]  a 给定数值
        */
        void Merge(Type a)
        {
            if(Min > a)
            {
                Min = a;
            }

            if(Max < a)
            {
                Max = a; 
            }
        }
        /*!
        * @brief      合并区间
        * @param[in]  other 给定区间
        */
        void Merge(const CInterval<Type>  & other)
        {
            Merge(other.Min);
            Merge(other.Max);
        }

        
        /*!
        * @brief      区间缩小，区间上下限都向区间中心缩小指定尺寸
        * @param [in] size 缩小尺度
        */
        void Deflate(Type size)
        {
            Min += size;
            Max -= size;
        }

        
        /*!
        * @brief      区间移动
        * @param [in] dmin 左端点移动尺度
        * @param [in] dmax 右端点移动尺度
        */
        void Shift(Type dmin, Type dmax)
        {
            Min += dmin;
            Max += dmax;
        }

        /*!
        * @brief      区间移动
        * @param [in] a 移动尺度
        */
        void Shift(Type a)
        {
            Min += a;
            Max += a;
        }

        /*!
        * @brief      区间左端点移动
        * @param [in] dmin 移动尺度
        */
        void ShiftMin(Type dmin)
        {
            Min += dmin;
        }

        /*!
        * @brief      区间右端点移动
        * @param [in] dmax 移动尺度
        */
        void ShiftMax(Type dmax)
        {
            Max += dmax;
        }

        
        /*!
        * @brief      移动上限使得区间为指定长度
        * @param [in] length 指定长度
        */
        void ShiftMaxToLength(Type length)
        {
            Max = Min + length;
        }

        
        /*!
        * @brief      移动下限使得区间为指定长度
        * @param [in] length 指定长度
        */
        void ShiftMinToLength(Type length)
        {
            Min = Max - length;
        }

        
        /*!
        * @brief      将整个区间压缩至下限点
        */
        void CompressToMin()
        {
            Max = Min;
        }

        
        /*!
        * @brief      将整个区间压缩至上限点
        */
        void CompressToMax()
        {
            Min = Max;
        }

        
        /*!
        * @brief      区间长度，上下限差的绝对值
        * @return     区间长度
        */
        Type Length() const
        {
            return abs(Max - Min);
        }

        
        /*!
        * @brief      区间尺寸，上限减去下限
        * @return     尺寸
        */
        Type Size() const
        {
            return Max - Min;
        }

        
        /*!
        * @brief      区间中点
        * @return     区间中点
        */
        Type Middle() const
        {
            return ((Type)0.5)*(Min + Max);
        }

        
        /*!
        * @brief      区间是否是非空区间
        *@return    是否为空区间
        * - true 是
        * - false 不是
        */
        bool NotEmpty() const
        {
            return Min <= Max;
        }

        
        /*!
        * @brief      两个区间的交集
        * @param [in] other 另一个区间
        * @return     交集
        */
        CInterval<Type> Intersection(const CInterval<Type>  & other) const
        {
            return CInterval<Type>(Min < other.Min ? other.Min : Min, //最小值中最大的
                                  Max < other.Max ? Max: other.Max); //最大值中最小的
        }

        
        /*!
        * @brief      两个非空区间是否相交
        * @param [in] other 另一个区间
        *@return    是否相交
        * - true 是
        * - false 不是
        */
        bool IsIntersect(const CInterval<Type>  & other) const
        {
            return Max >= other.Min && Min <= other.Max;
        }

        /*!
        * @brief      带容差的判断两个非空区间是否相交
        * @param [in] other 另一个区间
        * @param [in] epsilon 容差
        *@return    是否相交
        * - true 是
        * - false 不是
        */
        bool IsIntersect(const CInterval<Type>  & other, Type epsilon) const
        {
            return Max >= other.Min - epsilon && Min <= other.Max + epsilon;
        }

        /*!
        * @brief      带容差的判断两个区间是否相同
        * @param [in] other 另一个区间
        * @param [in] epsilon 容差
        *@return    是否相同
        * - true 是
        * - false 不是
        */
        bool Equals(const CInterval<Type>  & other, Type epsilon) const
        {
            return fabs(Min-other.Min) <= epsilon && fabs(Max-other.Max) <= epsilon;
        }

        /*!
        * @brief      值是否在区间内
        * @param [in] a 值
        *@return    是否在区间内
        * - true 是
        * - false 不是
        */
        bool Contain(Type a) const
        {
            assert(NotEmpty());
            return a >= Min && a <= Max;
        }

        /*!
        * @brief      带容差的判断值是否在区间内
        * @param [in] a 值
        * @param [in] epsilon 容差
        *@return    是否在区间内
        * - true 是
        * - false 不是
        */
        bool Contain(Type a, Type epsilon) const
        {
            assert(NotEmpty());
            return a >= Min - epsilon && a <= Max + epsilon;
        }

        /*!
        * @brief      区间是否在区间内
        * @param [in] other 区间
        * @return    是否在区间内
        * - true 是
        * - false 不是
        */
        bool Contain(const CInterval<Type>& other) const
        {
            assert(NotEmpty() && other.NotEmpty());
            return Contain(other.Min) && Contain(other.Max);
        }

        /*!
        * @brief      带容差的判断区间是否在区间内
        * @param [in] other 区间
        * @param [in] epsilon 容差
        *@return    是否在区间内
        * - true 是
        * - false 不是
        */
        bool Contain(const CInterval<Type>& other, Type epsilon) const
        {
            assert(NotEmpty() && other.NotEmpty());
            return Contain(other.Min, epsilon) && Contain(other.Max, epsilon);
        }

        /*!
        * @brief                     若数值在区间外, 将其拉回到临近的区间端点
        * @param[in]  value    数值
        * @return   是否被拉回
        * - true  数值被拉回到区间端点
        * - false  数值在区间内, 未被拉回  
        */
        bool Clamp(Type & value) const
        {
            bool clamped = false;
            if (value < Min)
            {
                value = Min;
                clamped = true;
            }
            if (value > Max)
            {
                value = Max;
                clamped = true;
            }

            return clamped;
        }

        /*!
        * @brief      判断值是否在区间的边界
        * @param [in] a 值
        * @param [in] epsilon 容差
        *@return    是否在区间边界
        * - true 是
        * - false 不是
        */
        bool IsBound(Type a, Type epsilon = 0) const
        {
            assert(NotEmpty());
            return fabs(a - Min) < epsilon || fabs(a - Max) < epsilon;
        }

        /*!
        * @brief      区间比例变换，上下限都乘以指定因子
        * @param [in] factor      指定因子
        * @return    不改变原区间，返回一个新区间
        */
        CInterval<Type> operator*(Type factor) const
        {
            CInterval<Type> r(factor*Min, factor*Max);
            r.Normalize();
            return r;
        }

        /*!
        * @brief      区间比例变换，上下限都乘以指定因子
        * @param [in] factor      指定因子
        * @param[in]  i 输入区间
        * @return    不改变原区间，返回一个新区间
        */
        friend CInterval<Type> operator*(Type factor, const CInterval<Type>& i)
        {
            CInterval<Type> r(factor*i.Min, factor*i.Max);
            r.Normalize();
            return r;
        }

        /*!
        *@brief 区间与数值相加
        *@param[in] value 数值
        *@return 结果区间
        */
        CInterval<Type> operator+(Type value) const
        {
            return CInterval<Type> (Min + value, Max + value);
        }

        /*!
        *@brief 区间与数值相减
        *@param[in] value 数值
        *@return 结果区间
        */
        CInterval<Type> operator-(Type value) const
        {
            return CInterval<Type> (Min - value, Max - value);
        }

        /*!
        *@brief 区间与数值相除
        *@param[in] value 数值
        *@return 结果区间
        */
        CInterval<Type> operator/(Type value) const
        {
            Type r1 = (Type)(1)/value;
            return operator*(r1);
        }

        /*!
        *@brief 区间乘法
        *@param[in] other 另一个区间
        *@return 运算结果
        */
        CInterval<Type> operator*(const CInterval<Type>& other) const
        {
            CInterval<Type> r;
            r.SetEmpty();
            r.Merge(Min*other.Min);
            r.Merge(Min*other.Max);
            r.Merge(Max*other.Min);
            r.Merge(Max*other.Max);

            return r;
        }

        /*!
        *@brief 区间加法
        *@param[in] other 另一个区间
        *@return 运算结果
        */
        CInterval<Type> operator+(const CInterval<Type>& other) const
        {
            return CInterval(Min + other.Min, Max + other.Max);
        }

        /*!
        *@brief 区间减法
        *@param[in] other 另一个区间
        *@return 运算结果
        */
        CInterval<Type> operator-(const CInterval<Type>& other) const
        {
            return CInterval(Min - other.Max, Max - other.Min);
        }

        /*!
        *@brief 区间除法
        *@param[in] other 另一个区间
        *@return 运算结果
        */
        CInterval<Type> operator/(const CInterval<Type>& other) const
        {
            CInterval<Type> r;
            r.SetEmpty();
            auto r1 = (Type)(1)/other.Min;
            auto r2 = (Type)(1)/other.Max;
            r.Merge(Min*r1);
            r.Merge(Min*r2);
            r.Merge(Max*r1);
            r.Merge(Max*r2);

            return r;
        }


        /*!
        * @brief      区间比例变换，this区间上下限都乘以指定因子
        * @param [in] factor      指定因子
        * @return     返回原区间的引用
        */
        CInterval <Type> & operator*=(Type factor)
        {
            Min *= factor;
            Max *= factor;
            return *this;
        }

        /*!
        * @brief      基于参考点对区间进行等比缩放
        * @param [in] factor    缩放因子
        * @param [in] ref     参考点       
        */
        void Scale(Type factor, Type ref = 0)
        {
            Min = factor * (Min - ref) + ref;
            Max = factor * (Max - ref) + ref;
        }

        /*!
        *@brief 判断区间是否相同
        *@param[in] v1 区间1
        *@param[in] v2 区间2
        *@return    是否相同
        * - true 是
        * - false 不是
        */
        friend bool operator ==(const CInterval<Type> & v1, const CInterval<Type> & v2)
        { 
            return v1.Min == v2.Min && v1.Max == v2.Max;
        }

        /*!
        *@brief 判断区间是否相同
        *@param[in] v1 区间1
        *@param[in] v2 区间2
        *@return    是否不相同
        * - true 是
        * - false 不是
        */
        friend bool operator !=(const CInterval<Type> & v1, const CInterval<Type> & v2)
        {
            return v1.Min != v2.Min || v1.Max != v2.Max;
        }
        
        /*!
        * @brief      区间集合减法，如果两个结果区间都不空，说明this区间包含other区间
        * @param [in] other  被减区间
        * @param [in] sub1         第一个结果区间
        * @param [in] sub2         第二个结果区间                          
        */
        void Substract(const CInterval<Type> & other, CInterval<Type> & sub1, CInterval<Type> & sub2) const
        {
            sub1.Set(Min, ggp::Min<Type>(Max, other.Min));
            sub2.Set(ggp::Max<Type>(Min, other.Max), Max);
        }

       
        /*!
        * @brief      转换区间成向量类型
        * @return     转换后的向量
        */
        CVec2<Type> & AsVector()
        {
            return *(CVec2<Type> *)(this);
        }


    public:
        Type Min, Max;
    };

    typedef CInterval<double> CIntervald;
    typedef CInterval<float>  CIntervalf;
    typedef CInterval<int>    CIntervali;

    //一个区间扣减多个区间
    template<typename Type, typename IterType, typename FunctorT>
    void IntervalSubstract(const CInterval<Type> & I, IterType begin, IterType end, vector<CInterval<Type>> & intervals, FunctorT getter)
    {
        const int kLeft = 0;
        const int kRight = 1;

        vector<std::pair<Type, int>> extreme_points;

        std::pair<Type, int> pt;

        for (auto it = begin; it != end; ++it)
        {
            auto & range = getter(it);
            extreme_points.push_back(std::make_pair(range.Min, kLeft));
            extreme_points.push_back(std::make_pair(range.Max, kRight));
        }

        std::sort(extreme_points.begin(), extreme_points.end(), 
        [&](std::pair<Type, int> & a, std::pair<Type, int> & b)->bool
        {
            if (a.first < b.first)
            {
                return true;
            }
            else if (a.first > b.first)
            {
                return false;
            }
            else if(a.second < b.second)
            {
                return true;
            }
            else
            {
                return false;
            }
        });

        int count = 0;

        CInterval<Type> range;
        range.Min = I.Min;



        for (auto it = extreme_points.begin(); it != extreme_points.end(); ++it)
        {
            int step;
            if (kLeft == it->second)
            {
                step = 1;
                ++count;
            }
            else
            {
                step = -1;
                --count;
            }

            if (count == 0)
            {
                range.Min = it->first;
            }
            else if ((count == 1) && (step == 1))
            {
                range.Max = it->first;
                if ((range.Min) < (range.Max))
                {
                    intervals.push_back(range);
                }
            }
        }

        range.Max = I.Max;

        if ((range.Min) < (range.Max))
        {
            intervals.push_back(range);
        }
    }

    /*! @} */

#endif  //__interval_h__
