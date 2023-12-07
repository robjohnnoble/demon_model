#ifndef MACROS_HPP
#define MACROS_HPP

#include <cmath>

template<typename T1, typename T2>
inline typename std::common_type<T1, T2>::type min(const T1& x, const T2& y) {
    return (x < y) ? x : y;
}
template<typename T1, typename T2>
inline typename std::common_type<T1, T2>::type max(const T1& x, const T2& y) {
    return (x > y) ? x : y;
}
template<typename T>
inline int sign(const T& x) {
    return (T(0) < x) - (x < T(0));
}
inline double roundit(double x) {
    return std::floor(x + 0.5);
}

#endif // MACROS_HPP