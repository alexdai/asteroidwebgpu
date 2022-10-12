//GMath.h
#ifndef GMATH_H
#define GMATH_H

#if defined(_WIN32) || defined(_MSC_VER)
#if !defined(GGP_STATIC)
    #ifdef GMATH_EXPORTS
        #define SYS_GMATH_API __declspec(dllexport)
    #else
        #define SYS_GMATH_API __declspec(dllimport)
    #endif
#else
    #define SYS_GMATH_API
#endif


#else
    #if !defined(GGP_STATIC)
        #ifdef GMATH_EXPORTS
            #define SYS_GMATH_API __attribute__((visibility("default")))
        #else
            #define SYS_GMATH_API
        #endif
    #else
        #define SYS_GMATH_API
    #endif
#endif

#endif
