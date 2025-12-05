// int *      signed, signed int   System dependent
// unsigned int      *      unsigned      System dependent
//__int8      1      char, signed char    -128 to 127
//__int16     2      short, short int, signed short int -32,768 to 32,767
//__int32     4      signed, signed int   -2,147,483,648 to 2,147,483,647
//__int64     8      none   -9,223,372,036,854,775,808 to 9,223,372,036,854,775,807
//char 1      signed char   -128 to 127
//unsigned char      1      none   0 to 255
//short       2      short int, signed short int -32,768 to 32,767
//unsigned short     2      unsigned short int   0 to 65,535
//long 4      long int, signed long int   -2,147,483,648 to 2,147,483,647
//unsigned long      4      unsigned long int    0 to 4,294,967,295
//enum *      none   Same as int
//float       4      none   3.4E +/- 38 (7 digits)
//double      8      none   1.7E +/- 308 (15 digits)
//long double 10     none   1.2E +/- 4932 (19 digits)

#ifndef IAEA_CONFIG
#define IAEA_CONFIG

#if (defined WIN32) || (defined WIN64)
#include <windows.h>
#endif
/* Without the above include file, gcc on Windows does not know about
   __int64
 */

#ifdef DOUBLE
typedef double IAEA_Float;
#else
typedef float  IAEA_Float;
#endif

typedef short IAEA_I16;
// typedef long  IAEA_I32; // RCN changed int to long to allow storage of EGS LATCH, Dec. 2006
typedef  int  IAEA_I32;    // Changed back on April 2011, following Daniel OBrien's comments
                           // It also corresponds to EGSnrc definition (see egs_config1.h file)
//typedef __int64 IAEA_I64;
#if (defined WIN32) || (defined WIN64)
typedef __int64 IAEA_I64;
#else
#if defined NO_LONG_LONG || defined LONG_IS_64
typedef long IAEA_I64;
#else
typedef long long IAEA_I64;
#endif
#endif

#ifdef __cplusplus
#define IAEA_EXTERN_C extern "C"
#else
#define IAEA_EXTERN_C extern
#endif

#if (defined WIN32) || (defined WIN64)

#ifdef BUILD_DLL
#define IAEA_EXPORT __declspec(dllexport)
#elif defined USE_DLL
#define IAEA_EXPORT __declspec(dllimport)
#else
#define IAEA_EXPORT
#endif
#define IAEA_LOCAL

#else

#ifdef HAVE_VISIBILITY
#define IAEA_EXPORT __attribute__ ((visibility ("default")))
#define IAEA_LOCAL  __attribute__ ((visibility ("hidden")))
#else
#define IAEA_EXPORT
#define IAEA_LOCAL
#endif

#endif

#endif


