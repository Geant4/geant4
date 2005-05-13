#ifndef HEPREP_TYPES_H
#define HEPREP_TYPES_H 1

namespace HEPREP {

#if defined(WIN32) && !defined(GNU_GCC)

// WIN32 and NOT GNU_GCC
typedef long int32;
typedef __int64 int64;
#define HEPREP_INT32_FORMAT "%d"
#define HEPREP_INT64_FORMAT "%ld"

#else // other than WIN32-MSVC
#if defined(_LP64) 

// 64 Bit Platforms
typedef int int32;
typedef long int64;
#define HEPREP_INT32_FORMAT "%d"
#define HEPREP_INT64_FORMAT "%ld"

#else 

// 32-Bit Platforms
typedef long int32;
typedef long long int64;
#define HEPREP_INT32_FORMAT "%ld"
#define HEPREP_INT64_FORMAT "%lld"

#endif // 32-Bit Platforms
#endif // other than WIN32-MSVC

} // namespace HEPREP

#endif
