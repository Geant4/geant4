// Copyright FreeHEP, 2005.
#ifndef CHEPREP_CONFIG_H
#define CHEPREP_CONFIG_H 1

/**
 * @author Mark Donszelmann
 * @version $Id: config.h,v 1.4 2005-06-02 21:28:45 duns Exp $
 */
namespace cheprep {

#if defined(WIN32) && !defined(GNU_GCC)

// WIN32 and NOT GNU_GCC
typedef __int64 int64;
typedef unsigned __int64 uint64;
#define CHEPREP_INT64_FORMAT "%ld"
#define CHEPREP_UINT64_FORMAT "%uld"

#else // other than WIN32-MSVC
#if defined(_LP64) 

// 64 Bit Platforms
typedef long int64;
typedef unsigned long uint64;
#define CHEPREP_INT64_FORMAT "%ld"
#define CHEPREP_UINT64_FORMAT "%uld"

#else 

// 32-Bit Platforms
typedef long long int64;
typedef unsigned long long uint64;
#define CHEPREP_INT64_FORMAT "%lld"
#define CHEPREP_UINT64_FORMAT "%ulld"

#endif // 32-Bit Platforms
#endif // other than WIN32-MSVC

} // namespace cheprep

#ifdef WIN32
#ifndef GNU_GCC
// Disable warning C4786: identifier was truncated to '255' characters in the debug information
    #pragma warning ( disable : 4786 )
// Disable warning C4250: inherits via dominance
    #pragma warning ( disable : 4250 )
#ifdef VC6
// FIX for KB 168440 - VC6
// Stream Operator << Cannot Handle __int64 Type
    #include<iostream>

    inline std::ostream& operator<<(std::ostream& os, __int64 i ) {
        char buf[20];
        sprintf(buf,"%I64d", i );
        os << buf;
        return os;
    }
#endif // VC6
#endif // GNU_GCC
#endif // WIN32


#endif  // CHEPREP_CONFIG_H
