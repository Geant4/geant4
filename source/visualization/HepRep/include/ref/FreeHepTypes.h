#ifndef FREEHEPTYPES_H
#define FREEHEPTYPES_H

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

#endif // FREEHEPTYPES_H

