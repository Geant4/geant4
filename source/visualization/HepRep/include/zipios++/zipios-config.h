// MD: defined for our purposes

/* Define to empty if the keyword does not work.  */
/* #undef const */

/* Define as __inline if that's what the C compiler calls it.  */
/* #undef inline */

/* Define if you have the ANSI C header files.  */
#define STDC_HEADERS 1

/* Define if the std compliant iostream library should be used (if present) */
#define USE_STD_IOSTREAM 1

/* Define if zlib has zError */
#define HAVE_ZERROR 1

/* Define if you have the <unistd.h> header file.  */
//#define HAVE_UNISTD_H 1

/* Name of package */
#define PACKAGE "zipios++"

/* Version number of package */
#define VERSION "0.1.5"

/* define if the compiler implements namespaces */
#define HAVE_NAMESPACES 1

/* define if the compiler supports Standard Template Library */
#define HAVE_STL 1

/* define if the compiler supports ISO C++ standard library */
#define HAVE_STD 1

/* define if the compiler has std compliant iostream library */
// MD defined for Geant4, Solaris and Linux 2.95.2 do not have this...
#if defined (WIN32) || __GNUC__ > 2
#define HAVE_STD_IOSTREAM 1
#endif

//#undef  S_ISSOCK
//#define S_ISSOCK(mode)	0

#ifdef WIN32
#ifndef GNU_GCC
#  ifdef _MSC_VER
// Disable class-browser warning about truncated template-names
#    pragma warning( disable : 4786 )
#  endif //_MSC_VER

// Needed for FilePath
#   define S_ISREG(mode)	(((mode) & _S_IFREG) == _S_IFREG)
#   define S_ISDIR(mode)	(((mode) & _S_IFDIR) == _S_IFDIR)
#   define S_ISCHR(mode)	(((mode) & _S_IFCHR) == _S_IFCHR)
#   define S_ISBLK(mode)	0
#   define S_ISSOCK(mode)	0
#   define S_ISFIFO(mode)	(((mode) & _S_IFIFO) == _S_IFIFO)
#endif
#endif

#include <cassert>
