/* config/scl_cf.h.  Generated automatically by configure.  */
/*
 * Check for headers
 */
#ifndef __SCL_CF_H__
#define __SCL_CF_H__

/* #undef DIRENT */
/* #undef SYSNDIR */
/* #undef SYSDIR */
/* #undef NDIR */

#define HAVE_DIRENT_H 1
/* #undef HAVE_SYS_NDIR_H */
/* #undef HAVE_SYS_DIR_H */
/* #undef HAVE_NDIR_H */
#define HAVE_STDARG_H 1

#define HAVE_STRING_H 1
#define HAVE_MEMORY_H 1

#ifndef WIN32
#define HAVE_UNISTD_H 1
#endif
/* #undef HAVE_SYSENT_H */
#define HAVE_SYS_STAT_H 1
/* #undef HAVE_STAT_H */

/* #undef NO_STDLIB_H */
/* #undef HAVE_VARARGS_H */
/* #undef HAVE_STROPTS_H */
/* #undef HAVE_SYSCONF_H */

/*
 * This section is for compile macros needed by
 * everything else.
 */

/*
 * Check for functions
 */
#define HAVE_MEMCPY 1
#define HAVE_SYSCONF 1
#define HAVE_MEMMOVE 1
#define HAVE_STRCHR 1
#define HAVE_ABS 1

#ifndef HAVE_STRCHR
#define strchr(s,c) index(s,c)
#endif /* HAVE_STRCHR */

/*
 * Check for tty/pty functions and structures
 */
/* #undef POSIX */

/*
 * Special hacks
 */
/* #undef CONVEX */
/* #undef SOLARIS */

#ifdef SOLARIS
#define __EXTENSIONS__
#endif /* SOLARIS */

#endif	/* __SCL_CF_H__ */
