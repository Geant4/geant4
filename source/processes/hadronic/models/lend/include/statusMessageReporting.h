/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef statusMessageReporting_h_included
#define statusMessageReporting_h_included

#include <stdio.h>
#include <stdarg.h>

#ifdef _WIN32
#define __func__ __FUNCTION__
#endif

#if defined __cplusplus
    extern "C" {
#endif

#define smr_unknownID 0
#define smr_tooManyIDs 1
#define smr_invalidID 2
#define smr_errnoID 3
#define smr_smrID 4

#define smr_maximumNumberOfRegisteredLibraries 128
#define smr_maximumFileNameSize 1024
#define smr_codeNULL 0
#define smr_codeMemoryAllocating 1
enum smr_status { smr_status_Ok = 0, smr_status_Info, smr_status_Warning, smr_status_Error };
typedef char *(*smr_userInterface)( void *userData );

typedef struct statusMessageReport {
    struct statusMessageReport *next;
    enum smr_status status;
    int libraryID;
    int code;
    int line;
    char fileName[smr_maximumFileNameSize+1];           /* Do not free this. */
    char function[smr_maximumFileNameSize+1];           /* Do not free this. */
    char *message;                                      /* User must free this when done. Should use smr_release. */
} statusMessageReport;

typedef struct statusMessageReporting {
    enum smr_status verbosity;
    statusMessageReport report;
} statusMessageReporting;

int smr_setup( void );
int smr_cleanup( void );

int smr_registerLibrary( char const *libraryName );
int smr_numberOfRegisteredLibraries( void );
char const *smr_getRegisteredLibrarysName( int ID );

statusMessageReporting *smr_new( statusMessageReporting *smr, enum smr_status verbosity );
int smr_initialize( statusMessageReporting *smr, enum smr_status verbosity );
statusMessageReporting *smr_clone( statusMessageReporting const *smr );
void smr_release( statusMessageReporting *smr );
void *smr_free( statusMessageReporting **smr );

int smr_setReportInfo(  statusMessageReporting *smr, void *userInterface, char const *file, int line, char const *function, int libraryID, int code, char const *fmt, ... );
int smr_vsetReportInfo( statusMessageReporting *smr, void *userInterface, char const *file, int line, char const *function, int libraryID, int code, char const *fmt, va_list *args );
int smr_setReportWarning(  statusMessageReporting *smr, void *userInterface, char const *file, int line, char const *function, int libraryID, int code, char const *fmt, ... );
int smr_vsetReportWarning( statusMessageReporting *smr, void *userInterface, char const *file, int line, char const *function, int libraryID, int code, char const *fmt, va_list *args );
int smr_setReportError(  statusMessageReporting *smr, void *userInterface, char const *file, int line, char const *function, int libraryID, int code, char const *fmt, ... );
int smr_vsetReportError( statusMessageReporting *smr, void *userInterface, char const *file, int line, char const *function, int libraryID, int code, char const *fmt, va_list *args );

enum smr_status smr_highestStatus( statusMessageReporting const *smr );
int smr_isOk( statusMessageReporting const *smr );
int smr_isInfo( statusMessageReporting const *smr );
int smr_isWarning( statusMessageReporting const *smr );
int smr_isError( statusMessageReporting const *smr );
int smr_isWarningOrError( statusMessageReporting const *smr );

int smr_isReportOk( statusMessageReport const *report );
int smr_isReportInfo( statusMessageReport const *report );
int smr_isReportWarning( statusMessageReport const *report );
int smr_isReportError( statusMessageReport const *report );
int smr_isReportWarningOrError( statusMessageReport const *report );

int smr_numberOfReports( statusMessageReporting const *smr );
statusMessageReport const *smr_firstReport( statusMessageReporting const *smr );
statusMessageReport const *smr_nextReport( statusMessageReport const *report );

enum smr_status smr_getVerbosity( statusMessageReporting const *smr );

int smr_getLibraryID( statusMessageReport const *report );
int smr_getCode( statusMessageReport const *report );
int smr_getLine( statusMessageReport const *report );
char const *smr_getFile( statusMessageReport const *report );
char const *smr_getFunction( statusMessageReport const *report );
char const *smr_getMessage( statusMessageReport const *report );
char *smr_copyMessage( statusMessageReport const *report );
char *smr_copyFullMessage( statusMessageReport const *report );
void smr_print( statusMessageReporting *smr, int clear );
void smr_write( statusMessageReporting *smr, FILE *f, int clear );
void smr_reportPrint( statusMessageReport const *report );
void smr_reportWrite( statusMessageReport const *report, FILE *f );

char const *smr_statusToString( enum smr_status status );

char *smr_allocateFormatMessage( char const *fmt, ... );
char *smr_vallocateFormatMessage( char const *fmt, va_list *args );

void *smr_malloc( statusMessageReporting *smr, size_t size, int zero, char const *forItem, char const *file, int line, char const *function );
void *smr_realloc( statusMessageReporting *smr, void *pOld, size_t size, char const *forItem, char const *file, int line, char const *function );
void *smr_freeMemory( void **p );
char *smr_allocateCopyString( statusMessageReporting *smr, char const *s, char const *forItem, char const *file, int line, char const *function );
char *smr_allocateCopyStringN( statusMessageReporting *smr, char const *s, size_t n, char const *forItem, char const *file, int line, char const *function );

#define smr_malloc2( smr, size, zero, forItem ) smr_malloc( smr, size, zero, forItem, __FILE__, __LINE__, __func__ )
#define smr_realloc2( smr, old, size, forItem ) smr_realloc( smr, old, size, forItem, __FILE__, __LINE__, __func__ )
#define smr_freeMemory2( p ) smr_freeMemory( (void **) &p )
#define smr_allocateCopyString2( smr, s, forItem ) smr_allocateCopyString( smr, s, forItem, __FILE__, __LINE__, __func__ )
#define smr_allocateCopyStringN2( smr, s, n, forItem ) smr_allocateCopyStringN( smr, s, n, forItem, __FILE__, __LINE__, __func__ )

#define smr_setReportInfo2( smr, libraryID, code, fmt, ... ) smr_setReportInfo( smr, NULL, __FILE__, __LINE__, __func__, libraryID, code, fmt, __VA_ARGS__ )
#define smr_setReportInfo2p( smr, libraryID, code, fmt ) smr_setReportInfo( smr, NULL, __FILE__, __LINE__, __func__, libraryID, code, fmt )
#define smr_vsetReportInfo2( smr, libraryID, code, fmt, args ) smr_vsetReportInfo( smr, NULL, __FILE__, __LINE__, __func__, libraryID, code, fmt, args )
#define smr_setReportWarning2( smr, libraryID, code, fmt, ... ) smr_setReportWarning( smr, NULL, __FILE__, __LINE__, __func__, libraryID, code, fmt, __VA_ARGS__ )
#define smr_setReportWarning2p( smr, libraryID, code, fmt ) smr_setReportWarning( smr, NULL, __FILE__, __LINE__, __func__, libraryID, code, fmt )
#define smr_vsetReportWarning2( smr, libraryID, code, fmt, args ) smr_vsetReportWarning( smr, NULL, __FILE__, __LINE__, __func__, libraryID, code, fmt, args )
#define smr_setReportError2( smr, libraryID, code, fmt, ... ) smr_setReportError( smr, NULL, __FILE__, __LINE__, __func__, libraryID, code, fmt, __VA_ARGS__ )
#define smr_setReportError2p( smr, libraryID, code, fmt ) smr_setReportError( smr, NULL, __FILE__, __LINE__, __func__, libraryID, code, fmt )
#define smr_vsetReportError2( smr, libraryID, code, fmt, args ) smr_vsetReportError( smr, NULL, __FILE__, __LINE__, __func__, libraryID, code, fmt, args )

#define smr_setReportInfo3( smr, userInterface, libraryID, code, fmt, ... ) smr_setReportInfo( smr, userInterface, __FILE__, __LINE__, __func__, libraryID, code, fmt, __VA_ARGS__ )
#define smr_setReportInfo3p( smr, userInterface, libraryID, code, fmt ) smr_setReportInfo( smr, userInterface, __FILE__, __LINE__, __func__, libraryID, code, fmt )
#define smr_vsetReportInfo3( smr, userInterface, libraryID, code, fmt, args ) smr_vsetReportInfo( smr, userInterface, __FILE__, __LINE__, __func__, libraryID, code, fmt, args )
#define smr_setReportWarning3( smr, userInterface, libraryID, code, fmt, ... ) smr_setReportWarning( smr, userInterface, __FILE__, __LINE__, __func__, libraryID, code, fmt, __VA_ARGS__ )
#define smr_setReportWarning3p( smr, userInterface, libraryID, code, fmt ) smr_setReportWarning( smr, userInterface, __FILE__, __LINE__, __func__, libraryID, code, fmt )
#define smr_vsetReportWarning3( smr, userInterface, libraryID, code, fmt, args ) smr_vsetReportWarning( smr, userInterface, __FILE__, __LINE__, __func__, libraryID, code, fmt, args )
#define smr_setReportError3( smr, userInterface, libraryID, code, fmt, ... ) smr_setReportError( smr, userInterface, __FILE__, __LINE__, __func__, libraryID, code, fmt, __VA_ARGS__ )
#define smr_setReportError3p( smr, userInterface, libraryID, code, fmt ) smr_setReportError( smr, userInterface, __FILE__, __LINE__, __func__, libraryID, code, fmt )
#define smr_vsetReportError3( smr, userInterface, libraryID, code, fmt, args ) smr_vsetReportError( smr, userInterface, __FILE__, __LINE__, __func__, libraryID, code, fmt, args )

#if defined __cplusplus
    }
#endif

#endif              /* End of statusMessageReporting_h_included. */
