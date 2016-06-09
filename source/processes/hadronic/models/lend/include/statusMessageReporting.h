/*
# <<BEGIN-copyright>>
# Copyright (c) 2010, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory 
# Written by Bret R. Beck, beck6@llnl.gov. 
# CODE-461393
# All rights reserved. 
#  
# This file is part of GIDI. For details, see nuclear.llnl.gov. 
# Please also read the "Additional BSD Notice" at nuclear.llnl.gov. 
# 
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met: 
#
#      1) Redistributions of source code must retain the above copyright notice, 
#         this list of conditions and the disclaimer below.
#      2) Redistributions in binary form must reproduce the above copyright notice, 
#         this list of conditions and the disclaimer (as noted below) in the 
#          documentation and/or other materials provided with the distribution.
#      3) Neither the name of the LLNS/LLNL nor the names of its contributors may be 
#         used to endorse or promote products derived from this software without 
#         specific prior written permission. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT 
# SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS 
# OR SERVICES;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
# AND ON  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
# <<END-copyright>>
*/
#ifndef statusMessageReporting_h_included
#define statusMessageReporting_h_included

#include <stdio.h>
#include <stdarg.h>

#if defined __cplusplus
    namespace GIDI {
    extern "C" {
#endif

#define smr_maximumPackageNameSize 256
enum smr_status { smr_status_Ok, smr_status_Info, smr_status_Error, smr_status_Fatal };
typedef int (*smr_userInterface)( void *userData, char **smr );

typedef struct statusMessageReporting {
    enum smr_status status;
    char packageName[smr_maximumPackageNameSize];       /* Do not free this. */
    int line;
    int code;
    char *message;                                      /* User must free this when done. Should use smr_release. */
} statusMessageReporting;

int smr_initialize( statusMessageReporting *smr );
int smr_release( statusMessageReporting *smr );
int smr_setMessageInfo(  statusMessageReporting *smr, void *userInterface, const char *file, int line, int code, const char *fmt, ... );
int smr_vsetMessageInfo( statusMessageReporting *smr, void *userInterface, const char *file, int line, int code, const char *fmt, va_list *args );
int smr_setMessageError( statusMessageReporting *smr, void *userInterface, const char *file, int line, int code, const char *fmt, ... );
int smr_vsetMessageError( statusMessageReporting *smr, void *userInterface, const char *file, int line, int code, const char *fmt, va_list *args );
char *smr_allocateFormatMessage( const char *fmt, ... );
char *smr_vallocateFormatMessage( const char *fmt, va_list *args );
int smr_isOk( statusMessageReporting *smr );
int smr_isInfo( statusMessageReporting *smr );
int smr_isError( statusMessageReporting *smr );
int smr_isFatal( statusMessageReporting *smr );
const char *smr_getMessage( statusMessageReporting *smr );
char *smr_getFullMessage( statusMessageReporting *smr );
void smr_print( statusMessageReporting *smr, FILE *f, int clear );

#if defined __cplusplus
    }
    }
#endif

#endif              /* End of statusMessageReporting_h_included. */
