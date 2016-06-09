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
#ifndef tpi_IDs_h_included
#define tpi_IDs_h_included

#if defined __cplusplus
    extern "C" {
#endif

#include <statusMessageReporting.h>
#include <gString.h>
#include <xData.h>

#if defined __cplusplus
    namespace GIDI {
#endif

typedef struct tpi_channelID_s tpi_channelID;
typedef struct tpi_spectralID_s tpi_spectralID;
typedef struct tpi_spectralIDQualifier_s tpi_spectralIDQualifier;

struct tpi_spectralIDQualifier_s {
    tpi_spectralIDQualifier *next;
    char *name;
    char *value;
};

struct tpi_spectralID_s {
    tpi_spectralID *next;
    char *name;
    tpi_spectralIDQualifier *qualifiers;
    char *revision;                         /* Points to revision value in qualifiers if present, otherwise NULL. */
    char *multiplicityStr;                  /* Points to multiplicity value in qualifiers if present, otherwise NULL. */
    int multiplicity;                       /* Value of multiplicityStr if it is an integer. */
    double level;                           /* Value of levelStr if it is a number. */
    char *levelStr;                         /* Points to level value in qualifiers if present, otherwise NULL. */
    tpi_channelID *decayChannel;
};

struct tpi_channelID_s {
    tpi_spectralID *spectralIDs;
};

tpi_channelID *tpi_channelID_allocate( statusMessageReporting *smr );
int tpi_channelID_initialize( statusMessageReporting *smr, tpi_channelID *channelID );
void *tpi_channelID_free( statusMessageReporting *smr, tpi_channelID *channelID );
int tpi_channelID_release( statusMessageReporting *smr, tpi_channelID *channelID );
tpi_channelID *tpi_channelID_parse( statusMessageReporting *smr, char const *str, char **EOP );
int tpi_channelID_toString( statusMessageReporting *smr, tpi_channelID *channelID, gString *gStr );
int tpi_channelID_toStringSanRevision( statusMessageReporting *smr, tpi_channelID *channelID, gString *gStr );
int tpi_channelID_toStringSans( statusMessageReporting *smr, tpi_channelID *channelID, gString *gStr, char *sans[] );

tpi_spectralID *tpi_spectralID_allocate( statusMessageReporting *smr );
int tpi_spectralID_initialize( statusMessageReporting *smr, tpi_spectralID *spectralID );
void *tpi_spectralID_free( statusMessageReporting *smr, tpi_spectralID *spectralID );
int tpi_spectralID_release( statusMessageReporting *smr, tpi_spectralID *spectralID );
tpi_spectralID *tpi_spectralID_parse( statusMessageReporting *smr, char const *str, char **EOP );
int tpi_spectralID_toString( statusMessageReporting *smr, tpi_spectralID *spectralID, gString *gStr );
int tpi_spectralID_toStringSanRevision( statusMessageReporting *smr, tpi_spectralID *spectralID, gString *gStr );
int tpi_spectralID_toStringSans( statusMessageReporting *smr, tpi_spectralID *spectralID, gString *gStr, char *sans[] );

char const *tpi_misc_firstNonWhiteSpace( statusMessageReporting *smr, char const *str );
char *tpi_misc_allocateAndCopyNCharacters( statusMessageReporting *smr, char const *str, int n );

#if defined __cplusplus
    }
    }
#endif

#endif              /* End of tpi_IDs_h_included. */
