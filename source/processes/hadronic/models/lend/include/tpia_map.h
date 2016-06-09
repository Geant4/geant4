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
#ifndef tpia_map_h_included
#define tpia_map_h_included

#if defined __cplusplus
    extern "C" {
#endif

#include <xData.h>

#if defined __cplusplus
    namespace GIDI {
#endif

enum tpia_map_status { tpia_map_status_Ok, tpia_map_status_memory, tpia_map_status_mapParsing,
    tpia_map_status_UnknownType };
enum tpia_mapEntry_type { tpia_mapEntry_type_target, tpia_mapEntry_type_path };

typedef struct tpia_map_s tpia_map;
typedef struct tpia_mapEntry_s tpia_mapEntry;
typedef struct tpia_map_smr_s tpia_map_smr;

struct tpia_map_smr_s {
    smr_userInterface smrUserInterface;
    tpia_map *map;
};

struct tpia_mapEntry_s {
    tpia_mapEntry *next;
    enum tpia_mapEntry_type type;
    char *schema;
    char *path;
    char *evaluation;
    char *projectile;
    char *targetName;
    tpia_map *map;
};

struct tpia_map_s {
    enum tpia_map_status status;
    tpia_map_smr smrUserInterface;
    char *path;
    char *mapFileName;
    int numberOfEntries;
    tpia_mapEntry *mapEntries;
};

tpia_map *tpia_map_create( statusMessageReporting *smr );
int tpia_map_initialize( statusMessageReporting *smr, tpia_map *map );
tpia_map *tpia_map_readFile( statusMessageReporting *smr, const char *basePath, const char *mapFileName );
void *tpia_map_free( statusMessageReporting *smr, tpia_map *map );
void tpia_map_release( statusMessageReporting *smr, tpia_map *map );
tpia_mapEntry *tpia_map_getFirstEntry( tpia_map *map );
tpia_mapEntry *tpia_map_getNextEntry( tpia_mapEntry *entry );
int tpia_map_addTarget( statusMessageReporting *smr, tpia_map *map, const char *method, const char *path, const char *evaluation, const char *projectile, const char *targetName );
int tpia_map_addPath( statusMessageReporting *smr, tpia_map *map, const char *path, const char *projectile );
char *tpia_map_findTarget( statusMessageReporting *smr, tpia_map *map, const char *evaluation, const char *projectile, const char *targetName );
tpia_map *tpia_map_findAllOfTarget( statusMessageReporting *smr, tpia_map *map, const char *projectile, const char *targetName );
char *tpia_map_getFullPath( statusMessageReporting *smr, tpia_map *map, const char *endPath );
int tpia_map_walkTree( statusMessageReporting *smr, tpia_map *map, int (*handler)( tpia_mapEntry *entry, int level, void *userData), void *userData );
char *tpia_map_toXMLString( statusMessageReporting *smr, tpia_map *map );
void tpia_map_simpleWrite( FILE *f, tpia_map *map );

#if defined __cplusplus
    }
    }
#endif

#endif          /* End of tpia_map_h_included. */
