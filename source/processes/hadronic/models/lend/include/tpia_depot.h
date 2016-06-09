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
#ifndef tpia_depot_h_included
#define tpia_depot_h_included

#if defined __cplusplus
    extern "C" {
#endif

#include <tpia_map.h>
#include <tpia_target.h>

#if defined __cplusplus
    namespace GIDI {
#endif

typedef struct tpia_targetEntry_s tpia_targetEntry;
typedef struct tpia_depot_s tpia_depot;

struct tpia_targetEntry_s {
    tpia_targetEntry *next;
    tpia_target *target;
};

struct tpia_depot_s {
    int status;
    char *projectileName;
    int numberOfTargets;
    tpia_targetEntry *targets;
    tpia_map *map;
};

tpia_depot *tpia_depot_create( statusMessageReporting *smr, const char *projectileName );
int tpia_depot_initialize( statusMessageReporting *smr, tpia_depot *depot, const char *projectileName );
tpia_depot *tpia_depot_free( tpia_depot *depot, int freeMap );
int tpia_depot_release( tpia_depot *depot, int freeMap );
int tpia_depot_setMap( statusMessageReporting *smr, tpia_depot *depot, tpia_map *map );
int tpia_depot_setMapFromFilename( statusMessageReporting *smr, tpia_depot *depot, const char *basePath, const char *mapFileName );
tpia_map *tpia_depot_releaseMap( tpia_depot *depot );
int tpia_depot_freeMap( tpia_depot *depot );
tpia_target *tpia_depot_addTarget( statusMessageReporting *smr, tpia_depot *depot, const char *evaluation, const char *targetName );
tpia_target *tpia_depot_addTargetFromMap( statusMessageReporting *smr, tpia_depot *depot, tpia_map *map, const char *evaluation, const char *targetName );

#if defined __cplusplus
    }
    }
#endif

#endif          /* End of tpia_depot_h_included. */
