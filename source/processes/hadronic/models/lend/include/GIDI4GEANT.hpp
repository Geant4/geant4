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
#ifndef GIDI4GEANT_h_included
#define GIDI4GEANT_h_included 1

#include <string>
#include <list>
#include <vector>
using namespace std;

#include <GIDI4GEANT_Misc.hpp>
#include <GIDI4GEANT_map.hpp>
#include <GIDI4GEANT_target.hpp>
#include <GIDI4GEANT_mass.hpp>

class GIDI4GEANT {

    private:
        int projectileID;
        list<GIDI4GEANT_map *> dataDirectories;
        vector<GIDI4GEANT_target *> targets;

        int init( int ip );

    public:
        GIDI4GEANT( int ip, string &dataDirectory );
        GIDI4GEANT( int ip, list<string> &dataDirectory );
        ~GIDI4GEANT( );

        int numberOfDataDirectories( void );
        int addDataDirectory( string &dataDirectory );
        int removeDataDirectory( string &dataDirectory );
        string getDataDirectoryAtIndex( int index );
        vector<string> *getDataDirectories( void );

        bool isThisDataAvailable( string &lib_name, int iZ, int iA, int iM = 0 );
        bool isThisDataAvailable( string &lib_name, string &targetName );

        char *dataFilename( string &lib_name, int iZ, int iA, int iM = 0 );
        char *dataFilename( string &lib_name, string &targetName );

        vector<string> *getNamesOfAvailableLibraries( int iZ, int iA, int iM = 0 );
        vector<string> *getNamesOfAvailableLibraries( string &targetName );

        vector<string> *getNamesOfAvailableTargets( void );

        GIDI4GEANT_target *readTarget( string &lib_name, int iZ, int iA, int iM = 0, bool bind = true );
        GIDI4GEANT_target *readTarget( string &lib_name, string &targetName, bool bind = true );

        GIDI4GEANT_target *getAlreadyReadTarget( int iZ, int iA, int iM = 0 );
        GIDI4GEANT_target *getAlreadyReadTarget( string &targetName );

        int freeTarget( int iZ, int iA, int iM = 0 );
        int freeTarget( string &targetSymbol );
        int freeTarget( GIDI4GEANT_target *target );

        vector<string> *getListOfReadTargetsNames( void );
};

#endif      // End of GIDI4GEANT_h_included
