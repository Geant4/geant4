//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
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
#ifndef G4GIDI_h_included
#define G4GIDI_h_included 1

#include <string>
#include <list>
#include <vector>

#include <G4GIDI_Misc.hh>
#include <G4GIDI_map.hh>
#include <G4GIDI_target.hh>
#include <G4GIDI_mass.hh>

class G4GIDI {

    private:
        int projectileID;
        std::list<G4GIDI_map *> dataDirectories;
        std::vector<G4GIDI_target *> targets;

        int init( int ip );

    public:
        G4GIDI( int ip, std::string &dataDirectory );
        G4GIDI( int ip, std::list<std::string> &dataDirectory );
        ~G4GIDI( );

        int numberOfDataDirectories( void );
        int addDataDirectory( std::string &dataDirectory );
        int removeDataDirectory( std::string &dataDirectory );
        std::string getDataDirectoryAtIndex( int index );
        std::vector<std::string> *getDataDirectories( void );

        bool isThisDataAvailable( std::string &lib_name, int iZ, int iA, int iM = 0 );
        bool isThisDataAvailable( std::string &lib_name, std::string &targetName );

        char *dataFilename( std::string &lib_name, int iZ, int iA, int iM = 0 );
        char *dataFilename( std::string &lib_name, std::string &targetName );

        std::vector<std::string> *getNamesOfAvailableLibraries( int iZ, int iA, int iM = 0 );
        std::vector<std::string> *getNamesOfAvailableLibraries( std::string &targetName );

        std::vector<std::string> *getNamesOfAvailableTargets( void );

        G4GIDI_target *readTarget( std::string &lib_name, int iZ, int iA, int iM = 0, bool bind = true );
        G4GIDI_target *readTarget( std::string &lib_name, std::string &targetName, bool bind = true );

        G4GIDI_target *getAlreadyReadTarget( int iZ, int iA, int iM = 0 );
        G4GIDI_target *getAlreadyReadTarget( std::string &targetName );

        int freeTarget( int iZ, int iA, int iM = 0 );
        int freeTarget( std::string &targetSymbol );
        int freeTarget( G4GIDI_target *target );

        std::vector<std::string> *getListOfReadTargetsNames( void );
};

#endif      // End of G4GIDI_h_included
