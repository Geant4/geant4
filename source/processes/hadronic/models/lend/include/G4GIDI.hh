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
# <<END-copyright>>
*/
#ifndef G4GIDI_h_included
#define G4GIDI_h_included 1

#include <string>
#include <list>
#include <vector>
//using namespace std;

#include <G4GIDI_Misc.hh>
#include <G4GIDI_map.hh>
#include <G4GIDI_target.hh>
#include <G4GIDI_mass.hh>

class G4GIDI {

    private:
        int projectileID;
        std::string projectile;
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
