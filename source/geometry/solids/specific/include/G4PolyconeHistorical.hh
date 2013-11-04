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
//
// $Id:$
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4PolyconeHistorical
//
// Class description:
//
//   Data structure for G4Polycone

// --------------------------------------------------------------------

#ifndef G4Polycone_historical_hh
#define G4Polycone_historical_hh

#include "G4Types.hh"

class G4PolyconeHistorical
{
  public:
    G4PolyconeHistorical();
    G4PolyconeHistorical( G4int z_planes );
    ~G4PolyconeHistorical();
    G4PolyconeHistorical( const G4PolyconeHistorical& source );
    G4PolyconeHistorical& operator=( const G4PolyconeHistorical& right );

    G4double Start_angle;
    G4double Opening_angle;
    G4int   Num_z_planes;
    G4double *Z_values;
    G4double *Rmin;
    G4double *Rmax;
};

#endif
