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
#ifndef G4UserLimitsForRD_h
#define G4UserLimitsForRD_h 1

#include "globals.hh"
#include "G4UserLimits.hh"

class G4Track;

////////////////////////////////////////////////////////////////////////////////
//
class G4UserLimitsForRD : public G4UserLimits
{
  // class description 
  //
  //   Derived class from G4UserLimits.  It is specific for
  //   RDM to select and deselect volumes in which RDM is applied
  //
  // class  description - end

public:
    G4UserLimitsForRD ( const G4String aType="RD") :
      G4UserLimits (aType)
  {  
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1)
      G4cerr <<"G4UserLimitsForRD constructor" <<endl;
#endif
  }
  ~G4UserLimitsForRD () {;}

  //  G4bool IsRDActive(const G4Track&);
};
#endif

