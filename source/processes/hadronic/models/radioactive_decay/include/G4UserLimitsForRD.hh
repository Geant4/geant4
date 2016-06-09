//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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

