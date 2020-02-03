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
// Author:  Michael Kelsey (SLAC)
// Date:    7 March 2013
//
// Description: pure virtual base class for final-state momentum
//              distributions in Bertini-style cascade
//
// 20130308  M. Kelsey -- Add access to name string for diagnostic utilities

#ifndef G4VMultiBodyMomDst_h
#define G4VMultiBodyMomDst_h 1

#include <vector>
#include "globals.hh"

class G4VMultiBodyMomDst {
public:
  G4VMultiBodyMomDst(const G4String& name, G4int verbose=0);
  virtual ~G4VMultiBodyMomDst() {;}
  
  virtual G4double GetMomentum(G4int ptype, const G4double& ekin) const = 0;

  virtual void setVerboseLevel(G4int verbose = 0) { verboseLevel = verbose; }
  virtual const G4String& GetName() const { return theName; }

protected:
  G4String theName;
  G4int verboseLevel;
};        

#endif	/* G4VMultiBodyMomDst_h */
