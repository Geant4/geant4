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
// $Id$
// Author:  Michael Kelsey (SLAC)
// Date:    22 April 2013
//
// Description: pure virtual base class for three-body final state angular
//              distributions in Bertini-style cascade, subclassed from
//		G4VTwoBodyAngDst (c.f. ROOT's TH1 -> TH2 -> TH3)
//

#ifndef G4VThreeBodyAngDst_h
#define G4VThreeBodyAngDst_h 1

#include "G4VTwoBodyAngDst.hh"

class G4VThreeBodyAngDst : public G4VTwoBodyAngDst {
public:
  G4VThreeBodyAngDst(const G4String& name, G4int verbose=0)
    : G4VTwoBodyAngDst(name, verbose) {;}
  virtual ~G4VThreeBodyAngDst() {;}

  // Three-body mode needs particle type and bullet energy
  virtual G4double GetCosTheta(G4int ptype, G4double ekin) const = 0;

  // Implement base-class interface to re-interpret 'pcm' as ptype
  virtual G4double GetCosTheta(const G4double& ekin, const G4double& pcm) const {
    return this->GetCosTheta((G4int)pcm, ekin);
  }
};        

#endif	/* G4VThreeBodyAngDst_h */
