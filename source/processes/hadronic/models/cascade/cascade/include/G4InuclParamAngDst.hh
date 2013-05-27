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
// Description: intermediate base class for INUCL parametrizations of
//		three-body angular distributions in Bertini-style cascade
//
// NOTE:  Coefficient arrays have fixed dimensions, validated by compiler

#ifndef G4InuclParamAngDst_h
#define G4InuclParamAngDst_h 1

#include "globals.hh"
#include "G4VThreeBodyAngDst.hh"


class G4InuclParamAngDst : public G4VThreeBodyAngDst {
public:
  // NOTE:  Array arguments must be STATIC, GLOBAL declarations
  G4InuclParamAngDst(const G4String& name, 
		     const G4double (&abnC)[2][4][4],
		     G4int verbose=0)
    : G4VThreeBodyAngDst(name, verbose), coeffAB(abnC) {;}

  virtual ~G4InuclParamAngDst() {;}
  
  virtual G4double GetCosTheta(G4int ptype, G4double ekin) const;

  // FIXME: Must re-declare base class interface with call-through
  //	    to avoid "hidden function" compiler warnings
  virtual G4double GetCosTheta(const G4double& ekin, const G4double& pcm) const {
    return G4VThreeBodyAngDst::GetCosTheta(ekin, pcm);
  }

protected:
  const G4double (&coeffAB)[2][4][4];	// (coeffs Ekin^0..3) * S^0..3
};        

#endif	/* G4InuclParamAngDst_h */
