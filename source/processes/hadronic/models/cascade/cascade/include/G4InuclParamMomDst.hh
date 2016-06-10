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
// $Id: G4InuclParamMomDst.hh 67796 2013-03-08 06:18:39Z mkelsey $
// Author:  Michael Kelsey (SLAC)
// Date:    7 March 2013
//
// Description: intermediate base class for INUCL parametrizations of
//		final-state momentum distributions in Bertini-style cascade
//
// NOTE:  Coefficient arrays have fixed dimensions, validated by compiler

#ifndef G4InuclParamMomDst_h
#define G4InuclParamMomDst_h 1

#include "globals.hh"
#include "G4VMultiBodyMomDst.hh"


class G4InuclParamMomDst : public G4VMultiBodyMomDst {
public:
  // NOTE:  Array arguments must be STATIC, GLOBAL declarations
  G4InuclParamMomDst(const G4String& name, 
		     const G4double (&pqprC)[2][4][4],
		     const G4double (&psC)[2][3],
		     G4int verbose=0)
    : G4VMultiBodyMomDst(name, verbose), coeffPR(pqprC), coeffPS(psC) {;}

  virtual ~G4InuclParamMomDst() {;}
  
  virtual G4double GetMomentum(G4int ptype, const G4double& ekin) const;

protected:
  const G4double (&coeffPR)[2][4][4];	// (coeffs Ekin^0..3) * S^0..3
  const G4double (&coeffPS)[2][3];	// PS = sum coeffs * Ekin^0..2 
};        

#endif	/* G4InuclParamMomDst_h */
