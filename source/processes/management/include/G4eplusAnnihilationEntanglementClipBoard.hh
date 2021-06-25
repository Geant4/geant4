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
//
// --------------------------------------------------------------------
// GEANT4 class header file
//
// Class Description:
// A clipboard for communicating between quantum entangled photons
// from e+ annihilation. Applies only to the two-gamma decay from the
// JP=0- state of e+e-, para-positronium. (This is the dominant process
// in materials such as human tissue, since the ortho-positronium JP=1-
// state, initially more likely, is "picked off" by a near-by electron
// and turned into para-positronium, because it has the much shorter
// meanlife.)
// Also applies to pi zero (JP=0-).
//
// ------------------ G4eplusAnnihilationEntanglementClipBoard ------------------
//
// Author: J.Allison, May 2017
//
// --------------------------------------------------------------------

// This is a concrete class based on G4VEntanglementClipBoard.  Information
// specific to this type of entanglement can be posted here.  See
// G4VEntanglementClipBoard.hh for more detail on how to use this class.

#ifndef G4eplusAnnihilationEntanglementClipBoard_hh
#define G4eplusAnnihilationEntanglementClipBoard_hh

#include "G4VEntanglementClipBoard.hh"

class G4eplusAnnihilationEntanglementClipBoard
: public G4VEntanglementClipBoard
{
public:
  G4eplusAnnihilationEntanglementClipBoard()
  : fComptonCosTheta1(0.)
  , fComptonPhi1(0.)
  {}
  ~G4eplusAnnihilationEntanglementClipBoard() {}

  // Cos(theta) and phi of the first Compton scattering of the first photon
  void SetComptonCosTheta1(G4double cosTheta1) {fComptonCosTheta1 = cosTheta1;}
  void SetComptonPhi1(G4double phi1) {fComptonPhi1 = phi1;}
  G4double GetComptonCosTheta1() const {return fComptonCosTheta1;}
  G4double GetComptonPhi1() const {return fComptonPhi1;}

protected:
  // Cos(theta) and phi of the first Compton scattering of the first photon
  G4double fComptonCosTheta1, fComptonPhi1;
};

#endif
