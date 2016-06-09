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
// $Id: G4QHadronBuilder.hh,v 1.2 2006/12/12 11:02:22 mkossov Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class header file
//
//                 History: 
//     Created by Mikhail Kossov, October 2006
//     Simple service class building hadron out of two partons; 
//     For comparison mirror member functions are taken from G4 class:
//     G4HadronBuilder
// -----------------------------------------------------------------------------
//

#ifndef G4QHadronBuilder_h
#define G4QHadronBuilder_h 1

#include "globals.hh"
#include <vector>
#include "G4QHadron.hh"
#include "G4QParton.hh"

class G4QHadronBuilder
{
public:
  G4QHadronBuilder(); 
  G4QHadron* Build(G4QParton* black, G4QParton* white);
  G4QHadron* BuildLowSpin(G4QParton* black, G4QParton* white);
  G4QHadron* BuildHighSpin(G4QParton* black, G4QParton* white);
private:
  enum Spin {SpinZero=1, SpinHalf=2, SpinOne=3, SpinThreeHalf=4};
  G4QHadron* Meson(G4QParton* black, G4QParton* white, Spin spin);
  G4QHadron* Baryon(G4QParton* black,G4QParton* white, Spin spin);
		// Body
  G4double mesonSpinMix;
  G4double baryonSpinMix;
  std::vector<G4double> scalarMesonMixings;
  std::vector<G4double> vectorMesonMixings;
};

// G4QHcreate type is an interface to Build/BuildLowSpin/BuildHighSpin member functions
typedef G4QHadron* (G4QHadronBuilder::*G4QHcreate) (G4QParton*,G4QParton*);
#endif
