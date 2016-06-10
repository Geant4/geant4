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
/// \file processes/phonon/src/G4PhononPolarization.cc
/// \brief implementation of the G4PhononPolarization enum
//
// $Id: G4PhononPolarization.cc 75725 2013-11-05 16:52:30Z mkelsey $
//

#include "G4PhononPolarization.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhononLong.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"


G4int G4PhononPolarization::Get(const G4ParticleDefinition* aPD) {
  if (aPD == G4PhononLong::Definition())      return Long;
  if (aPD == G4PhononTransSlow::Definition()) return TransSlow;
  if (aPD == G4PhononTransFast::Definition()) return TransFast;
  return UNKNOWN;
}

G4ParticleDefinition* G4PhononPolarization::Get(G4int pol) {
  switch (pol) {
  case Long:      return G4PhononLong::Definition(); break;
  case TransSlow: return G4PhononTransSlow::Definition(); break;
  case TransFast: return G4PhononTransFast::Definition(); break;
  default: ;
  }

  return 0;
}
