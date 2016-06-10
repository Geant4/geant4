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
/// \file processes/phonon/include/G4PhononPolarization.hh
/// \brief Definition of the G4PhononPolarization enum
//
// $Id: G4PhononPolarization.hh 75725 2013-11-05 16:52:30Z mkelsey $
//
#ifndef G4PhononPolarization_h
#define G4PhononPolarization_h 1

#include "globals.hh"

class G4ParticleDefinition;


namespace G4PhononPolarization {
  enum { Long=0, TransSlow=1, TransFast=2, UNKNOWN=-1 };

  G4int Get(const G4ParticleDefinition* aPD);
  G4ParticleDefinition* Get(G4int pol);
}

#endif	/* G4PhononPolarization_h */
