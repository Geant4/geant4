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
// IMPORTANT: This is a trick that allows to use
// G4 random engine even WITHIN the FLUKA INTERNAL functions.
// The idea is to provide a function with the same signature as the FLUKA one,
// but that relies on G4 random engine.
//
// The FLUKA flrnlp object file is replaced 
// by the object generated after compilation of this file
// (see FlukaInterface GNUmakefile).
//
// Author: G.Hugo, 01 August 2022
//
// ***************************************************************************
#ifdef G4_USE_FLUKA

#include "flrnlp.hh"
// G4
#include "Randomize.hh"


void flrnlp_(G4double rndvec[], const G4int& nvect) {

	for (G4int randomNumberIndex = 0; randomNumberIndex < nvect; ++randomNumberIndex) {
		rndvec[randomNumberIndex] = G4UniformRand();
	}

}


#endif // G4_USE_FLUKA
