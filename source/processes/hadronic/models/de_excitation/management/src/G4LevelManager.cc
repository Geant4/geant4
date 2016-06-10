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
// $Id: G4LevelManager.cc 88375 2015-02-16 17:31:21Z vnivanch $
//
// -------------------------------------------------------------------
//
//      GEANT4 source file 
//
//      File name:     G4LevelManager
//
//      Author:        V.Ivanchenko
// 
//      Creation date: 4 January 2012
//
//      Modifications: 
//  13.02.2015 Design change for gamma de-excitation 
//      
// -------------------------------------------------------------------

#include "G4LevelManager.hh"

G4LevelManager::G4LevelManager(const std::vector<G4float>& energies,
			       const std::vector<G4float>& lifetime,
			       const std::vector<G4float>& lifetimegamma,
			       const std::vector<G4int>& spin,
			       const std::vector<const G4NucLevel*>& levels)
  : fLevelEnergy(energies),fLifeTime(lifetime),
    fLifeTimeGamma(lifetimegamma),fSpin(spin),fLevels(levels)
{ 
  nTransitions = fLevelEnergy.size() - 1; 
  //G4cout << "New G4LevelManager N= " << nTransitions << " " << fLevelEnergy.size() 
  //	 << " <" << this << ">" << G4endl;
}

G4LevelManager::~G4LevelManager()
{
  for(size_t i=0; i<=nTransitions; ++i) { delete fLevels[i]; }
}
