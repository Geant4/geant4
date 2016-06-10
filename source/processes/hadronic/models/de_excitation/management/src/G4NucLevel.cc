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
// $Id: G4NucLevel.cc 88375 2015-02-16 17:31:21Z vnivanch $
//
// -------------------------------------------------------------------
//
//      GEANT4 header file 
//
//      File name:     G4NucLevel
//
//      Author:        V.Ivanchenko
// 
//      Creation date: 4 January 2012
//
//      Modifications:
//  13.02.2015 Design change for gamma de-excitation 
//      
// -------------------------------------------------------------------

#include "G4NucLevel.hh"

G4NucLevel::G4NucLevel(const std::vector<G4float>&  eTransition,
		       const std::vector<G4float>&  wLevelGamma,
		       const std::vector<G4float>&  wLevelGammaE,
		       const std::vector<G4float>&  wGamma,
		       const std::vector<G4int>&  vTrans,
		       const std::vector<const std::vector<G4float>*>& wShell)
  :fFinalEnergy(eTransition),fGammaCumProbability(wLevelGamma),
   fGammaECumProbability(wLevelGammaE),fGammaProbability(wGamma),
   fTrans(vTrans),fShellProbability(wShell)
{
  length = eTransition.size();
}

G4NucLevel::~G4NucLevel()
{
  //G4cout << "G4NucLevel: length= " << length << " " 
  //	 << fGammaProbability.size() << G4endl;
  for(size_t i=0; i<length; ++i) {
    delete fShellProbability[i];
  }
}
