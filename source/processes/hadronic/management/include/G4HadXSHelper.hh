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
// Geant4 header G4HadXSHelper
//
// Author V.Ivanchenko 18.05.2022
//
// Utilities used at initialisation of hadronic physics
//

#ifndef G4HadXSHelper_h
#define G4HadXSHelper_h 1

#include "globals.hh"
#include "G4HadXSTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadronicProcess.hh"
#include <vector>

class G4HadXSHelper
{
public:

  // find energy of cross section maximum for all couples
  static std::vector<G4double>* 
  FindCrossSectionMax(G4HadronicProcess*, const G4ParticleDefinition*,
                      const G4double tmin, const G4double tmax);

  // fill structure describing more than one peak in cross sections
  static std::vector<G4TwoPeaksHadXS*>*
  FillPeaksStructure(G4HadronicProcess*, const G4ParticleDefinition*,
                     const G4double tmin, const G4double tmax);
};

#endif


