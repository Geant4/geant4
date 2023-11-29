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
// Geant4 header G4EmUtility
//
// Author V.Ivanchenko 14.03.2022
//
// Utilities used at initialisation of EM physics
//

#ifndef G4EmUtility_h
#define G4EmUtility_h 1

#include <vector>
#include "globals.hh"
#include "G4Region.hh"
#include "G4PhysicsTable.hh"
#include "G4EmTableType.hh"
#include "G4ParticleDefinition.hh"
#include "G4VDiscreteProcess.hh"
#include "G4LossTableBuilder.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4DataVector.hh"
#include "G4VEmModel.hh"

class G4EmUtility
{
public:

  // find G4Region pointer by name, by default no verbosity
  static const G4Region* FindRegion(const G4String& regionName,
                                    const G4int verbose = 0);

  // sample random G4Element for the case if cross section is
  // proportional to the number of electrons in an atom 
  static const G4Element* SampleRandomElement(const G4Material*);

  // sample random G4Isotope
  static const G4Isotope* SampleRandomIsotope(const G4Element*);

  // find energy of cross section maximum for all couples
  static std::vector<G4double>* FindCrossSectionMax(G4PhysicsTable*);
  static std::vector<G4double>* 
  FindCrossSectionMax(G4VDiscreteProcess*, const G4ParticleDefinition*);

  // fill structure describing more than one peak in cross sections
  static std::vector<G4TwoPeaksXS*>*
  FillPeaksStructure(G4PhysicsTable*, G4LossTableBuilder*);

  // model initialisation
  static void InitialiseElementSelectors(G4VEmModel*,
                                         const G4ParticleDefinition*,
                                         const G4DataVector& cuts,
                                         const G4double emin,
                                         const G4double emax);
};

#endif


