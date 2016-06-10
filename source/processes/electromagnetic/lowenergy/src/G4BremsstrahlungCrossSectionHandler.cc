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
// $Id: G4BremsstrahlungCrossSectionHandler.cc 66241 2012-12-13 18:34:42Z gunter $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4BremsstrahlungCrossSectionHandler
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 25 September 2001
//
// Modifications:
// 
// 10.10.2001 MGP Revision to improve code quality and consistency with design
// 21.01.2003 VI  cut per region
// 03.03.2009 LP  Added public method to make a easier migration of
//                G4LowEnergyBremsstrahlung to G4LivermoreBremsstrahlungModel
//
// 15 Jul 2009   Nicolas A. Karakatsanis
//
//                           - BuildCrossSectionForMaterials method was revised in order to calculate the 
//                             logarithmic values of the loaded data. 
//                             It retrieves the data values from the G4EMLOW data files but, then, calculates the
//                             respective log values and loads them to seperate data structures.
//                             The EM data sets, initialized this way, contain both non-log and log values.
//                             These initialized data sets can enhance the computing performance of data interpolation
//                             operations
//
//
//
//
// -------------------------------------------------------------------

#include "G4BremsstrahlungCrossSectionHandler.hh"
#include "G4eBremsstrahlungSpectrum.hh"
#include "G4DataVector.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4Material.hh"
#include "G4ProductionCutsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4BremsstrahlungCrossSectionHandler::G4BremsstrahlungCrossSectionHandler(const G4VEnergySpectrum* spec,
									 G4VDataSetAlgorithm* )
  : theBR(spec)
{
  interp = new G4SemiLogInterpolation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4BremsstrahlungCrossSectionHandler::~G4BremsstrahlungCrossSectionHandler()
{
  delete interp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4VEMDataSet*>*
G4BremsstrahlungCrossSectionHandler::BuildCrossSectionsForMaterials(const G4DataVector& energyVector,
								    const G4DataVector* energyCuts)
{
  std::vector<G4VEMDataSet*>* set = new std::vector<G4VEMDataSet*>;

  G4DataVector* energies;
  G4DataVector* cs;

  G4DataVector* log_energies;
  G4DataVector* log_cs;

  G4int nOfBins = energyVector.size();

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  for (size_t mLocal=0; mLocal<numOfCouples; mLocal++) {

    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(mLocal);
    const G4Material* material= couple->GetMaterial();
    const G4ElementVector* elementVector = material->GetElementVector();
    const G4double* nAtomsPerVolume = material->GetVecNbOfAtomsPerVolume();
    G4int nElements = material->GetNumberOfElements();

    G4double tcut  = (*energyCuts)[mLocal];

    G4VDataSetAlgorithm* algo = interp->Clone();
    G4VEMDataSet* setForMat = new G4CompositeEMDataSet(algo,1.,1.);

    for (G4int i=0; i<nElements; i++) {

      G4int Z = (G4int) ((*elementVector)[i]->GetZ());

      energies = new G4DataVector;
      cs       = new G4DataVector;

      log_energies = new G4DataVector;
      log_cs       = new G4DataVector;

      G4double density = nAtomsPerVolume[i];

      for (G4int bin=0; bin<nOfBins; bin++) {

        G4double e = energyVector[bin];
        energies->push_back(e);
        if (e==0.) e=1e-300;
        log_energies->push_back(std::log10(e));
        G4double value = 0.0;

        if(e > tcut) {
          G4double elemCs = FindValue(Z, e);

          value  = theBR->Probability(Z, tcut, e, e);

          value *= elemCs*density;
	}
        cs->push_back(value);

        if (value==0.) value=1e-300;
        log_cs->push_back(std::log10(value));
      }
      G4VDataSetAlgorithm* algol = interp->Clone();

      //G4VEMDataSet* elSet = new G4EMDataSet(i,energies,cs,algol,1.,1.);

      G4VEMDataSet* elSet = new G4EMDataSet(i,energies,cs,log_energies,log_cs,algol,1.,1.);

      setForMat->AddComponent(elSet);
    }
    set->push_back(setForMat);
  }

  return set;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BremsstrahlungCrossSectionHandler::GetCrossSectionAboveThresholdForElement(G4double energy,
										      G4double cutEnergy,
										      G4int Z)
{
  G4double value = 0.;
  if(energy > cutEnergy)
    {
      G4double elemCs = FindValue(Z, energy);
      value  = theBR->Probability(Z,cutEnergy, energy, energy);	  
      value *= elemCs;
    }
  return value;
}
