//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eIonisationCrossSectionHandler
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 25 Sept 2001
//
// Modifications: 
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eIonisationCrossSectionHandler.hh"
#include "G4VEnergySpectrum.hh"
#include "G4DataVector.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eIonisationCrossSectionHandler::G4eIonisationCrossSectionHandler(
    const G4VEnergySpectrum* spec, G4VDataSetAlgorithm* alg, 
    G4double emin, G4double emax, G4int nbin)
 :  G4VCrossSectionHandler(),
    theParam(spec)
{
  G4VCrossSectionHandler::Initialise(alg, emin, emax, nbin);
  interp = new G4SemiLogInterpolation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eIonisationCrossSectionHandler::~G4eIonisationCrossSectionHandler() 
{
  delete interp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4std::vector<G4VEMDataSet*>* 
       G4eIonisationCrossSectionHandler::BuildCrossSectionsForMaterials(
                                const G4DataVector& energyVector, 
				const G4DataVector* energyCuts)
{

  G4std::vector<G4VEMDataSet*>* set = new G4std::vector<G4VEMDataSet*>;

  G4DataVector* energies;
  G4DataVector* cs;
  G4int nBins = energyVector.size();

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
    G4Exception("G4VCrossSectionHandler::G4VCrossSectionHandler - no MaterialTable found)");

  G4int nMaterials = G4Material::GetNumberOfMaterials();

  for (G4int m=0; m<nMaterials; m++) {

    const G4Material* material = (*materialTable)[m];
    const G4ElementVector* elementVector = material->GetElementVector();
    const G4double* nAtomsPerVolume = material->GetVecNbOfAtomsPerVolume();
    G4int nElements = material->GetNumberOfElements();

    G4double tcut  = (*energyCuts)[m];

    G4VEMDataSet* setForMat = new G4CompositeEMDataSet(interp,1.,1.);

    for (G4int i=0; i<nElements; i++) {
 
      G4int Z = (G4int) (*elementVector)[i]->GetZ();
      G4int nShells = NumberOfComponents(Z);
      energies = new G4DataVector;
      cs       = new G4DataVector;
      G4double density = nAtomsPerVolume[i];

      for (G4int bin=0; bin<nBins; bin++) {

        G4double e = energyVector[bin];
        energies->push_back(e);
        G4double value = 0.0;

        if(e > tcut) {
          for (G4int n=0; n<nShells; n++) {
            G4double cross = FindValue(Z, e, n);
            G4double p = theParam->Probability(Z, tcut, e, e, n);
        
            value += cross * p * density;
	  }
	}
        cs->push_back(value);
      } 
      G4VEMDataSet* elSet = new G4EMDataSet(i,energies,cs,interp,1.,1.);
      setForMat->AddComponent(elSet);
    }
    set->push_back(setForMat);
  }
 
  return set;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




