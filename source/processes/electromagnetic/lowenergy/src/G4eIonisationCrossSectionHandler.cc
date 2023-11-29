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
// 10 Oct 2001  M.G. Pia     Revision to improve code quality and consistency with design
// 19 Jul 2002   VI          Create composite data set for material
// 21 Jan 2003  V.Ivanchenko Cut per region
// 28 Jan 2009  L.Pandola    Added public method to make a easier migration of 
//                           G4LowEnergyIonisation to G4LivermoreIonisationModel
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
// -------------------------------------------------------------------

#include "G4eIonisationCrossSectionHandler.hh"
#include "G4SystemOfUnits.hh"
#include "G4VEnergySpectrum.hh"
#include "G4DataVector.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LinLogLogInterpolation.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4Material.hh"
#include "G4ProductionCutsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 G4eIonisationCrossSectionHandler::G4eIonisationCrossSectionHandler(
    const G4VEnergySpectrum* spec, G4VDataSetAlgorithm* alg,
    G4double emin, G4double emax, G4int nbin)
 :  G4VCrossSectionHandler(),
    theParam(spec),verbose(0)
{
  G4VCrossSectionHandler::Initialise(alg, emin, emax, nbin);
  interp = new G4LinLogLogInterpolation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eIonisationCrossSectionHandler::~G4eIonisationCrossSectionHandler()
{
  delete interp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4VEMDataSet*>* G4eIonisationCrossSectionHandler::BuildCrossSectionsForMaterials(
                        const G4DataVector& energyVector,
                        const G4DataVector* energyCuts)
{
  std::vector<G4VEMDataSet*>* set = new std::vector<G4VEMDataSet*>;

  G4DataVector* energies;
  G4DataVector* cs;

  G4DataVector* log_energies;
  G4DataVector* log_cs;

  std::size_t nOfBins = energyVector.size();

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();

  for (G4int mLocal=0; mLocal<numOfCouples; ++mLocal) {

    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(mLocal);
    const G4Material* material= couple->GetMaterial();
    const G4ElementVector* elementVector = material->GetElementVector();
    const G4double* nAtomsPerVolume = material->GetAtomicNumDensityVector();
    G4int nElements = (G4int)material->GetNumberOfElements();

    if(verbose > 0) 
      {
	G4cout << "eIonisation CS for " << mLocal << "th material "
	       << material->GetName()
	       << "  eEl= " << nElements << G4endl;
      }

    G4double tcut  = (*energyCuts)[mLocal];

    G4VDataSetAlgorithm* algo = interp->Clone();
    G4VEMDataSet* setForMat = new G4CompositeEMDataSet(algo,1.,1.);

    for (G4int i=0; i<nElements; ++i) {

      G4int Z = (G4int) (*elementVector)[i]->GetZ();
      G4int nShells = NumberOfComponents(Z);

      energies = new G4DataVector;
      cs       = new G4DataVector;

      log_energies = new G4DataVector;
      log_cs       = new G4DataVector;

      G4double density = nAtomsPerVolume[i];

      for (std::size_t bin=0; bin<nOfBins; ++bin) {

        G4double e = energyVector[bin];
        energies->push_back(e);
        log_energies->push_back(std::log10(e));
        G4double value = 0.0;
        G4double log_value = -300;

        if(e > tcut) {
          for (G4int n=0; n<nShells; n++) {
            G4double cross = FindValue(Z, e, n);
            G4double p = theParam->Probability(Z, tcut, e, e, n);
            value += cross * p * density;

	    if(verbose>0 && mLocal == 0 && e>=1. && e<=0.) 
	    {
              G4cout << "G4eIonCrossSH: e(MeV)= " << e/MeV
                     << " n= " << n
                     << " cross= " << cross
                     << " p= " << p
                     << " value= " << value
                     << " tcut(MeV)= " << tcut/MeV
                     << " rho= " << density
                     << " Z= " << Z
                     << G4endl;
	      }
	  }
          if (value == 0.) value = 1e-300;
          log_value = std::log10(value);
	}
        cs->push_back(value);
        log_cs->push_back(log_value);
      }
      G4VDataSetAlgorithm* algoLocal = interp->Clone();
      G4VEMDataSet* elSet = new G4EMDataSet(i,energies,cs,log_energies,log_cs,algoLocal,1.,1.);

      setForMat->AddComponent(elSet);
    }
    set->push_back(setForMat);
  }

  return set;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eIonisationCrossSectionHandler::GetCrossSectionAboveThresholdForElement(G4double energy,
										   G4double cutEnergy,
										   G4int Z)
{
  G4int nShells = NumberOfComponents(Z);
  G4double value = 0.;
  if(energy > cutEnergy) 
    {
      for (G4int n=0; n<nShells; ++n) {
	G4double cross = FindValue(Z, energy, n);
	G4double p = theParam->Probability(Z, cutEnergy, energy, energy, n);
	value += cross * p;
      }
    }
  return value;
}
