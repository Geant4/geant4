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
// File name:     G4RDeIonisationCrossSectionHandler
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 25 Sept 2001
//
// Modifications: 
// 10 Oct 2001  M.G. Pia     Revision to improve code quality and consistency with design
// 19 Jul 2002   VI          Create composite data set for material
// 21 Jan 2003  V.Ivanchenko Cut per region
//
// -------------------------------------------------------------------

#include "G4RDeIonisationCrossSectionHandler.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RDVEnergySpectrum.hh"
#include "G4DataVector.hh"
#include "G4RDCompositeEMDataSet.hh"
#include "G4RDVDataSetAlgorithm.hh"
#include "G4RDLinLogLogInterpolation.hh"
#include "G4RDSemiLogInterpolation.hh"
#include "G4RDVEMDataSet.hh"
#include "G4RDEMDataSet.hh"
#include "G4Material.hh"
#include "G4ProductionCutsTable.hh"


G4RDeIonisationCrossSectionHandler::G4RDeIonisationCrossSectionHandler(
    const G4RDVEnergySpectrum* spec, G4RDVDataSetAlgorithm* alg,
    G4double emin, G4double emax, G4int nbin)
 :  G4RDVCrossSectionHandler(),
    theParam(spec)
{
  G4RDVCrossSectionHandler::Initialise(alg, emin, emax, nbin);
  interp = new G4RDLinLogLogInterpolation();
}


G4RDeIonisationCrossSectionHandler::~G4RDeIonisationCrossSectionHandler()
{
  delete interp;
}


std::vector<G4RDVEMDataSet*>* G4RDeIonisationCrossSectionHandler::BuildCrossSectionsForMaterials(
                        const G4DataVector& energyVector,
                        const G4DataVector* energyCuts)
{
  G4int verbose = 0;
  std::vector<G4RDVEMDataSet*>* set = new std::vector<G4RDVEMDataSet*>;

  G4DataVector* energies;
  G4DataVector* cs;
  G4int nOfBins = energyVector.size();

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  for (size_t m=0; m<numOfCouples; m++) {

    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(m);
    const G4Material* material= couple->GetMaterial();
    const G4ElementVector* elementVector = material->GetElementVector();
    const G4double* nAtomsPerVolume = material->GetAtomicNumDensityVector();
    G4int nElements = material->GetNumberOfElements();

    if(verbose > 0) {
      G4cout << "eIonisation CS for " << m << "th material "
             << material->GetName()
             << "  eEl= " << nElements << G4endl;
    }

    G4double tcut  = (*energyCuts)[m];

    G4RDVDataSetAlgorithm* algo = interp->Clone();
    G4RDVEMDataSet* setForMat = new G4RDCompositeEMDataSet(algo,1.,1.);

    for (G4int i=0; i<nElements; i++) {

      G4int Z = (G4int) (*elementVector)[i]->GetZ();
      G4int nShells = NumberOfComponents(Z);
      energies = new G4DataVector;
      cs       = new G4DataVector;
      G4double density = nAtomsPerVolume[i];

      for (G4int bin=0; bin<nOfBins; bin++) {

        G4double e = energyVector[bin];
        energies->push_back(e);
        G4double value = 0.0;

        if(e > tcut) {
          for (G4int n=0; n<nShells; n++) {
            G4double cross = FindValue(Z, e, n);
            G4double p = theParam->Probability(Z, tcut, e, e, n);
            value += cross * p * density;

            if(verbose>0 && m == 0 && e>=1. && e<=0.) {
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
	}
        cs->push_back(value);
      }
      G4RDVDataSetAlgorithm* algo = interp->Clone();
      G4RDVEMDataSet* elSet = new G4RDEMDataSet(i,energies,cs,algo,1.,1.);
      setForMat->AddComponent(elSet);
    }
    set->push_back(setForMat);
  }

  return set;
}


