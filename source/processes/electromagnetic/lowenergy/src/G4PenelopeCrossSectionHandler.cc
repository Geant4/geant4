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
// GEANT4 Class file
//
//
// File name:     G4PenelopeCrossSectionHandler
//
// Author:        Luciano Pandola
// 
// Creation date: 04 Jul 2003
//
// -------------------------------------------------------------------
#include "G4PenelopeCrossSectionHandler.hh"
#include "G4PenelopeIonisation.hh"
#include "G4DataVector.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LinLogLogInterpolation.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4Material.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleDefinition.hh"

G4PenelopeCrossSectionHandler::G4PenelopeCrossSectionHandler(G4PenelopeIonisation* theProcess,
				const G4ParticleDefinition& aParticleType,
				G4VDataSetAlgorithm* alg,
				G4double emin, 
				G4double emax, 
				G4int nbin)   
 :  G4VCrossSectionHandler(),
    process(theProcess),
    particle(aParticleType)
{
  G4VCrossSectionHandler::Initialise(alg, emin, emax, nbin);
  interp = new G4LinLogLogInterpolation();
}


G4PenelopeCrossSectionHandler::~G4PenelopeCrossSectionHandler()
{
  delete interp;
}


std::vector<G4VEMDataSet*>* G4PenelopeCrossSectionHandler::BuildCrossSectionsForMaterials(
                        const G4DataVector& energyVector,
                        const G4DataVector* energyCuts)
{
  //G4int verbose = 0;
  std::vector<G4VEMDataSet*>* set = new std::vector<G4VEMDataSet*>;

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
    G4double electronVolumeDensity = 
      material->GetTotNbOfElectPerVolume();  //electron density
    G4int nElements = material->GetNumberOfElements();

    /*
    if(verbose > 0) {
      G4cout << "Penelope CS for " << m << "th material "
             << material->GetName()
             << "  eEl= " << nElements << G4endl;
    }
    */

    G4double tcut  = (*energyCuts)[m];

    G4VDataSetAlgorithm* algo = interp->Clone();
    G4VEMDataSet* setForMat = new G4CompositeEMDataSet(algo,1.,1.);

    for (G4int i=0; i<nElements; i++) {

      G4int Z = (G4int) (*elementVector)[i]->GetZ();
      energies = new G4DataVector;
      cs       = new G4DataVector;
      G4double density = nAtomsPerVolume[i];

      for (G4int bin=0; bin<nOfBins; bin++) {

        G4double e = energyVector[bin];
        energies->push_back(e);
        G4double value = 0.0;

        if(e > tcut) {
            G4double cross = FindValue(Z, e);
            G4double p = process->CalculateCrossSectionsRatio(e,tcut,Z,electronVolumeDensity,
							      particle);
            value += cross * p * density;

	    /*
            if(verbose>0 && m == 0 && e>=1. && e<=0.) {
              G4cout << "G4PenIonCrossSH: e(MeV)= " << e/MeV
                     << " cross= " << cross
                     << " p= " << p
                     << " value= " << value
                     << " tcut(MeV)= " << tcut/MeV
                     << " rho= " << density
                     << " Z= " << Z
                     << G4endl;
	    }
	   */

	}
        cs->push_back(value);
      }

      G4VDataSetAlgorithm* algo = interp->Clone();
      G4VEMDataSet* elSet = new G4EMDataSet(i,energies,cs,algo,1.,1.);
      setForMat->AddComponent(elSet);
    }
    set->push_back(setForMat);
  }
  return set;
}


