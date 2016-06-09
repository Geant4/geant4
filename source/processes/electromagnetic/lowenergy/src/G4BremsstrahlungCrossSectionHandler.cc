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
// $Id: G4BremsstrahlungCrossSectionHandler.cc,v 1.9 2006/06/29 19:38:42 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
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
// 10.10.2001 MGP Revision to improve code quality and consistency with design
// 21.01.2003 VI  cut per region
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

G4BremsstrahlungCrossSectionHandler::G4BremsstrahlungCrossSectionHandler(const G4VEnergySpectrum* spec,
									 G4VDataSetAlgorithm* )
  : theBR(spec)
{
  interp = new G4SemiLogInterpolation();
}


G4BremsstrahlungCrossSectionHandler::~G4BremsstrahlungCrossSectionHandler()
{
  delete interp;
}


std::vector<G4VEMDataSet*>*
G4BremsstrahlungCrossSectionHandler::BuildCrossSectionsForMaterials(const G4DataVector& energyVector,
								    const G4DataVector* energyCuts)
{
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
    const G4double* nAtomsPerVolume = material->GetVecNbOfAtomsPerVolume();
    G4int nElements = material->GetNumberOfElements();

    G4double tcut  = (*energyCuts)[m];

    G4VDataSetAlgorithm* algo = interp->Clone();
    G4VEMDataSet* setForMat = new G4CompositeEMDataSet(algo,1.,1.);

    for (G4int i=0; i<nElements; i++) {

      G4int Z = (G4int) ((*elementVector)[i]->GetZ());
      energies = new G4DataVector;
      cs       = new G4DataVector;
      G4double density = nAtomsPerVolume[i];

      for (G4int bin=0; bin<nOfBins; bin++) {

        G4double e = energyVector[bin];
        energies->push_back(e);
        G4double value = 0.0;

        if(e > tcut) {
          G4double elemCs = FindValue(Z, e);

          value  = theBR->Probability(Z, tcut, e, e);

          value *= elemCs*density;
	}
        cs->push_back(value);
      }
      G4VDataSetAlgorithm* algol = interp->Clone();
      G4VEMDataSet* elSet = new G4EMDataSet(i,energies,cs,algol,1.,1.);
      setForMat->AddComponent(elSet);
    }
    set->push_back(setForMat);
  }

  return set;
}
