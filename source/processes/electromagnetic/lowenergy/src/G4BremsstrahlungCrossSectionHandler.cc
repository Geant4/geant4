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
// $Id: G4BremsstrahlungCrossSectionHandler.cc,v 1.3 2001-10-10 11:59:35 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4MaterialTable.hh"

G4BremsstrahlungCrossSectionHandler::G4BremsstrahlungCrossSectionHandler(const G4VEnergySpectrum* spec, 
									 G4VDataSetAlgorithm* alg)
  : theBR(spec)
{
  interp = new G4SemiLogInterpolation();
}


G4BremsstrahlungCrossSectionHandler::~G4BremsstrahlungCrossSectionHandler() 
{
  delete interp;
}


G4std::vector<G4VEMDataSet*>* 
G4BremsstrahlungCrossSectionHandler::BuildCrossSectionsForMaterials(const G4DataVector& energyVector, 
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
 
      G4int Z = (G4int) ((*elementVector)[i]->GetZ());
      energies = new G4DataVector;
      cs       = new G4DataVector;
      G4double density = nAtomsPerVolume[i];

      for (G4int bin=0; bin<nBins; bin++) {

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
      G4VEMDataSet* elSet = new G4EMDataSet(i,energies,cs,interp,1.,1.);
      setForMat->AddComponent(elSet);
    }
    set->push_back(setForMat);
  }
 
  return set;
}
