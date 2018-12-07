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
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 1  Aug 2001   MGP        Created
// 19 Jul 2002   VI         Create composite data set for material
// 24 Apr 2003   VI         Cut per region mfpt
//
// -------------------------------------------------------------------

#include "G4RDCrossSectionHandler.hh"
#include "G4RDVDataSetAlgorithm.hh"
#include "G4RDVEMDataSet.hh"
#include "G4RDEMDataSet.hh"
#include "G4RDCompositeEMDataSet.hh"
#include "G4RDShellEMDataSet.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "Randomize.hh"
#include <map>
#include <vector>

#include "G4RDLogLogInterpolation.hh"

G4RDCrossSectionHandler::G4RDCrossSectionHandler()
{ }

G4RDCrossSectionHandler::~G4RDCrossSectionHandler()
{ }

std::vector<G4RDVEMDataSet*>*
G4RDCrossSectionHandler::BuildCrossSectionsForMaterials(const G4DataVector& energyVector,
						      const G4DataVector*)
{
  G4DataVector* energies;
  G4DataVector* data;

  std::vector<G4RDVEMDataSet*>* matCrossSections = new std::vector<G4RDVEMDataSet*>;

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  size_t nOfBins = energyVector.size();
  const G4RDVDataSetAlgorithm* interpolationAlgo = CreateInterpolation();

  for (size_t m=0; m<numOfCouples; m++)
    {
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(m);
      const G4Material* material= couple->GetMaterial();
      G4int nElements = material->GetNumberOfElements();
      const G4ElementVector* elementVector = material->GetElementVector();
      const G4double* nAtomsPerVolume = material->GetAtomicNumDensityVector();

      G4RDVDataSetAlgorithm* algo = interpolationAlgo->Clone();

      G4RDVEMDataSet* setForMat = new G4RDCompositeEMDataSet(algo,1.,1.);

      for (G4int i=0; i<nElements; i++) {
 
        G4int Z = (G4int) (*elementVector)[i]->GetZ();
        G4double density = nAtomsPerVolume[i];

        energies = new G4DataVector;
        data = new G4DataVector;


        for (size_t bin=0; bin<nOfBins; bin++)
	  {
	    G4double e = energyVector[bin];
	    energies->push_back(e);
	    G4double cross = density*FindValue(Z,e);
	    data->push_back(cross);
	  }

        G4RDVDataSetAlgorithm* algo1 = interpolationAlgo->Clone();
        G4RDVEMDataSet* elSet = new G4RDEMDataSet(i,energies,data,algo1,1.,1.);
        setForMat->AddComponent(elSet);
      }

      matCrossSections->push_back(setForMat);
    }
  return matCrossSections;
}

