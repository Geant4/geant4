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
// $Id: G4CrossSectionHandler.cc,v 1.18 2006/06/29 19:38:48 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
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

#include "G4CrossSectionHandler.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4ShellEMDataSet.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "Randomize.hh"
#include <map>
#include <vector>

#include "G4LogLogInterpolation.hh"

G4CrossSectionHandler::G4CrossSectionHandler()
{ }

G4CrossSectionHandler::~G4CrossSectionHandler()
{ }

std::vector<G4VEMDataSet*>*
G4CrossSectionHandler::BuildCrossSectionsForMaterials(const G4DataVector& energyVector,
						      const G4DataVector*)
{
  G4DataVector* energies;
  G4DataVector* data;

  std::vector<G4VEMDataSet*>* matCrossSections = new std::vector<G4VEMDataSet*>;

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  size_t nOfBins = energyVector.size();
  const G4VDataSetAlgorithm* interpolationAlgo = CreateInterpolation();

  for (size_t m=0; m<numOfCouples; m++)
    {
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(m);
      const G4Material* material= couple->GetMaterial();
      G4int nElements = material->GetNumberOfElements();
      const G4ElementVector* elementVector = material->GetElementVector();
      const G4double* nAtomsPerVolume = material->GetAtomicNumDensityVector();

      G4VDataSetAlgorithm* algo = interpolationAlgo->Clone();

      G4VEMDataSet* setForMat = new G4CompositeEMDataSet(algo,1.,1.);

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

        G4VDataSetAlgorithm* algo1 = interpolationAlgo->Clone();
        G4VEMDataSet* elSet = new G4EMDataSet(i,energies,data,algo1,1.,1.);
        setForMat->AddComponent(elSet);
      }

      matCrossSections->push_back(setForMat);
    }
  return matCrossSections;
}

