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
// $Id: G4CrossSectionHandler.cc,v 1.14 2002-07-19 17:32:48 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 1  Aug 2001   MGP        Created
// 19 Jul 2002   VI         Create composite data set for material
//
// -------------------------------------------------------------------

#include "G4CrossSectionHandler.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4ShellEMDataSet.hh"
#include "G4MaterialTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "Randomize.hh" 
#include "g4std/map"
#include "g4std/vector"
#include "g4std/fstream"
#include "g4std/strstream"

#include "G4LogLogInterpolation.hh"

G4CrossSectionHandler::G4CrossSectionHandler()
{ }

G4CrossSectionHandler::~G4CrossSectionHandler()
{ }

G4std::vector<G4VEMDataSet*>* 
G4CrossSectionHandler::BuildCrossSectionsForMaterials(const G4DataVector& energyVector,
						      const G4DataVector*)
{
  G4DataVector* energies;
  G4DataVector* data;

  G4std::vector<G4VEMDataSet*>* matCrossSections = new G4std::vector<G4VEMDataSet*>;

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  G4int nMaterials = G4Material::GetNumberOfMaterials();
  
  size_t nOfBins = energyVector.size();
  const G4VDataSetAlgorithm* interpolationAlgo = CreateInterpolation();

  for (G4int m=0; m<nMaterials; m++)
    {
      const G4Material* material= (*materialTable)[m];
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

