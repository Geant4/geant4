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
// $Id: G4CrossSectionHandler.cc,v 1.11 2001-10-05 18:24:19 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 1 Aug 2001   MGP        Created
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

G4CrossSectionHandler::G4CrossSectionHandler()
{ }

G4CrossSectionHandler::~G4CrossSectionHandler()
{ }

G4std::vector<G4VEMDataSet*>* G4CrossSectionHandler::BuildCrossSectionsForMaterials(const G4DataVector& energyVector,
										    const G4DataVector* energyCuts)
{
  G4DataVector* energies;
  G4DataVector* data;

  G4std::vector<G4VEMDataSet*>* matCrossSections = new G4std::vector<G4VEMDataSet*>;

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
    G4Exception("G4CrossSectionHandler::G4CrossSectionHandler - no MaterialTable found)");

  G4int nMaterials = G4Material::GetNumberOfMaterials();
  
  size_t nOfBins = energyVector.size();

  for (G4int m=0; m<nMaterials; m++)
    {
      const G4Material* material= (*materialTable)[m];
      energies = new G4DataVector;
      data = new G4DataVector;
      for (size_t bin=0; bin<nOfBins; bin++)
	{
	  G4double e = energyVector[bin];
	  energies->push_back(e);
	  G4double materialCrossSection = ValueForMaterial(material,e);
	  if (materialCrossSection > 0.)
	    {
	      data->push_back(materialCrossSection);
	    }
	  else
	    {
	      data->push_back(DBL_MAX);
	    }
	}
      const G4VDataSetAlgorithm* interpolationAlgo = GetInterpolation();
      G4VEMDataSet* dataSet = new G4EMDataSet(m,energies,data,interpolationAlgo,1.,1.); 
      matCrossSections->push_back(dataSet);
    }
  return matCrossSections;
}

