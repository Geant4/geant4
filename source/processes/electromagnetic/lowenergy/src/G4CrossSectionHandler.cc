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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CrossSectionHandler::G4CrossSectionHandler()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CrossSectionHandler::~G4CrossSectionHandler()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4VEMDataSet*>*
G4CrossSectionHandler::BuildCrossSectionsForMaterials(const G4DataVector& energyVector,
						      const G4DataVector*)
{
  G4DataVector* energies;
  G4DataVector* data;

  G4DataVector* log_energies;
  G4DataVector* log_data;

  std::vector<G4VEMDataSet*>* matCrossSections = new std::vector<G4VEMDataSet*>;

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();

  std::size_t nOfBins = energyVector.size();
  const G4VDataSetAlgorithm* interpolationAlgo = CreateInterpolation();

  for (G4int mLocal=0; mLocal<numOfCouples; ++mLocal)
    {
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(mLocal);
      const G4Material* material= couple->GetMaterial();
      G4int nElements = (G4int)material->GetNumberOfElements();
      const G4ElementVector* elementVector = material->GetElementVector();
      const G4double* nAtomsPerVolume = material->GetAtomicNumDensityVector();

      G4VDataSetAlgorithm* algo = interpolationAlgo->Clone();

      G4VEMDataSet* setForMat = new G4CompositeEMDataSet(algo,1.,1.);

      for (G4int i=0; i<nElements; ++i) {
 
        G4int Z = (G4int) (*elementVector)[i]->GetZ();
        G4double density = nAtomsPerVolume[i];

        energies = new G4DataVector;
        data = new G4DataVector;

        log_energies = new G4DataVector;
        log_data = new G4DataVector;

        for (std::size_t bin=0; bin<nOfBins; ++bin)
	  {
	    G4double e = energyVector[bin];
	    energies->push_back(e);
            if (e==0.) e=1e-300;
            log_energies->push_back(std::log10(e));
	    G4double cross = density*FindValue(Z,e);
	    data->push_back(cross);
            if (cross==0.) cross=1e-300;
            log_data->push_back(std::log10(cross));
	  }

        G4VDataSetAlgorithm* algo1 = interpolationAlgo->Clone();
        G4VEMDataSet* elSet = new G4EMDataSet(i,energies,data,log_energies,log_data,algo1,1.,1.);

        setForMat->AddComponent(elSet);
      }

      matCrossSections->push_back(setForMat);
    }
  delete interpolationAlgo;
  return matCrossSections;
}

