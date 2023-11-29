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
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 21 Sept 2001 Elena Guardincerri     Created
// 25 Mar  2002  V.Ivanchenko          Change AverageNOfIonisations int->double
//
// -------------------------------------------------------------------

#include "G4ShellVacancy.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Material.hh"
#include "G4Poisson.hh"
#include "G4VEMDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4ShellVacancy::G4ShellVacancy()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4ShellVacancy::~G4ShellVacancy()
{
  std::size_t size = xsis.size();
  for (std::size_t k =0; k<size; ++k)
    {
      delete xsis[k];
      xsis[k] = nullptr;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4ShellVacancy::AddXsiTable(G4VEMDataSet* p)
{
  xsis.push_back(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4int> G4ShellVacancy::GenerateNumberOfIonisations(const G4MaterialCutsCouple*
								 couple,
								 G4double
								 incidentEnergy,
								 G4double eLoss) const
{
  std::vector<G4int> numberOfIonisations;
  const G4Material* material = couple->GetMaterial();
  G4int numberOfElements = (G4int)material->GetNumberOfElements();

  for (G4int i = 0; i<numberOfElements; ++i)
    {
      G4double averageNumberOfIonisations = AverageNOfIonisations(couple,
	  	   					          i,
							          incidentEnergy,
							          eLoss);
      G4int ionisations = 0;
      if(averageNumberOfIonisations > 0.0) {
        ionisations = (G4int) G4Poisson(averageNumberOfIonisations);
      }

      numberOfIonisations.push_back(ionisations);

    }
  return numberOfIonisations;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ShellVacancy::AverageNOfIonisations(const G4MaterialCutsCouple* couple,
	  				             G4int index,
					             G4double energy,
					             G4double eLoss) const

{
  G4double averageEnergy = energy - eLoss/2.;
  std::size_t indexInMaterialTable = couple->GetIndex();

  G4VEMDataSet* aSetOfXsi = xsis[indexInMaterialTable];

  G4double aXsi = aSetOfXsi->FindValue(averageEnergy,index);

  return aXsi * eLoss;
}
