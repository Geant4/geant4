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
// $Id: G4ShellVacancy.cc
// GEANT4 tag $Name: 
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
#include "G4Material.hh"
#include "G4Poisson.hh"
#include "G4VEMDataSet.hh"

G4ShellVacancy::G4ShellVacancy()

{ }

G4ShellVacancy::~G4ShellVacancy()

{ 
  G4int size = xsis.size();
  for (G4int k =0; k<size; k++)
    {
      delete xsis[k];
      xsis[k] = 0;
    } 
}

void G4ShellVacancy::AddXsiTable(G4VEMDataSet* set)

{
  xsis.push_back(set);
}

G4std::vector<G4int> G4ShellVacancy::GenerateNumberOfIonisations(const G4Material* 
								 material, 
								 G4double 
								 incidentEnergy, 
								 G4double eLoss) const

{ 
  G4std::vector<G4int> numberOfIonisations; 

  G4int numberOfElements = material->GetNumberOfElements();

  for (G4int i = 0; i<numberOfElements; i++)
    {
      G4double averageNumberOfIonisations = AverageNOfIonisations(material,
	  	   					          i,
							          incidentEnergy,
							          eLoss);
      G4int ionisations = 0;
      if(averageNumberOfIonisations > 0.0) {
        ionisations = (G4int) G4Poisson(averageNumberOfIonisations);
      }

      numberOfIonisations.push_back(ionisations);
    
    }
  return  numberOfIonisations;

}

G4double G4ShellVacancy::AverageNOfIonisations(const G4Material* material,
	  				             G4int index, 
					             G4double energy,
					             G4double eLoss) const

{
  G4int indexOfElementInMaterial= -1;
  
  G4double averageEnergy = energy - eLoss/2.;
  
  size_t indexInMaterialTable = material->GetIndex();

  G4VEMDataSet* aSetOfXsi = xsis[indexInMaterialTable];

  G4double aXsi = aSetOfXsi->FindValue(averageEnergy,index); 
        
  return aXsi * eLoss;  
}
