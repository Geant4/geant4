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
// History:
// -----------
//  01 Oct 2011   A.M., S.I. - 1st implementation
// 
// Class description
// ----------------
//  Computation of K, L & M shell ECPSSR ionisation cross sections for protons and alphas
//  Based on the work of A. Taborda et al. 
//  X-Ray Spectrom. 2011, 40, 127-134
// ---------------------------------------------------------------------------------------

#include <fstream>
#include <iomanip>

#include "G4ecpssrFormFactorKxsModel.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include "G4EMDataSet.hh"
#include "G4LogLogInterpolation.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ecpssrFormFactorKxsModel::G4ecpssrFormFactorKxsModel()
{ 
  interpolation = new G4LogLogInterpolation();

  for (G4int i=3; i<93; i++) 
  {
      protonDataSetMap[i] = new G4EMDataSet(i,interpolation);
      protonDataSetMap[i]->LoadData("pixe/ecpssr/proton/k-i01m001c01-");
  }

  for (G4int i=3; i<93; i++) 
  {
      alphaDataSetMap[i] = new G4EMDataSet(i,interpolation);
      alphaDataSetMap[i]->LoadData("pixe/ecpssr/alpha/k-i02m004c02-");
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ecpssrFormFactorKxsModel::~G4ecpssrFormFactorKxsModel()
{ 
  protonDataSetMap.clear();
  alphaDataSetMap.clear();
  delete interpolation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ecpssrFormFactorKxsModel::CalculateCrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)
{
  G4Proton* aProton = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();  
  G4double sigma = 0;

  if (energyIncident > 0.1*MeV && energyIncident < 100.*MeV && zTarget < 93 && zTarget > 2) {

    if (massIncident == aProton->GetPDGMass())
      {      	
	sigma = protonDataSetMap[zTarget]->FindValue(energyIncident/MeV);  
	if (sigma !=0 && energyIncident > protonDataSetMap[zTarget]->GetEnergies(0).back()*MeV) return 0.;
      }
    else if (massIncident == aAlpha->GetPDGMass())
      {
	sigma = alphaDataSetMap[zTarget]->FindValue(energyIncident/MeV); 
	if (sigma !=0 && energyIncident > alphaDataSetMap[zTarget]->GetEnergies(0).back()*MeV) return 0.;
      }
    else
      { 
	sigma = 0.;
      }
  }
  
  // sigma is in internal units: it has been converted from 
  // the input file in barns bt the EmDataset
  return sigma;
}
