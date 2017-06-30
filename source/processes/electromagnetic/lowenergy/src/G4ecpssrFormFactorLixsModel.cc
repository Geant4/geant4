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

#include "G4ecpssrFormFactorLixsModel.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include "G4EMDataSet.hh"
#include "G4LinInterpolation.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ecpssrFormFactorLixsModel::G4ecpssrFormFactorLixsModel()
{ 
  interpolation = new G4LinInterpolation();

  for (G4int i=11; i<93; i++) 
  {
      protonL1DataSetMap[i] = new G4EMDataSet(i,interpolation);
      protonL1DataSetMap[i]->LoadData("pixe/ecpssr/proton/l1-i01m001c01-");

      protonL2DataSetMap[i] = new G4EMDataSet(i,interpolation);
      protonL2DataSetMap[i]->LoadData("pixe/ecpssr/proton/l2-i01m001c01-");

      protonL3DataSetMap[i] = new G4EMDataSet(i,interpolation);
      protonL3DataSetMap[i]->LoadData("pixe/ecpssr/proton/l3-i01m001c01-");
  }

  for (G4int i=11; i<93; i++) 
  {
      alphaL1DataSetMap[i] = new G4EMDataSet(i,interpolation);
      alphaL1DataSetMap[i]->LoadData("pixe/ecpssr/alpha/l1-i02m004c02-");
	     
      alphaL2DataSetMap[i] = new G4EMDataSet(i,interpolation);
      alphaL2DataSetMap[i]->LoadData("pixe/ecpssr/alpha/l2-i02m004c02-");
	     
      alphaL3DataSetMap[i] = new G4EMDataSet(i,interpolation);
      alphaL3DataSetMap[i]->LoadData("pixe/ecpssr/alpha/l3-i02m004c02-");
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ecpssrFormFactorLixsModel::~G4ecpssrFormFactorLixsModel()
{ 
  protonL1DataSetMap.clear();
  alphaL1DataSetMap.clear();
  
  protonL2DataSetMap.clear();
  alphaL2DataSetMap.clear();
  
  protonL3DataSetMap.clear();
  alphaL3DataSetMap.clear();
  
  delete interpolation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ecpssrFormFactorLixsModel::CalculateL1CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)
{
  G4Proton* aProton = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();  
  G4double sigma = 0;

  if (energyIncident > 0.1*MeV && energyIncident < 100.*MeV && zTarget < 93 && zTarget > 10) {

    if (massIncident == aProton->GetPDGMass())
      {      
	sigma = protonL1DataSetMap[zTarget]->FindValue(energyIncident/MeV);  
        if (sigma !=0 && energyIncident > protonL1DataSetMap[zTarget]->GetEnergies(0).back()*MeV) return 0.;
      }
    else if (massIncident == aAlpha->GetPDGMass())
      {
        sigma = alphaL1DataSetMap[zTarget]->FindValue(energyIncident/MeV); 
        if (sigma !=0 && energyIncident > alphaL1DataSetMap[zTarget]->GetEnergies(0).back()*MeV) return 0.;
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ecpssrFormFactorLixsModel::CalculateL2CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)
{
  G4Proton* aProton = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();  
  G4double sigma = 0;

  if (energyIncident > 0.1*MeV && energyIncident < 100.*MeV && zTarget < 93 && zTarget > 10) {

    if (massIncident == aProton->GetPDGMass())
      {      
	sigma = protonL2DataSetMap[zTarget]->FindValue(energyIncident/MeV);  
        if (sigma !=0 && energyIncident > protonL2DataSetMap[zTarget]->GetEnergies(0).back()*MeV) return 0.;
      }
    else if (massIncident == aAlpha->GetPDGMass())
      {
        sigma = alphaL2DataSetMap[zTarget]->FindValue(energyIncident/MeV); 
        if (sigma !=0 && energyIncident > alphaL2DataSetMap[zTarget]->GetEnergies(0).back()*MeV) return 0.;
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ecpssrFormFactorLixsModel::CalculateL3CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)
{
  G4Proton* aProton = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();  
  G4double sigma = 0;

  if (energyIncident > 0.1*MeV && energyIncident < 100.*MeV && zTarget < 93 && zTarget > 10) {

    if (massIncident == aProton->GetPDGMass())
      {      
	sigma = protonL3DataSetMap[zTarget]->FindValue(energyIncident/MeV);  
        if (sigma !=0 && energyIncident > protonL3DataSetMap[zTarget]->GetEnergies(0).back()*MeV) return 0.;
      }
    else if (massIncident == aAlpha->GetPDGMass())
      {
        sigma = alphaL3DataSetMap[zTarget]->FindValue(energyIncident/MeV); 
        if (sigma !=0 && energyIncident > alphaL3DataSetMap[zTarget]->GetEnergies(0).back()*MeV) return 0.;
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
