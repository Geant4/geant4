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
//  EXRS2012 proceedings
// ---------------------------------------------------------------------------------------

#include <fstream>
#include <iomanip>

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "G4EMDataSet.hh"
#include "G4LinInterpolation.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"

#include "G4ecpssrFormFactorMixsModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ecpssrFormFactorMixsModel::G4ecpssrFormFactorMixsModel()
{ 
  interpolation = new G4LinInterpolation();

  for (G4int i=29; i<93; i++) 
  {
      protonM1DataSetMap[i] = new G4EMDataSet(i,interpolation);
      protonM1DataSetMap[i]->LoadData("pixe/ecpssr/proton/m1-i01m001c01-");
	    
      protonM2DataSetMap[i] = new G4EMDataSet(i,interpolation);
      protonM2DataSetMap[i]->LoadData("pixe/ecpssr/proton/m2-i01m001c01-");
	    
      protonM3DataSetMap[i] = new G4EMDataSet(i,interpolation);
      protonM3DataSetMap[i]->LoadData("pixe/ecpssr/proton/m3-i01m001c01-");
	    
      protonM4DataSetMap[i] = new G4EMDataSet(i,interpolation);
      protonM4DataSetMap[i]->LoadData("pixe/ecpssr/proton/m4-i01m001c01-");
	    
      protonM5DataSetMap[i] = new G4EMDataSet(i,interpolation);
      protonM5DataSetMap[i]->LoadData("pixe/ecpssr/proton/m5-i01m001c01-");
  }

  protonMiXsVector.push_back(protonM1DataSetMap);
  protonMiXsVector.push_back(protonM2DataSetMap);
  protonMiXsVector.push_back(protonM3DataSetMap);
  protonMiXsVector.push_back(protonM4DataSetMap);
  protonMiXsVector.push_back(protonM5DataSetMap);


  for (G4int i=29; i<93; i++) 
  {
      alphaM1DataSetMap[i] = new G4EMDataSet(i,interpolation);
      alphaM1DataSetMap[i]->LoadData("pixe/ecpssr/alpha/m1-i02m004c02-");
	     
      alphaM2DataSetMap[i] = new G4EMDataSet(i,interpolation);
      alphaM2DataSetMap[i]->LoadData("pixe/ecpssr/alpha/m2-i02m004c02-");
	     
      alphaM3DataSetMap[i] = new G4EMDataSet(i,interpolation);
      alphaM3DataSetMap[i]->LoadData("pixe/ecpssr/alpha/m3-i02m004c02-");
	     
      alphaM4DataSetMap[i] = new G4EMDataSet(i,interpolation);
      alphaM4DataSetMap[i]->LoadData("pixe/ecpssr/alpha/m4-i02m004c02-");
	     
      alphaM5DataSetMap[i] = new G4EMDataSet(i,interpolation);
      alphaM5DataSetMap[i]->LoadData("pixe/ecpssr/alpha/m5-i02m004c02-");

  }

  alphaMiXsVector.push_back(alphaM1DataSetMap);
  alphaMiXsVector.push_back(alphaM2DataSetMap);
  alphaMiXsVector.push_back(alphaM3DataSetMap);
  alphaMiXsVector.push_back(alphaM4DataSetMap);
  alphaMiXsVector.push_back(alphaM5DataSetMap);



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ecpssrFormFactorMixsModel::~G4ecpssrFormFactorMixsModel()
{ 
  protonM1DataSetMap.clear();
  alphaM1DataSetMap.clear();
  
  protonM2DataSetMap.clear();
  alphaM2DataSetMap.clear();
  
  protonM3DataSetMap.clear();
  alphaM3DataSetMap.clear();
  
  protonM4DataSetMap.clear();
  alphaM4DataSetMap.clear();
  
  protonM5DataSetMap.clear();
  alphaM5DataSetMap.clear();
  
  delete interpolation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ecpssrFormFactorMixsModel::CalculateMiCrossSection(G4int zTarget,G4double massIncident, G4double energyIncident, G4int mShellId)
{
  G4Proton* aProton = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();  
  G4double sigma = 0;
  G4int mShellIndex = mShellId -1;

  if (energyIncident > 0.1*MeV && energyIncident < 100*MeV && zTarget < 93 && zTarget > 28) {

    if (massIncident == aProton->GetPDGMass())
      {      
	sigma = protonMiXsVector[mShellIndex][zTarget]->FindValue(energyIncident/MeV);  
        if (sigma !=0 && energyIncident > protonMiXsVector[mShellIndex][zTarget]->GetEnergies(0).back()*MeV) return 0.;
      }
    else if (massIncident == aAlpha->GetPDGMass())
      {
        sigma = alphaMiXsVector[mShellIndex][zTarget]->FindValue(energyIncident/MeV); 
        if (sigma !=0 && energyIncident > alphaMiXsVector[mShellIndex][zTarget]->GetEnergies(0).back()*MeV) return 0.;
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

G4double G4ecpssrFormFactorMixsModel::CalculateM1CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)
{

  // mShellId
  return  CalculateMiCrossSection (zTarget, massIncident, energyIncident, 1); 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ecpssrFormFactorMixsModel::CalculateM2CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)
{

  // mShellId
  return  CalculateMiCrossSection (zTarget, massIncident, energyIncident, 2); 

  /*

  G4Proton* aProton = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();  
  G4double sigma = 0;

  if (energyIncident > 0.1*MeV && energyIncident < 10*MeV && zTarget < 93 && zTarget > 61) {

    if (massIncident == aProton->GetPDGMass())
      {      
	sigma = protonM2DataSetMap[zTarget]->FindValue(energyIncident/MeV);  
        if (sigma !=0 && energyIncident > protonM2DataSetMap[zTarget]->GetEnergies(0).back()*MeV) return 0.;
      }
    else if (massIncident == aAlpha->GetPDGMass())
      {
        sigma = alphaM2DataSetMap[zTarget]->FindValue(energyIncident/MeV); 
        if (sigma !=0 && energyIncident > alphaM2DataSetMap[zTarget]->GetEnergies(0).back()*MeV) return 0.;
      }
    else
      { 
	sigma = 0.;
      }
  }
  
  // sigma is in internal units: it has been converted from 
  // the input file in barns bt the EmDataset
  return sigma;
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ecpssrFormFactorMixsModel::CalculateM3CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)
{

  return  CalculateMiCrossSection (zTarget, massIncident, energyIncident, 3); 
  /*


  G4Proton* aProton = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();  
  G4double sigma = 0;

  if (energyIncident > 0.1*MeV && energyIncident < 10*MeV && zTarget < 93 && zTarget > 61) {

    if (massIncident == aProton->GetPDGMass())
      {      
	sigma = protonM3DataSetMap[zTarget]->FindValue(energyIncident/MeV);  
        if (sigma !=0 && energyIncident > protonM3DataSetMap[zTarget]->GetEnergies(0).back()*MeV) return 0.;
      }
    else if (massIncident == aAlpha->GetPDGMass())
      {
        sigma = alphaM3DataSetMap[zTarget]->FindValue(energyIncident/MeV); 
        if (sigma !=0 && energyIncident > alphaM3DataSetMap[zTarget]->GetEnergies(0).back()*MeV) return 0.;
      }
    else
      { 
	sigma = 0.;
      }
  }
  
  // sigma is in internal units: it has been converted from 
  // the input file in barns bt the EmDataset
  return sigma;
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ecpssrFormFactorMixsModel::CalculateM4CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)
{

  return  CalculateMiCrossSection (zTarget, massIncident, energyIncident, 4); 
  /*
  G4Proton* aProton = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();  
  G4double sigma = 0;

  if (energyIncident > 0.1*MeV && energyIncident < 10*MeV && zTarget < 93 && zTarget > 61) {

    if (massIncident == aProton->GetPDGMass())
      {      
	sigma = protonM3DataSetMap[zTarget]->FindValue(energyIncident/MeV);  
        if (sigma !=0 && energyIncident > protonM3DataSetMap[zTarget]->GetEnergies(0).back()*MeV) return 0.;
      }
    else if (massIncident == aAlpha->GetPDGMass())
      {
        sigma = alphaM3DataSetMap[zTarget]->FindValue(energyIncident/MeV); 
        if (sigma !=0 && energyIncident > alphaM3DataSetMap[zTarget]->GetEnergies(0).back()*MeV) return 0.;
      }
    else
      { 
	sigma = 0.;
      }
  }
  
  // sigma is in internal units: it has been converted from 
  // the input file in barns bt the EmDataset
  return sigma;
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ecpssrFormFactorMixsModel::CalculateM5CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)
{

  return  CalculateMiCrossSection (zTarget, massIncident, energyIncident, 5); 
  /*
  G4Proton* aProton = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();  
  G4double sigma = 0;

  if (energyIncident > 0.1*MeV && energyIncident < 10*MeV && zTarget < 93 && zTarget > 61) {

    if (massIncident == aProton->GetPDGMass())
      {      
	sigma = protonM3DataSetMap[zTarget]->FindValue(energyIncident/MeV);  
        if (sigma !=0 && energyIncident > protonM3DataSetMap[zTarget]->GetEnergies(0).back()*MeV) return 0.;
      }
    else if (massIncident == aAlpha->GetPDGMass())
      {
        sigma = alphaM3DataSetMap[zTarget]->FindValue(energyIncident/MeV); 
        if (sigma !=0 && energyIncident > alphaM3DataSetMap[zTarget]->GetEnergies(0).back()*MeV) return 0.;
      }
    else
      { 
	sigma = 0.;
      }
  }
  
  // sigma is in internal units: it has been converted from 
  // the input file in barns bt the EmDataset
  return sigma;
  */
}
