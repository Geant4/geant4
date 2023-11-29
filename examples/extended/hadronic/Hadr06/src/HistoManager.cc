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
/// \file HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);     //enable inactivation of histograms
  
  // Define histograms start values
  const G4int kMaxHisto = 25;
  const G4String id[] = {"0","1","2","3","4","5","6","7","8","9",
                         "10","11","12","13","14","15","16","17","18","19",
			 "20","21","22","23","24" };
  const G4String title[] = 
       { "dummy",                                                        //0
         "total energy deposit",                                         //1
         "Edep (MeV/mm) along absorber",                                 //2
         "total kinetic energy flow",                                    //3
         "energy spectrum of gamma at creation",                         //4
         "energy spectrum of e+- at creation",                           //5
         "energy spectrum of neutrons at creation",                      //6
         "energy spectrum of protons at creation",                       //7
         "energy spectrum of deuterons at creation",                     //8
         "energy spectrum of alphas at creation",                        //9
         "energy spectrum of all others ions at creation",               //10
         "energy spectrum of all others baryons at creation",            //11
         "energy spectrum of all others mesons at creation",             //12
         "energy spectrum of all others leptons (neutrinos) at creation", //13
         "energy spectrum of emerging gamma",                            //14
         "energy spectrum of emerging e+-",                              //15
         "energy spectrum of emerging neutrons",                         //16
         "energy spectrum of emerging protons",                          //17
         "energy spectrum of emerging deuterons",                        //18
         "energy spectrum of emerging alphas",                           //19
         "energy spectrum of all others emerging ions",                  //20
         "energy spectrum of all others emerging baryons",               //21
         "energy spectrum of all others emerging mesons",                //22
         "energy spectrum of all others emerging leptons (neutrinos)",   //23
         "total energy released : edep + eflow"                          //24 
       };

  // Default values (to be reset via /analysis/h1/set command) 
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
