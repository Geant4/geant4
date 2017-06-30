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
// $Id: HistoManager.cc 67909 2013-03-12 18:51:09Z vnivanch $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("Activation")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);     //enable inactivation of histograms
  
  // Define histograms start values
  ///const G4int kMaxHisto1 = 24, kMaxHisto2 = 44;
  const G4String id[] = {"0","1","2","3","4","5","6","7","8","9",
                         "10","11","12","13",
                         "14","15","16","17","18","19","20","21","22","23",
                         "24","25","26","27","28","29","30","31","32","33",
                         "34","35","36","37","38","39","40","41","42","43" };

  const G4String title[] = 
       { "dummy",                                                        //0
         "total energy deposit",                                         //1
         "Edep (MeV/mm) along beam directiom",                           //2
         "total kinetic energy emerging",                                //3
         "energy spectrum of emerging gamma",                            //4
         "energy spectrum of emerging e+-",                              //5
         "energy spectrum of emerging neutrons",                         //6
         "energy spectrum of emerging protons",                          //7
         "energy spectrum of emerging deuterons",                        //8
         "energy spectrum of emerging alphas",                           //9
         "energy spectrum of all others emerging ions",                  //10
         "energy spectrum of all others emerging baryons",               //11
         "energy spectrum of all others emerging mesons",                //12
         "energy spectrum of all others emerging leptons (neutrinos)",   //13
         "dN/dt (becquerel) of emerging gamma",                          //14
         "dN/dt (becquerel) of emerging e+- ",                           //15
         "dN/dt (becquerel) of emerging neutrons",                       //16
         "dN/dt (becquerel) of emerging protons",                        //17
         "dN/dt (becquerel) of emerging deuterons",                      //18
         "dN/dt (becquerel) of emerging alphas",                         //19
         "dN/dt (becquerel) of all others emerging ions",                //20
         "dN/dt (becquerel) of all others emerging baryons",             //21
         "dN/dt (becquerel) of all others emerging mesons",              //22
         "dN/dt (becquerel) of all others emerging leptons (neutrinos)"  //23
       };

  // Default values (to be reset via /analysis/h1/set command) 
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto1; k++) {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }
  
  // population of species
  for (G4int k=kMaxHisto1; k<kMaxHisto2; k++) {
    G4int ih = analysisManager->CreateH1(id[k], " population",nbins,vmin,vmax);
    analysisManager->SetH1Activation(ih, false);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
