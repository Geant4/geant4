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
/// \file electromagnetic/TestEm5/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
// $Id: HistoManager.cc 104417 2017-05-30 08:30:48Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("testem5")
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
  analysisManager->SetActivation(true);    // enable inactivation of histograms

  // Define histograms start values
  const G4int kMaxHisto = 62;
  const G4String id[] = { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                         "10","11","12","13","14","15","16","17","18","19",
                         "20","21","22","23","24","25","26","27","28","29",
                         "30","31","32","33","34","35","36","37","38","39",
                         "40","41","42","43","44","45","46","47","48","49",
                         "50","51","52","53","54","55","56","57","58","59",
                         "60","61"
                        };
                        
  const G4String title[] =
                { "dummy",                                                //0
                  "energy deposit in absorber: dN/dE",                    //1
                  "energy of charged secondaries at creation",            //2
                  "energy of neutral secondaries at creation",            //3
                  "energy of charged at creation (log scale)",            //4
                  "energy of neutral at creation (log scale)",            //5
                  "x_vertex of charged secondaries (all)",                //6
                  "x_vertex of charged secondaries (not absorbed)",       //7
                  "dummy","dummy",                                        //8-9
                  "(transmit, charged) : kinetic energy at exit: dN/dE",  //10
                  "(transmit, charged) : ener fluence: dE(MeV)/dOmega",   //11
                  "(transmit, charged) : space angle: dN/dOmega",         //12
                  "(transmit, charged) : projected angle at exit",        //13
                  "(transmit, charged) : projected position at exit",     //14
                  "(transmit, charged) : radius at exit",                 //15
                  "dummy",                                                //16
                  "dummy",                                                //17
                  "dummy",                                                //18
                  "dummy",                                                //19
                  "(transmit, neutral) : kinetic energy at exit",         //20
                  "(transmit, neutral) : ener fluence: dE(MeV)/dOmega",   //21
                  "(transmit, neutral) : space angle: dN/dOmega",         //22
                  "(transmit, neutral) : projected angle at exit",        //23
                  "dummy","dummy","dummy","dummy","dummy","dummy",       //24-29
                  "(reflect , charged) : kinetic energy at exit",         //30
                  "(reflect , charged) : ener fluence: dE(MeV)/dOmega",   //31
                  "(reflect , charged) : space angle: dN/dOmega",         //32
                  "(reflect , charged) : projected angle at exit",        //33
                  "dummy","dummy","dummy","dummy","dummy","dummy",       //34-39
                  "(reflect , neutral) : kinetic energy at exit",         //40
                  "(reflect , neutral) : ener fluence: dE(MeV)/dOmega",   //41
                  "(reflect , neutral) : space angle: dN/dOmega",         //42
                  "(reflect , neutral) : projected angle at exit",        //43
                  "dummy",                                                //44
                  "dummy",                                                //45
                  "dummy",                                                //46
                  "dummy",                                                //47
                  "dummy",                                                //48
                  "dummy",                                                //49
                  "energy of Auger e- at creation",                       //50
                  "energy of fluorescence gamma at creation",             //51
                  "energy of Auger e- at creation (log scale)",           //52
                  "energy of fluorescence gamma at creation (log scale)", //53
                  "energy of PIXE Auger e- at creation",                  //54
                  "energy of PIXE gamma at creation",                     //55
                  "energy of PIXE Auger e- at creation (log scale)",      //56
                  "energy of PIXE gamma at creation (log scale)",         //57
                  "energy of DNA Auger e- at creation",                   //58
                  "energy of DNA gamma at creation",                      //59
                  "energy of DNA Auger e- at creation (log scale)",       //60
                  "energy of DNA gamma at creation (log scale)"           //61
                 };

  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1("h"+id[k], title[k], nbins,vmin,vmax);
    analysisManager->SetH1Activation(ih, false);
  }
}
