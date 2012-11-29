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
/// \file medical/fanoCavity2/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
// $Id$
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "HistoManager.hh"
#include "HistoMessenger.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
{
  fileName[0]  = "fanocavity2";
  factoryOn = false;
  fNbHist   = 0;

  // histograms
  for (G4int k=0; k<MaxHisto; k++) {
    fHistId[k] = 0;
    fHistPt[k] = 0;    
    fExist[k] = false;
    fUnit[k]  = 1.0;
    fWidth[k] = 1.0;
    fAscii[k] = false;    
  }

  fHistoMessenger = new HistoMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete fHistoMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{
  // if no histos, do nothing
  if (fNbHist == 0) return;
  
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  G4String extension = analysisManager->GetFileType();
  fileName[1] = fileName[0] + "." + extension;
  
  // Open an output file
  //
  G4bool fileOpen = analysisManager->OpenFile(fileName[0]);
  if (!fileOpen) {
    G4cout << "\n---> HistoManager::book(): cannot open " << fileName[1] 
           << G4endl;
    return;
  }  

  // create selected histograms
  //
  analysisManager->SetFirstHistoId(1);
    
  for (G4int k=0; k<MaxHisto; k++) {
    if (fExist[k]) {
      fHistId[k] = analysisManager->CreateH1( fLabel[k], fTitle[k],
                                              fNbins[k], fVmin[k], fVmax[k]);
      fHistPt[k] = analysisManager->GetH1(fHistId[k]);
      factoryOn = true;
    }
  }

  if (factoryOn)  
    G4cout << "\n----> Histogram file is opened in " << fileName[1] << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{
  if (factoryOn) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
    analysisManager->Write();
    analysisManager->CloseFile();
    saveAscii();                    // Write fAscii file, if any
    G4cout << "\n----> Histograms are saved in " << fileName[1] << G4endl;
      
    delete G4AnalysisManager::Instance();
    factoryOn = false;
  }         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double e, G4double weight)
{
  if (ih > MaxHisto) {
    G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
           << "does not fExist; e= " << e << " w= " << weight << G4endl;
    return;
  }

  if (fHistPt[ih]) fHistPt[ih]->fill(e/fUnit[ih], weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetHisto(G4int ih,
               G4int nbins, G4double vmin, G4double vmax, const G4String& unit)
{
  if (ih > MaxHisto) {
    G4cout << "---> warning from HistoManager::SetHisto() : histo " << ih
           << "does not fExist" << G4endl;
    return;
  }
  
  const G4String id[] = {"0","1","2","3","4","5","6","7","8","9" ,"10", "11"};
  
  const G4String title[] =
                { "dummy",                                                //0
                  "emission point of e-",                                //1
                  "energy spectrum of e-",                                 //2
                  "theta distribution of e-",                                //3
                  "emission point of e- hitting cavity",                //4
                  "energy spectrum of e- when entering in cavity",        //5
                  "theta distribution of e- before enter in cavity",        //6
                  "theta distribution of e- at first step in cavity",        //7
                  "track segment of e- in cavity",                        //8
                  "step size of e- in wall",                                //9
                  "step size of e- in cavity",                                //10
                  "energy deposit in cavity per track"                        //11
                 };  

  G4String titl = title[ih];
  fUnit[ih] = 1.;

  if (unit != "none") {
    titl = title[ih] + " (" + unit + ")";
    fUnit[ih] = G4UnitDefinition::GetValueOf(unit);
  }

  fExist[ih] = true;
  fLabel[ih] = id[ih];
  fTitle[ih] = titl;
  fNbins[ih] = nbins;
  fVmin[ih]  = vmin;
  fVmax[ih]  = vmax;
  fWidth[ih] = fUnit[ih]*(vmax-vmin)/nbins;
  
  fNbHist++;

  G4cout << "----> SetHisto " << ih << ": " << titl << ";  "
         << nbins << " bins from "
         << vmin << " " << unit << " to " << vmax << " " << unit << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double fac)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::Normalize() : histo " << ih
           << "  fac= " << fac << G4endl;
    return;
  }

  if (fHistPt[ih]) fHistPt[ih]->scale(fac);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::PrintHisto(G4int ih)
{
 if (ih < MaxHisto) { fAscii[ih] = true; fAscii[0] = true; }
 else
    G4cout << "---> warning from HistoManager::PrintHisto() : histo " << ih
           << "does not exist" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <fstream>

void HistoManager::saveAscii()
{
 if (!fAscii[0]) return;
 
 G4String name = fileName[0] + ".ascii";
 std::ofstream File(name, std::ios::out);
 if (!File) {
   G4cout 
       << "\n---> HistoManager::saveAscii(): cannot open " << name << G4endl;
    return;
  }  

 File.setf( std::ios::scientific, std::ios::floatfield );
  
 //write selected histograms
 for (G4int ih=0; ih<MaxHisto; ih++) {
    if (fHistPt[ih] && fAscii[ih]) {
          
      File << "\n  1D histogram " << ih << ": " << fTitle[ih] 
           << "\n \n \t     X \t\t     Y" << G4endl;
     
      for (G4int iBin=0; iBin<fNbins[ih]; iBin++) {
         File << "  " << iBin << "\t" 
              << fHistPt[ih]->axis().bin_center(iBin) << "\t"
              << fHistPt[ih]->bin_height(iBin) 
              << G4endl;
      } 
    }
 } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

