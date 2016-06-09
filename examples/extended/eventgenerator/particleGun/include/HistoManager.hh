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
/// \file eventgenerator/particleGun/include/HistoManager.hh
/// \brief Definition of the HistoManager class
//
//
// $Id$
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include "g4root.hh"
////#include "g4xml.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoMessenger;

const G4int MaxHisto = 9;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager
{
  public:

    HistoManager();
   ~HistoManager();

    void SetFileName   (const G4String& name) { fileName[0] = name;};
    void book();
    void save();
    void SetHisto (G4int,G4int,G4double,G4double,const G4String& unit="none");  
    void FillHisto(G4int id, G4double e, G4double weight = 1.0);
    void Normalize(G4int id, G4double fac);    
    void PrintHisto  (G4int);
    
    G4bool    HistoExist  (G4int id) {return fExist[id];}
    G4double  GetHistoUnit(G4int id) {return fUnit[id];}
    G4double  GetBinWidth (G4int id) {return fWidth[id];}

  private:

    G4String         fileName[2];
    G4bool           factoryOn;

    G4int            fNbHist;
    G4int            fHistId[MaxHisto];
    G4AnaH1*         fHistPt[MaxHisto];
    G4bool           fExist[MaxHisto];
    G4String         fLabel[MaxHisto];
    G4String         fTitle[MaxHisto];
    G4int            fNbins[MaxHisto];
    G4double         fVmin [MaxHisto];
    G4double         fVmax [MaxHisto];
    G4double         fUnit [MaxHisto];
    G4double         fWidth[MaxHisto];
    G4bool           fAscii[MaxHisto];
        
    HistoMessenger*  fHistoMessenger;
    
  private:
    void saveAscii();            
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

