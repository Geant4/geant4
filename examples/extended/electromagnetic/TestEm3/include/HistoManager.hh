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
// $Id: HistoManager.hh,v 1.2 2004/06/21 10:52:55 maire Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HistoManager_h
#define HistoManager_h 1

#ifdef G4ANALYSIS_USE

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef USE_AIDA
namespace AIDA {
 class ITree;
 class IHistogramFactory;
 class IHistogram1D;
} 
#endif

class HistoMessenger;

#include "DetectorConstruction.hh"
  const G4int MaxHisto = MaxAbsor;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager
{
  public:
  
    HistoManager();
   ~HistoManager();
   
    void SetFileName (G4String name) { fileName = name;};
    void SetFactory  ();
    void SaveFactory ();    
    void SetHisto (G4int, G4int, G4double, G4double, G4String unit="none");
    void RemoveHisto (G4int);
    
#ifdef USE_AIDA    
    AIDA::ITree*             GetTree()              {return tree;}
    AIDA::IHistogramFactory* GetHistogramFactory()  {return hf;}        
    AIDA::IHistogram1D*      GetHisto(G4int id)     {return histo[id];}
#endif
    
    G4double                 GetHistoUnit(G4int id) {return Unit[id];}
    G4double                 GetBinWidth (G4int id) {return Width[id];}
    
  private:
  
    G4String                 fileName;
    
#ifdef USE_AIDA    
    AIDA::ITree*             tree;
    AIDA::IHistogramFactory* hf;    
    AIDA::IHistogram1D*      histo[MaxHisto];
#endif
    
    G4bool                   exist[MaxHisto];
    G4String                 Label[MaxHisto];
    G4String                 Title[MaxHisto];
    G4int                    Nbins[MaxHisto];
    G4double                 Vmin [MaxHisto];
    G4double                 Vmax [MaxHisto];        
    G4double                 Unit [MaxHisto];
    G4double                 Width[MaxHisto];
    G4bool                   factoryOn;        
    HistoMessenger*          histoMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif     //G4ANALYSIS_USE
#endif

