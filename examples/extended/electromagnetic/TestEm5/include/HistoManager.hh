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
// $Id: HistoManager.hh,v 1.3 2003/11/03 12:58:52 maire Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HistoManager_h
#define HistoManager_h 1

#ifdef G4ANALYSIS_USE

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace AIDA {
 class ITree;
 class IHistogramFactory;
 class IHistogram1D;
} 

class HistoMessenger;

  const G4int MaxHisto = 17;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager
{
  public:
    HistoManager();
   ~HistoManager();
   
    void SetFactory(G4String);
    void SetHisto  (G4int, G4int, G4double, G4double, G4String unit="none");
    
    AIDA::ITree*             GetTree()              {return tree;}
    AIDA::IHistogramFactory* GetHistogramFactory()  {return hf;}        
    AIDA::IHistogram1D*      GetHisto(G4int id)     {return histo[id];}
    G4double                 GetHistoUnit(G4int id) {return histoUnit[id];}
    G4double                 GetBinWidth (G4int id) {return binWidth[id];}
  private:
    AIDA::ITree*             tree;
    AIDA::IHistogramFactory* hf;    
    AIDA::IHistogram1D*      histo[MaxHisto];
    G4double                 histoUnit[MaxHisto];
    G4double                 binWidth[MaxHisto];
    G4bool                   factoryOn;        
    HistoMessenger*          histoMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
#endif

