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
// $Id: HistoManager.hh,v 1.2 2003-10-28 18:29:54 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//#ifdef G4ANALYSIS_USE

namespace AIDA {
 class ITree;
 class IHistogramFactory;
 class IHistogram1D;
} 

//#endif

class HistoMessenger;

const G4int MaxHisto = 17;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager
{
  public:
    HistoManager();
   ~HistoManager();

#ifdef G4ANALYSIS_USE
   
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

    HistoMessenger*          histoMessenger;

#endif

    G4bool                   factoryOn;        
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

