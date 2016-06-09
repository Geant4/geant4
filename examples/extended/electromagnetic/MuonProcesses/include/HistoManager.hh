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
// $Id: HistoManager.hh,v 1.2 2004/06/30 15:48:56 maire Exp $
// GEANT4 tag $Name: geant4-06-02-patch-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace AIDA {
 class ITree;
 class IHistogramFactory;
 class IHistogram1D;
}

class HistoMessenger;

const G4int MaxHisto = 5;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager
{
  public:

    HistoManager();
   ~HistoManager();

    void SetFileName (const G4String& name) { fileName = name;};
    void SetFileType (const G4String& name) { fileType = name;};
    void book();
    void save();
    void SetHisto (G4int,G4int,G4double,G4double,const G4String& unit="none");  
    void FillHisto(G4int id, G4double e, G4double weight = 1.0);
    void RemoveHisto (G4int);
    
    AIDA::IHistogramFactory* GetHistogramFactory()  {return hf;}
    AIDA::IHistogram1D*      GetHisto(G4int id)     {return histo[id];}
                   
    G4bool    HistoExist  (G4int id) {return exist[id];} 
    G4String  GetTitle    (G4int id) {return Title[id];}
    G4int     GetNbins    (G4int id) {return Nbins[id];}
    G4double  GetVmin     (G4int id) {return Vmin[id];}
    G4double  GetVmax     (G4int id) {return Vmax[id];}
    G4double  GetHistoUnit(G4int id) {return Unit[id];}
    G4double  GetBinWidth (G4int id) {return Width[id];}

  private:

    G4String                 fileName;
    G4String                 fileType;
    AIDA::ITree*             tree;
    AIDA::IHistogramFactory* hf;
    AIDA::IHistogram1D*      histo[MaxHisto];
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

#endif

