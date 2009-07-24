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
// $Id: HistoManager.hh,v 1.3 2009-07-24 13:02:02 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace AIDA {
 class IAnalysisFactory;
 class ITree;
 class IHistogram1D;
 class ITuple;
} 

class HistoMessenger;

  const G4int MaxHisto = 6; const G4int MaxNtupl = 2;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager
{
  public:
  
    HistoManager();
   ~HistoManager();
   
    void SetFileName   (const G4String& name) { fileName[0] = name;};
    void SetFileType   (const G4String& name) { fileType    = name;};
    void SetFileOption (const G4String& name) { fileOption  = name;};        
    void book();
    void save();
    void SetHisto (G4int,G4int,G4double,G4double,const G4String& unit="none");  
    void FillHisto(G4int id, G4double bin, G4double weight = 1.0);
    void Normalize(G4int id, G4double fac);    
    void RemoveHisto(G4int);   
    void PrintHisto (G4int);
    
    void SetNtuple (G4int);      
    void FillNtuple(G4int id, G4int column, G4double value);
    void AddRowNtuple(G4int id);
                
    G4bool    HistoExist  (G4int id) {return exist[id];}    
    G4double  GetHistoUnit(G4int id) {return Unit[id];}
    G4double  GetBinWidth (G4int id) {return Width[id];}    
    G4String  GetFileName ()         {return fileName[0];}
        
  private:
  
    G4String                 fileName[2];
    G4String                 fileType;
    G4String                 fileOption;
    G4bool                   factoryOn;                    
    AIDA::IAnalysisFactory*  af;        
    AIDA::ITree*             tree;
    
    AIDA::IHistogram1D*      histo[MaxHisto];    
    G4bool                   exist[MaxHisto];
    G4String                 Label[MaxHisto];
    G4String                 Title[MaxHisto];
    G4int                    Nbins[MaxHisto];
    G4double                 Vmin [MaxHisto];
    G4double                 Vmax [MaxHisto];        
    G4double                 Unit [MaxHisto];
    G4double                 Width[MaxHisto];
    G4bool                   ascii[MaxHisto];
        
    AIDA::ITuple*      	        ntupl[MaxNtupl];    
    G4bool                    existNt[MaxNtupl];
    G4String                  LabelNt[MaxNtupl];
    G4String                  TitleNt[MaxNtupl];
    G4String                 ColumnNt[MaxNtupl];    

    HistoMessenger*          histoMessenger;
    
  private:
    void saveAscii();         
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

