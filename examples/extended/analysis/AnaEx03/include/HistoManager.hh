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
/// \file analysis/AnaEx03/include/HistoManager.hh
/// \brief Definition of the HistoManager class
//
// $Id: HistoManager.hh 98059 2016-07-01 16:23:33Z gcosmo $
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

const G4int kMaxHisto = 4;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager
{
  public:
    HistoManager();
   ~HistoManager();
   
    void Book();
    void Save();

    void FillHisto(G4int id, G4double bin, G4double weight = 1.0);
    void Normalize(G4int id, G4double fac);    

    void FillNtuple(G4double EnergyAbs, G4double EnergyGap,
                    G4double TrackLAbs, G4double TrackLGap);
    
    void PrintStatistic();
        
  private:
    AIDA::IAnalysisFactory*  fAF;        
    AIDA::ITree*             fTree;
    
    AIDA::IHistogram1D*      fHisto[kMaxHisto];            
    AIDA::ITuple*            fNtuple1;    
    AIDA::ITuple*            fNtuple2;    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

