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
// $Id: RunAction.hh,v 1.8 2004/01/21 17:29:26 maire Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "DetectorConstruction.hh"

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class RunActionMessenger;

#ifdef G4ANALYSIS_USE
namespace AIDA {
 class ITree;
 class IHistogramFactory;
 class IHistogram1D;
}
#endif

class RunAction : public G4UserRunAction
{
  public:

    RunAction(DetectorConstruction*);
   ~RunAction();

    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);

    void fillPerEvent(G4int,G4double,G4double,G4double);

#ifdef G4ANALYSIS_USE
    AIDA::IHistogram1D* GetHisto(G4int id) {return histo[id];}
    G4double GetHistoUnit(G4int id) {return histoUnit[id];}
#endif

    void SetFileName(const G4String& s) {fileName = s;};
    void SetHisto (G4int, G4int, G4double, G4double, G4String);

    void PrintDedxTables();
    
  private:

    G4double sumEAbs [MaxAbsor], sum2EAbs [MaxAbsor]; 
    G4double sumLAbs [MaxAbsor], sum2LAbs [MaxAbsor];
    G4double sumEleav[MaxAbsor], sum2Eleav[MaxAbsor];           

    DetectorConstruction* Detector;    
    RunActionMessenger*   runMessenger;
            
    G4String fileName;
    G4String hid  [MaxAbsor],  htitle[MaxAbsor];
    G4int    hbins[MaxAbsor];
    G4double hmin [MaxAbsor],  hmax  [MaxAbsor];
    G4double histoUnit[MaxAbsor];
    
#ifdef G4ANALYSIS_USE    
    AIDA::ITree* tree;
    AIDA::IHistogramFactory* hf;
    AIDA::IHistogram1D* histo[MaxAbsor];
#endif      
             
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

