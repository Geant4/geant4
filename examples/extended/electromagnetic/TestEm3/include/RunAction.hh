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
// $Id: RunAction.hh,v 1.7 2004-01-16 08:20:08 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "DetectorConstruction.hh"

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>

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

  void SetHisto (G4int, G4int, G4double, G4double, G4String);

  void SetFileName(const G4String& s) {filename = s;};

  void PrintDedxTables();
    
private:
    
  std::vector<G4double> sumEAbs; 
  std::vector<G4double> sum2EAbs; 
  std::vector<G4double> sumLAbs; 
  std::vector<G4double> sum2LAbs;
  std::vector<G4double> sumEleav; 
  std::vector<G4double> sum2Eleav;       

  DetectorConstruction* Detector;    
  RunActionMessenger*   runMessenger;        
  G4String              filename;
  std::vector<G4String> hid;    
  std::vector<G4String> htitle;    
  std::vector<G4int>    hbins;    
  std::vector<G4double> hmin;    
  std::vector<G4double> hmax;    
  std::vector<G4double> histoUnit;
  G4int  nmax;    
    
#ifdef G4ANALYSIS_USE    
  AIDA::ITree* tree;
  AIDA::IHistogramFactory* hf;
  std::vector<AIDA::IHistogram1D*> histo;
#endif      
             
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

