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
// $Id: Em3RunAction.hh,v 1.11 2002-12-12 11:19:37 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Em3RunAction_h
#define Em3RunAction_h 1

#include "Em3DetectorConstruction.hh"

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class Em3RunActionMessenger;

#ifndef G4NOHIST
namespace AIDA {
 class ITree;
 class IHistogramFactory;
 class IHistogram1D;
} 
#endif

class Em3RunAction : public G4UserRunAction
{
  public:
  
    Em3RunAction(Em3DetectorConstruction*);
   ~Em3RunAction();

    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);
    
    void fillPerEvent(G4int,G4double,G4double);
    
#ifndef G4NOHIST
    AIDA::IHistogram1D* GetHisto(G4int id) {return histo[id];}    
#endif

    void SetHisto (G4int id, G4int nbins, G4double vmin, G4double vmax);
    
    void PrintDedxTables();
            
  private:
    void bookHisto();
    
  private:
    
    G4double sumEAbs[MaxAbsor], sum2EAbs[MaxAbsor];
    G4double sumLAbs[MaxAbsor], sum2LAbs[MaxAbsor];

    Em3DetectorConstruction* Detector;    
    Em3RunActionMessenger*   runMessenger;        
    
#ifndef G4NOHIST    
    AIDA::ITree* tree;
    AIDA::IHistogramFactory* hf;    
    AIDA::IHistogram1D* histo[MaxAbsor];
#endif               
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void Em3RunAction::fillPerEvent(G4int kAbs, G4double EAbs, G4double LAbs)
{      
  //accumulate statistic
  //
  sumEAbs[kAbs] += EAbs; sum2EAbs[kAbs] += EAbs*EAbs;
  sumLAbs[kAbs] += LAbs; sum2LAbs[kAbs] += LAbs*LAbs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

