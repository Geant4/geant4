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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef fluoTestEventAction_h
#define fluoTestEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#ifdef G4ANALYSIS_USE
#include "fluoTestAnalysisManager.hh"
#endif


class fluoTestEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class fluoTestEventAction : public G4UserEventAction
{
  public:
#ifdef G4ANALYSIS_USE
  fluoTestEventAction(fluoTestAnalysisManager* analysisMgr);
#else
    fluoTestEventAction();
#endif   
 virtual ~fluoTestEventAction();

  public:
    virtual void   BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
    void SetDrawFlag   (G4String val)  {drawFlag = val;};
    void SetPrintModulo(G4int    val)  {printModulo = val;};
    
  private:
    G4int                       siCollID;
    G4int                       sampleCollID;
    G4int                       hpgeCollID;
    G4String                    drawFlag;
    G4int                       printModulo;                         
    fluoTestEventActionMessenger*  eventMessenger;
#ifdef G4ANALYSIS_USE
    fluoTestAnalysisManager* analysisManager;
#endif
};

#endif

    
