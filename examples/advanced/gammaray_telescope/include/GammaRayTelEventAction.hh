// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelEventAction.hh,v 1.2 2000-11-15 20:27:39 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelEventAction  ------
//           by R.Giannitrapani, F. Longo & G.Santin (13 nov 2000)
//
// ************************************************************


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef GammaRayTelEventAction_h
#define GammaRayTelEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#ifdef G4HIS_USE_AIDA
#include "GammaRayTelHistogram.hh"
#endif
//class GammaRayTelEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelEventAction : public G4UserEventAction
{
  public:
#ifdef G4HIS_USE_AIDA
  GammaRayTelEventAction(GammaRayTelHistogram* histoMgr);
#else
  GammaRayTelEventAction();
#endif
  virtual ~GammaRayTelEventAction();

  public:
  virtual void   BeginOfEventAction(const G4Event*);
  virtual void   EndOfEventAction(const G4Event*);
  
  void SetDrawFlag   (G4String val)  {drawFlag = val;};
  void SetPrintModulo(G4int    val)  {printModulo = val;};

  private:
    G4int       trackerCollID;                
    G4String    drawFlag;
    G4int       printModulo;                         
#ifdef G4HIS_USE_AIDA
    GammaRayTelHistogram* histoManager;
#endif
};

#endif

    




