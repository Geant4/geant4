// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// ********************************************************************

#ifndef LISAStackingAction_h
#define LISAStackingAction_h 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class LISAStackingActionMessenger;
class G4Track;

class LISAStackingAction : public G4UserStackingAction {

  public:
    LISAStackingAction();
    virtual ~LISAStackingAction();

    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    virtual void NewStage();
    virtual void PrepareNewEvent();

  public:
    inline void SetPrimarySurvey (const G4bool val) {PrimarySurveyFlag  = val;}
    inline void SetParticleSurvey(const G4bool val) {ParticleSurveyFlag = val;}

  private:
    LISAStackingActionMessenger* stackingMessenger;
    G4bool PrimarySurveyFlag;
    G4bool ParticleSurveyFlag;


};

#endif

