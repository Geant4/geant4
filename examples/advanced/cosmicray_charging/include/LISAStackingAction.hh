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

