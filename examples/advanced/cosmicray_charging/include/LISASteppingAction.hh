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

#ifndef LISASteppingAction_h
#define LISASteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class LISAEventAction;
class LISASteppingActionMessenger;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class LISASteppingAction : public G4UserSteppingAction {

  public:
    LISASteppingAction();
    virtual ~LISASteppingAction();

    virtual void UserSteppingAction(const G4Step*);

  public:
    inline void Discharge(void)          {charge_in[0]=0; charge_out[0]=0;
                                         charge_in[1]=0; charge_out[1]=0;};

    inline G4int GetChargeIn(const G4int tm)       {return charge_in[tm];}
    inline G4int GetChargeOut(const G4int tm)      {return charge_out[tm];}


    inline void SetFlagSpectrum(const G4bool val)  {FlagSpectrum  = val;}

  private:
    G4int charge_in[2];
    G4int charge_out[2];

    LISASteppingActionMessenger* steppingMessenger;
    G4bool FlagSpectrum;



};

#endif
