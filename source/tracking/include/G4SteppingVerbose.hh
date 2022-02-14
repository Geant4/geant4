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
// G4SteppingVerbose
//
// Class description:
//
// This class manages the verbose outputs in G4SteppingManager. 

// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//---------------------------------------------------------------------
#ifndef G4SteppingVerose_hh
#define G4SteppingVerose_hh 1

#include "G4VSteppingVerbose.hh"

class G4SteppingVerbose : public G4VSteppingVerbose
{
  public:

    // Constructor/Destructor

    G4SteppingVerbose();
    virtual ~G4SteppingVerbose();

    virtual G4VSteppingVerbose* Clone()
    { return new G4SteppingVerbose; }

    // Methods to be invoked in the SteppingManager

    virtual void NewStep();
    virtual void AtRestDoItInvoked();
    virtual void AlongStepDoItAllDone();
    virtual void PostStepDoItAllDone();
    virtual void AlongStepDoItOneByOne();
    virtual void PostStepDoItOneByOne();
    virtual void StepInfo();
    virtual void TrackingStarted();
    virtual void DPSLStarted();
    virtual void DPSLUserLimit();
    virtual void DPSLPostStep();
    virtual void DPSLAlongStep();
    virtual void VerboseTrack();
    virtual void VerboseParticleChange();
    virtual void ShowStep() const;

  private:
    static G4int useBestUnitPrecision;

  public:
    static void UseBestUnit(G4int prec = 4);
    static G4int BestUnitPrecision();

};

#endif
