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

#ifndef ExN07StackingAction_H
#define ExN07StackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class G4Track;

class ExN07StackingAction : public G4UserStackingAction
{
  public:
    ExN07StackingAction();
    virtual ~ExN07StackingAction();

  public:
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    virtual void PrepareNewEvent();

  private:
    static G4int nGamma[6];
    static G4int nElectron[6];
    static G4int nPositron[6];
    static G4double eMinGamma[6];
    static G4double eMinElectron[6];
    static G4double eMinPositron[6];

  public:
    static G4int GetNGamma(G4int i) { return nGamma[i]; }
    static G4int GetNElectron(G4int i) { return nElectron[i]; }
    static G4int GetNPositron(G4int i) { return nPositron[i]; }
    static G4double GetEMinGamma(G4int i) { return eMinGamma[i]; }
    static G4double GetEMinElectron(G4int i) { return eMinElectron[i]; }
    static G4double GetEMinPositron(G4int i) { return eMinPositron[i]; }
};

#endif

