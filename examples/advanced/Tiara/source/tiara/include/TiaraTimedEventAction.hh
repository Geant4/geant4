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
// $Id: TiaraTimedEventAction.hh,v 1.1.1.1 2003/06/12 13:08:24 dressel Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// ----------------------------------------------------------------------
// Class TiaraTimedEventAction
//
// Class description:
//
// Event action used for visualization.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------


#ifndef TiaraTimedEventAction_h
#define TiaraTimedEventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"

class TiaraCellScorerStore;

class TiaraTimedEventAction : public G4UserEventAction
{
public:
  TiaraTimedEventAction(G4int time);
  virtual ~TiaraTimedEventAction();

  virtual void   BeginOfEventAction(const G4Event*);
  virtual void   EndOfEventAction(const G4Event*);
  void SetScorerStore(TiaraCellScorerStore *scorerStore);
  G4int GetTotalProcessedTime() const;
  void SetRnadomNumFilename(const G4String &fname);

private:
  TiaraTimedEventAction(const TiaraTimedEventAction &);
  TiaraTimedEventAction &operator=(const TiaraTimedEventAction &);
  TiaraCellScorerStore *fScorerStore;
  G4int fEvStartTime;
  G4int fCurrentRunTime;
  G4int fMaxRunTime;
  G4int fTimeFromPreviousRuns;
  G4String fFilename;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



#endif

    
