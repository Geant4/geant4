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
// $Id: TiaraVisEventAction.hh,v 1.1.1.1 2003/06/12 13:08:25 dressel Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// ----------------------------------------------------------------------
// Class TiaraVisEventAction
//
// Class description:
//
// Event action used for visualization.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------


#ifndef TiaraVisEventAction_h
#define TiaraVisEventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"

class G4Trajectory;
class TiaraCellScorerStore;

class TiaraVisEventAction : public G4UserEventAction
{
public:
  TiaraVisEventAction();
  virtual ~TiaraVisEventAction();

  virtual void   BeginOfEventAction(const G4Event*);
  virtual void   EndOfEventAction(const G4Event*);
  virtual void Clear();
  void DrawTrajectory(G4Trajectory &trj) const;
  void SetScorerStore(TiaraCellScorerStore *scorerStore);
private:
  TiaraVisEventAction(const TiaraVisEventAction &);
  TiaraVisEventAction &operator=(const TiaraVisEventAction &);
  TiaraCellScorerStore *fScorerStore;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



#endif

    
