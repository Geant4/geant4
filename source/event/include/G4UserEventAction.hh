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
// G4UserEventAction
//
// Class description:
//
// This is the base class of one of the user's optional action classes.
// The two methods BeginOfEventAction() and EndOfEventAction() are invoked
// at the beginning and the end of one event processing. These methods are
// invoked by G4EventManager.
// Be aware that BeginOfEventAction() is invoked when a G4Event object is
// sent to G4EventManager. Thus the primary vertexes/particles have already
// been made by the primary generator. In case the user wants to do something
// before generating primaries (i.e., store random number status), do it in
// the G4VUserPrimaryGeneratorAction concrete class

// Author: Makoto Asai (SLAC)
// Adding MergeSubEvent - Sep/11/2023 Makoto Asai (JLab)
// --------------------------------------------------------------------
#ifndef G4UserEventAction_hh
#define G4UserEventAction_hh 1

class G4EventManager;
class G4Event;

class G4UserEventAction 
{
  public:

    G4UserEventAction();
    virtual ~G4UserEventAction() = default;
    virtual void SetEventManager(G4EventManager* value);

    virtual void BeginOfEventAction(const G4Event* anEvent);
    virtual void EndOfEventAction(const G4Event* anEvent);
      // Two virtual method the user can override.

    virtual void MergeSubEvent(G4Event* masterEvent, const G4Event* subEvent);
      // A virtual method to merge the results of a sub-event into the master
      // event. The ownership of "subEvent" and its contents blong to the
      // worker thread. 
      // Merging trajectories and scores are taken care by G4Event and
      // G4ScoringManager so the user does not need to take care of them.
      // But merging hits collections and UserEventInformation must be taken 
      // care by this method.
      // This method is invoked only for the case of sub-event parallelism.

  protected:

      G4EventManager* fpEventManager = nullptr; // not owned
};

#endif
