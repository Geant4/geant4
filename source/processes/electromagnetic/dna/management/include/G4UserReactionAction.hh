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
// $Id: G4UserReactionAction.hh 64057 2012-10-30 15:04:49Z gcosmo $
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// Author: Mathieu Karamitros (kara@cenbg.in2p3.fr)
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef G4VUSERITACTION_H
#define G4VUSERITACTION_H

#include "globals.hh"
#include "G4Track.hh"
#include "G4TrackFastVector.hh"

/**
  * G4UserReactionAction is used by G4ITStepManager.
  * - StartProcessing called before processing
  * - TimeStepAction called at every global step
  * - UserReactionAction called when a reaction occurs
  * - EndProcessing called after processing
  */

class G4UserReactionAction
{
public:
   G4UserReactionAction();
   G4UserReactionAction(const G4UserReactionAction& );
   virtual ~G4UserReactionAction();

   virtual void StartProcessing(){;}

   /** In this method, the user can use :
    * G4ITStepManager::Instance()->GetGlobalTime(), to know the current simulation time
    * G4ITStepManager::Instance()->GetMinTime(), to know the selected minimum time
    * WARNING : The call of this method happens before the call of DoIT methods
    */
   virtual void TimeStepAction(){;}

   /**
    * This method enables to kill products right after they are generated
    */
   virtual void UserReactionAction(const G4Track& /*trackA*/,const G4Track& /*trackB*/,
                                   const G4TrackFastVector& /*products*/,
                                   int /*nbProducts*/){;}
   virtual void EndProcessing(){;}

private:
    G4UserReactionAction& operator=(const G4UserReactionAction& );
};

#endif // G4VUSERITACTION_H
