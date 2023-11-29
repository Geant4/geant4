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
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// Author: Mathieu Karamitros

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, so do not hesitate to
// send us your feedback!
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial.
// If you use Geant4-DNA chemistry and you publish papers about your software,
// in addition to the general paper on Geant4-DNA:
//
// The Geant4-DNA project,
// S. Incerti et al., Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we ask that you please cite the following reference papers on chemistry:
//
// Diffusion-controlled reactions modelling in Geant4-DNA,
// M. Karamitros et al., 2014
// Modeling Radiation Chemistry in the Geant4 Toolkit,
// M. Karamitros et al., Prog. Nucl. Sci. Tec. 2 (2011) 503-508


#ifndef G4VUSERITACTION_H
#define G4VUSERITACTION_H

#include "globals.hh"
#include "G4Track.hh"
#include <vector>

/**
 * G4UserTimeStepAction is used by G4Scheduler.
 * - StartProcessing called before processing
 * - TimeStepAction called at every global step
 * - UserReactionAction called when a reaction occurs
 * - EndProcessing called after processing
 */

class G4UserTimeStepAction
{
public:
  G4UserTimeStepAction();
  G4UserTimeStepAction(const G4UserTimeStepAction& );
  virtual ~G4UserTimeStepAction();

  virtual void StartProcessing(){;}
  virtual void NewStage(){;}

  /** In this method, the user can use :
   * G4Scheduler::Instance()->GetGlobalTime(), to know the current simulation time
   * G4Scheduler::Instance()->GetTimeStep(), to know the selected minimum time
   * WARNING : The call of this method happens before the call of DoIT methods
   */
  virtual void UserPreTimeStepAction(){;}
  virtual void UserPostTimeStepAction(){;}

  /**
   * Inform about a reaction
   */
  virtual void UserReactionAction(const G4Track& /*trackA*/,
      const G4Track& /*trackB*/,
      const std::vector<G4Track*>* /*products*/){;}
  virtual void EndProcessing(){;}

protected:
  void SetMinimumTimeSteps(std::map<double,double>*);
  void AddTimeStep(double /*startingTime*/, double /*timeStep*/);

private:
  G4UserTimeStepAction& operator=(const G4UserTimeStepAction& );
};

#endif // G4VUSERITACTION_H
