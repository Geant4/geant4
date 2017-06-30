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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// and papers
// M. Batmunkh et al. J Radiat Res Appl Sci 8 (2015) 498-507
// O. Belov et al. Physica Medica 32 (2016) 1510-1520 
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// -------------------------------------------------------------------
// November 2016
// -------------------------------------------------------------------
//
// $Id$
//
/// \file TimeStepAction.hh
/// \brief Definition of the TimeStepAction class

#ifndef ITACTION_H
#define ITACTION_H

#include "G4UserTimeStepAction.hh"

class TimeStepAction : public G4UserTimeStepAction
{
public:
  TimeStepAction();
  virtual ~TimeStepAction();
  TimeStepAction(const TimeStepAction& other);
  TimeStepAction& operator=(const TimeStepAction& other);

  virtual void StartProcessing(){;}

  /** In this method, the user can use :
   * G4ITTimeStepper::Instance()->GetGlobalTime(),
   *    to know the current simulation time
   * G4ITTimeStepper::Instance()->GetTimeStep(),
   *    to know the selected minimum time
   * WARNING :
   *    The call of this method happens before the call of DoIT methods
   */
  virtual void UserPreTimeStepAction(){;}
  virtual void UserPostTimeStepAction();

  /**
   * Inform about a reaction
   */
  virtual void UserReactionAction(const G4Track& /*trackA*/,
                                  const G4Track& /*trackB*/,
                                  const std::vector<G4Track*>* /*products*/);

  virtual void EndProcessing(){;}
};

#endif // ITACTION_H
