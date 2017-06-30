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
/// \file ITTrackingInteractivity.hh
/// \brief Definition of the ITTrackingInteractivity class

#ifndef ITTRACKINGINTERACTIVITY_HH
#define ITTRACKINGINTERACTIVITY_HH

#include "G4ITTrackingInteractivity.hh"
#include <vector>

class G4VTrajectory;

/*
 * This class should be modified only by advanced users
 */
class ITTrackingInteractivity : public G4ITTrackingInteractivity
{
  G4UserTrackingAction* fpUserTrackingAction;
  G4UserSteppingAction* fpUserSteppingAction;
  int fStoreTrajectory;
  std::vector<G4VTrajectory*> fTrajectories;

public:
  ITTrackingInteractivity();
  virtual ~ITTrackingInteractivity();

  virtual void Initialize();
  virtual void StartTracking(G4Track*);
  virtual void AppendStep(G4Track* track, G4Step* step);
  virtual void EndTracking(G4Track*);
  virtual void Finalize();

  void SetUserAction(G4UserTrackingAction*);
  inline G4UserTrackingAction* GetUserTrackingAction();

  void SetUserAction(G4UserSteppingAction*);
  inline G4UserSteppingAction* GetUserSteppingAction();
};

inline
void ITTrackingInteractivity::SetUserAction(G4UserTrackingAction* trackAct)
{
  fpUserTrackingAction = trackAct;
}

inline
void ITTrackingInteractivity::SetUserAction(G4UserSteppingAction* stepAct)
{
  fpUserSteppingAction = stepAct;
}

inline G4UserSteppingAction*
ITTrackingInteractivity::GetUserSteppingAction()
{
  return fpUserSteppingAction;
}

inline G4UserTrackingAction*
ITTrackingInteractivity::GetUserTrackingAction()
{
  return fpUserTrackingAction;
}

#endif // ITTRACKINGINTERACTIVITY_HH
