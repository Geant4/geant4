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
// $Id: G4UserTimeStepAction.cc 100802 2016-11-02 14:55:27Z gcosmo $
//
// Author: Mathieu Karamitros
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include <G4VScheduler.hh>
#include "G4UserTimeStepAction.hh"

G4UserTimeStepAction::G4UserTimeStepAction()
{;}

G4UserTimeStepAction::~G4UserTimeStepAction()
{;}

G4UserTimeStepAction::G4UserTimeStepAction(const G4UserTimeStepAction& /*other*/){;}

G4UserTimeStepAction& G4UserTimeStepAction::operator=(const G4UserTimeStepAction& /*rhs*/)
{
//    if (this == &rhs) return *this;
    return *this;
}

void G4UserTimeStepAction::SetMinimumTimeSteps(std::map<double, double>* timeSteps)
{
	G4VScheduler::Instance()-> SetTimeSteps(timeSteps);
}

void G4UserTimeStepAction::AddTimeStep(double startingTime, double timeStep)
{
  G4VScheduler::Instance()-> AddTimeStep(startingTime,timeStep);
}
