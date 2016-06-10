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
// $Id: G4VITTimeStepComputer.cc 81769 2014-06-05 08:30:21Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4VITTimeStepComputer.hh"

G4ThreadLocal G4double G4VITTimeStepComputer::fCurrentGlobalTime = -1;
G4ThreadLocal G4double G4VITTimeStepComputer::fUserMinTimeStep   = -1;

G4VITTimeStepComputer::G4VITTimeStepComputer()
{
    fVerbose = 0;
//    fReactants = 0;
    fReactants.reset();
    fSampledMinTimeStep = 0 ;
    fpReactionTable      = 0;
}

G4VITTimeStepComputer::~G4VITTimeStepComputer()
{;}

G4VITTimeStepComputer& G4VITTimeStepComputer::operator=(const G4VITTimeStepComputer& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

G4VITTimeStepComputer::G4VITTimeStepComputer(const G4VITTimeStepComputer& right)
{
    fVerbose            = right.fVerbose ;
    fpReactionTable      = right.fpReactionTable ;
//    fReactants          = 0;
    fReactants          .reset();
    fSampledMinTimeStep = 0 ;
}

void G4VITTimeStepComputer::SetTimes(const G4double& currentGlobalTime,
                            const G4double& userMinStepTime)
{
    fCurrentGlobalTime = currentGlobalTime ;
    fUserMinTimeStep   = userMinStepTime ;
}
