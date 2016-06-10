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
// $Id: G4VITRestDiscreteProcess.cc 93883 2015-11-03 08:25:04Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------


#include "G4VITRestDiscreteProcess.hh"
G4VITRestDiscreteProcess::G4VITRestDiscreteProcess()
                   :G4VITProcess("No Name Discrete Process")
{
  G4Exception("G4VITRestDiscreteProcess::G4VITRestDiscreteProcess","Illegal operation",
	      JustWarning,"default constructor is called");
}

G4VITRestDiscreteProcess::G4VITRestDiscreteProcess(const G4String& aName , G4ProcessType aType)
                  : G4VITProcess(aName, aType)
{
  enableAlongStepDoIt  = false;
}

G4VITRestDiscreteProcess::~G4VITRestDiscreteProcess()
{
}

G4VITRestDiscreteProcess::G4VITRestDiscreteProcess(const G4VITRestDiscreteProcess& right)
                  : G4VITProcess(right)
{
}
