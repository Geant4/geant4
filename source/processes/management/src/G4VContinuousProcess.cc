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
// $Id: G4VContinuousProcess.cc,v 1.7 2010-12-22 09:14:54 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
// --------------------------------------------------------------
//   New Physics scheme           8 Jan. 1997  H.Kurahige
// ------------------------------------------------------------

#include "G4VContinuousProcess.hh"
G4VContinuousProcess::G4VContinuousProcess()
  :G4VProcess("No Name Continuous Process"),
   valueGPILSelection(CandidateForSelection) 
{
  G4Exception("G4VContinuousProcess::G4VContinuousProcess()","ProcMan102",
	      JustWarning,"Default constructor is called");
}

G4VContinuousProcess::G4VContinuousProcess(const G4String& aName , G4ProcessType aType)
  : G4VProcess(aName, aType),
    valueGPILSelection(CandidateForSelection) 
{
  enableAtRestDoIt = false;
  enablePostStepDoIt = false;
}

G4VContinuousProcess::~G4VContinuousProcess()
{
}

G4VContinuousProcess::G4VContinuousProcess(G4VContinuousProcess& right)
                  : G4VProcess(right),
		    valueGPILSelection(right.valueGPILSelection)
{
}











