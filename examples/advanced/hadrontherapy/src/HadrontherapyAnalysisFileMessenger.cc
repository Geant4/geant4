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
// See more at: http://workgroup.lngs.infn.it/geant4lns/
//
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------
#include "HadrontherapyAnalysisFileMessenger.hh"
#include "HadrontherapyAnalysisManager.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

#ifdef ANALYSIS_USE
//
/////////////////////////////////////////////////////////////////////////////
HadrontherapyAnalysisFileMessenger::HadrontherapyAnalysisFileMessenger(HadrontherapyAnalysisManager* amgr)
:AnalysisManager(amgr)
{ 
  FileNameCmd = new G4UIcmdWithAString("/analysis/setAnalysisFile",this);
  FileNameCmd->SetGuidance("Set the .root filename for the root-output");
  FileNameCmd->SetDefaultValue("default.root");
  FileNameCmd->SetParameterName("choice",true); ///<doc did not say what second boolean really does
  FileNameCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyAnalysisFileMessenger::~HadrontherapyAnalysisFileMessenger()
{
  delete FileNameCmd;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisFileMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if (command == FileNameCmd)
    {
	AnalysisManager->SetAnalysisFileName(newValue);
	AnalysisManager->flush(); //< fills matrix, writes it into file and books a new file and histograms etc.
    //AnalysisManager->book();
	}
}

#endif
