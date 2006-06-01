//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ElectronNuclearProcess.cc,v 1.1 2006-06-01 15:32:51 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4ElectronNuclearProcess.hh" 
#include "G4Electron.hh"

G4ElectronNuclearProcess::
G4ElectronNuclearProcess(const G4String& processName)
  : G4HadronInelasticProcess( processName, G4Electron::Electron() )
{ 
  G4CrossSectionDataStore * theStore = GetCrossSectionDataStore();
  theStore->AddDataSet(&theData);
} 
    
G4ElectronNuclearProcess::~G4ElectronNuclearProcess()
{
}

