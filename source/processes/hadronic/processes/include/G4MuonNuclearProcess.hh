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
// $Id: $
//
// Author:      D.H. Wright
// Date:        2 February 2011
//
// Description: muon interacts with nucleus by exchange of virtual
//              gamma, which then interacts hadronically    


#ifndef G4MuonNuclearProcess_h
#define G4MuonNuclearProcess_h 1
 
#include "globals.hh"
#include "G4HadronicProcess.hh"

class G4MuonNuclearProcess : public G4HadronicProcess
{
public:

  G4MuonNuclearProcess(const G4String& processName ="muonNuclear");
  
  virtual ~G4MuonNuclearProcess();
 
  virtual G4bool IsApplicable(const G4ParticleDefinition& aParticleType);

  virtual void ProcessDescription(std::ostream& outFile) const;
};
#endif
