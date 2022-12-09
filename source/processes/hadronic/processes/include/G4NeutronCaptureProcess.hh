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
// G4 Neutron capture process -- header file
// F.W. Jones, TRIUMF, 03-DEC-96
//  
// For further comments see G4NeutronCaptureProcess.cc.
//
// 27-MAR-97 FWJ: first version for Alpha release
// 14-APR-97 FWJ: cross section data class name changed
//
// 19-MAY-98 FWJ: variant G4HadronCapture process for
// G4CrossSectionDataSet/DataStore class design.
// 29-JUN-98 FWJ: default data set G4HadronCrossSections
// 01-SEP-2008 V.Ivanchenko: use methods from the base class
//

// Class Description
// Process for capture of neutral hadrons; 
// to be used in your physics list in case you need this physics.
// Class Description - End


#ifndef G4NeutronCaptureProcess_h
#define G4NeutronCaptureProcess_h 1
 
#include "globals.hh"
#include "G4HadronicProcess.hh"

class G4NeutronCaptureProcess final : public G4HadronicProcess
{
public:

  explicit G4NeutronCaptureProcess(const G4String& processName ="nCapture");
  
  ~G4NeutronCaptureProcess() final = default;
 
  G4bool IsApplicable(const G4ParticleDefinition&) final;

  void ProcessDescription(std::ostream& outFile) const final;
};
#endif
