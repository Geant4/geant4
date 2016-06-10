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
// $Id: G4PhotoNuclearProcess.hh 69176 2013-04-22 06:51:19Z gcosmo $
//
// Class Description
// Process for photon nuclear inelastic scattering; 
// to be used in your physics list in case you need this physics.
//

// Hadronic Process: Ion Inelastic Process
// J.P. Wellisch, CERN, Apr. 14 2000
// Last modified: 03-Apr-1997

#ifndef G4PhotoNuclearProcess_h
#define G4PhotoNuclearProcess_h 1
 
#include "G4HadronInelasticProcess.hh"
#include "G4PhotoNuclearCrossSection.hh"


class G4PhotoNuclearProcess : public G4HadronInelasticProcess
{
  public:
    
    G4PhotoNuclearProcess(const G4String& processName = "photonNuclear")
      : G4HadronInelasticProcess(processName, G4Gamma::Gamma() )
    {
      G4CrossSectionDataStore* theStore = GetCrossSectionDataStore();
      theStore->AddDataSet(new G4PhotoNuclearCrossSection);
    }

  ~G4PhotoNuclearProcess() {}

  virtual void ProcessDescription(std::ostream& outFile) const;
};

#endif

