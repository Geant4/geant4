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
//
// G4 Hadron Elastic Scattering Process -- header file
// F.W. Jones, TRIUMF, 04-JUN-96
//  
// For further comments see G4HadronElasticProcess.cc
//
// Modified by FWJ 03-DEC-96: uses G4LCrossSectionData
// Modified by FWJ    FEB-97: adapted to new tracking design
// Modified by FWJ    MAR-97: newer new tracking design
//
// Modified by FWJ 27-MAR-97: first version for Alpha release
// 14-APR-97 FWJ: cross section data class name changed
//
// changing design to use Model-scheme: HPW, 20th June 1997
//
// 14-APR-98 FWJ: variant G4HadronElastic process for
// G4CrossSectionDataSet/DataStore class design.
// 29-JUN-98 FWJ: default data set G4HadronCrossSections
// 01-SEP-2008 V.Ivanchenko: use methods from the base class
//

// Class Description
// Process for hadron nuclear elastic scattering; 
// to be used in your physics list in case you need this physics.
// Class Description - End


#ifndef G4HadronElasticProcess_h
#define G4HadronElasticProcess_h 1
 
#include "globals.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronElasticDataSet.hh"

class G4HadronElasticProcess : public G4HadronicProcess
{
public:

   G4HadronElasticProcess(const G4String& processName = "HadronElastic");

   virtual ~G4HadronElasticProcess();
 
   virtual G4bool IsApplicable(const G4ParticleDefinition& aParticleType);

};
#endif
