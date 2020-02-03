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
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronElasticPhysics
//
// Author: 23 November 2006 V. Ivanchenko
//
// Modified:
// 29.07.2010 V.Ivanchenko rename this class from G4HadronHElasticPhysics to
//                         G4HadronElasticPhysics, old version of the class
//                         is renamed to G4HadronElasticPhysics93
// 03.06.2011 V.Ivanchenko change design - added access to pointers of the
//                         neutron process and model
//
//----------------------------------------------------------------------------
//

#ifndef G4HadronElasticPhysics_h
#define G4HadronElasticPhysics_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4HadronElastic.hh"
#include "G4HadronicProcess.hh"

class G4VCrossSectionDataSet;
class G4ParticleDefinition;

class G4HadronElasticPhysics : public G4VPhysicsConstructor
{
public: 

  explicit G4HadronElasticPhysics(G4int ver = 0, 
				  const G4String& nam = "hElasticWEL_CHIPS_XS"); 

  virtual ~G4HadronElasticPhysics();

  // This method will be invoked in the Construct() method. 
  // each particle type will be instantiated
  void ConstructParticle() override;
 
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type 
  void ConstructProcess() override;

  G4HadronicProcess* GetElasticProcess(const G4ParticleDefinition* part) const;

  G4HadronElastic* GetElasticModel(const G4ParticleDefinition* part) const;

  G4HadronicProcess* GetNeutronProcess() const;

  G4HadronElastic* GetNeutronModel() const;

  void AddXSection(const G4ParticleDefinition*, G4VCrossSectionDataSet*) const; 

private:

  // copy constructor and hide assignment operator
  G4HadronElasticPhysics(G4HadronElasticPhysics &) = delete;
  G4HadronElasticPhysics & operator=(const G4HadronElasticPhysics &right) = delete;

protected:
  G4int verbose;
};
#endif








