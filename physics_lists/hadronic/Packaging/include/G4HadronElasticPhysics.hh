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
// $Id: G4HadronElasticPhysics.hh,v 1.2 2006-05-02 08:00:48 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronElasticPhysics
//
// Author: 11 April 2006 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4HadronElasticPhysics_h
#define G4HadronElasticPhysics_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4UHadronElasticProcess.hh"
#include "G4HadronElastic.hh"
#include <vector>

class G4HadronElasticPhysics : public G4VPhysicsConstructor
{
public: 
  G4HadronElasticPhysics(const G4String& name = "elastic",
			 G4int ver = 1, G4bool hp = false);
  virtual ~G4HadronElasticPhysics();

public: 
  // This method will be invoked in the Construct() method. 
  // each particle type will be instantiated
  virtual void ConstructParticle();
 
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type 
  virtual void ConstructProcess();

  void SetLimitSWave(G4double);

  void SetLimitIonKineticEnergy(G4double);

private:

  std::vector<G4HadronicProcess*> p_list;
  G4HadronicInteraction* model;

  G4double pLimit;
  G4double edepLimit;

  G4int    verbose;
  G4bool   hpFlag;
  G4bool   elasticFlag;
  G4bool   wasActivated;
};

inline void G4HadronElasticPhysics::SetLimitSWave(G4double val)
{
  pLimit = val;
}

inline void G4HadronElasticPhysics::SetLimitIonKineticEnergy(G4double val)
{
  edepLimit = val;
}


#endif








