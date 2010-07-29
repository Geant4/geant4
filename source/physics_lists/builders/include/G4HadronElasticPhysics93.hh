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
// $Id: G4HadronElasticPhysics93.hh,v 1.1 2010-07-29 10:52:13 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronElasticPhysics93
//
// Author: 11 April 2006 V. Ivanchenko
//
// Modified:
// 05.07.2006 V.Ivanchenko fix problem of initialisation of HP
// 23.11.2006 V.Ivanchenko remove variables
//
//----------------------------------------------------------------------------
//

#ifndef G4HadronElasticPhysics93_h
#define G4HadronElasticPhysics93_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4UHadronElasticProcess.hh"
#include "G4HadronElastic.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"

class G4HadronElasticPhysics93 : public G4VPhysicsConstructor
{
public: 
   G4HadronElasticPhysics93(G4int ver = 1);
   G4HadronElasticPhysics93(const G4String& name,
			 G4int ver = 0, G4bool hp = false,
                         G4bool glauber = false);
  virtual ~G4HadronElasticPhysics93();

public: 
  // This method will be invoked in the Construct() method. 
  // each particle type will be instantiated
  virtual void ConstructParticle();
 
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type 
  virtual void ConstructProcess();

private:

  G4HadronicInteraction* model;
  G4HadronicInteraction* neutronModel;
  G4HadronicInteraction* neutronHPModel;

  G4String mname;

  G4int    verbose;
  G4bool   hpFlag;
  G4bool   glFlag;
  G4bool   wasActivated;
};


#endif








