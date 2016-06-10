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
// $Id: G4HadronHElasticPhysics.hh 90757 2015-06-09 07:45:14Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronHElasticPhysics
//
// Author: 23 November 2006 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4HadronHElasticPhysics_h
#define G4HadronHElasticPhysics_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"

class G4DiffElasticRatio;

class G4HadronHElasticPhysics : public G4VPhysicsConstructor
{
public: 

  G4HadronHElasticPhysics(G4int ver = 0, G4bool diffraction=false); 

  virtual ~G4HadronHElasticPhysics();

  // This method will be invoked in the Construct() method. 
  // each particle type will be instantiated
  virtual void ConstructParticle();
 
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type 
  virtual void ConstructProcess();

  inline void SetDiffraction(G4bool val) {fDiffraction = val;};

private:

  G4HadronHElasticPhysics(G4HadronHElasticPhysics &);
  G4HadronHElasticPhysics & operator=(const G4HadronHElasticPhysics &right);

  G4int  verbose;
  G4bool fDiffraction;

  static G4ThreadLocal G4DiffElasticRatio* diffRatio;

};


#endif








