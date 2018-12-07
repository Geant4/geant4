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
// ClassName:   G4HadronElasticPhysicsLEND
//
// Author: 02 June 2011 Tatsumi Koi
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4HadronElasticPhysicsLEND_h
#define G4HadronElasticPhysicsLEND_h 1

#include "G4HadronElasticPhysics.hh"

class G4HadronElasticPhysicsLEND : public G4HadronElasticPhysics
{
public: 

  explicit G4HadronElasticPhysicsLEND(G4int ver = 1, 
				      const G4String& eval=""); 

  virtual ~G4HadronElasticPhysicsLEND();
 
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type 
  void ConstructProcess() final;

private:

  G4HadronElasticPhysicsLEND(G4HadronElasticPhysicsLEND &) = delete;
  G4HadronElasticPhysicsLEND & operator=
  (const G4HadronElasticPhysicsLEND &right) = delete;

  G4String evaluation;
};


#endif








