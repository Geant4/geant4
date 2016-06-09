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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4MiscCHIPSBuilder
//
// Author: 2010, G.Folger
//
// created by modifying G4QInelasticCHIPSBuilder to apply only to MISC particles
// -----------------------------------------------------------------------------
#ifndef G4MiscCHIPSBuilder_h
#define G4MiscCHIPSBuilder_h 1

#include "globals.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4QInelastic.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

class G4MiscCHIPSBuilder
{
 public: 
  G4MiscCHIPSBuilder(G4int verbose=0);
  virtual ~G4MiscCHIPSBuilder();

 public: 
  void Build();

 protected:
  // the particle table has the complete List of existing particle types
  G4ParticleTable* theParticleTable;
  G4ParticleTable::G4PTblDicIterator* theParticleIterator;

 private:
  G4int  verbose;
  G4bool wasActivated;
  G4QInelastic* inelastic;
};

#endif
