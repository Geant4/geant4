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
// Header:    G4IonPhysicsPHP
//
// Author:    A.Ribon  24-May-2016
//
// Ion physics with ParticleHP, used below 200 MeV/n for d, t, He3, alpha;
// Binary Cascade used above to 190 MeV/n and then FTFP at higher energies.
// This is as G4IonPhysics, except that ParticleHP is used below 200 MeV/n
// for d, t, He3, alpha. 
// 
// Modified:
//
//---------------------------------------------------------------------------

#ifndef G4IonPhysicsPHP_h
#define G4IonPhysicsPHP_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class G4HadronicInteraction;
class G4VCrossSectionDataSet;
class G4FTFBuilder;
class G4ParticleHPInelasticData;

class G4IonPhysicsPHP : public G4VPhysicsConstructor {
  public:

    G4IonPhysicsPHP( G4int ver = 0 );
    G4IonPhysicsPHP( const G4String& nname );
    virtual ~G4IonPhysicsPHP();

    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
    void ConstructParticle();
    void ConstructProcess();

  private:

    void AddProcess( const G4String&, G4ParticleDefinition*, 
                     G4ParticleHPInelasticData*, G4HadronicInteraction*, 
                     G4HadronicInteraction*, G4HadronicInteraction*,
		     G4VCrossSectionDataSet*);

    static G4ThreadLocal G4FTFBuilder*             theBuilder;

    G4int  verbose;
};

#endif
