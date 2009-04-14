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
// ClassName:   G4MiscQGSCBuilder
//
// Author: 2009 M. Kossov (on the basis of the G4MiscLHEPBuilder)
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef G4MiscQGSCBuilder_h
#define G4MiscQGSCBuilder_h 1

#include "globals.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4TheoFSGenerator.hh"
#include "G4StringChipsParticleLevelInterface.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4QuasiElasticChannel.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

//class G4HadronicProcess;
//class G4TheoFSGenerator;
//class G4StringChipsParticleLevelInterface;
//class G4ExcitedStringDecay;
//class G4HadronProcessStore;
//class G4QuasiElasticChannel;

class G4MiscQGSCBuilder
//class G4MiscQGSCBuilder : public G4VPhysicsConstructor
{
  public: 
    G4MiscQGSCBuilder(G4int verbose);
    virtual ~G4MiscQGSCBuilder();

  public: 
    void Build();

 protected:
  // the particle table has the complete List of existing particle types
  G4ParticleTable* theParticleTable;
  G4ParticleTable::G4PTblDicIterator* theParticleIterator;

 private:
    G4TheoFSGenerator * theModel;
    G4StringChipsParticleLevelInterface * theCascade;
    G4QGSModel< G4QGSParticipants > * theQGSCModel;
    G4ExcitedStringDecay * theQGSCDecay;
    G4QuasiElasticChannel * theQuasiElastic;

    G4int  verbose;
    G4bool wasActivated;
};
// 2009 by M. Kossov

#endif
