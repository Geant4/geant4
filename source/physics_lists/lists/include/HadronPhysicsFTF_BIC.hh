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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsFTF_BIC
//
// Author: 2007 Gunter Folger
//
// Modified:
// 19.06.2008 G.Folger: change default for QE to NOT use Chips QE
//
//----------------------------------------------------------------------------
//
#ifndef HadronPhysicsFTF_BIC_h
#define HadronPhysicsFTF_BIC_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4PionBuilder.hh"
#include "G4KaonBuilder.hh"
#include "G4BinaryPionBuilder.hh"
#include "G4BertiniKaonBuilder.hh"
#include "G4FTFBinaryPionBuilder.hh"
#include "G4FTFBinaryKaonBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4LEPProtonBuilder.hh"
#include "G4FTFBinaryProtonBuilder.hh"
#include "G4BinaryProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4LEPNeutronBuilder.hh"
#include "G4FTFBinaryNeutronBuilder.hh"
#include "G4BinaryNeutronBuilder.hh"

#include "G4HyperonFTFPBuilder.hh"
#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"

class HadronPhysicsFTF_BIC : public G4VPhysicsConstructor
{
  public: 
    HadronPhysicsFTF_BIC(G4int verbose =1);
    HadronPhysicsFTF_BIC(const G4String& name,G4bool quasiElastic=false);
    virtual ~HadronPhysicsFTF_BIC();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    void CreateModels();
    G4NeutronBuilder * theNeutrons;
    G4LEPNeutronBuilder * theLEPNeutron;
    G4FTFBinaryNeutronBuilder * theFTFBinaryNeutron;
    G4BinaryNeutronBuilder * theBinaryNeutron;
    
    G4PionBuilder * thePion;
    G4KaonBuilder * theKaon;
    G4BinaryPionBuilder * theBICPion;
    G4BertiniKaonBuilder * theBertiniKaon;
    G4FTFBinaryPionBuilder * theFTFBinaryPion;
    G4FTFBinaryKaonBuilder * theFTFBinaryKaon;
    
    G4ProtonBuilder * thePro;
    G4FTFBinaryProtonBuilder * theFTFBinaryPro; 
    G4BinaryProtonBuilder * theBinaryPro;
    
    G4HyperonFTFPBuilder * theHyperon;
    
    G4AntiBarionBuilder * theAntiBaryon;
    G4FTFPAntiBarionBuilder * theFTFPAntiBaryon;

    G4bool QuasiElastic;
};

#endif

