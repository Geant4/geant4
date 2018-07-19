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
// $Id: G4HadronPhysicsQGSP_BERT.hh 105736 2017-08-16 13:01:11Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsQGSP_BERT
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 23.11.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko: remove stopping
// 25.04.2007 G.Folger: Add quasielastic option, use this by default
// 10.12.2007 G.Folger: Add projectilediffrative option for proton/neutron, off by default
// 31.10.2012 A.Ribon: Use G4MiscBuilder
// 19.03.2013 A.Ribon: Replace LEP with FTFP
//
//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsQGSP_BERT_h
#define G4HadronPhysicsQGSP_BERT_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4Cache.hh"

class G4ComponentGGHadronNucleusXsc;
class G4VCrossSectionDataSet;

class G4HadronPhysicsQGSP_BERT : public G4VPhysicsConstructor
{
  public: 
    G4HadronPhysicsQGSP_BERT(G4int verbose =1);
    G4HadronPhysicsQGSP_BERT(const G4String& name, G4bool quasiElastic=true);
    virtual ~G4HadronPhysicsQGSP_BERT();

    virtual void ConstructParticle() override;
    virtual void ConstructProcess() override;
    virtual void TerminateWorker() override;

  protected:
    G4bool QuasiElasticFTF;
    G4bool QuasiElasticQGS;
    void CreateModels();
    virtual void Neutron();
    virtual void Proton();
    virtual void Pion();
    virtual void Kaon() { /*Done together w/ Pion*/ }
    virtual void Others();
    virtual void DumpBanner() {}
    //This contains extra configurataion specific to this PL
    virtual void ExtraConfiguration();

    G4double minQGSP_proton;
    G4double minQGSP_neutron;
    G4double minQGSP_pik;
    G4double minFTFP_proton;
    G4double minFTFP_neutron;
    G4double minFTFP_pik;
    G4double maxFTFP_proton;
    G4double maxFTFP_neutron;
    G4double maxFTFP_pik;
    G4double minBERT_proton;
    G4double minBERT_neutron;
    G4double minBERT_pik;
    G4double maxBERT_proton;
    G4double maxBERT_neutron;
    G4double maxBERT_pik;

    //Thread-private data write them here to delete them
    G4VectorCache<G4VCrossSectionDataSet*> xs_ds;
    G4Cache<G4ComponentGGHadronNucleusXsc*> xs_k;
};

#endif

