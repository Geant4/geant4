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
// ClassName: HadronPhysicsCRMC_FTFP_BERT
//
// Author:    2018 Alberto Ribon
//
// This is a variant of HadronPhysicsFTFP_BERT whereby CRMC is used 
// for modeling final-state for pion- , kaon- , proton- and neutron-nuclear
// inelastic interactions at very high energies.
// For other hadron projectile types (e.g. hyperons, antinucleons and
// antihyperons) the usual FTFP_BERT approach is used at all energies.
// The inelastic hadronic cross sections are, for all hadron projectiles
// and energies, the usual ones (exactly as in FTFP_BERT).
//
// Modified:
//----------------------------------------------------------------------------
//
#ifndef HadronPhysicsCRMC_FTFP_BERT_h
#define HadronPhysicsCRMC_FTFP_BERT_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4Cache.hh"

class G4ComponentGGHadronNucleusXsc;
class G4VCrossSectionDataSet;


class HadronPhysicsCRMC_FTFP_BERT : public G4VPhysicsConstructor {
  public:
    HadronPhysicsCRMC_FTFP_BERT( G4int verbose = 1 );
    HadronPhysicsCRMC_FTFP_BERT( const G4String& name );
    virtual ~HadronPhysicsCRMC_FTFP_BERT();
    virtual void ConstructParticle() override;
    virtual void ConstructProcess() override;
  protected:
    void CreateModels();
    virtual void Neutron();
    virtual void Proton();
    virtual void Pion();
    virtual void Kaon();
    virtual void Others();
    virtual void DumpBanner();
    // This contains extra configurataion specific to this PL
    virtual void ExtraConfiguration();
    G4double minCRMC;
    G4double minFTFP;
    G4double maxFTFP;
    G4double minBERT;
    G4double maxBERT;
    // Thread-private data write them here to delete them
    G4VectorCache< G4VCrossSectionDataSet* > xs_ds;
    G4Cache< G4ComponentGGHadronNucleusXsc* > xs_k;
};

#endif

