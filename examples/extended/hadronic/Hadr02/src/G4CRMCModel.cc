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
// *                                                                  *
// * Parts of this code which have been  developed by Abdel-Waged     *
// * et al under contract (31-465) to the King Abdul-Aziz City for    *
// * Science and Technology (KACST), the National Centre of           *
// * Mathematics and Physics (NCMP), Saudi Arabia.                    *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr02/src/G4UrQMD1_3Model.cc
/// \brief Implementation of the G4CRMCModel class
//
//
//---------------------------------------------------------------------------
//
// ClassName: CRMCNeutronBuilder
//
// Author:    2018 Alberto Ribon
//
// Modified:
//
//----------------------------------------------------------------------------
//  
#ifdef G4_USE_CRMC

#include "G4CRMCModel.hh"
#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4IonTable.hh"
#include "G4CollisionOutput.hh"
#include "G4V3DNucleus.hh"
#include "G4Track.hh"
#include "G4Nucleus.hh"
#include "G4LorentzRotation.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Version.hh"
#include "G4AntiHe3.hh"
#include "G4AntiDeuteron.hh"
#include "G4AntiTriton.hh"
#include "G4AntiAlpha.hh"
#include <fstream>
#include <string>
#include "G4HadronicParameters.hh"


G4CRMCModel::G4CRMCModel( const int model ) : G4HadronicInteraction( "CRMC" ),  verbose( 0 ), 
  fModel( model ) {
  if (verbose > 3) {
    G4cout << " >>> G4CRMCModel default constructor" << G4endl;
  }
  WelcomeMessage();
  CurrentEvent = 0;

  G4int ranseed = 1234567; 
  G4cout << "\n seed:  " << ranseed << G4endl;

  const char* crmc_param = "crmc.param";  // CRMC default, see CRMCoptions.cc in the CRMC package
  
  fInterface = new CRMCinterface;
  fInterface->init( model );

  // Open FORTRAN IO at first call.
  // Notice that the energy unit must be GeV.
  fInterface->crmc_init( G4HadronicParameters::Instance()->GetMaxEnergy()/GeV, 
                         ranseed, model, 0, 0, crmc_param, "", 0 );

  fParticleTable = G4ParticleTable::GetParticleTable();
  fIonTable = fParticleTable->GetIonTable();
}


G4CRMCModel::~G4CRMCModel() {}


/*
const std::pair< G4double, G4double > G4CRMCModel::GetFatalEnergyCheckLevels() const {
  // The default in the base class, G4HadronicInteraction, is 2% and 1 GeV .
  // Note that when both relative and absolute need fail the final-state is rejected
  // and the interaction is re-sampled.
  return std::pair< G4double, G4double >( 10.0*perCent, DBL_MAX );
}
*/


G4HadFinalState* G4CRMCModel::ApplyYourself( const G4HadProjectile &theProjectile, 
                                             G4Nucleus &theTarget ) {

  // Clean up data vectors
  gCRMC_data.Clean();

  // The original track will always be discontinued and secondaries followed.
  fFinalState.Clear();
  fFinalState.SetStatusChange( stopAndKill );

  // Get relevant information about the projectile and target (A, Z, energy/nuc, momentum, etc)
  const G4ParticleDefinition* definitionP = theProjectile.GetDefinition();
  G4int AP = G4lrint( definitionP->GetBaryonNumber() );
  G4int ZP = G4lrint( definitionP->GetPDGCharge() );
  G4int AT = theTarget.GetA_asInt();
  G4int ZT = theTarget.GetZ_asInt();
  G4int idProj = definitionP->GetPDGEncoding();
  G4int idTarg = ZT*10000 + AT*10;
  G4ThreeVector pBefore = theProjectile.Get4Momentum().vect();
  G4double E = theProjectile.GetTotalEnergy();        
  G4double totalEbefore = E*AP + theTarget.AtomicMass( AT, ZT ) + theTarget.GetEnergyDeposit();

  // Note: because of the rotation of the projectile along the z-axis (made by
  //       the calling method G4HadronicProcess::PostStepDoIt), it is equivalent
  //       to take the momentum component along the z-axis or the whole momentum. 
  G4double pProj = theProjectile.Get4Momentum().vect().mag() / GeV;  // Energy unit must be GeV

  // Note: from my understanding of the CRMC interface, it seems that in the case
  //       of nucleus projectile, the energy per nucleon (instead of the energy of
  //       the whole projectile) should be provided.
  //       We consider the absolute value of the baryon number to cover also the
  //       case of anti-nuclei.
  if ( std::abs( AP ) > 1 ) pProj /= static_cast< G4double >( std::abs( AP ) );  // Energy per nucleon

  G4double pTarg = 0.0;

  fInterface->crmc_set( 1,         // fNCollision,
                        pProj,     // fCfg.fProjectileMomentum,
                        pTarg,     // fCfg.fTargetMomentum,
                        idProj,    // fCfg.fProjectileId,
                        idTarg );  // fCfg.fTargetId);

  // Sample 1 interaction
  fInterface->crmc_generate( 0,    // fCfg.fTypoaut,
                             1,    // iColl+1,
                             gCRMC_data.fNParticles,
                             gCRMC_data.fImpactParameter,
                             gCRMC_data.fPartId[0],
                             gCRMC_data.fPartPx[0],
                             gCRMC_data.fPartPy[0],
                             gCRMC_data.fPartPz[0],
                             gCRMC_data.fPartEnergy[0],
                             gCRMC_data.fPartMass[0],
                             gCRMC_data.fPartStatus[0] );

  // Save secondary particles for output
  for ( G4int i = 0; i < gCRMC_data.fNParticles; i++ ) { 
    // Keep only final state particles 
    // (-9 is the beam, 2 is a particle which decayed and 1 is final)
    if ( gCRMC_data.fPartStatus[i] != 1 ) continue;
    G4ParticleDefinition* pdef = GetParticleDefinition( gCRMC_data.fPartId[i] );
    G4DynamicParticle* part = new G4DynamicParticle( pdef, 
                                                     G4ThreeVector( gCRMC_data.fPartPx[i]*GeV,
                                                                    gCRMC_data.fPartPy[i]*GeV,
                                                                    gCRMC_data.fPartPz[i]*GeV ) );
    fFinalState.AddSecondary( part );
  }

  if ( verbose >= 3 ) {
    G4double totalEafter = 0.0;
    G4ThreeVector totalPafter;
    G4double charge = 0.0;
    G4int baryon = 0;
    G4int nSecondaries = fFinalState.GetNumberOfSecondaries();
    for ( G4int j = 0; j < nSecondaries; j++ ) {
      totalEafter += fFinalState.GetSecondary(j)->GetParticle()->GetTotalEnergy();
      totalPafter += fFinalState.GetSecondary(j)->GetParticle()->GetMomentum();
      G4ParticleDefinition* pd = fFinalState.GetSecondary(j)->GetParticle()->GetDefinition();
      charge += pd->GetPDGCharge();
      baryon += pd->GetBaryonNumber();
    }
    G4cout << "----------------------------------------" << G4endl
           << "Total energy before collision   = " << totalEbefore << " MeV" << G4endl
           << "Total energy after collision    = " << totalEafter  << " MeV" << G4endl
           << "Total momentum before collision = " << pBefore      << " MeV/c" << G4endl
           << "Total momentum after collision  = " << totalPafter  << " MeV/c" << G4endl;
    if ( verbose >= 4 ) {
      G4cout << "Total charge before collision   = " << (ZP + ZT)*eplus << G4endl
             << "Total charge after collision    = " << charge <<G4endl
             << "Total baryon number before collision = " << (AP + AT) << G4endl
             << "Total baryon number after collision  = "<< baryon << G4endl;
    }
    G4cout << "----------------------------------------" << G4endl;
  }

  return &fFinalState;
}


void G4CRMCModel::WelcomeMessage () const {
  G4cout << G4endl
         << " *****************************************************************" << G4endl;
  if ( fModel == 0 ) {
    G4cout << "                              CRMC - EPOS LHC " << G4endl;
  } else if ( fModel == 1 ) {
    G4cout << "                              CRMC - EPOS 1.99" << G4endl; 
  } else if ( fModel == 6 ) {
    G4cout << "                              CRMC - SIBYLL 2.3c" << G4endl; 
  } else if ( fModel == 12 ) {
    G4cout << "                              CRMC - DPMJET 3" << G4endl; 
  }
  G4cout << " *****************************************************************" << G4endl
         << G4endl;
  return;
}


G4ParticleDefinition* G4CRMCModel::GetParticleDefinition( long particle_id ) {
  // CRMC ion definition : id = crmc_ion_coef_0 + crmc_ion_coef_z*Z + crmc_ion_coef_a*A
  const G4int crmc_ion_coef_0 = 1000000000;
  const G4int crmc_ion_coef_z = 10000;
  const G4int crmc_ion_coef_a = 10;
  G4ParticleDefinition* pdef = fParticleTable->FindParticle( particle_id );
  if ( ! pdef  &&  particle_id > crmc_ion_coef_0 ) {
    int Z = ( particle_id - crmc_ion_coef_0 ) / crmc_ion_coef_z;
    int A = ( particle_id - crmc_ion_coef_0 - crmc_ion_coef_z*Z ) / crmc_ion_coef_a;
    pdef = fIonTable->GetIon( Z, A );
  }
  return pdef;
}

#endif //G4_USE_CRMC

