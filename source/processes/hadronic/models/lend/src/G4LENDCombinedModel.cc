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
#include "G4LENDCombinedModel.hh"
#include "G4LENDCombinedCrossSection.hh"
#include "G4LENDElastic.hh"
#include "G4LENDInelastic.hh"
#include "G4LENDCapture.hh"
#include "G4LENDFission.hh"
#include "G4DynamicParticle.hh"

G4LENDCombinedModel::G4LENDCombinedModel( G4ParticleDefinition* pd )
:G4LENDModel( "LENDCombinedModel" ) {
   proj = pd;
   crossSection = new G4LENDCombinedCrossSection( pd );
   elastic = new G4LENDElastic( pd );
   inelastic = new G4LENDInelastic( pd );
   capture = new G4LENDCapture( pd );
   fission = new G4LENDFission( pd );
   channels[0] = elastic;
   channels[1] = inelastic;
   channels[2] = capture;
   channels[3] = fission;
}

void G4LENDCombinedModel::BuildPhysicsTable(const G4ParticleDefinition& projectile) {
   crossSection->BuildPhysicsTable( projectile );
   create_used_target_map();
}

G4bool G4LENDCombinedModel::HasData( const G4DynamicParticle* , G4int iZ, G4int iA , G4int iM,
                                     const G4Isotope* , const G4Element* , const G4Material* ) {

   G4bool result = false;
   if ( get_target_from_map ( lend_manager->GetNucleusEncoding ( iZ, iA , iM ) ) != NULL ) result=true;
   return result;
}


G4HadFinalState * G4LENDCombinedModel::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTarg ){

   G4LENDModel* channel = NULL;

   G4int iZ = aTarg.GetZ_asInt();
   G4int iA = aTarg.GetA_asInt();
   //To pass kinetic energy, need to generate dynamic particle  
   G4DynamicParticle* dp = new G4DynamicParticle( proj , G4ThreeVector(0.,0.,1.) , aTrack.GetKineticEnergy() );
   G4int ichannel = crossSection->SelectChannel( dp , iZ , iA , aTarg.GetIsotope(), NULL , aTrack.GetMaterial() );
   delete dp;
   //ichannel=1;
   channel = channels[ichannel];
   return channel->ApplyYourself(aTrack,aTarg);
}
