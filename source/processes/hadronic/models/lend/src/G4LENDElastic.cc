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

#include "G4LENDElastic.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Nucleus.hh"
#include "G4ParticleTable.hh"

G4HadFinalState * G4LENDElastic::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTarg )
{

   G4double temp = aTrack.GetMaterial()->GetTemperature();

   //G4int iZ = int ( aTarg.GetZ() );
   //G4int iA = int ( aTarg.GetN() );
   //migrate to integer A and Z (GetN_asInt returns number of neutrons in the nucleus since this) 
   G4int iZ = aTarg.GetZ_asInt();
   G4int iA = aTarg.GetA_asInt();

   G4double ke = aTrack.GetKineticEnergy();

   //G4HadFinalState* theResult = new G4HadFinalState();
   G4HadFinalState* theResult = &theParticleChange;
   theResult->Clear();

   G4GIDI_target* aTarget = usedTarget_map.find( lend_manager->GetNucleusEncoding( iZ , iA ) )->second->GetTarget();
   G4double aMu = aTarget->getElasticFinalState( ke*MeV, temp, NULL, NULL );

   G4double phi = twopi*G4UniformRand();
   G4double theta = std::acos( aMu );
   //G4double sinth = std::sin( theta );

   G4ReactionProduct theNeutron( const_cast<G4ParticleDefinition *>( aTrack.GetDefinition() ) );
   theNeutron.SetMomentum( aTrack.Get4Momentum().vect() );
   theNeutron.SetKineticEnergy( ke );

//G4cout << "iZ " << iZ << " iA " << iA  << G4endl;

   G4ReactionProduct theTarget( G4ParticleTable::GetParticleTable()->FindIon( iZ , iA , 0 , iZ ) );

   G4double mass = G4ParticleTable::GetParticleTable()->FindIon( iZ , iA , 0 , iZ )->GetPDGMass();

// add Thermal motion 
   G4double kT = k_Boltzmann*temp;
   G4ThreeVector v ( G4RandGauss::shoot() * std::sqrt( kT*mass ) 
                   , G4RandGauss::shoot() * std::sqrt( kT*mass ) 
                   , G4RandGauss::shoot() * std::sqrt( kT*mass ) );
   theTarget.SetMomentum( v );

     G4ThreeVector the3Neutron = theNeutron.GetMomentum();
     G4double nEnergy = theNeutron.GetTotalEnergy();
     G4ThreeVector the3Target = theTarget.GetMomentum();
     G4double tEnergy = theTarget.GetTotalEnergy();
     G4ReactionProduct theCMS;
     G4double totE = nEnergy+tEnergy;
     G4ThreeVector the3CMS = the3Target+the3Neutron;
     theCMS.SetMomentum(the3CMS);
     G4double cmsMom = std::sqrt(the3CMS*the3CMS);
     G4double sqrts = std::sqrt((totE-cmsMom)*(totE+cmsMom));
     theCMS.SetMass(sqrts);
     theCMS.SetTotalEnergy(totE);

       theNeutron.Lorentz(theNeutron, theCMS);
       theTarget.Lorentz(theTarget, theCMS);
       G4double en = theNeutron.GetTotalMomentum(); // already in CMS.
       G4ThreeVector cms3Mom=theNeutron.GetMomentum(); // for neutron direction in CMS
       G4double cms_theta=cms3Mom.theta();
       G4double cms_phi=cms3Mom.phi();
       G4ThreeVector tempVector;
       tempVector.setX( std::cos(theta)*std::sin(cms_theta)*std::cos(cms_phi)
                       +std::sin(theta)*std::cos(phi)*std::cos(cms_theta)*std::cos(cms_phi)
                       -std::sin(theta)*std::sin(phi)*std::sin(cms_phi) );
       tempVector.setY( std::cos(theta)*std::sin(cms_theta)*std::sin(cms_phi)
                       +std::sin(theta)*std::cos(phi)*std::cos(cms_theta)*std::sin(cms_phi)
                       +std::sin(theta)*std::sin(phi)*std::cos(cms_phi) );
       tempVector.setZ( std::cos(theta)*std::cos(cms_theta)
                       -std::sin(theta)*std::cos(phi)*std::sin(cms_theta) );
       tempVector *= en;
       theNeutron.SetMomentum(tempVector);
       theTarget.SetMomentum(-tempVector);
       G4double tP = theTarget.GetTotalMomentum();
       G4double tM = theTarget.GetMass();
       theTarget.SetTotalEnergy(std::sqrt((tP+tM)*(tP+tM)-2.*tP*tM));


       theNeutron.Lorentz(theNeutron, -1.*theCMS);

//110913 Add Protection for very low energy (1e-6eV) scattering 
      if ( theNeutron.GetKineticEnergy() <= 0 )
      {
         theNeutron.SetTotalEnergy ( theNeutron.GetMass() * ( 1 + std::pow( 10 , -15.65 ) ) );
      }

      theTarget.Lorentz(theTarget, -1.*theCMS);
      if ( theTarget.GetKineticEnergy() < 0 )
      {
         theTarget.SetTotalEnergy ( theTarget.GetMass() * ( 1 + std::pow( 10 , -15.65 ) ) );
      }
//110913 END

       theTarget.Lorentz(theTarget, -1.*theCMS);

     theResult->SetEnergyChange(theNeutron.GetKineticEnergy());
     theResult->SetMomentumChange(theNeutron.GetMomentum().unit());
     G4DynamicParticle* theRecoil = new G4DynamicParticle;

//     theRecoil->SetDefinition( ionTable->GetIon( iZ , iA ) ); 
       theRecoil->SetDefinition( G4ParticleTable::GetParticleTable()->FindIon( iZ, iA , 0, iZ ));
     theRecoil->SetMomentum( theTarget.GetMomentum() );

     theResult->AddSecondary( theRecoil );

   return theResult; 

}

