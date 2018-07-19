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
// 080505 Fixed and changed sampling method of impact parameter by T. Koi 
// 080602 Fix memory leaks by T. Koi 
// 080612 Delete unnecessary dependency and unused functions
//        Change criterion of reaction by T. Koi
// 081107 Add UnUseGEM (then use the default channel of G4Evaporation)
//            UseFrag (chage criterion of a inelastic reaction)
//        Fix bug in nucleon projectiles  by T. Koi    
// 090122 Be8 -> Alpha + Alpha 
// 090331 Change member shenXS and genspaXS object to pointer 
// 091119 Fix for incidence of neutral particles 
//
#include "G4QMDReaction.hh"
#include "G4QMDNucleus.hh"
#include "G4QMDGroundStateNucleus.hh"
#include "G4Pow.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"

#include "G4CrossSectionDataSetRegistry.hh"

G4QMDReaction::G4QMDReaction()
: G4HadronicInteraction("QMDModel")
, system ( NULL )
, deltaT ( 1 ) // in fsec (c=1)
, maxTime ( 100 ) // will have maxTime-th time step
, envelopF ( 1.05 ) // 10% for Peripheral reactions
, gem ( true )
, frag ( false )
{

   //090331
   shenXS = new G4IonsShenCrossSection();
   //genspaXS = new G4GeneralSpaceNNCrossSection();
   piNucXS = (G4PiNuclearCrossSection*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4PiNuclearCrossSection::Default_Name());
   meanField = new G4QMDMeanField();
   collision = new G4QMDCollision();

   excitationHandler = new G4ExcitationHandler;
   excitationHandler->SetDeexChannelsType( fCombined );
   //fEvaporation - 8 default channels
   //fCombined    - 8 default + 60 GEM
   //fGEM         - 2 default + 66 GEM
   evaporation = new G4Evaporation;
   excitationHandler->SetEvaporation( evaporation );
   setEvaporationCh();

   coulomb_collision_gamma_proj = 0.0;
   coulomb_collision_rx_proj = 0.0;
   coulomb_collision_rz_proj = 0.0;
   coulomb_collision_px_proj = 0.0;
   coulomb_collision_pz_proj = 0.0;

   coulomb_collision_gamma_targ = 0.0;
   coulomb_collision_rx_targ = 0.0;
   coulomb_collision_rz_targ = 0.0;
   coulomb_collision_px_targ = 0.0;
   coulomb_collision_pz_targ = 0.0;

}



G4QMDReaction::~G4QMDReaction()
{
   delete evaporation; 
   delete excitationHandler;
   delete collision;
   delete meanField;
}



G4HadFinalState* G4QMDReaction::ApplyYourself( const G4HadProjectile & projectile , G4Nucleus & target )
{
   //G4cout << "G4QMDReaction::ApplyYourself" << G4endl;

   theParticleChange.Clear();

   system = new G4QMDSystem;

   G4int proj_Z = 0;
   G4int proj_A = 0;
   const G4ParticleDefinition* proj_pd = ( const G4ParticleDefinition* ) projectile.GetDefinition();
   if ( proj_pd->GetParticleType() == "nucleus" )
   {
      proj_Z = proj_pd->GetAtomicNumber();
      proj_A = proj_pd->GetAtomicMass();
   }
   else
   {
      proj_Z = (int)( proj_pd->GetPDGCharge()/eplus );
      proj_A = 1;
   }
   //G4int targ_Z = int ( target.GetZ() + 0.5 );
   //G4int targ_A = int ( target.GetN() + 0.5 );
   //migrate to integer A and Z (GetN_asInt returns number of neutrons in the nucleus since this) 
   G4int targ_Z = target.GetZ_asInt();
   G4int targ_A = target.GetA_asInt();
   const G4ParticleDefinition* targ_pd = G4IonTable::GetIonTable()->GetIon( targ_Z , targ_A , 0.0 );


   //G4NistManager* nistMan = G4NistManager::Instance();
//   G4Element* G4NistManager::FindOrBuildElement( targ_Z );

   const G4DynamicParticle* proj_dp = new G4DynamicParticle ( proj_pd , projectile.Get4Momentum() );
   //const G4Element* targ_ele =  nistMan->FindOrBuildElement( targ_Z ); 
   //G4double aTemp = projectile.GetMaterial()->GetTemperature();

     //090331
  
   G4VCrossSectionDataSet* theXS = shenXS;

   if ( proj_pd->GetParticleType() == "meson" ) theXS = piNucXS;

   //G4double xs_0 = genspaXS->GetCrossSection ( proj_dp , targ_ele , aTemp );
   //G4double xs_0 = theXS->GetCrossSection ( proj_dp , targ_ele , aTemp );
   //110822 
   G4double xs_0 = theXS->GetIsoCrossSection ( proj_dp , targ_Z , targ_A );

     G4double bmax_0 = std::sqrt( xs_0 / pi );
     //std::cout << "bmax_0 in fm (fermi) " <<  bmax_0/fermi << std::endl;

     //delete proj_dp; 

   G4bool elastic = true;
   
   std::vector< G4QMDNucleus* > nucleuses; // Secondary nuceluses 
   G4ThreeVector boostToReac; // ReactionSystem (CM or NN); 
   G4ThreeVector boostBackToLAB; // Reaction System to LAB; 

   G4LorentzVector targ4p( G4ThreeVector( 0.0 ) , targ_pd->GetPDGMass()/GeV );
   G4ThreeVector boostLABtoCM = targ4p.findBoostToCM( proj_dp->Get4Momentum()/GeV ); // CM of target and proj; 

   G4double p1 = proj_dp->GetMomentum().mag()/GeV/proj_A; 
   G4double m1 = proj_dp->GetDefinition()->GetPDGMass()/GeV/proj_A;
   G4double e1 = std::sqrt( p1*p1 + m1*m1 ); 
   G4double e2 = targ_pd->GetPDGMass()/GeV/targ_A;
   G4double beta_nn = -p1 / ( e1+e2 );

   G4ThreeVector boostLABtoNN ( 0. , 0. , beta_nn ); // CM of NN; 

   G4double beta_nncm = ( - boostLABtoCM.beta() + boostLABtoNN.beta() ) / ( 1 - boostLABtoCM.beta() * boostLABtoNN.beta() ) ;  

   //std::cout << targ4p << std::endl; 
   //std::cout << proj_dp->Get4Momentum()<< std::endl; 
   //std::cout << beta_nncm << std::endl; 
   G4ThreeVector boostNNtoCM( 0. , 0. , beta_nncm ); // 
   G4ThreeVector boostCMtoNN( 0. , 0. , -beta_nncm ); // 

   boostToReac = boostLABtoNN; 
   boostBackToLAB = -boostLABtoNN; 

   delete proj_dp; 

   G4int icounter = 0;
   G4int icounter_max = 1024;
   while ( elastic ) // Loop checking, 11.03.2015, T. Koi
   {
      icounter++;
      if ( icounter > icounter_max ) { 
	 G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
         break;
      }

// impact parameter 
      //G4double bmax = 1.05*(bmax_0/fermi);  // 10% for Peripheral reactions
      G4double bmax = envelopF*(bmax_0/fermi);
      G4double b = bmax * std::sqrt ( G4UniformRand() );
//071112
      //G4double b = 0;
      //G4double b = bmax;
      //G4double b = bmax/1.05 * 0.7 * G4UniformRand();

      //G4cout << "G4QMDRESULT bmax_0 = " << bmax_0/fermi << " fm, bmax = " << bmax << " fm , b = " << b  << " fm " << G4endl; 

      G4double plab = projectile.GetTotalMomentum()/GeV;
      G4double elab = ( projectile.GetKineticEnergy() + proj_pd->GetPDGMass() + targ_pd->GetPDGMass() )/GeV;

      calcOffSetOfCollision( b , proj_pd , targ_pd , plab , elab , bmax , boostCMtoNN );

// Projectile
      G4LorentzVector proj4pLAB = projectile.Get4Momentum()/GeV;

      G4QMDGroundStateNucleus* proj(NULL); 
      if ( projectile.GetDefinition()->GetParticleType() == "nucleus" 
        || projectile.GetDefinition()->GetParticleName() == "proton"
        || projectile.GetDefinition()->GetParticleName() == "neutron" )
      {

         proj_Z = proj_pd->GetAtomicNumber();
         proj_A = proj_pd->GetAtomicMass();

         proj = new G4QMDGroundStateNucleus( proj_Z , proj_A );
         //proj->ShowParticipants();


         meanField->SetSystem ( proj );
         proj->SetTotalPotential( meanField->GetTotalPotential() );
         proj->CalEnergyAndAngularMomentumInCM();

      }

// Target
      //G4int iz = int ( target.GetZ() );
      //G4int ia = int ( target.GetN() );
      //migrate to integer A and Z (GetN_asInt returns number of neutrons in the nucleus since this) 
      G4int iz = int ( target.GetZ_asInt() );
      G4int ia = int ( target.GetA_asInt() );

      G4QMDGroundStateNucleus* targ = new G4QMDGroundStateNucleus( iz , ia );

      meanField->SetSystem (targ );
      targ->SetTotalPotential( meanField->GetTotalPotential() );
      targ->CalEnergyAndAngularMomentumInCM();
   
      //G4LorentzVector targ4p( G4ThreeVector( 0.0 ) , targ->GetNuclearMass()/GeV );
// Boost Vector to CM
      //boostToCM = targ4p.findBoostToCM( proj4pLAB );

//    Target 
      for ( G4int i = 0 ; i < targ->GetTotalNumberOfParticipant() ; i++ )
      {

         G4ThreeVector p0 = targ->GetParticipant( i )->GetMomentum();
         G4ThreeVector r0 = targ->GetParticipant( i )->GetPosition();

         G4ThreeVector p ( p0.x() + coulomb_collision_px_targ 
                         , p0.y() 
                         , p0.z() * coulomb_collision_gamma_targ + coulomb_collision_pz_targ ); 

         G4ThreeVector r ( r0.x() + coulomb_collision_rx_targ 
                         , r0.y() 
                         , r0.z() / coulomb_collision_gamma_targ + coulomb_collision_rz_targ ); 
     
         system->SetParticipant( new G4QMDParticipant( targ->GetParticipant( i )->GetDefinition() , p , r ) );
         system->GetParticipant( i )->SetTarget();

      }

      G4LorentzVector proj4pCM = CLHEP::boostOf ( proj4pLAB , boostToReac );
      G4LorentzVector targ4pCM = CLHEP::boostOf ( targ4p , boostToReac );

//    Projectile
      if ( proj != NULL )
      {

//    projectile is nucleus

         for ( G4int i = 0 ; i < proj->GetTotalNumberOfParticipant() ; i++ )
         {

            G4ThreeVector p0 = proj->GetParticipant( i )->GetMomentum();
            G4ThreeVector r0 = proj->GetParticipant( i )->GetPosition();

            G4ThreeVector p ( p0.x() + coulomb_collision_px_proj 
                            , p0.y() 
                            , p0.z() * coulomb_collision_gamma_proj + coulomb_collision_pz_proj ); 

            G4ThreeVector r ( r0.x() + coulomb_collision_rx_proj 
                            , r0.y() 
                            , r0.z() / coulomb_collision_gamma_proj + coulomb_collision_rz_proj ); 
     
            system->SetParticipant( new G4QMDParticipant( proj->GetParticipant( i )->GetDefinition() , p  , r ) );
            system->GetParticipant ( i + targ->GetTotalNumberOfParticipant() )->SetProjectile();
         }

      }
      else
      {

//       projectile is particle

         // avoid multiple set in "elastic" loop
         if ( system->GetTotalNumberOfParticipant() == targ->GetTotalNumberOfParticipant() )
         {

            G4int i = targ->GetTotalNumberOfParticipant(); 
      
            G4ThreeVector p0( 0 ); 
            G4ThreeVector r0( 0 );

            G4ThreeVector p ( p0.x() + coulomb_collision_px_proj 
                            , p0.y() 
                            , p0.z() * coulomb_collision_gamma_proj + coulomb_collision_pz_proj ); 

            G4ThreeVector r ( r0.x() + coulomb_collision_rx_proj 
                            , r0.y() 
                            , r0.z() / coulomb_collision_gamma_proj + coulomb_collision_rz_proj ); 

            system->SetParticipant( new G4QMDParticipant( (G4ParticleDefinition*)projectile.GetDefinition() , p , r ) );
            // This is not important becase only 1 projectile particle.
            system->GetParticipant ( i )->SetProjectile();
         }

      }
      //system->ShowParticipants();

      delete targ;
      delete proj;

   meanField->SetSystem ( system );
   collision->SetMeanField ( meanField );

// Time Evolution 
   //std::cout << "Start time evolution " << std::endl;
   //system->ShowParticipants();
   for ( G4int i = 0 ; i < maxTime ; i++ )
   {
      //G4cout << " do Paropagate " << i << " th time step. " << G4endl;
      meanField->DoPropagation( deltaT );
      //system->ShowParticipants();
      collision->CalKinematicsOfBinaryCollisions( deltaT );

      if ( i / 10 * 10 == i ) 
      {
         //G4cout << i << " th time step. " << G4endl;
         //system->ShowParticipants();
      } 
      //system->ShowParticipants();
   }
   //system->ShowParticipants();


   //std::cout << "Doing Cluster Judgment " << std::endl;

   nucleuses = meanField->DoClusterJudgment();

// Elastic Judgment  

   G4int numberOfSecondary = int ( nucleuses.size() ) + system->GetTotalNumberOfParticipant(); 

   G4int sec_a_Z = 0;
   G4int sec_a_A = 0;
   const G4ParticleDefinition* sec_a_pd = NULL;
   G4int sec_b_Z = 0;
   G4int sec_b_A = 0;
   const G4ParticleDefinition* sec_b_pd = NULL;

   if ( numberOfSecondary == 2 )
   {

      G4bool elasticLike_system = false;
      if ( nucleuses.size() == 2 ) 
      {

         sec_a_Z = nucleuses[0]->GetAtomicNumber();
         sec_a_A = nucleuses[0]->GetMassNumber();
         sec_b_Z = nucleuses[1]->GetAtomicNumber();
         sec_b_A = nucleuses[1]->GetMassNumber();

         if ( ( sec_a_Z == proj_Z && sec_a_A == proj_A && sec_b_Z == targ_Z && sec_b_A == targ_A )
           || ( sec_a_Z == targ_Z && sec_a_A == targ_A && sec_b_Z == proj_Z && sec_b_A == proj_A ) )
         {
            elasticLike_system = true;
         } 

      }
      else if ( nucleuses.size() == 1 ) 
      {

         sec_a_Z = nucleuses[0]->GetAtomicNumber();
         sec_a_A = nucleuses[0]->GetMassNumber();
         sec_b_pd = system->GetParticipant( 0 )->GetDefinition();

         if ( ( sec_a_Z == proj_Z && sec_a_A == proj_A && sec_b_pd == targ_pd )
           || ( sec_a_Z == targ_Z && sec_a_A == targ_A && sec_b_pd == proj_pd ) )
         {
            elasticLike_system = true;
         } 

      }  
      else
      {

         sec_a_pd = system->GetParticipant( 0 )->GetDefinition();
         sec_b_pd = system->GetParticipant( 1 )->GetDefinition();
 
         if ( ( sec_a_pd == proj_pd && sec_b_pd == targ_pd ) 
           || ( sec_a_pd == targ_pd && sec_b_pd == proj_pd ) ) 
         {
            elasticLike_system = true;
         } 

      } 

      if ( elasticLike_system == true )
      {

         G4bool elasticLike_energy = true;
//    Cal ExcitationEnergy 
         for ( G4int i = 0 ; i < int ( nucleuses.size() ) ; i++ )
         { 

            //meanField->SetSystem( nucleuses[i] );
            meanField->SetNucleus( nucleuses[i] );
            //nucleuses[i]->SetTotalPotential( meanField->GetTotalPotential() );
            //nucleuses[i]->CalEnergyAndAngularMomentumInCM();

            if ( nucleuses[i]->GetExcitationEnergy()*GeV > 1.0*MeV ) elasticLike_energy = false;  

         } 

//    Check Collision 
         G4bool withCollision = true;
         if ( system->GetNOCollision() == 0 ) withCollision = false;

//    Final judegement for Inelasitc or Elastic;
//
//       ElasticLike without Collision 
         //if ( elasticLike_energy == true && withCollision == false ) elastic = true;  // ielst = 0
//       ElasticLike with Collision 
         //if ( elasticLike_energy == true && withCollision == true ) elastic = true;   // ielst = 1 
//       InelasticLike without Collision 
         //if ( elasticLike_energy == false ) elastic = false;                          // ielst = 2                
         if ( frag == true )
            if ( elasticLike_energy == false ) elastic = false;
//       InelasticLike with Collision 
         if ( elasticLike_energy == false && withCollision == true ) elastic = false; // ielst = 3

      }

      }
      else
      {

//       numberOfSecondary != 2 
         elastic = false;

      }

//071115
      //G4cout << elastic << G4endl;
      // if elastic is true try again from sampling of impact parameter 

      if ( elastic == true )
      {
         // delete this nucleues
         for ( std::vector< G4QMDNucleus* >::iterator
               it = nucleuses.begin() ; it != nucleuses.end() ; it++ )
         {
            delete *it;
         }
         nucleuses.clear();
      }

   } 


// Statical Decay Phase

   for ( std::vector< G4QMDNucleus* >::iterator it
       = nucleuses.begin() ; it != nucleuses.end() ; it++ )
   {

/*
      G4cout << "G4QMDRESULT "
             << (*it)->GetAtomicNumber() 
             << " " 
             << (*it)->GetMassNumber() 
             << " " 
             << (*it)->Get4Momentum() 
             << " " 
             << (*it)->Get4Momentum().vect() 
             << " " 
             << (*it)->Get4Momentum().restMass() 
             << " " 
             << (*it)->GetNuclearMass()/GeV 
             << G4endl;
*/

      meanField->SetNucleus ( *it );

      if ( (*it)->GetAtomicNumber() == 0  // neutron cluster
        || (*it)->GetAtomicNumber() == (*it)->GetMassNumber() ) // proton cluster
      {
         // push back system 
         for ( G4int i = 0 ; i < (*it)->GetTotalNumberOfParticipant() ; i++ )
         {
            G4QMDParticipant* aP = new G4QMDParticipant( ( (*it)->GetParticipant( i ) )->GetDefinition() , ( (*it)->GetParticipant( i ) )->GetMomentum() , ( (*it)->GetParticipant( i ) )->GetPosition() );  
            system->SetParticipant ( aP );  
         } 
         continue;  
      }

      G4double nucleus_e = std::sqrt ( G4Pow::GetInstance()->powN ( (*it)->GetNuclearMass()/GeV , 2 ) + G4Pow::GetInstance()->powN ( (*it)->Get4Momentum().vect().mag() , 2 ) );
      G4LorentzVector nucleus_p4CM ( (*it)->Get4Momentum().vect() , nucleus_e ); 

//      std::cout << "G4QMDRESULT nucleus deltaQ " << deltaQ << std::endl;

      G4int ia = (*it)->GetMassNumber();
      G4int iz = (*it)->GetAtomicNumber();

      G4LorentzVector lv ( G4ThreeVector( 0.0 ) , (*it)->GetExcitationEnergy()*GeV + G4IonTable::GetIonTable()->GetIonMass( iz , ia ) );

      G4Fragment* aFragment = new G4Fragment( ia , iz , lv );

      G4ReactionProductVector* rv;
      rv = excitationHandler->BreakItUp( *aFragment );
      G4bool notBreak = true;
      for ( G4ReactionProductVector::iterator itt
          = rv->begin() ; itt != rv->end() ; itt++ )
      {

          notBreak = false;
          // Secondary from this nucleus (*it) 
          const G4ParticleDefinition* pd = (*itt)->GetDefinition();

          G4LorentzVector p4 ( (*itt)->GetMomentum()/GeV , (*itt)->GetTotalEnergy()/GeV );  //in nucleus(*it) rest system
          G4LorentzVector p4_CM = CLHEP::boostOf( p4 , -nucleus_p4CM.findBoostToCM() );  // Back to CM
          G4LorentzVector p4_LAB = CLHEP::boostOf( p4_CM , boostBackToLAB ); // Back to LAB  


//090122
          //theParticleChange.AddSecondary( dp ); 
          if ( !( pd->GetAtomicNumber() == 4 && pd->GetAtomicMass() == 8 ) )
          {
             G4DynamicParticle* dp = new G4DynamicParticle( pd , p4_LAB*GeV );  
             theParticleChange.AddSecondary( dp ); 
          }
          else
          {
             //Be8 -> Alpha + Alpha + Q
             G4ThreeVector randomized_direction( G4UniformRand() , G4UniformRand() , G4UniformRand() );
             randomized_direction = randomized_direction.unit();
             G4double q_decay = (*itt)->GetMass() - 2*G4Alpha::Alpha()->GetPDGMass();
             G4double p_decay = std::sqrt ( G4Pow::GetInstance()->powN(G4Alpha::Alpha()->GetPDGMass()+q_decay/2,2) - G4Pow::GetInstance()->powN(G4Alpha::Alpha()->GetPDGMass() , 2 ) ); 
             G4LorentzVector p4_a1 ( p_decay*randomized_direction , G4Alpha::Alpha()->GetPDGMass()+q_decay/2 );  //in Be8 rest system
             
             G4LorentzVector p4_a1_Be8 = CLHEP::boostOf ( p4_a1/GeV , -p4.findBoostToCM() );
             G4LorentzVector p4_a1_CM = CLHEP::boostOf ( p4_a1_Be8 , -nucleus_p4CM.findBoostToCM() );
             G4LorentzVector p4_a1_LAB = CLHEP::boostOf ( p4_a1_CM , boostBackToLAB );

             G4LorentzVector p4_a2 ( -p_decay*randomized_direction , G4Alpha::Alpha()->GetPDGMass()+q_decay/2 );  //in Be8 rest system
             
             G4LorentzVector p4_a2_Be8 = CLHEP::boostOf ( p4_a2/GeV , -p4.findBoostToCM() );
             G4LorentzVector p4_a2_CM = CLHEP::boostOf ( p4_a2_Be8 , -nucleus_p4CM.findBoostToCM() );
             G4LorentzVector p4_a2_LAB = CLHEP::boostOf ( p4_a2_CM , boostBackToLAB );
             
             G4DynamicParticle* dp1 = new G4DynamicParticle( G4Alpha::Alpha() , p4_a1_LAB*GeV );  
             G4DynamicParticle* dp2 = new G4DynamicParticle( G4Alpha::Alpha() , p4_a2_LAB*GeV );  
             theParticleChange.AddSecondary( dp1 ); 
             theParticleChange.AddSecondary( dp2 ); 
          }
//090122

/*
          G4cout
                << "Regist Secondary "
                << (*itt)->GetDefinition()->GetParticleName()
                << " "
                << (*itt)->GetMomentum()/GeV
                << " "
                << (*itt)->GetKineticEnergy()/GeV
                << " "
                << (*itt)->GetMass()/GeV
                << " "
                << (*itt)->GetTotalEnergy()/GeV
                << " "
                << (*itt)->GetTotalEnergy()/GeV * (*itt)->GetTotalEnergy()/GeV
                 - (*itt)->GetMomentum()/GeV * (*itt)->GetMomentum()/GeV
                << " "
                << nucleus_p4CM.findBoostToCM() 
                << " "
                << p4
                << " "
                << p4_CM
                << " "
                << p4_LAB
                << G4endl;
*/

      }
      if ( notBreak == true )
      {

         const G4ParticleDefinition* pd = G4IonTable::GetIonTable()->GetIon( (*it)->GetAtomicNumber() , (*it)->GetMassNumber(), (*it)->GetExcitationEnergy()*GeV );
         G4LorentzVector p4_CM = nucleus_p4CM;
         G4LorentzVector p4_LAB = CLHEP::boostOf( p4_CM , boostBackToLAB ); // Back to LAB  
         G4DynamicParticle* dp = new G4DynamicParticle( pd , p4_LAB*GeV );  
         theParticleChange.AddSecondary( dp ); 

      }

      for ( G4ReactionProductVector::iterator itt
          = rv->begin() ; itt != rv->end() ; itt++ )
      {
          delete *itt;
      }
      delete rv;

      delete aFragment;

   }



   for ( G4int i = 0 ; i < system->GetTotalNumberOfParticipant() ; i++ )
   {

      // Secondary particles 

      const G4ParticleDefinition* pd = system->GetParticipant( i )->GetDefinition();
      G4LorentzVector p4_CM = system->GetParticipant( i )->Get4Momentum();
      G4LorentzVector p4_LAB = CLHEP::boostOf( p4_CM , boostBackToLAB );
      G4DynamicParticle* dp = new G4DynamicParticle( pd , p4_LAB*GeV );  
      theParticleChange.AddSecondary( dp ); 

/*
      G4cout << "G4QMDRESULT "
      << "r" << i << " " << system->GetParticipant ( i ) -> GetPosition() << " "
      << "p" << i << " " << system->GetParticipant ( i ) -> Get4Momentum()
      << G4endl;
*/

   }

   for ( std::vector< G4QMDNucleus* >::iterator it
       = nucleuses.begin() ; it != nucleuses.end() ; it++ )
   {
      delete *it;  // delete nulceuse 
   }
   nucleuses.clear();

   system->Clear();
   delete system; 

   theParticleChange.SetStatusChange( stopAndKill );

   return &theParticleChange;

}



void G4QMDReaction::calcOffSetOfCollision( G4double b , 
const G4ParticleDefinition* pd_proj ,
const G4ParticleDefinition* pd_targ ,
G4double ptot , G4double etot , G4double bmax , G4ThreeVector boostToCM )
{

   G4double mass_proj = pd_proj->GetPDGMass()/GeV;
   G4double mass_targ = pd_targ->GetPDGMass()/GeV;
  
   G4double stot = std::sqrt ( etot*etot - ptot*ptot );

   G4double pstt = std::sqrt ( ( stot*stot - ( mass_proj + mass_targ ) * ( mass_proj + mass_targ ) 
                  ) * ( stot*stot - ( mass_proj - mass_targ ) * ( mass_proj - mass_targ ) ) ) 
                 / ( 2.0 * stot );

   G4double pzcc = pstt;
   G4double eccm = stot - ( mass_proj + mass_targ );
   
   G4int zp = 1;
   G4int ap = 1;
   if ( pd_proj->GetParticleType() == "nucleus" )
   {
      zp = pd_proj->GetAtomicNumber();
      ap = pd_proj->GetAtomicMass();
   }
   else 
   {
      // proton, neutron, mesons
      zp = int ( pd_proj->GetPDGCharge()/eplus + 0.5 );  
      // ap = 1;
   }
   

   G4int zt = pd_targ->GetAtomicNumber();
   G4int at = pd_targ->GetAtomicMass();


   //G4double rmax0 = 8.0;  // T.K dicide parameter value  // for low energy
   G4double rmax0 = bmax + 4.0;
   G4double rmax = std::sqrt( rmax0*rmax0 + b*b );

   G4double ccoul = 0.001439767;
   G4double pcca = 1.0 - double ( zp * zt ) * ccoul / eccm / rmax - ( b / rmax )*( b / rmax );

   G4double pccf = std::sqrt( pcca );

   //Fix for neutral particles
   G4double aas1 = 0.0;
   G4double bbs = 0.0;

   if ( zp != 0 )
   {
      G4double aas = 2.0 * eccm * b / double ( zp * zt ) / ccoul;
      bbs = 1.0 / std::sqrt ( 1.0 + aas*aas );
      aas1 = ( 1.0 + aas * b / rmax ) * bbs;
   }

   G4double cost = 0.0;
   G4double sint = 0.0;
   G4double thet1 = 0.0;
   G4double thet2 = 0.0;
   if ( 1.0 - aas1*aas1 <= 0 || 1.0 - bbs*bbs <= 0.0 )   
   {
      cost = 1.0;
      sint = 0.0;
   } 
   else 
   {
      G4double aat1 = aas1 / std::sqrt ( 1.0 - aas1*aas1 );
      G4double aat2 = bbs / std::sqrt ( 1.0 - bbs*bbs );

      thet1 = std::atan ( aat1 );
      thet2 = std::atan ( aat2 );

//    TK enter to else block  
      G4double theta = thet1 - thet2;
      cost = std::cos( theta );
      sint = std::sin( theta );
   }

   G4double rzpr = -rmax * cost * ( mass_targ ) / ( mass_proj + mass_targ );
   G4double rzta =  rmax * cost * ( mass_proj ) / ( mass_proj + mass_targ );

   G4double rxpr = rmax / 2.0 * sint;

   G4double rxta = -rxpr;


   G4double pzpc = pzcc * (  cost * pccf + sint * b / rmax ); 
   G4double pxpr = pzcc * ( -sint * pccf + cost * b / rmax ); 

   G4double pztc = - pzpc;
   G4double pxta = - pxpr;

   G4double epc = std::sqrt ( pzpc*pzpc + pxpr*pxpr + mass_proj*mass_proj );
   G4double etc = std::sqrt ( pztc*pztc + pxta*pxta + mass_targ*mass_targ );

   G4double pzpr = pzpc;
   G4double pzta = pztc;
   G4double epr = epc;
   G4double eta = etc;

// CM -> NN
   G4double gammacm = boostToCM.gamma();
   //G4double betacm = -boostToCM.beta();
   G4double betacm = boostToCM.z();
   pzpr = pzpc + betacm * gammacm * ( gammacm / ( 1. + gammacm ) * pzpc * betacm + epc );
   pzta = pztc + betacm * gammacm * ( gammacm / ( 1. + gammacm ) * pztc * betacm + etc );
   epr = gammacm * ( epc + betacm * pzpc );
   eta = gammacm * ( etc + betacm * pztc );

   //G4double betpr = pzpr / epr;
   //G4double betta = pzta / eta;

   G4double gammpr = epr / ( mass_proj );
   G4double gammta = eta / ( mass_targ );
      
   pzta = pzta / double ( at );
   pxta = pxta / double ( at );
      
   pzpr = pzpr / double ( ap );
   pxpr = pxpr / double ( ap );

   G4double zeroz = 0.0; 

   rzpr = rzpr -zeroz;
   rzta = rzta -zeroz;

   // Set results 
   coulomb_collision_gamma_proj = gammpr;
   coulomb_collision_rx_proj = rxpr;
   coulomb_collision_rz_proj = rzpr;
   coulomb_collision_px_proj = pxpr;
   coulomb_collision_pz_proj = pzpr;

   coulomb_collision_gamma_targ = gammta;
   coulomb_collision_rx_targ = rxta;
   coulomb_collision_rz_targ = rzta;
   coulomb_collision_px_targ = pxta;
   coulomb_collision_pz_targ = pzta;

}



void G4QMDReaction::setEvaporationCh()
{

   if ( gem == true ) 
      evaporation->SetGEMChannel();
   else
      evaporation->SetDefaultChannel();

}

void G4QMDReaction::ModelDescription(std::ostream& outFile) const
{
   outFile << "Lorentz covarianted Quantum Molecular Dynamics model for nucleus (particle) vs nucleus reactions\n";
}
