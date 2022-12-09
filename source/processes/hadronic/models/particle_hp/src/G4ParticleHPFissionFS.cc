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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 12-Apr-06 fix in delayed neutron and photon emission without FS data by T. Koi
// 07-Sep-11 M. Kelsey -- Follow change to G4HadFinalState interface
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//

#include "G4Exp.hh"
#include "G4ParticleHPFissionFS.hh"
#include "G4PhysicalConstants.hh"
#include "G4Nucleus.hh"
#include "G4DynamicParticleVector.hh"
#include "G4ParticleHPFissionERelease.hh"
#include "G4IonTable.hh"
#include "G4PhysicsModelCatalog.hh"


 G4ParticleHPFissionFS::G4ParticleHPFissionFS()
 {
   secID = G4PhysicsModelCatalog::GetModelID( "model_NeutronHPFission" );
   hasXsec = false;
   produceFissionFragments = false;
 }

 void G4ParticleHPFissionFS::Init (G4double A, G4double Z, G4int M,
                                   G4String& dirName, G4String& aFSType,
                                   G4ParticleDefinition* projectile )
 {
    theFS.Init(A, Z, M, dirName, aFSType, projectile);
    theFC.Init(A, Z, M, dirName, aFSType, projectile);
    theSC.Init(A, Z, M, dirName, aFSType, projectile);
    theTC.Init(A, Z, M, dirName, aFSType, projectile);
    theLC.Init(A, Z, M, dirName, aFSType, projectile);

    theFF.Init(A, Z, M, dirName, aFSType, projectile);
    if ( G4ParticleHPManager::GetInstance()->GetProduceFissionFragments()
      && theFF.HasFSData() ) 
    {
       G4cout << "Fission fragment production is now activated in HP package for " 
       << "Z = " << (G4int)Z
       << ", A = " << (G4int)A
       << G4endl;
       G4cout << "As currently modeled this option precludes production of delayed neutrons from fission fragments." << G4endl;
       produceFissionFragments = true; 
    }
 }

 G4HadFinalState * G4ParticleHPFissionFS::ApplyYourself(const G4HadProjectile & theTrack)
 {  
   // Because it may change by UI command 
   produceFissionFragments=G4ParticleHPManager::GetInstance()->GetProduceFissionFragments(); 

   // prepare neutron
   if ( theResult.Get() == NULL ) theResult.Put( new G4HadFinalState );
   theResult.Get()->Clear();
   G4double eKinetic = theTrack.GetKineticEnergy();
   const G4HadProjectile *incidentParticle = &theTrack;
   G4ReactionProduct theNeutron( const_cast<G4ParticleDefinition *>(incidentParticle->GetDefinition()) );
   theNeutron.SetMomentum( incidentParticle->Get4Momentum().vect() );
   theNeutron.SetKineticEnergy( eKinetic );

   // prepare target
   G4Nucleus aNucleus;
   G4ReactionProduct theTarget; 
   G4double targetMass = theFS.GetMass();
   G4ThreeVector neuVelo = (1./incidentParticle->GetDefinition()->GetPDGMass())*theNeutron.GetMomentum();
   theTarget = aNucleus.GetBiasedThermalNucleus( targetMass, neuVelo, theTrack.GetMaterial()->GetTemperature());
   theTarget.SetDefinition( G4IonTable::GetIonTable()->GetIon( G4int(theBaseZ), G4int(theBaseA) , 0.0 ) );  //TESTPHP

   // set neutron and target in the FS classes 
   theFS.SetNeutronRP(theNeutron);
   theFS.SetTarget(theTarget);
   theFC.SetNeutronRP(theNeutron);
   theFC.SetTarget(theTarget);
   theSC.SetNeutronRP(theNeutron);
   theSC.SetTarget(theTarget);
   theTC.SetNeutronRP(theNeutron);
   theTC.SetTarget(theTarget);
   theLC.SetNeutronRP(theNeutron);
   theLC.SetTarget(theTarget);
   theFF.SetNeutronRP(theNeutron);
   theFF.SetTarget(theTarget);

   // boost to target rest system and decide on channel.
   theNeutron.Lorentz(theNeutron, -1*theTarget);

   // dice the photons

   G4DynamicParticleVector * thePhotons;    
   thePhotons = theFS.GetPhotons();

   // select the FS in charge

   eKinetic = theNeutron.GetKineticEnergy();    
   G4double xSec[4];
   xSec[0] = theFC.GetXsec(eKinetic);
   xSec[1] = xSec[0]+theSC.GetXsec(eKinetic);
   xSec[2] = xSec[1]+theTC.GetXsec(eKinetic);
   xSec[3] = xSec[2]+theLC.GetXsec(eKinetic);
   G4int it;
   unsigned int i=0;
   G4double random = G4UniformRand();
   if(xSec[3]==0) 
   {
     it=-1;
   }
   else
   {
     for(i=0; i<4; i++)
     {
       it =i;
       if(random<xSec[i]/xSec[3]) break;
     }
   }

   // dice neutron multiplicities, energies and momenta in Lab. @@
   // no energy conservation on an event-to-event basis. we rely on the data to be ok. @@
   // also for mean, we rely on the consistancy of the data. @@

   G4int Prompt=0, delayed=0, all=0;
   G4DynamicParticleVector * theNeutrons = 0;
   switch(it) // check logic, and ask, if partials can be assumed to correspond to individual particles @@@
   {
     case 0:
       theFS.SampleNeutronMult(all, Prompt, delayed, eKinetic, 0);
       if(Prompt==0&&delayed==0) Prompt=all;
       theNeutrons = theFC.ApplyYourself(Prompt); // delayed always in FS 
       // take 'U' into account explicitly (see 5.4) in the sampling of energy @@@@
       break;
     case 1:
       theFS.SampleNeutronMult(all, Prompt, delayed, eKinetic, 1);
       if(Prompt==0&&delayed==0) Prompt=all;
       theNeutrons = theSC.ApplyYourself(Prompt); // delayed always in FS, off done in FSFissionFS
       break;
     case 2:
       theFS.SampleNeutronMult(all, Prompt, delayed, eKinetic, 2);
       if(Prompt==0&&delayed==0) Prompt=all;
       theNeutrons = theTC.ApplyYourself(Prompt); // delayed always in FS
       break;
     case 3:
       theFS.SampleNeutronMult(all, Prompt, delayed, eKinetic, 3);
       if(Prompt==0&&delayed==0) Prompt=all;
       theNeutrons = theLC.ApplyYourself(Prompt); // delayed always in FS
       break;
     default:
       break;
   }

   // dice delayed neutrons and photons, and fallback 
   // for Prompt in case channel had no FS data; add all paricles to FS.

   if ( produceFissionFragments ) delayed=0;

   G4double * theDecayConstants;

   if( theNeutrons != 0)
   {
     theDecayConstants = new G4double[delayed];
     for(i=0; i<theNeutrons->size(); ++i)
     {
       theResult.Get()->AddSecondary(theNeutrons->operator[](i), secID);
     }
     delete theNeutrons;  

     G4DynamicParticleVector * theDelayed = 0;
     theDelayed = theFS.ApplyYourself(0, delayed, theDecayConstants);
     for(i=0; i<theDelayed->size(); i++)
     {
       G4double time = -G4Log(G4UniformRand())/theDecayConstants[i];
       time += theTrack.GetGlobalTime();
       theResult.Get()->AddSecondary(theDelayed->operator[](i), secID);
       theResult.Get()->GetSecondary(theResult.Get()->GetNumberOfSecondaries()-1)->SetTime(time);
     }
     delete theDelayed;                  
   }
   else
   {
     theFS.SampleNeutronMult(all, Prompt, delayed, eKinetic, 0);
     theDecayConstants = new G4double[delayed];
     if(Prompt==0&&delayed==0) Prompt=all;
     theNeutrons = theFS.ApplyYourself(Prompt, delayed, theDecayConstants);
     G4int i0;
     for(i0=0; i0<Prompt; ++i0)
     {
       theResult.Get()->AddSecondary(theNeutrons->operator[](i0), secID);
     }

     for(i0=Prompt; i0<Prompt+delayed; ++i0)
     {
       // Protect against the very rare case of division by zero
       G4double time = 0.0;
       if ( theDecayConstants[i0-Prompt] > 1.0e-30 ) {
         time = -G4Log(G4UniformRand())/theDecayConstants[i0-Prompt];
       } else {
         G4ExceptionDescription ed;
         ed << " theDecayConstants[i0-Prompt]=" << theDecayConstants[i0-Prompt]
            << "   -> cannot sample the time : set it to 0.0 !" << G4endl;
         G4Exception( "G4ParticleHPFissionFS::ApplyYourself ", "HAD_FISSIONHP_001", JustWarning, ed );
       }

       time += theTrack.GetGlobalTime();        
       theResult.Get()->AddSecondary(theNeutrons->operator[](i0), secID);
       theResult.Get()->GetSecondary(theResult.Get()->GetNumberOfSecondaries()-1)->SetTime(time);
     }
     delete theNeutrons;   
   }
   delete [] theDecayConstants;

   std::size_t nPhotons = 0;
   if(thePhotons!=0)
   {
     nPhotons = thePhotons->size();
     for(i=0; i<nPhotons; ++i)
     {
       theResult.Get()->AddSecondary(thePhotons->operator[](i), secID);
     }
     delete thePhotons; 
   }

   // finally deal with local energy depositions.

   G4ParticleHPFissionERelease * theERelease = theFS.GetEnergyRelease();
   G4double eDepByFragments = theERelease->GetFragmentKinetic();
   //theResult.SetLocalEnergyDeposit(eDepByFragments);
   if ( !produceFissionFragments ) theResult.Get()->SetLocalEnergyDeposit(eDepByFragments);
   // clean up the primary neutron
   theResult.Get()->SetStatusChange(stopAndKill);

   if ( produceFissionFragments )
   {
      G4int fragA_Z=0;
      G4int fragA_A=0;
      G4int fragA_M=0;
      // System is traget rest!
      theFF.GetAFissionFragment(eKinetic,fragA_Z,fragA_A,fragA_M);
      G4int fragB_Z=(G4int)theBaseZ-fragA_Z;
      G4int fragB_A=(G4int)theBaseA-fragA_A-Prompt;

      G4IonTable* pt = G4IonTable::GetIonTable();
      //Excitation energy is not taken into account 
      G4ParticleDefinition* pdA = pt->GetIon( fragA_Z , fragA_A , 0.0 );
      G4ParticleDefinition* pdB = pt->GetIon( fragB_Z , fragB_A , 0.0 );

      //Isotropic Distribution 
      G4double phi = twopi*G4UniformRand();
      // Bug #1745 DHW G4double theta = pi*G4UniformRand();
      G4double costheta = 2.*G4UniformRand()-1.;
      G4double theta = std::acos(costheta); 
      G4double sinth = std::sin(theta);
      G4ThreeVector direction(sinth*std::cos(phi), sinth*std::sin(phi), costheta);

      // Just use ENDF value for this 
      G4double ER = eDepByFragments; 
      G4double ma = pdA->GetPDGMass();
      G4double mb = pdB->GetPDGMass();
      G4double EA = ER / ( 1 + ma/mb);
      G4double EB = ER - EA;
      G4DynamicParticle* dpA = new G4DynamicParticle( pdA , direction , EA);
      G4DynamicParticle* dpB = new G4DynamicParticle( pdB , -direction , EB);
      theResult.Get()->AddSecondary(dpA, secID);
      theResult.Get()->AddSecondary(dpB, secID);
   }

   return theResult.Get();
 }
