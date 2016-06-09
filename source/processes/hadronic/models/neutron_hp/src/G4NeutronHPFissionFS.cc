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
//
#include "G4NeutronHPFissionFS.hh"
#include "G4Nucleus.hh"
#include "G4DynamicParticleVector.hh"
#include "G4NeutronHPFissionERelease.hh"
 void G4NeutronHPFissionFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
 {
    theFS.Init(A, Z, dirName, aFSType);
    theFC.Init(A, Z, dirName, aFSType);
    theSC.Init(A, Z, dirName, aFSType);
    theTC.Init(A, Z, dirName, aFSType);
    theLC.Init(A, Z, dirName, aFSType);
 }
 G4HadFinalState * G4NeutronHPFissionFS::ApplyYourself(const G4HadProjectile & theTrack)
 {  
// prepare neutron
   theResult.Clear();
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

// set neutron and target in the FS classes 
  theFS.SetNeutron(theNeutron);
  theFS.SetTarget(theTarget);
  theFC.SetNeutron(theNeutron);
  theFC.SetTarget(theTarget);
  theSC.SetNeutron(theNeutron);
  theSC.SetTarget(theTarget);
  theTC.SetNeutron(theNeutron);
  theTC.SetTarget(theTarget);
  theLC.SetNeutron(theNeutron);
  theLC.SetTarget(theTarget);

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
       // take 'U' into account explicitely (see 5.4) in the sampling of energy @@@@
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

   G4double * theDecayConstants;

   if(theNeutrons != 0)
   {
     theDecayConstants = new G4double[delayed];
     G4int nPhotons = 0;
     if(thePhotons!=0) nPhotons = thePhotons->size();
     for(i=0; i<theNeutrons->size(); i++)
     {
       theResult.AddSecondary(theNeutrons->operator[](i));
     }
     delete theNeutrons;  

     G4DynamicParticleVector * theDelayed = 0;
     theDelayed = theFS.ApplyYourself(0, delayed, theDecayConstants);
     for(i=0; i<theDelayed->size(); i++)
     {
       G4double time = -std::log(G4UniformRand())/theDecayConstants[i];
       time += theTrack.GetGlobalTime();
       G4HadSecondary * track = new G4HadSecondary(theDelayed->operator[](i));
       track->SetTime(time);
       theResult.AddSecondary(track);
     }
     delete theDelayed;                  
   }
   else
   {
//    cout << " all = "<<all<<G4endl;
     theFS.SampleNeutronMult(all, Prompt, delayed, eKinetic, 0);
     theDecayConstants = new G4double[delayed];
     if(Prompt==0&&delayed==0) Prompt=all;
     theNeutrons = theFS.ApplyYourself(Prompt, delayed, theDecayConstants);
     G4int nPhotons = 0;
     if(thePhotons!=0) nPhotons = thePhotons->size();
     G4int i0;
     for(i0=0; i0<Prompt; i0++)
     {
       theResult.AddSecondary(theNeutrons->operator[](i0));
     }
     for(i0=Prompt; i0<Prompt+delayed; i0++)
     {
       G4double time = -std::log(G4UniformRand())/theDecayConstants[i0-Prompt];
       time += theTrack.GetGlobalTime();        
       //G4HadSecondary * track = new G4HadSecondary(theNeutrons->operator[](i)); this line will be delete
       G4HadSecondary * track = new G4HadSecondary( theNeutrons->operator[]( i0 ) );
       track->SetTime(time);
       theResult.AddSecondary(track);
     }
     delete theNeutrons;   
   }
   delete [] theDecayConstants;
//    cout << "all delayed "<<delayed<<G4endl; 
   unsigned int nPhotons = 0;
   if(thePhotons!=0)
   {
     nPhotons = thePhotons->size();
     for(i=0; i<nPhotons; i++)
     {
       theResult.AddSecondary(thePhotons->operator[](i));
     }
     delete thePhotons; 
   }

// finally deal with local energy depositions.
//    G4cout <<"Number of secondaries = "<<theResult.GetNumberOfSecondaries()<< G4endl;
//    G4cout <<"Number of photons = "<<nPhotons<<G4endl;
//    G4cout <<"Number of Prompt = "<<Prompt<<G4endl;
//    G4cout <<"Number of delayed = "<<delayed<<G4endl;

   G4NeutronHPFissionERelease * theERelease = theFS.GetEnergyRelease();
   G4double eDepByFragments = theERelease->GetFragmentKinetic();
   theResult.SetLocalEnergyDeposit(eDepByFragments);
//    cout << "local energy deposit" << eDepByFragments<<G4endl;
// clean up the primary neutron
   theResult.SetStatusChange(stopAndKill);
   return &theResult;
 }
