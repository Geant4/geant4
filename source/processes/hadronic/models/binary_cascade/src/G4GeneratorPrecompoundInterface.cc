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
// $Id: G4GeneratorPrecompoundInterface.cc 80152 2014-04-03 14:04:18Z gcosmo $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class file
//
//      History: first implementation
//      HPW, 10DEC 98, the decay part originally written by Gunter Folger
//                in his FTF-test-program.
//
//      M.Kelsey, 28 Jul 2011 -- Replace loop to decay input secondaries
//		with new utility class, simplify cleanup loops
// -----------------------------------------------------------------------------

#include <algorithm>
#include <vector>

#include "G4GeneratorPrecompoundInterface.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticleVector.hh"
#include "G4KineticTrackVector.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"

#include "G4V3DNucleus.hh"
#include "G4Nucleon.hh"

#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "G4AntiDeuteron.hh"
#include "G4AntiTriton.hh"
#include "G4AntiHe3.hh"
#include "G4AntiAlpha.hh"

#include "G4FragmentVector.hh"
#include "G4ReactionProduct.hh"
#include "G4ReactionProductVector.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4DecayKineticTracks.hh"
#include "G4HadronicInteractionRegistry.hh"

G4GeneratorPrecompoundInterface::G4GeneratorPrecompoundInterface(G4VPreCompoundModel* preModel)
: CaptureThreshold(10*MeV)
{
   proton = G4Proton::Proton();
   neutron = G4Neutron::Neutron();

   deuteron=G4Deuteron::Deuteron();
   triton  =G4Triton::Triton();
   He3     =G4He3::He3();
   He4     =G4Alpha::Alpha();

   ANTIproton=G4AntiProton::AntiProton();
   ANTIneutron=G4AntiNeutron::AntiNeutron();

   ANTIdeuteron=G4AntiDeuteron::AntiDeuteron();
   ANTItriton  =G4AntiTriton::AntiTriton();
   ANTIHe3     =G4AntiHe3::AntiHe3();
   ANTIHe4     =G4AntiAlpha::AntiAlpha();

   if(preModel) { SetDeExcitation(preModel); }
   else  {
      G4HadronicInteraction* hadi =
            G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
      G4VPreCompoundModel* pre = static_cast<G4VPreCompoundModel*>(hadi);
      if(!pre) { pre = new G4PreCompoundModel(); }
      SetDeExcitation(pre);
   }
}

G4GeneratorPrecompoundInterface::~G4GeneratorPrecompoundInterface()
{
}

//---------------------------------------------------------------------
// choose to calculate excitation energy from energy balance
#define exactExcitationEnergy
//#define debugPrecoInt
//#define G4GPI_debug_excitation

G4ReactionProductVector* G4GeneratorPrecompoundInterface::
Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus)
{
   G4ReactionProductVector * theTotalResult = new G4ReactionProductVector;

   // decay the strong resonances
   G4DecayKineticTracks decay(theSecondaries);

   // prepare the fragment
   G4int anA=theNucleus->GetMassNumber();
   G4int aZ=theNucleus->GetCharge();
   G4int numberOfEx = 0;
   G4int numberOfCh = 0;
   G4int numberOfHoles = 0;
   G4double exEnergy = 0.0;
   G4double R = theNucleus->GetNuclearRadius();
   G4ThreeVector exciton3Momentum(0.,0.,0.);
   G4ThreeVector captured3Momentum(0.,0.,0.);
   G4ThreeVector wounded3Momentum(0.,0.,0.);

   // loop over secondaries
#ifdef exactExcitationEnergy
   G4LorentzVector secondary4Momemtum(0,0,0,0);
#endif
   G4KineticTrackVector::iterator iter;
   for(iter=theSecondaries->begin(); iter !=theSecondaries->end(); ++iter)
   {
      const G4ParticleDefinition* part = (*iter)->GetDefinition();
      G4double e = (*iter)->Get4Momentum().e();
      G4double mass = (*iter)->Get4Momentum().mag();
      G4ThreeVector mom = (*iter)->Get4Momentum().vect();
      if((part != proton && part != neutron) ||
            (e > mass + CaptureThreshold) ||
            ((*iter)->GetPosition().mag() > R)) {
         G4ReactionProduct * theNew = new G4ReactionProduct(part);
         theNew->SetMomentum(mom);
         theNew->SetTotalEnergy(e);
         theTotalResult->push_back(theNew);
#ifdef exactExcitationEnergy
         secondary4Momemtum += (*iter)->Get4Momentum();
#endif
      } else {
         // within the nucleus, neutron or proton
         // now calculate  A, Z of the fragment, momentum, number of exciton states
         ++anA;
         ++numberOfEx;
         G4int Z = G4int(part->GetPDGCharge()/eplus + 0.1);
         aZ += Z;
         numberOfCh += Z;
         captured3Momentum += mom;
         exEnergy += (e - mass);
      }
      delete (*iter);
   }
   delete theSecondaries;

   // loop over wounded nucleus
   G4Nucleon * theCurrentNucleon =
         theNucleus->StartLoop() ? theNucleus->GetNextNucleon() : 0;
   while(theCurrentNucleon) {
      if(theCurrentNucleon->AreYouHit()) {
         ++numberOfHoles;
         ++numberOfEx;
         --anA;
         aZ -= G4int(theCurrentNucleon->GetDefinition()->GetPDGCharge()/eplus + 0.1);
         wounded3Momentum += theCurrentNucleon->Get4Momentum().vect();
         //G4cout << "hit nucleon " << theCurrentNucleon->Get4Momentum() << G4endl;
         exEnergy += theCurrentNucleon->GetBindingEnergy();
      }
      theCurrentNucleon = theNucleus->GetNextNucleon();
   }
   exciton3Momentum = captured3Momentum - wounded3Momentum;

   if(anA == 0) return theTotalResult;

   if(anA >= aZ)
   {
      G4double fMass =  G4NucleiProperties::GetNuclearMass(anA, aZ);

#ifdef exactExcitationEnergy
      // recalculate exEnergy from Energy balance....
      const G4HadProjectile * primary = GetPrimaryProjectile();
      G4double Einitial= primary->Get4Momentum().e()
            											            + G4NucleiProperties::GetNuclearMass(theNucleus->GetMassNumber(),
            											                  theNucleus->GetCharge());
      // Uzhi        G4double Efinal = fMass + secondary4Momemtum.e();
      G4double Efinal = std::sqrt(exciton3Momentum.mag2() + fMass*fMass)
      + secondary4Momemtum.e();
      if ( (Einitial - Efinal) > 0 ) {
         // G4cout << "G4GPI::Propagate() : positive exact excitation Energy "
         //        << (Einitial - Efinal)/MeV << " MeV, exciton estimate "
         //        << exEnergy/MeV << " MeV" << G4endl;

         //          exEnergy=Einitial - Efinal;
         G4LorentzVector PrimMom=primary->Get4Momentum(); PrimMom.setE(Einitial);

         exEnergy=(PrimMom - secondary4Momemtum).mag() - fMass;
      }
      else {
         //  G4cout << "G4GeneratorPrecompoundInterface::Propagate() : "
         //         << "negative exact excitation Energy "
         //         << (Einitial - Efinal)/MeV
         //         << " MeV, setting  excitation to 0 MeV" << G4endl;
         exEnergy=0.;
      }
#endif

      if(exEnergy < 0.) exEnergy=0.;   // Uzhi 11 Dec. 2012

      fMass += exEnergy;

      G4ThreeVector balance=primary->Get4Momentum().vect() -
            secondary4Momemtum.vect() - exciton3Momentum;

#ifdef G4GPI_debug_excitation
      G4cout << "momentum balance" << balance
            << " value " << balance.mag()                   <<G4endl
            << "primary         "<< primary->Get4Momentum() <<G4endl
            << "secondary       "<< secondary4Momemtum      <<G4endl
            << "captured        "<< captured3Momentum       <<G4endl
            << "wounded         "<< wounded3Momentum        <<G4endl
            << "exciton         "<< exciton3Momentum        <<G4endl
            << "second + exciton"
            << secondary4Momemtum.vect() + exciton3Momentum << G4endl;
#endif
      //#ifdef exactExcitationEnergy
      //        G4LorentzVector exciton4Momentum(exciton3Momentum, fMass);
      //        G4LorentzVector exciton4Momentum(exciton3Momentum,
      //                std::sqrt(exciton3Momentum.mag2() + fMass*fMass));
      //#else
      G4LorentzVector exciton4Momentum(exciton3Momentum,
            std::sqrt(exciton3Momentum.mag2() + fMass*fMass));
      //#endif
      //G4cout<<"exciton4Momentum "<<exciton4Momentum<<G4endl;
      // Need to de-excite the remnant nucleus only if excitation energy > 0.
      G4Fragment anInitialState(anA, aZ, exciton4Momentum);
      anInitialState.SetNumberOfParticles(numberOfEx-numberOfHoles);
      anInitialState.SetNumberOfCharged(numberOfCh);
      anInitialState.SetNumberOfHoles(numberOfHoles);

      G4ReactionProductVector * aPrecoResult =
            theDeExcitation->DeExcite(anInitialState);
      // fill pre-compound part into the result, and return
      theTotalResult->insert(theTotalResult->end(),aPrecoResult->begin(),
            aPrecoResult->end() );
      delete aPrecoResult;
   }

   return theTotalResult;
}

G4HadFinalState* G4GeneratorPrecompoundInterface::
ApplyYourself(const G4HadProjectile &, G4Nucleus & )
{
   G4cout << "G4GeneratorPrecompoundInterface: ApplyYourself interface called stand-allone."
         << G4endl;
   G4cout << "This class is only a mediator between generator and precompound"<<G4endl;
   G4cout << "Please remove from your physics list."<<G4endl;
   throw G4HadronicException(__FILE__, __LINE__, "SEVERE: G4GeneratorPrecompoundInterface model interface called stand-allone.");
   return new G4HadFinalState;
}
void G4GeneratorPrecompoundInterface::PropagateModelDescription(std::ostream& outFile) const
{
   outFile << "G4GeneratorPrecompoundInterface interfaces a high\n"
         << "energy model through the wounded nucleus to precompound de-excition.\n"
         << "Low energy protons and neutron present among secondaries produced by \n"
         << "the high energy generator and within the nucleus are captured. The wounded\n"
         << "nucleus and the captured particles form an excited nuclear fragment. This\n"
         << "fragment is passed to the Geant4 pre-compound model for de-excitation.\n"
         << "Nuclear de-excitation:\n";
   // preco

}

// Uzhi Nov. 2012 ------------------------------------------------
G4ReactionProductVector* G4GeneratorPrecompoundInterface::
PropagateNuclNucl(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus,
      G4V3DNucleus* theProjectileNucleus)
{
#ifdef debugPrecoInt
   G4cout<<"G4GeneratorPrecompoundInterface::PropagateNuclNucl "<<G4endl;
#endif

   G4ReactionProductVector * theTotalResult = new G4ReactionProductVector;

   // prepare the target residual
   G4int anA=theNucleus->GetMassNumber();
   G4int aZ=theNucleus->GetCharge();
   G4int numberOfEx = 0;
   G4int numberOfCh = 0;
   G4int numberOfHoles = 0;
   G4double exEnergy = 0.0;
   G4double R = theNucleus->GetNuclearRadius();
   G4LorentzVector Target4Momentum(0,0,0,0);

#ifdef debugPrecoInt
   G4cout<<"Target A Z "<<anA<<" "<<aZ<<G4endl;
#endif

   // loop over wounded target nucleus
   G4Nucleon * theCurrentNucleon =
         theNucleus->StartLoop() ? theNucleus->GetNextNucleon() : 0;
   while(theCurrentNucleon) {
      if(theCurrentNucleon->AreYouHit()) {
         ++numberOfHoles;
         ++numberOfEx;
         --anA;
         aZ -= G4int(theCurrentNucleon->GetDefinition()->GetPDGCharge()/
               eplus + 0.1);
         exEnergy += theCurrentNucleon->GetBindingEnergy();
         Target4Momentum -=theCurrentNucleon->Get4Momentum();
      }
      theCurrentNucleon = theNucleus->GetNextNucleon();
   }

#ifdef debugPrecoInt
   G4cout<<"Residual Target A Z E* 4mom "<<anA<<" "<<aZ<<" "<<exEnergy<<" "
         <<Target4Momentum<<G4endl;
#endif

   // prepare the projectile residual
#ifdef debugPrecoInt
   G4cout<<"Primary BaryonNumber "
         <<GetPrimaryProjectile()->GetDefinition()->GetBaryonNumber()<<G4endl;
#endif

   G4bool ProjectileIsAntiNucleus=
         GetPrimaryProjectile()->GetDefinition()->GetBaryonNumber() < -1;

   G4ThreeVector bst = GetPrimaryProjectile()->Get4Momentum().boostVector();

   G4int anAb=theProjectileNucleus->GetMassNumber();
   G4int aZb=theProjectileNucleus->GetCharge();
   G4int numberOfExB = 0;
   G4int numberOfChB = 0;
   G4int numberOfHolesB = 0;
   G4double exEnergyB = 0.0;
   G4double Rb = theProjectileNucleus->GetNuclearRadius();
   G4LorentzVector Projectile4Momentum(0,0,0,0);

#ifdef debugPrecoInt
   G4cout<<"Projectile A Z "<<anAb<<" "<<aZb<<G4endl;
#endif

   // loop over wounded projectile nucleus
   theCurrentNucleon =
         theProjectileNucleus->StartLoop() ? theProjectileNucleus->GetNextNucleon() : 0;
   while(theCurrentNucleon) {
      if(theCurrentNucleon->AreYouHit()) {
         ++numberOfHolesB;
         ++numberOfExB;
         --anAb;
         if(!ProjectileIsAntiNucleus) {
            aZb -= G4int(theCurrentNucleon->GetDefinition()->GetPDGCharge()/
                  eplus + 0.1);
         } else {
            aZb += G4int(theCurrentNucleon->GetDefinition()->GetPDGCharge()/
                  eplus - 0.1);
         }
         exEnergyB += theCurrentNucleon->GetBindingEnergy();
         Projectile4Momentum -=theCurrentNucleon->Get4Momentum();
      }
      theCurrentNucleon = theProjectileNucleus->GetNextNucleon();
   }

   G4bool ExistTargetRemnant   =  G4double (numberOfHoles)        <
         0.3* G4double (numberOfHoles + anA);
   G4bool ExistProjectileRemnant= G4double (numberOfHolesB)       <
         0.3*G4double (numberOfHolesB + anAb);

#ifdef debugPrecoInt
   G4cout<<"Projectile residual A Z E* 4mom "<<anAb<<" "<<aZb<<" "<<exEnergyB<<" "
         <<Projectile4Momentum<<G4endl;
   G4cout<<" ExistTargetRemnant ExistProjectileRemnant "
         <<ExistTargetRemnant<<" "<< ExistProjectileRemnant<<G4endl;
#endif
   //-----------------------------------------------------------------------------
   // decay the strong resonances
   G4DecayKineticTracks decay(theSecondaries);

#ifdef debugPrecoInt
   G4LorentzVector secondary4Momemtum(0,0,0,0);
   G4int SecondrNum(0);
#endif

   // loop over secondaries
   G4KineticTrackVector::iterator iter;
   for(iter=theSecondaries->begin(); iter !=theSecondaries->end(); ++iter)
   {
      const G4ParticleDefinition* part = (*iter)->GetDefinition();
      G4LorentzVector aTrack4Momentum=(*iter)->Get4Momentum();

      if( part != proton && part != neutron &&
            (part != ANTIproton  && ProjectileIsAntiNucleus) &&
            (part != ANTIneutron && ProjectileIsAntiNucleus)   )
      {
         G4ReactionProduct * theNew = new G4ReactionProduct(part);
         theNew->SetMomentum(aTrack4Momentum.vect());
         theNew->SetTotalEnergy(aTrack4Momentum.e());
         theTotalResult->push_back(theNew);
#ifdef debugPrecoInt
         SecondrNum++;
         secondary4Momemtum += (*iter)->Get4Momentum();
         G4cout<<"Secondary  "<<SecondrNum<<" "
               <<theNew->GetDefinition()->GetParticleName()<<" "
               <<secondary4Momemtum<<G4endl;
#endif
         delete (*iter);
         continue;
      }

      G4bool CanBeCapturedByTarget = false;
      if( part == proton || part == neutron)
      {
         CanBeCapturedByTarget = ExistTargetRemnant    &&
               (CaptureThreshold >
   (aTrack4Momentum       + Target4Momentum).mag() -
   aTrack4Momentum.mag() - Target4Momentum.mag())   &&
   ((*iter)->GetPosition().mag() < R);
      }
      // ---------------------------
      G4LorentzVector Position((*iter)->GetPosition(),
            (*iter)->GetFormationTime());
      Position.boost(bst);

      G4bool CanBeCapturedByProjectile = false;

      if( !ProjectileIsAntiNucleus &&
            ( part == proton || part == neutron))
      {
         CanBeCapturedByProjectile = ExistProjectileRemnant &&
               (CaptureThreshold >
         (aTrack4Momentum       + Projectile4Momentum).mag() -
         aTrack4Momentum.mag() - Projectile4Momentum.mag())    &&
         (Position.vect().mag() < Rb);
      }

      if( ProjectileIsAntiNucleus &&
            ( part == ANTIproton || part == ANTIneutron))
      {
         CanBeCapturedByProjectile = ExistProjectileRemnant &&
               (CaptureThreshold >
         (aTrack4Momentum       + Projectile4Momentum).mag() -
         aTrack4Momentum.mag() - Projectile4Momentum.mag())    &&
         (Position.vect().mag() < Rb);
      }

      if(CanBeCapturedByTarget && CanBeCapturedByProjectile)
      {
         if(G4UniformRand() < 0.5)
         { CanBeCapturedByTarget = true; CanBeCapturedByProjectile = false;}
         else
         { CanBeCapturedByTarget = false; CanBeCapturedByProjectile = true;}
      }

      if(CanBeCapturedByTarget)
      {
         // within the target nucleus, neutron or proton
         // now calculate  A, Z of the fragment, momentum,
         // number of exciton states
#ifdef debugPrecoInt
         G4cout<<"Track is CapturedByTarget "<<" "
               <<aTrack4Momentum<<" "<<aTrack4Momentum.mag()<<G4endl;
#endif
         ++anA;
         ++numberOfEx;
         G4int Z = G4int(part->GetPDGCharge()/eplus + 0.1);
         aZ += Z;
         numberOfCh += Z;
         Target4Momentum +=aTrack4Momentum;
         delete (*iter);
      } else if(CanBeCapturedByProjectile)
      {
         // within the projectile nucleus, neutron or proton
         // now calculate  A, Z of the fragment, momentum,
         // number of exciton states
#ifdef debugPrecoInt
         G4cout<<"Track is CapturedByProjectile"<<" "
               <<aTrack4Momentum<<" "<<aTrack4Momentum.mag()<<G4endl;
#endif
         ++anAb;
         ++numberOfExB;
         G4int Z = G4int(part->GetPDGCharge()/eplus + 0.1);
         if( ProjectileIsAntiNucleus ) Z=-Z;
         aZb += Z;
         numberOfChB += Z;
         Projectile4Momentum +=aTrack4Momentum;
         delete (*iter);
      } else
      { // the track is not captured
         G4ReactionProduct * theNew = new G4ReactionProduct(part);
         theNew->SetMomentum(aTrack4Momentum.vect());
         theNew->SetTotalEnergy(aTrack4Momentum.e());
         theTotalResult->push_back(theNew);

#ifdef debugPrecoInt
         SecondrNum++;
         secondary4Momemtum += (*iter)->Get4Momentum();
         G4cout<<"Secondary  "<<SecondrNum<<" "
               <<theNew->GetDefinition()->GetParticleName()<<" "
               <<secondary4Momemtum<<G4endl;
#endif
         delete (*iter);
         continue;
      }
   }
   delete theSecondaries;
   //-----------------------------------------------------

#ifdef debugPrecoInt
   G4cout<<"Final target residual A Z E* 4mom "<<anA<<" "<<aZ<<" "
         <<exEnergy<<" "<<Target4Momentum<<G4endl;
#endif

   if(0!=anA )
   {
      G4double fMass =  G4NucleiProperties::GetNuclearMass(anA, aZ);

      if((anA == theNucleus->GetMassNumber()) && (exEnergy <= 0.))
      {Target4Momentum.setE(fMass);}

      G4double RemnMass=Target4Momentum.mag();
      if(RemnMass < fMass)
      {
         RemnMass=fMass + exEnergy;
         Target4Momentum.setE(std::sqrt(Target4Momentum.vect().mag2() +
               RemnMass*RemnMass));
      } else
      { exEnergy=RemnMass-fMass;}

      if( exEnergy < 0.) exEnergy=0.;

      // Need to de-excite the remnant nucleus
      G4Fragment anInitialState(anA, aZ, Target4Momentum);
      anInitialState.SetNumberOfParticles(numberOfEx-numberOfHoles);
      anInitialState.SetNumberOfCharged(numberOfCh);
      anInitialState.SetNumberOfHoles(numberOfHoles);

      G4ReactionProductVector * aPrecoResult =
            theDeExcitation->DeExcite(anInitialState);

      // fill pre-compound part into the result, and return
      for(unsigned int ll=0; ll<aPrecoResult->size(); ++ll)
      {
         theTotalResult->push_back(aPrecoResult->operator[](ll));
#ifdef debugPrecoInt
         G4cout<<"Tr frag "<<aPrecoResult->operator[](ll)->GetDefinition()->GetParticleName()
    										               <<" "<<aPrecoResult->operator[](ll)->GetMomentum()<<G4endl;
#endif
      }
      delete aPrecoResult;
   }

   //-----------------------------------------------------
   if((anAb == theProjectileNucleus->GetMassNumber())&& (exEnergyB <= 0.))
   {Projectile4Momentum = GetPrimaryProjectile()->Get4Momentum();}

#ifdef debugPrecoInt
   G4cout<<"Final projectile residual A Z E* Pmom "<<anAb<<" "<<aZb<<" "
         <<exEnergyB<<" "<<Projectile4Momentum<<G4endl;
#endif

   if(0!=anAb)
   {
      G4double fMass =  G4NucleiProperties::GetNuclearMass(anAb, aZb);
      G4double RemnMass=Projectile4Momentum.mag();

      if(RemnMass < fMass)
      {
         RemnMass=fMass + exEnergyB;
         Projectile4Momentum.setE(std::sqrt(Projectile4Momentum.vect().mag2() +
               RemnMass*RemnMass));
      } else
      { exEnergyB=RemnMass-fMass;}

      if( exEnergyB < 0.) exEnergyB=0.;

      // Need to de-excite the remnant nucleus
      G4Fragment anInitialState(anAb, aZb, Projectile4Momentum);
      anInitialState.SetNumberOfParticles(numberOfExB-numberOfHolesB);
      anInitialState.SetNumberOfCharged(numberOfChB);
      anInitialState.SetNumberOfHoles(numberOfHolesB);

      G4ReactionProductVector * aPrecoResult =
            theDeExcitation->DeExcite(anInitialState);

      // fill pre-compound part into the result, and return
      for(unsigned int ll=0; ll<aPrecoResult->size(); ++ll)
      {
         if(ProjectileIsAntiNucleus)
         {

#ifdef debugPrecoInt
            G4cout<<"aPrecoRes  "<<aPrecoResult->operator[](ll)->GetDefinition()->GetParticleName()
        										            <<" "<<aPrecoResult->operator[](ll)->GetMomentum()
        										            <<" "<<aPrecoResult->operator[](ll)->GetTotalEnergy()
        										            <<" "<<aPrecoResult->operator[](ll)->GetMass()<<G4endl;
#endif

            const G4ParticleDefinition * aFragment=aPrecoResult->operator[](ll)->GetDefinition();
            const G4ParticleDefinition * LastFragment=aFragment;
            if     (aFragment == proton)  {LastFragment=G4AntiProton::AntiProtonDefinition();}
            else if(aFragment == neutron) {LastFragment=G4AntiNeutron::AntiNeutronDefinition();}
            else if(aFragment == deuteron){LastFragment=G4AntiDeuteron::AntiDeuteronDefinition();}
            else if(aFragment == triton)  {LastFragment=G4AntiTriton::AntiTritonDefinition();}
            else if(aFragment == He3)     {LastFragment=G4AntiHe3::AntiHe3Definition();}
            else if(aFragment == He4)     {LastFragment=G4AntiAlpha::AntiAlphaDefinition();}
            else {}

            aPrecoResult->operator[](ll)->SetDefinitionAndUpdateE(LastFragment);
         }

#ifdef debugPrecoInt
         G4cout<<"aPrecoResA "<<aPrecoResult->operator[](ll)->GetDefinition()->GetParticleName()
        										            <<" "<<aPrecoResult->operator[](ll)->GetMomentum()
        										            <<" "<<aPrecoResult->operator[](ll)->GetTotalEnergy()
        										            <<" "<<aPrecoResult->operator[](ll)->GetMass()<<G4endl;
#endif
         theTotalResult->push_back(aPrecoResult->operator[](ll));
      }

      delete aPrecoResult;
   }

   return theTotalResult;
}

// Uzhi Nov. 2012 ------------------------------------------------

