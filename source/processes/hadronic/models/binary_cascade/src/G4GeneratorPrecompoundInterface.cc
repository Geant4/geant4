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
// -----------------------------------------------------------------------------
//      GEANT 4 class file
//
//      History: first implementation
//      HPW, 10DEC 98, the decay part originally written by Gunter Folger
//                in his FTF-test-program.
//
//      M.Kelsey, 28 Jul 2011 -- Replace loop to decay input secondaries
//		with new utility class, simplify cleanup loops
//
//      A.Ribon, 27 Oct 2021 -- Extended the method PropagateNuclNucl
//               to deal with projectile hypernuclei and anti-hypernuclei
//
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
#include "G4Lambda.hh"

#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"

#include "G4V3DNucleus.hh"
#include "G4Nucleon.hh"

#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "G4AntiLambda.hh"
#include "G4AntiDeuteron.hh"
#include "G4AntiTriton.hh"
#include "G4AntiHe3.hh"
#include "G4AntiAlpha.hh"

#include "G4HyperTriton.hh"
#include "G4HyperH4.hh"
#include "G4HyperAlpha.hh"
#include "G4HyperHe5.hh"
#include "G4DoubleHyperH4.hh"
#include "G4DoubleHyperDoubleNeutron.hh"

#include "G4AntiHyperTriton.hh"
#include "G4AntiHyperH4.hh"
#include "G4AntiHyperAlpha.hh"
#include "G4AntiHyperHe5.hh"
#include "G4AntiDoubleHyperH4.hh"
#include "G4AntiDoubleHyperDoubleNeutron.hh"

#include "G4FragmentVector.hh"
#include "G4ReactionProduct.hh"
#include "G4ReactionProductVector.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4DecayKineticTracks.hh"
#include "G4HadronicInteractionRegistry.hh"

#include "G4PhysicsModelCatalog.hh"
#include "G4HyperNucleiProperties.hh"
//---------------------------------------------------------------------
#include "Randomize.hh"
#include "G4Log.hh"

//#define debugPrecoInt

G4GeneratorPrecompoundInterface::G4GeneratorPrecompoundInterface(G4VPreCompoundModel* preModel)
  : CaptureThreshold(70*MeV), DeltaM(5.0*MeV), DeltaR(0.0), secID(-1)
{
   proton = G4Proton::Proton();
   neutron = G4Neutron::Neutron();
   lambda = G4Lambda::Lambda();

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

   secID = G4PhysicsModelCatalog::GetModelID("model_PRECO");
}

G4GeneratorPrecompoundInterface::~G4GeneratorPrecompoundInterface()
{
}

G4ReactionProductVector* G4GeneratorPrecompoundInterface::
Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus)
{
   #ifdef debugPrecoInt
      G4cout<<G4endl<<"G4GeneratorPrecompoundInterface::Propagate"<<G4endl;
      G4cout<<"Target A and Z "<<theNucleus->GetMassNumber()<<" "<<theNucleus->GetCharge()<<G4endl;
      G4cout<<"Directly produced particles number "<<theSecondaries->size()<<G4endl;
   #endif

   G4ReactionProductVector * theTotalResult = new G4ReactionProductVector;

   // decay the strong resonances
   G4DecayKineticTracks decay(theSecondaries);
   #ifdef debugPrecoInt
      G4cout<<"Final stable particles number "<<theSecondaries->size()<<G4endl;
   #endif

   // prepare the fragment (it is assumed that target nuclei are never hypernuclei)
   G4int anA=theNucleus->GetMassNumber();
   G4int aZ=theNucleus->GetCharge();
// G4double TargetNucleusMass =  G4NucleiProperties::GetNuclearMass(anA, aZ);

   G4int numberOfEx = 0;
   G4int numberOfCh = 0;
   G4int numberOfHoles = 0;

   G4double R = theNucleus->GetNuclearRadius();

   G4LorentzVector captured4Momentum(0.,0.,0.,0.);
   G4LorentzVector Residual4Momentum(0.,0.,0.,0.);  // TargetNucleusMass is not need at the moment 
   G4LorentzVector Secondary4Momentum(0.,0.,0.,0.);

   // loop over secondaries
   G4KineticTrackVector::iterator iter;
   for(iter=theSecondaries->begin(); iter !=theSecondaries->end(); ++iter)
   {
      const G4ParticleDefinition* part = (*iter)->GetDefinition();
      G4double e = (*iter)->Get4Momentum().e();
      G4double mass = (*iter)->Get4Momentum().mag();
      G4ThreeVector mom = (*iter)->Get4Momentum().vect();
      if((part != proton && part != neutron) ||
            ((*iter)->GetPosition().mag() > R)) {
         G4ReactionProduct * theNew = new G4ReactionProduct(part);
         theNew->SetMomentum(mom);
         theNew->SetTotalEnergy(e);
	 theNew->SetCreatorModelID((*iter)->GetCreatorModelID());
	 theNew->SetParentResonanceDef((*iter)->GetParentResonanceDef());
	 theNew->SetParentResonanceID((*iter)->GetParentResonanceID());
         theTotalResult->push_back(theNew);
         Secondary4Momentum += (*iter)->Get4Momentum();
         #ifdef debugPrecoInt
            G4cout<<"Secondary 4Mom "<<part->GetParticleName()<<" "<<(*iter)->Get4Momentum()<<" "
                  <<(*iter)->Get4Momentum().mag()<<G4endl;
         #endif
      } else {
         if( e-mass > -CaptureThreshold*G4Log( G4UniformRand()) ) {
            G4ReactionProduct * theNew = new G4ReactionProduct(part);
            theNew->SetMomentum(mom);
            theNew->SetTotalEnergy(e);
	    theNew->SetCreatorModelID((*iter)->GetCreatorModelID());
	    theNew->SetParentResonanceDef((*iter)->GetParentResonanceDef());
	    theNew->SetParentResonanceID((*iter)->GetParentResonanceID());	    
            theTotalResult->push_back(theNew);
            Secondary4Momentum += (*iter)->Get4Momentum();
            #ifdef debugPrecoInt
               G4cout<<"Secondary 4Mom "<<part->GetParticleName()<<" "<<(*iter)->Get4Momentum()<<" "
                     <<(*iter)->Get4Momentum().mag()<<G4endl;
            #endif
         } else {
            // within the nucleus, neutron or proton
            // now calculate  A, Z of the fragment, momentum, number of exciton states
            ++anA;
            ++numberOfEx;
            G4int Z = G4int(part->GetPDGCharge()/eplus + 0.1);
            aZ += Z;
            numberOfCh += Z;
            captured4Momentum += (*iter)->Get4Momentum();
            #ifdef debugPrecoInt
               G4cout<<"Captured  4Mom "<<part->GetParticleName()<<(*iter)->Get4Momentum()<<G4endl;
            #endif
         }
      }
      delete (*iter);
   }
   delete theSecondaries;

   // loop over wounded nucleus
   G4Nucleon * theCurrentNucleon =
         theNucleus->StartLoop() ? theNucleus->GetNextNucleon() : 0;
   while(theCurrentNucleon)    /* Loop checking, 31.08.2015, G.Folger */
	{
      if(theCurrentNucleon->AreYouHit()) {
         ++numberOfHoles;
         ++numberOfEx;
         --anA;
         aZ -= G4int(theCurrentNucleon->GetDefinition()->GetPDGCharge()/eplus + 0.1);

         Residual4Momentum -= theCurrentNucleon->Get4Momentum();
      }
      theCurrentNucleon = theNucleus->GetNextNucleon();
   }

   #ifdef debugPrecoInt
      G4cout<<G4endl;
      G4cout<<"Secondary 4Mom "<<Secondary4Momentum<<G4endl;
      G4cout<<"Captured  4Mom "<<captured4Momentum<<G4endl;
      G4cout<<"Sec + Captured "<<Secondary4Momentum+captured4Momentum<<G4endl;
      G4cout<<"Residual4Mom   "<<Residual4Momentum<<G4endl;
      G4cout<<"Sum 4 momenta  "
            <<Secondary4Momentum + captured4Momentum + Residual4Momentum <<G4endl;
   #endif

   // Check that we use QGS model; loop over wounded nucleus
   G4bool QGSM(false);
   theCurrentNucleon = theNucleus->StartLoop() ? theNucleus->GetNextNucleon() : 0;
   while(theCurrentNucleon)   /* Loop checking, 31.08.2015, G.Folger */
	{
      if(theCurrentNucleon->AreYouHit()) 
      {
       if(theCurrentNucleon->Get4Momentum().mag() < 
          theCurrentNucleon->GetDefinition()->GetPDGMass()) QGSM=true;
      }
      theCurrentNucleon = theNucleus->GetNextNucleon();
   }

   #ifdef debugPrecoInt
      if(!QGSM){
        G4cout<<G4endl;
        G4cout<<"Residual A and Z "<<anA<<" "<<aZ<<G4endl;
        G4cout<<"Residual  4Mom "<<Residual4Momentum<<G4endl;
        if(numberOfEx == 0) 
        {G4cout<<"Residual  4Mom = 0 means that there were not wounded and captured nucleons"<<G4endl;}
      }
   #endif

   if(anA == 0) return theTotalResult;

   G4LorentzVector exciton4Momentum(0.,0.,0.,0.);
   if(anA >= aZ)
   {
    if(!QGSM)
    {          // FTF model was used  
      G4double fMass =  G4NucleiProperties::GetNuclearMass(anA, aZ);

//      G4LorentzVector exciton4Momentum = Residual4Momentum + captured4Momentum;
      exciton4Momentum = Residual4Momentum + captured4Momentum;
//exciton4Momentum.setE(std::sqrt(exciton4Momentum.vect().mag2()+sqr(fMass)));
      G4double ActualMass = exciton4Momentum.mag();      
      if(ActualMass <= fMass ) {
        exciton4Momentum.setE(std::sqrt(exciton4Momentum.vect().mag2()+sqr(fMass)));
      }

      #ifdef debugPrecoInt
         G4double exEnergy = 0.0;
         if(ActualMass <= fMass ) {exEnergy = 0.;}
         else                     {exEnergy = ActualMass - fMass;}
         G4cout<<"Ground state residual Mass "<<fMass<<" E* "<<exEnergy<<G4endl;
      #endif
    }
    else
    {          // QGS model was used 
     G4double InitialTargetMass = 
              G4NucleiProperties::GetNuclearMass(theNucleus->GetMassNumber(), theNucleus->GetCharge());

     exciton4Momentum = 
              GetPrimaryProjectile()->Get4Momentum() + G4LorentzVector(0.,0.,0.,InitialTargetMass)
             -Secondary4Momentum;

     G4double fMass =  G4NucleiProperties::GetNuclearMass(anA, aZ);
     G4double ActualMass = exciton4Momentum.mag();

     #ifdef debugPrecoInt
        G4cout<<G4endl;
        G4cout<<"Residual A and Z "<<anA<<" "<<aZ<<G4endl;
        G4cout<<"Residual4Momentum "<<exciton4Momentum<<G4endl;
        G4cout<<"ResidualMass, GroundStateMass and E* "<<ActualMass<<" "<<fMass<<" "
              <<ActualMass - fMass<<G4endl;
     #endif

     if(ActualMass - fMass < 0.)
     {
      G4double ResE = std::sqrt(exciton4Momentum.vect().mag2() + sqr(fMass+10*MeV));
      exciton4Momentum.setE(ResE);
      #ifdef debugPrecoInt
         G4cout<<"ActualMass - fMass < 0. "<<ActualMass<<" "<<fMass<<" "<<ActualMass - fMass<<G4endl;
      #endif
     }
    }

    // Need to de-excite the remnant nucleus only if excitation energy > 0.
    G4Fragment anInitialState(anA, aZ, exciton4Momentum);
    anInitialState.SetNumberOfParticles(numberOfEx-numberOfHoles);
    anInitialState.SetNumberOfCharged(numberOfCh);
    anInitialState.SetNumberOfHoles(numberOfHoles);
    anInitialState.SetCreatorModelID(secID);

    G4ReactionProductVector * aPrecoResult =
      theDeExcitation->DeExcite(anInitialState);
    // fill pre-compound part into the result, and return
    #ifdef debugPrecoInt
    G4cout<<"Target fragment number "<<aPrecoResult->size()<<G4endl;
    #endif
    for(unsigned int ll=0; ll<aPrecoResult->size(); ++ll)
    {
      theTotalResult->push_back(aPrecoResult->operator[](ll));
      #ifdef debugPrecoInt
      G4cout<<"Fragment "<<ll<<" "
	    <<aPrecoResult->operator[](ll)->GetDefinition()->GetParticleName()<<" "
	    <<aPrecoResult->operator[](ll)->GetMomentum()<<" "
	    <<aPrecoResult->operator[](ll)->GetTotalEnergy()<<" "
	    <<aPrecoResult->operator[](ll)->GetDefinition()->GetPDGMass()<<G4endl;
       #endif
    }
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
         << "energy model through the wounded nucleus to precompound de-excitation.\n"
         << "Low energy protons and neutron present among secondaries produced by \n"
         << "the high energy generator and within the nucleus are captured. The wounded\n"
         << "nucleus and the captured particles form an excited nuclear fragment. This\n"
         << "fragment is passed to the Geant4 pre-compound model for de-excitation.\n"
         << "Nuclear de-excitation:\n";
   // preco

}


G4ReactionProductVector* G4GeneratorPrecompoundInterface::
PropagateNuclNucl(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus,
      G4V3DNucleus* theProjectileNucleus)
{
#ifdef debugPrecoInt
   G4cout<<G4endl<<"G4GeneratorPrecompoundInterface::PropagateNuclNucl "<<G4endl;
   G4cout<<"Projectile A and Z (and numberOfLambdas) "<<theProjectileNucleus->GetMassNumber()<<" "
	 <<theProjectileNucleus->GetCharge()<<" ("
	 <<theProjectileNucleus->GetNumberOfLambdas()<<")"<<G4endl;
   G4cout<<"Target     A and Z "<<theNucleus->GetMassNumber()<<" "
         <<theNucleus->GetCharge()<<" ("
	 <<theNucleus->GetNumberOfLambdas()<<")"<<G4endl;
   G4cout<<"Directly produced particles number "<<theSecondaries->size()<<G4endl;
   G4cout<<"Projectile 4Mom and mass "<<GetPrimaryProjectile()->Get4Momentum()<<" "
                                      <<GetPrimaryProjectile()->Get4Momentum().mag()<<G4endl<<G4endl;
#endif

   // prepare the target residual (assumed to be never a hypernucleus)
   G4int anA=theNucleus->GetMassNumber();
   G4int aZ=theNucleus->GetCharge();
   //G4int aL=theNucleus->GetNumberOfLambdas();  // Should be 0
   G4int numberOfEx = 0;
   G4int numberOfCh = 0;
   G4int numberOfHoles = 0;
   G4double exEnergy = 0.0;
   G4double R = theNucleus->GetNuclearRadius();
   G4LorentzVector Target4Momentum(0.,0.,0.,0.);

   // loop over the wounded target nucleus
   G4Nucleon * theCurrentNucleon =
         theNucleus->StartLoop() ? theNucleus->GetNextNucleon() : 0;
   while(theCurrentNucleon)   /* Loop checking, 31.08.2015, G.Folger */
	{
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
   G4cout<<"Residual Target A Z (numberOfLambdas) E* 4mom "<<anA<<" "<<aZ<<" (0"//<<aL
         <<") "<<exEnergy<<" "<<Target4Momentum<<G4endl;
#endif

   // prepare the projectile residual - which can be a hypernucleus or anti-hypernucleus

   G4bool ProjectileIsAntiNucleus=
         GetPrimaryProjectile()->GetDefinition()->GetBaryonNumber() < -1;

   G4ThreeVector bst = GetPrimaryProjectile()->Get4Momentum().boostVector();

   G4int anAb=theProjectileNucleus->GetMassNumber();
   G4int aZb=theProjectileNucleus->GetCharge();
   G4int aLb=theProjectileNucleus->GetNumberOfLambdas();  // Non negative number of (anti-)lambdas in (anti-)nucleus
   G4int numberOfExB = 0;
   G4int numberOfChB = 0;
   G4int numberOfHolesB = 0;
   G4double exEnergyB = 0.0;
   G4double Rb = theProjectileNucleus->GetNuclearRadius();
   G4LorentzVector Projectile4Momentum(0.,0.,0.,0.);

   // loop over the wounded projectile nucleus or anti-nucleus
   theCurrentNucleon =
         theProjectileNucleus->StartLoop() ? theProjectileNucleus->GetNextNucleon() : 0;
   while(theCurrentNucleon)    /* Loop checking, 31.08.2015, G.Folger */
	{
      if(theCurrentNucleon->AreYouHit()) {
         ++numberOfHolesB;
         ++numberOfExB;
         --anAb;
         if(!ProjectileIsAntiNucleus) {
            aZb -= G4int(theCurrentNucleon->GetDefinition()->GetPDGCharge()/
                  eplus + 0.1);
	    if (theCurrentNucleon->GetParticleType()==G4Lambda::Definition()) --aLb;
         } else {
            aZb += G4int(theCurrentNucleon->GetDefinition()->GetPDGCharge()/
                  eplus - 0.1);
	    if (theCurrentNucleon->GetParticleType()==G4AntiLambda::Definition()) --aLb;
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
   G4cout<<"Projectile residual A Z (numberOfLambdas) E* 4mom "<<anAb<<" "<<aZb<<" ("<<aLb
	 <<") "<<exEnergyB<<" "<<Projectile4Momentum<<G4endl;
   G4cout<<" ExistTargetRemnant ExistProjectileRemnant "
         <<ExistTargetRemnant<<" "<< ExistProjectileRemnant<<G4endl;
#endif
   //-----------------------------------------------------------------------------
   // decay the strong resonances
   G4ReactionProductVector * theTotalResult = new G4ReactionProductVector;
   G4DecayKineticTracks decay(theSecondaries);

   MakeCoalescence(theSecondaries);

   #ifdef debugPrecoInt
      G4cout<<"Secondary stable particles number "<<theSecondaries->size()<<G4endl;
   #endif

#ifdef debugPrecoInt
   G4LorentzVector secondary4Momemtum(0,0,0,0);
   G4int SecondrNum(0);
#endif

   // Loop over secondaries.
   // We are assuming that only protons and neutrons - for nuclei -
   // and only antiprotons and antineutrons - for antinuclei - can be absorbed,
   // not instead lambdas (or hyperons more generally) - for nuclei - or anti-lambdas
   // (or anti-hyperons more generally) - for antinuclei. This is a simplification,
   // to be eventually reviewed later on, in particular when generic hypernuclei and
   // anti-hypernuclei are introduced, instead of the few light hypernuclei and
   // anti-hypernuclei which currently exist in Geant4.
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
	 theNew->SetCreatorModelID((*iter)->GetCreatorModelID());
	 theNew->SetParentResonanceDef((*iter)->GetParentResonanceDef());
	 theNew->SetParentResonanceID((*iter)->GetParentResonanceID());	 
         theTotalResult->push_back(theNew);
#ifdef debugPrecoInt
         SecondrNum++;
         secondary4Momemtum += (*iter)->Get4Momentum();
         G4cout<<"Secondary  "<<SecondrNum<<" "
               <<theNew->GetDefinition()->GetParticleName()<<" "
               <<theNew->GetMomentum()<<" "<<theNew->GetTotalEnergy()<<G4endl;

#endif
         delete (*iter);
         continue;
      }

      G4bool CanBeCapturedByTarget = false;
      if( part == proton || part == neutron)
      {
         CanBeCapturedByTarget = ExistTargetRemnant    &&
               (-CaptureThreshold*G4Log( G4UniformRand()) >
              (aTrack4Momentum       + Target4Momentum).mag() -
               aTrack4Momentum.mag() - Target4Momentum.mag())   &&
                   ((*iter)->GetPosition().mag() < R);
      }
      // ---------------------------
      G4LorentzVector Position((*iter)->GetPosition(), (*iter)->GetFormationTime());
      Position.boost(bst);

      G4bool CanBeCapturedByProjectile = false;

      if( !ProjectileIsAntiNucleus &&
            ( part == proton || part == neutron))
      {
         CanBeCapturedByProjectile = ExistProjectileRemnant &&
               (-CaptureThreshold*G4Log( G4UniformRand()) >
              (aTrack4Momentum       + Projectile4Momentum).mag() -
               aTrack4Momentum.mag() - Projectile4Momentum.mag())    &&
                   (Position.vect().mag() < Rb);
      }

      if( ProjectileIsAntiNucleus &&
            ( part == ANTIproton || part == ANTIneutron))
      {
         CanBeCapturedByProjectile = ExistProjectileRemnant &&
               (-CaptureThreshold*G4Log( G4UniformRand()) >
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
         G4cout<<"Track is CapturedByTarget "<<" "<<part->GetParticleName()<<" "
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
         G4cout<<"Track is CapturedByProjectile"<<" "<<part->GetParticleName()<<" "
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
	 theNew->SetCreatorModelID((*iter)->GetCreatorModelID());
	 theNew->SetParentResonanceDef((*iter)->GetParentResonanceDef());
	 theNew->SetParentResonanceID((*iter)->GetParentResonanceID());	 
         theTotalResult->push_back(theNew);

#ifdef debugPrecoInt
         SecondrNum++;
         secondary4Momemtum += (*iter)->Get4Momentum();
/*
         G4cout<<"Secondary  "<<SecondrNum<<" "
               <<theNew->GetDefinition()->GetParticleName()<<" "
               <<secondary4Momemtum<<G4endl;
*/
#endif
         delete (*iter);
         continue;
      }
   }
   delete theSecondaries;
   //-----------------------------------------------------

   #ifdef debugPrecoInt
   G4cout<<"Final target residual A Z (numberOfLambdas) E* 4mom "<<anA<<" "<<aZ<<" (0"//<<aL
	    <<") "<<exEnergy<<" "<<Target4Momentum<<G4endl;
   #endif

   if(0!=anA )
   {
      // We assume that the target residual is never a hypernucleus
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
      anInitialState.SetCreatorModelID(secID);

      G4ReactionProductVector * aPrecoResult =
            theDeExcitation->DeExcite(anInitialState);

      #ifdef debugPrecoInt
         G4cout<<"Target fragment number "<<aPrecoResult->size()<<G4endl;
      #endif

      // fill pre-compound part into the result, and return
      for(unsigned int ll=0; ll<aPrecoResult->size(); ++ll)
      {
         theTotalResult->push_back(aPrecoResult->operator[](ll));
         #ifdef debugPrecoInt
            G4cout<<"Target fragment "<<ll<<" "
                  <<aPrecoResult->operator[](ll)->GetDefinition()->GetParticleName()<<" "
                  <<aPrecoResult->operator[](ll)->GetMomentum()<<" "
                  <<aPrecoResult->operator[](ll)->GetTotalEnergy()<<" "
                  <<aPrecoResult->operator[](ll)->GetMass()<<G4endl;
         #endif
      }
      delete aPrecoResult;
   }

   //-----------------------------------------------------
   if((anAb == theProjectileNucleus->GetMassNumber())&& (exEnergyB <= 0.))
   {Projectile4Momentum = GetPrimaryProjectile()->Get4Momentum();}

   #ifdef debugPrecoInt
   G4cout<<"Final projectile residual A Z (numberOfLambdas) E* Pmom Pmag2 "<<anAb<<" "<<aZb<<" ("
	 <<aLb<<") "<<exEnergyB<<" "<<Projectile4Momentum<<" "
                            <<Projectile4Momentum.mag2()<<G4endl;
   #endif

   if(0!=anAb)
   {
      // The projectile residual can be a hypernucleus or anti-hypernucleus
      G4double fMass = 0.0;
      if ( aLb > 0 ) {
	fMass = G4HyperNucleiProperties::GetNuclearMass(anAb, aZb, aLb);
      } else {
	fMass = G4NucleiProperties::GetNuclearMass(anAb, aZb);
      }        
      G4double RemnMass=Projectile4Momentum.mag();

      if(RemnMass < fMass)
      {
         RemnMass=fMass + exEnergyB;
         Projectile4Momentum.setE(std::sqrt(Projectile4Momentum.vect().mag2() +
               RemnMass*RemnMass));
      } else
      { exEnergyB=RemnMass-fMass;}

      if( exEnergyB < 0.) exEnergyB=0.;

      G4ThreeVector bstToCM =Projectile4Momentum.findBoostToCM();
      Projectile4Momentum.boost(bstToCM);

      // Need to de-excite the remnant nucleus
      G4Fragment anInitialState(anAb, aZb, aLb, Projectile4Momentum);
      anInitialState.SetNumberOfParticles(numberOfExB-numberOfHolesB);
      anInitialState.SetNumberOfCharged(numberOfChB);
      anInitialState.SetNumberOfHoles(numberOfHolesB);
      anInitialState.SetCreatorModelID(secID);

      G4ReactionProductVector * aPrecoResult =
            theDeExcitation->DeExcite(anInitialState);

      #ifdef debugPrecoInt
         G4cout<<"Projectile fragment number "<<aPrecoResult->size()<<G4endl;
      #endif

      // fill pre-compound part into the result, and return
      for(unsigned int ll=0; ll<aPrecoResult->size(); ++ll)
      {
         G4LorentzVector tmp=G4LorentzVector(aPrecoResult->operator[](ll)->GetMomentum(),
                                             aPrecoResult->operator[](ll)->GetTotalEnergy());
         tmp.boost(-bstToCM); // Transformation to the system of original remnant
         aPrecoResult->operator[](ll)->SetMomentum(tmp.vect());
         aPrecoResult->operator[](ll)->SetTotalEnergy(tmp.e());

         if(ProjectileIsAntiNucleus)
         {
            const G4ParticleDefinition * aFragment=aPrecoResult->operator[](ll)->GetDefinition();
            const G4ParticleDefinition * LastFragment=aFragment;
            if     (aFragment == proton)  {LastFragment=G4AntiProton::AntiProtonDefinition();}
            else if(aFragment == neutron) {LastFragment=G4AntiNeutron::AntiNeutronDefinition();}
            else if(aFragment == lambda)  {LastFragment=G4AntiLambda::AntiLambdaDefinition();}
            else if(aFragment == deuteron){LastFragment=G4AntiDeuteron::AntiDeuteronDefinition();}
            else if(aFragment == triton)  {LastFragment=G4AntiTriton::AntiTritonDefinition();}
            else if(aFragment == He3)     {LastFragment=G4AntiHe3::AntiHe3Definition();}
            else if(aFragment == He4)     {LastFragment=G4AntiAlpha::AntiAlphaDefinition();}
            else {}

            if (aLb > 0) {  // Anti-hypernucleus
	      if        (aFragment == G4HyperTriton::Definition()) {
		LastFragment=G4AntiHyperTriton::Definition();
	      } else if (aFragment == G4HyperH4::Definition()) {
		LastFragment=G4AntiHyperH4::Definition();
	      } else if (aFragment == G4HyperAlpha::Definition()) {
		LastFragment=G4AntiHyperAlpha::Definition();
	      } else if (aFragment == G4HyperHe5::Definition()) {
		LastFragment=G4AntiHyperHe5::Definition();
	      } else if (aFragment == G4DoubleHyperH4::Definition()) {
		LastFragment=G4AntiDoubleHyperH4::Definition();
	      } else if (aFragment == G4DoubleHyperDoubleNeutron::Definition()) {
		LastFragment=G4AntiDoubleHyperDoubleNeutron::Definition();
	      }
	    }
	    
            aPrecoResult->operator[](ll)->SetDefinitionAndUpdateE(LastFragment);
         }

         #ifdef debugPrecoInt
            G4cout<<"Projectile fragment "<<ll<<" "
                  <<aPrecoResult->operator[](ll)->GetDefinition()->GetParticleName()<<" "
                  <<aPrecoResult->operator[](ll)->GetMomentum()<<" "
                  <<aPrecoResult->operator[](ll)->GetTotalEnergy()<<" "
                  <<aPrecoResult->operator[](ll)->GetMass()<<G4endl;
         #endif

         theTotalResult->push_back(aPrecoResult->operator[](ll));
      }

      delete aPrecoResult;
   }

   return theTotalResult;
}


void G4GeneratorPrecompoundInterface::MakeCoalescence(G4KineticTrackVector *tracks) {
  // This method replaces pairs of proton-neutron - in the case of nuclei - or
  // antiproton-antineutron - in the case of anti-nuclei - which are close in
  // momentum, with, respectively, deuterons and anti-deuterons.
  // Note that in the case of hypernuclei or anti-hypernuclei, lambdas or anti-lambdas
  // are not considered for coalescence because hyper-deuteron or anti-hyper-deuteron
  // are assumed not to exist.
  
  if (!tracks) return;

  G4double MassCut = deuteron->GetPDGMass() + DeltaM;   // In MeV

  for ( std::size_t i = 0; i < tracks->size(); ++i ) {  // search for protons

    G4KineticTrack* trackP = (*tracks)[i];
    if ( ! trackP ) continue;
    if (trackP->GetDefinition() != proton) continue;

    G4LorentzVector Prot4Mom = trackP->Get4Momentum();
    G4LorentzVector ProtSPposition = G4LorentzVector(trackP->GetPosition(), trackP->GetFormationTime());

    for ( std::size_t j = 0; j < tracks->size(); ++j ) {  // search for neutron
                                                
      G4KineticTrack* trackN = (*tracks)[j];
      if (! trackN ) continue;
      if (trackN->GetDefinition() != neutron) continue;

      G4LorentzVector Neut4Mom = trackN->Get4Momentum();
      G4LorentzVector NeutSPposition = G4LorentzVector( trackN->GetPosition(), trackN->GetFormationTime()*hbarc/fermi);
      G4double EffMass = (Prot4Mom + Neut4Mom).mag();

      if ( EffMass <= MassCut ) {  // && (EffDistance <= SpaceCut)) { // Create deuteron
        G4KineticTrack* aDeuteron = 
          new G4KineticTrack( deuteron,
                              (trackP->GetFormationTime() +  trackN->GetFormationTime())/2.0,
                              (trackP->GetPosition()      +  trackN->GetPosition()     )/2.0,
                              ( Prot4Mom                  +  Neut4Mom                        ));
        aDeuteron->SetCreatorModelID(secID);
        tracks->push_back(aDeuteron);
        delete trackP; delete trackN;
        (*tracks)[i] = nullptr; (*tracks)[j] = nullptr;
        break;
      }
    }
  }

  // Find and remove null pointers created by decays above
  for ( G4int jj = (G4int)tracks->size()-1; jj >= 0; --jj ) {
    if ( ! (*tracks)[jj] ) tracks->erase(tracks->begin()+jj);
  }
}
