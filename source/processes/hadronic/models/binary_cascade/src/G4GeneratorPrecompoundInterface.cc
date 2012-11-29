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
// $Id$
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
#include "G4V3DNucleus.hh"
#include "G4Nucleon.hh"
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
        G4ParticleDefinition* part = (*iter)->GetDefinition();
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

    if(anA>0 && aZ>0) {
        G4double fMass =  G4NucleiProperties::GetNuclearMass(anA, aZ);
#ifdef exactExcitationEnergy
        // recalculate exEnergy from Energy balance....
        const G4HadProjectile * primary = GetPrimaryProjectile();
        G4double Einitial= primary->Get4Momentum().e()
            		        + G4NucleiProperties::GetNuclearMass(theNucleus->GetMassNumber(),theNucleus->GetCharge());
        G4double Efinal = fMass + secondary4Momemtum.e();
        if ( (Einitial - Efinal) > 0 ) {
            // G4cout << "G4GPI::Propagate() : positive exact excitation Energy "
            //  << (Einitial - Efinal)/MeV << " MeV, exciton estimate " << exEnergy/MeV << " MeV" << G4endl;
            exEnergy=Einitial - Efinal;
        }
        else {
            //  G4cout << "G4GeneratorPrecompoundInterface::Propagate() : negative exact excitation Energy "
            //  << (Einitial - Efinal)/MeV << " MeV, setting  excitation to 0 MeV" << G4endl;
            exEnergy=0.;
        }
#endif
        fMass += exEnergy;
        G4ThreeVector balance=primary->Get4Momentum().vect() - secondary4Momemtum.vect() - exciton3Momentum;
        #ifdef G4GPI_debug_excitation
        G4cout << "momentum balance init/final  " << balance << " value " << balance.mag() << G4endl
                << "primary / secondaries "<< primary->Get4Momentum() << " / "
                << secondary4Momemtum << " captured/wounded: " << captured3Momentum << " / " <<  wounded3Momentum
                << "  exciton " << exciton3Momentum << G4endl
                << secondary4Momemtum.vect() + exciton3Momentum << G4endl;
        #endif
#ifdef exactExcitationEnergy
        G4LorentzVector exciton4Momentum(exciton3Momentum, fMass);
#else
        G4LorentzVector exciton4Momentum(exciton3Momentum,
                std::sqrt(exciton3Momentum.mag2() + fMass*fMass));
#endif
        if ( exEnergy > 0.0 ) {  // Need to de-excite the remnant nucleus only if excitation energy > 0.
            G4Fragment anInitialState(anA, aZ, exciton4Momentum);
            anInitialState.SetNumberOfParticles(numberOfEx-numberOfHoles);
            anInitialState.SetNumberOfCharged(numberOfCh);
            anInitialState.SetNumberOfHoles(numberOfHoles);

            G4ReactionProductVector * aPrecoResult = theDeExcitation->DeExcite(anInitialState);
            // fill pre-compound part into the result, and return
            theTotalResult->insert(theTotalResult->end(),aPrecoResult->begin(),aPrecoResult->end() );
            delete aPrecoResult;

        } else {  // No/negative excitation energy, we only need to create the remnant nucleus
            //  energy is not conserved, ignore exciton momentum, i.e. remnant nucleus will be at rest
            G4ParticleDefinition* theKindOfFragment = 0;
            if (anA == 1 && aZ == 0) {
                theKindOfFragment = G4Neutron::NeutronDefinition();
            } else if (anA == 1 && aZ == 1) {
                theKindOfFragment = G4Proton::ProtonDefinition();
            } else if (anA == 2 && aZ == 1) {
                theKindOfFragment = G4Deuteron::DeuteronDefinition();
            } else if (anA == 3 && aZ == 1) {
                theKindOfFragment = G4Triton::TritonDefinition();
            } else if (anA == 3 && aZ == 2) {
                theKindOfFragment = G4He3::He3Definition();
            } else if (anA == 4 && aZ == 2) {
                theKindOfFragment = G4Alpha::AlphaDefinition();;
            } else {
                theKindOfFragment =
                        G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(aZ,anA,0.0);
            }
            if (theKindOfFragment != 0) {
                G4ReactionProduct * theNew = new G4ReactionProduct(theKindOfFragment);
                theNew->SetMomentum(G4ThreeVector(0.,0.,0.));
                theNew->SetTotalEnergy(fMass);
                //theNew->SetFormationTime(??0.??);
                theTotalResult->push_back(theNew);
            }
        }
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
