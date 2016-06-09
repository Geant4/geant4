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
// $Id: G4PreCompoundModel.cc,v 1.17 2008/12/09 14:09:59 ahoward Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// by V. Lara
//
//J. M. Quesada (Apr.08). Several changes. Soft cut-off switched off. 
//(May. 08). Protection against non-physical preeq. transitional regime has 
// been set
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option
// JMQ (06 September 2008) Also external choices have been added for:
//                      - superimposed Coulomb barrier (useSICB=true) 
//                      - "never go back"  hipothesis (useNGB=true) 
//                      - soft cutoff from preeq. to equlibrium (useSCO=true)
//                      - CEM transition probabilities (useCEMtr=true)  


#include "G4PreCompoundModel.hh"
#include "G4PreCompoundEmission.hh"
#include "G4PreCompoundTransitions.hh"
#include "G4GNASHTransitions.hh"
#include "G4ParticleDefinition.hh"


#ifdef PRECOMPOUND_TEST
G4Fragment G4PreCompoundModel::theInitialFragmentForTest;
std::vector<G4String*> G4PreCompoundModel::theCreatorModels;
#endif

const G4PreCompoundModel & G4PreCompoundModel::operator=(const G4PreCompoundModel &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4PreCompoundModel::operator= meant to not be accessable");
  return *this;
}


G4bool G4PreCompoundModel::operator==(const G4PreCompoundModel &) const
{
  return false;
}

G4bool G4PreCompoundModel::operator!=(const G4PreCompoundModel &) const
{
  return true;
}



// Additional Declarations

G4HadFinalState * G4PreCompoundModel::ApplyYourself(const G4HadProjectile & thePrimary,
                                                    G4Nucleus & theNucleus)
{  
  // prepare fragment
  G4Fragment anInitialState;
  // This si for GNASH transitions
  anInitialState.SetParticleDefinition(const_cast<G4ParticleDefinition *>(thePrimary.GetDefinition()));

  G4int anA=static_cast<G4int>(theNucleus.GetN());
  anA += thePrimary.GetDefinition()->GetBaryonNumber();

  anInitialState.SetA(anA);
  
  G4int aZ=static_cast<G4int>(theNucleus.GetZ());
  aZ += static_cast<G4int>(thePrimary.GetDefinition()->GetPDGCharge());

  anInitialState.SetZ(aZ);
  
  // Assume the projectile is a nucleon
  
  // Number of Excited Particles
  anInitialState.SetNumberOfParticles(1+thePrimary.GetDefinition()->GetBaryonNumber());
  
  // Number of Charged Excited Particles
  // JMQ/AH modify number of charged particles with probability of the Z/A ratio of the nucleus:
  //  if(G4UniformRand() <= aZ/anA) BUG! - integer arithmetic
  if(G4UniformRand() <= (static_cast<G4double>(aZ))/(static_cast<G4double>(anA))) 
      anInitialState.SetNumberOfCharged(static_cast<G4int>(thePrimary.GetDefinition()->GetPDGCharge()+.01) + 1);
  else
      anInitialState.SetNumberOfCharged(static_cast<G4int>(thePrimary.GetDefinition()->GetPDGCharge()+.01));
    
//AH     anInitialState.SetNumberOfCharged(static_cast<G4int>(thePrimary.GetDefinition()->GetPDGCharge()+.01) + 
//AH 				    static_cast<G4int>(0.5+G4UniformRand()));

  // Number of Holes 
  anInitialState.SetNumberOfHoles(1);
  
  // pre-compound nucleus energy.
  G4double anEnergy = 0;
  G4double nucleusMass =  G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(static_cast<G4int>(theNucleus.GetZ()),
                                                                                         static_cast<G4int>(theNucleus.GetN()));
  anEnergy =  nucleusMass + thePrimary.GetTotalEnergy();
  
  // Momentum
  G4ThreeVector p = thePrimary.Get4Momentum().vect();

  // 4-momentum
  G4LorentzVector momentum(p, anEnergy);
  anInitialState.SetMomentum(momentum);
  
#ifdef PRECOMPOUND_TEST
  G4PreCompoundModel::theInitialFragmentForTest = anInitialState;
#endif
  
  // call excitation handler
  const G4Fragment aFragment(anInitialState);
  G4ReactionProductVector * result = DeExcite(aFragment);

#ifdef PRECOMPOUND_TEST
  for (std::vector<G4String*>::iterator icm = theCreatorModels.begin(); 
       icm != theCreatorModels.end(); ++icm )
    {
      delete (*icm);
    }
  theCreatorModels.clear();
#endif
  // fill particle change
  theResult.Clear();
  theResult.SetStatusChange(stopAndKill);
  for(G4ReactionProductVector::iterator i= result->begin(); i != result->end(); ++i)
    {
      G4DynamicParticle * aNew = 
	new G4DynamicParticle((*i)->GetDefinition(),
			      (*i)->GetTotalEnergy(),
			      (*i)->GetMomentum());
#ifdef PRECOMPOUND_TEST
      theCreatorModels.push_back(new G4String((*i)->GetCreatorModel()));
#endif
      delete (*i);
      theResult.AddSecondary(aNew);
    }
  delete result;
  
  //return the filled particle change
  return &theResult;
}


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

G4ReactionProductVector* G4PreCompoundModel::DeExcite(const G4Fragment & theInitialState) const
{
  
  G4ReactionProductVector * Result = new G4ReactionProductVector;
  
  // Copy of the initial state 
  G4Fragment aFragment(theInitialState);

  if (aFragment.GetExcitationEnergy() < 10*eV)
    {
      // Perform Equilibrium Emission
      PerformEquilibriumEmission(aFragment,Result);
      return Result;
    }
  
  if (aFragment.GetA() < 5) {
    G4ReactionProduct * theRP = new G4ReactionProduct(G4ParticleTable::GetParticleTable()->
						      GetIon(static_cast<G4int>(aFragment.GetZ()),
							     static_cast<G4int>(aFragment.GetA()),
							     aFragment.GetExcitationEnergy()));
    theRP->SetMomentum(aFragment.GetMomentum().vect());
    theRP->SetTotalEnergy(aFragment.GetMomentum().e());	  
    Result->push_back(theRP);
    return Result;
  }
  
  G4PreCompoundEmission aEmission;
  if (useHETCEmission) aEmission.SetHETCModel();
  aEmission.SetUp(theInitialState);
  
  //for cross section options
  
  if (OPTxs!= 0 && OPTxs!=1 && OPTxs !=2 && OPTxs !=3 && OPTxs !=4  ) {
    std::ostringstream errOs;
    errOs << "BAD CROSS SECTION OPTION in G4PreCompoundModel.cc !!"  <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, errOs.str());}
  else aEmission.SetOPTxs(OPTxs); 
  
  //for the choice of superimposed Coulomb Barrier for inverse cross sections
  
   aEmission.UseSICB(useSICB);
 
  
  //----------
  
  G4VPreCompoundTransitions * aTransition = 0;
  if (useGNASHTransition) 
    {
      aTransition = new G4GNASHTransitions;
    }
  else 
    {
      aTransition = new G4PreCompoundTransitions;
      // for the choice of "never go back" hypothesis and CEM transition probabilities  
      if (useNGB) aTransition->UseNGB(useNGB);
      if (useCEMtr) aTransition->UseCEMtr(useCEMtr);
    }
  
  // Main loop. It is performed until equilibrium deexcitation.
  //G4int fragment=0;
  
  for (;;) {
    
    //fragment++;
    //G4cout<<"-------------------------------------------------------------------"<<G4endl;
    //G4cout<<"Fragment number .. "<<fragment<<G4endl;
    
    // Initialize fragment according with the nucleus parameters
    aEmission.Initialize(aFragment);
    
    
    
    G4double g = (6.0/pi2)*aFragment.GetA()*
      G4PreCompoundParameters::GetAddress()->GetLevelDensity();
    
    
    
    
    G4int EquilibriumExcitonNumber = static_cast<G4int>(std::sqrt(2.0*g*aFragment.GetExcitationEnergy())+ 0.5);
//   
//    G4cout<<"Neq="<<EquilibriumExcitonNumber<<G4endl;
//
// J. M. Quesada (Jan. 08)  equilibrium hole number could be used as preeq.- evap. delimiter (IAEA report)
//    G4int EquilibriumHoleNumber = static_cast<G4int>(0.2*std::sqrt(g*aFragment.GetExcitationEnergy())+ 0.5);
    
// Loop for transitions, it is performed while there are preequilibrium transitions.
    G4bool ThereIsTransition = false;
    
    //        G4cout<<"----------------------------------------"<<G4endl;
    //        G4double NP=aFragment.GetNumberOfParticles();
    //        G4double NH=aFragment.GetNumberOfHoles();
    //        G4double NE=aFragment.GetNumberOfExcitons();
    //        G4cout<<" Ex. Energy="<<aFragment.GetExcitationEnergy()<<G4endl;
    //        G4cout<<"N. excitons="<<NE<<"  N. Part="<<NP<<"N. Holes ="<<NH<<G4endl;
        
    
    //G4int transition=0;
    do 
      {
        //transition++;
        //G4cout<<"transition number .."<<transition<<G4endl;
        //G4cout<<" n ="<<aFragment.GetNumberOfExcitons()<<G4endl;
        G4bool go_ahead = false;
        // soft cutoff criterium as an "ad-hoc" solution to force increase in  evaporation  
        //       G4double test = static_cast<G4double>(aFragment.GetNumberOfHoles());
        G4double test = static_cast<G4double>(aFragment.GetNumberOfExcitons());
        
        
        if (test < EquilibriumExcitonNumber) go_ahead=true;
        //J. M. Quesada (Apr. 08): soft-cutoff switched off by default
        if (useSCO) {
          if (test < EquilibriumExcitonNumber)
            //  if (test < EquilibriumHoleNumber)
            {
              test /= static_cast<G4double>(EquilibriumExcitonNumber); 
              //     test /= static_cast<G4double>(EquilibriumHoleNumber);
              test -= 1.0;
              test = test*test;
              test /= 0.32;
              test = 1.0 - std::exp(-test);
              go_ahead = (G4UniformRand() < test);
              
            }
        } 
        
        //JMQ: WARNING:  CalculateProbability MUST be called prior to Get methods !! (O values would be returned otherwise)
        G4double TotalTransitionProbability = aTransition->CalculateProbability(aFragment);
        G4double P1=aTransition->GetTransitionProb1();
        G4double P2=aTransition->GetTransitionProb2();
        G4double P3=aTransition->GetTransitionProb3();
        //       G4cout<<"P1="<<P1<<" P2="<<P2<<"  P3="<<P3<<G4endl;
        
        
        //J.M. Quesada (May. 08). Physical criterium (lamdas)  PREVAILS over approximation (critical exciton number)
        if(P1<=(P2+P3)) go_ahead=false;
        
        if (go_ahead &&  aFragment.GetA() > 4)
	  {				
	    G4double TotalEmissionProbability = aEmission.GetTotalProbability(aFragment);
            //
            //  G4cout<<"TotalEmissionProbability="<<TotalEmissionProbability<<G4endl;
            //
	    // Check if number of excitons is greater than 0
	    // else perform equilibrium emission
	    if (aFragment.GetNumberOfExcitons() <= 0) 
	      {
		// Perform Equilibrium Emission
#ifdef debug // ------------- debug -----------------------------------------
		CheckConservation(theInitialState,aFragment,Result);
#endif // ------------------- debug -----------------------------------------
		PerformEquilibriumEmission(aFragment,Result);
		delete aTransition;
		return Result;
	      }
	    
	    //	    G4PreCompoundTransitions aTransition(aFragment);
            
            //J.M.Quesada (May 08) this has already been done in order to decide what to do (preeq-eq) 
	    // Sum of transition probabilities
            //	    G4double TotalTransitionProbability = aTransition->CalculateProbability(aFragment);
            
	    // Sum of all probabilities
	    G4double TotalProbability = TotalEmissionProbability + TotalTransitionProbability;
            
	    // Select subprocess
	    if (G4UniformRand() > TotalEmissionProbability/TotalProbability) 
	      {
		// It will be transition to state with a new number of excitons
		ThereIsTransition = true;		
		// Perform the transition
		aFragment = aTransition->PerformTransition(aFragment);
	      } 
	    else 
	      {
		// It will be fragment emission
		ThereIsTransition = false;
		Result->push_back(aEmission.PerformEmission(aFragment));
	      }
	  } 
	else 
	  {
	    // Perform Equilibrium Emission
#ifdef debug
	    CheckConservation(theInitialState,aFragment,Result);
#endif
	    PerformEquilibriumEmission(aFragment,Result);
	    delete aTransition;
	    return Result;
	  }
      } while (ThereIsTransition);   // end of do loop
  } // end of for (;;) loop
}




void G4PreCompoundModel::PerformEquilibriumEmission(const G4Fragment & aFragment,
						    G4ReactionProductVector * Result) const 
{
  G4ReactionProductVector * theEquilibriumResult;

  theEquilibriumResult = GetExcitationHandler()->BreakItUp(aFragment);
  
  Result->insert(Result->end(),theEquilibriumResult->begin(), theEquilibriumResult->end());

  delete theEquilibriumResult;
  return;
}


#ifdef debug
void G4PreCompoundModel::CheckConservation(const G4Fragment & theInitialState,
					   const G4Fragment & aFragment,
					   G4ReactionProductVector * Result) const
{
  G4double ProductsEnergy = aFragment.GetMomentum().e();
  G4ThreeVector ProductsMomentum = aFragment.GetMomentum();
  G4int ProductsA = static_cast<G4int>(aFragment.GetA());
  G4int ProductsZ = static_cast<G4int>(aFragment.GetZ());
  for (G4ReactionProductVector::iterator h = Result->begin(); 
       h != Result->end(); ++h) 
    {
      ProductsEnergy += (*h)->GetTotalEnergy();
      ProductsMomentum += (*h)->GetMomentum();
      ProductsA += static_cast<G4int>((*h)->GetDefinition()->GetBaryonNumber());
      ProductsZ += static_cast<G4int>((*h)->GetDefinition()->GetPDGCharge());
    }

  if (ProductsA != theInitialState.GetA()) 
    {
      G4cout << "!!!!!!!!!! Baryonic Number Conservation Violation !!!!!!!!!!\n"
	     << "G4PreCompoundModel.cc: Barionic Number Conservation test for just preequilibrium fragments\n" 
	     << "Initial A = " << theInitialState.GetA() 
	     << "   Fragments A = " << ProductsA << "   Diference --> " 
	     << theInitialState.GetA() - ProductsA << '\n';
    }
  if (ProductsZ != theInitialState.GetZ()) 
    {
      G4cout << "!!!!!!!!!! Charge Conservation Violation !!!!!!!!!!\n"
	     << "G4PreCompoundModel.cc: Charge Conservation test for just preequilibrium fragments\n" 
	     << "Initial Z = " << theInitialState.GetZ() 
	     << "   Fragments Z = " << ProductsZ << "   Diference --> " 
	     << theInitialState.GetZ() - ProductsZ << '\n';
    }
  if (std::abs(ProductsEnergy-theInitialState.GetMomentum().e()) > 1.0*keV) 
    {
      G4cout << "!!!!!!!!!! Energy Conservation Violation !!!!!!!!!!\n" 
	     << "G4PreCompoundModel.cc: Energy Conservation test for just preequilibrium fragments\n"  
	     << "Initial E = " << theInitialState.GetMomentum().e()/MeV << " MeV"
	     << "   Fragments E = " << ProductsEnergy/MeV  << " MeV   Diference --> " 
	     << (theInitialState.GetMomentum().e() - ProductsEnergy)/MeV << " MeV\n";
    } 
  if (std::abs(ProductsMomentum.x()-theInitialState.GetMomentum().x()) > 1.0*keV || 
      std::abs(ProductsMomentum.y()-theInitialState.GetMomentum().y()) > 1.0*keV ||
      std::abs(ProductsMomentum.z()-theInitialState.GetMomentum().z()) > 1.0*keV) 
    {
      G4cout << "!!!!!!!!!! Momentum Conservation Violation !!!!!!!!!!\n"
	     << "G4PreCompoundModel.cc: Momentum Conservation test for just preequilibrium fragments\n" 
	     << "Initial P = " << theInitialState.GetMomentum().vect() << " MeV"
	     << "   Fragments P = " << ProductsMomentum  << " MeV   Diference --> " 
	     << theInitialState.GetMomentum().vect() - ProductsMomentum << " MeV\n";
    }
  return;
}

#endif



