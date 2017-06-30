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
// $Id: G4PreCompoundModel.cc 104867 2017-06-23 14:00:18Z gcosmo $
//
// by V. Lara
//
// Modified:
// 01.04.2008 J.M.Quesada Several changes. Soft cut-off switched off. 
// 01.05.2008 J.M.Quesada Protection against non-physical preeq. 
//                        transitional regime has been set
// 03.09.2008 J.M.Quesada for external choice of inverse cross section option
// 06.09.2008 J.M.Quesada Also external choices have been added for:
//                      - superimposed Coulomb barrier (useSICB=true) 
//                      - "never go back"  hipothesis (useNGB=true) 
//                      - soft cutoff from preeq. to equlibrium (useSCO=true)
//                      - CEM transition probabilities (useCEMtr=true)  
// 20.08.2010 V.Ivanchenko Cleanup of the code: 
//                      - integer Z and A;
//                      - emission and transition classes created at 
//                        initialisation
//                      - options are set at initialisation
//                      - do not use copy-constructors for G4Fragment  
// 03.01.2012 V.Ivanchenko Added pointer to G4ExcitationHandler to the 
//                         constructor

#include "G4PreCompoundModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PreCompoundEmission.hh"
#include "G4PreCompoundTransitions.hh"
#include "G4GNASHTransitions.hh"
#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

#include "G4NucleiProperties.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "Randomize.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4LorentzVector.hh"
#include "G4Exp.hh"

////////////////////////////////////////////////////////////////////////////////

G4PreCompoundModel::G4PreCompoundModel(G4ExcitationHandler* ptr) 
  : G4VPreCompoundModel(ptr,"PRECO"),theEmission(nullptr),theTransition(nullptr),
    useSCO(false),isInitialised(false),minZ(3),minA(5) 
{
  //G4cout << "### NEW PrecompoundModel " << this << G4endl;
  if(!ptr) { SetExcitationHandler(new G4ExcitationHandler()); }

  proton = G4Proton::Proton();
  neutron = G4Neutron::Neutron();
  fLevelDensity = fLimitEnergy = 0.0;
}

////////////////////////////////////////////////////////////////////////////////

G4PreCompoundModel::~G4PreCompoundModel() 
{
  delete theEmission;
  delete theTransition;
  delete GetExcitationHandler();
}

////////////////////////////////////////////////////////////////////////////////

void G4PreCompoundModel::BuildPhysicsTable(const G4ParticleDefinition&) 
{
  InitialiseModel();
}

////////////////////////////////////////////////////////////////////////////////

void G4PreCompoundModel::InitialiseModel() 
{
  if(isInitialised) { return; }
  isInitialised = true;

  //G4cout << "G4PreCompoundModel::InitialiseModel() started" << G4endl;

  G4DeexPrecoParameters* param = 
    G4NuclearLevelData::GetInstance()->GetParameters();
  
  // 12/pi2 factor is used in real computation
  fLevelDensity = param->GetLevelDensity()*12.0/CLHEP::pi2;
  fLimitEnergy  = param->GetPrecoLowEnergy();

  useSCO = param->UseSoftCutoff();

  minZ = param->GetMinZForPreco();
  minA = param->GetMinAForPreco();

  theEmission = new G4PreCompoundEmission();
  if(param->UseHETC()) { theEmission->SetHETCModel(); }
  else { theEmission->SetDefaultModel(); }
  theEmission->SetOPTxs(param->GetPrecoModelType());

  if(param->UseGNASH()) { theTransition = new G4GNASHTransitions; }
  else { theTransition = new G4PreCompoundTransitions(); }
  theTransition->UseNGB(param->NeverGoBack());
  theTransition->UseCEMtr(param->UseCEM());

  GetExcitationHandler()->Initialise();
}

////////////////////////////////////////////////////////////////////////////////

G4HadFinalState* 
G4PreCompoundModel::ApplyYourself(const G4HadProjectile & thePrimary,
				  G4Nucleus & theNucleus)
{  
  const G4ParticleDefinition* primary = thePrimary.GetDefinition();
  if(primary != neutron && primary != proton) {
    G4ExceptionDescription ed;
    ed << "G4PreCompoundModel is used for ";
    if(primary) { ed << primary->GetParticleName(); }
    G4Exception("G4PreCompoundModel::ApplyYourself()","had0033",FatalException,
                ed,"");
    return 0;
  }

  G4int Zp = 0;
  G4int Ap = 1;
  if(primary == proton) { Zp = 1; }

  G4double timePrimary=thePrimary.GetGlobalTime();

  G4int A = theNucleus.GetA_asInt();
  G4int Z = theNucleus.GetZ_asInt();
   
  //G4cout << "### G4PreCompoundModel::ApplyYourself: A= " << A << " Z= " << Z
  //	 << " Ap= " << Ap << " Zp= " << Zp << G4endl; 
  // 4-Momentum
  G4LorentzVector p = thePrimary.Get4Momentum();
  G4double mass = G4NucleiProperties::GetNuclearMass(A, Z);
  p += G4LorentzVector(0.0,0.0,0.0,mass);
  //G4cout << "Primary 4-mom " << p << "  mass= " << mass << G4endl;

  // prepare fragment
  G4Fragment anInitialState(A + Ap, Z + Zp, p);
  anInitialState.SetNumberOfExcitedParticle(2, 1);
  anInitialState.SetNumberOfHoles(1,0);
  anInitialState.SetCreationTime(thePrimary.GetGlobalTime());
  
  // call excitation handler
  G4ReactionProductVector * result = DeExcite(anInitialState);

  // fill particle change
  theResult.Clear();
  theResult.SetStatusChange(stopAndKill);
  for(G4ReactionProductVector::iterator i= result->begin(); 
      i != result->end(); ++i)
    {
      G4DynamicParticle * aNewDP =
	       new G4DynamicParticle((*i)->GetDefinition(),
			      (*i)->GetTotalEnergy(),
			      (*i)->GetMomentum());
      G4HadSecondary aNew =  G4HadSecondary(aNewDP);
      G4double time=(*i)->GetFormationTime();
      if(time < 0.0) { time = 0.0; }
      aNew.SetTime(timePrimary + time);
      aNew.SetCreatorModelType((*i)->GetCreatorModel());
      delete (*i);
      theResult.AddSecondary(aNew);
    }
  delete result;
  
  //return the filled particle change
  return &theResult;
}

////////////////////////////////////////////////////////////////////////////////

G4ReactionProductVector* G4PreCompoundModel::DeExcite(G4Fragment& aFragment)
{
  if(!isInitialised) { InitialiseModel(); }

  G4ReactionProductVector * Result = new G4ReactionProductVector;
  G4double Eex = aFragment.GetExcitationEnergy();
  G4int Z = aFragment.GetZ_asInt(); 
  G4int A = aFragment.GetA_asInt(); 

  //G4cout << "### G4PreCompoundModel::DeExcite" << G4endl;
  //G4cout << aFragment << G4endl;
 
  // Perform Equilibrium Emission 
  if ((Z < minZ && A < minA) || Eex < fLimitEnergy*A) {
    PerformEquilibriumEmission(aFragment, Result);
    return Result;
  }
  
  // main loop  
  G4int count = 0;
  static const G4int countmax = 1000;
  for (;;) {
    //G4cout << "### PreCompound loop over fragment" << G4endl;
    //G4cout << aFragment << G4endl;
    G4int EquilibriumExcitonNumber = 
      G4lrint(std::sqrt(aFragment.GetExcitationEnergy()
			*aFragment.GetA_asInt()*fLevelDensity));
    //   
    //    G4cout<<"Neq="<<EquilibriumExcitonNumber<<G4endl;
    //
    // J. M. Quesada (Jan. 08)  equilibrium hole number could be used as preeq.
    // evap. delimiter (IAEA report)
    
    // Loop for transitions, it is performed while there are 
    // preequilibrium transitions.
    G4bool ThereIsTransition = false;
    
    //        G4cout<<"----------------------------------------"<<G4endl;
    //        G4double NP=aFragment.GetNumberOfParticles();
    //        G4double NH=aFragment.GetNumberOfHoles();
    //        G4double NE=aFragment.GetNumberOfExcitons();
    //        G4cout<<" Ex. Energy="<<aFragment.GetExcitationEnergy()<<G4endl;
    //   G4cout<<"N. excitons="<<NE<<"  N. Part="<<NP<<"N. Holes ="<<NH<<G4endl;
    do {
      ++count;
      //G4cout<<"transition number .."<<count
      //      <<" n ="<<aFragment.GetNumberOfExcitons()<<G4endl;
      G4bool go_ahead = false;
      // soft cutoff criterium as an "ad-hoc" solution to force 
      // increase in  evaporation  
      G4int test = aFragment.GetNumberOfExcitons();
      if (test <= EquilibriumExcitonNumber) { go_ahead=true; }

      //J. M. Quesada (Apr. 08): soft-cutoff switched off by default
      if (useSCO && go_ahead)
	{
	  G4double x = G4double(test)/G4double(EquilibriumExcitonNumber) - 1;
	  if( G4UniformRand() < 1.0 -  G4Exp(-x*x/0.32) ) { go_ahead = false; }
	} 
        
      // JMQ: WARNING:  CalculateProbability MUST be called prior to Get!! 
      // (O values would be returned otherwise)
      G4double TotalTransitionProbability = 
	theTransition->CalculateProbability(aFragment);
      G4double P1 = theTransition->GetTransitionProb1();
      G4double P2 = theTransition->GetTransitionProb2();
      G4double P3 = theTransition->GetTransitionProb3();
      //G4cout<<"#0 P1="<<P1<<" P2="<<P2<<"  P3="<<P3<<G4endl;
      
      //J.M. Quesada (May 2008) Physical criterium (lamdas)  PREVAILS over 
      //                        approximation (critical exciton number)
      //V.Ivanchenko (May 2011) added check on number of nucleons
      //                        to send a fragment to FermiBreakUp
      if(!go_ahead || P1 <= P2+P3 || 
	 aFragment.GetZ_asInt() < minZ || aFragment.GetA_asInt() < minA ||
	 (aFragment.GetExcitationEnergy() <= fLimitEnergy*aFragment.GetA_asInt()))
	{
	  //G4cout<<"#4 EquilibriumEmission"<<G4endl; 
	  PerformEquilibriumEmission(aFragment,Result);
	  return Result;
	}
      else 
	{
	  //
	  // Check if number of excitons is greater than 0
	  // else perform equilibrium emission
	  if (aFragment.GetNumberOfExcitons() <= 0) 
	    {
	      PerformEquilibriumEmission(aFragment,Result);
	      return Result;
	    }
	    
	  G4double TotalEmissionProbability = 
	    theEmission->GetTotalProbability(aFragment);
	  //
	  //G4cout<<"#1 TotalEmissionProbability="<<TotalEmissionProbability
	  // <<" Nex= " <<aFragment.GetNumberOfExcitons()<<G4endl;
	  //J.M.Quesada (May 08) this has already been done in order to decide  
	  //                     what to do (preeq-eq) 
	  // Sum of all probabilities
	  G4double TotalProbability = TotalEmissionProbability 
	    + TotalTransitionProbability;
            
	  // Select subprocess
	  if (TotalProbability*G4UniformRand() > TotalEmissionProbability) 
	    {
	      //G4cout<<"#2 Transition"<<G4endl; 
	      // It will be transition to state with a new number of excitons
	      ThereIsTransition = true;		
	      // Perform the transition
	      theTransition->PerformTransition(aFragment);
	    } 
	  else 
	    {
	      //G4cout<<"#3 Emission"<<G4endl; 
	      // It will be fragment emission
	      ThereIsTransition = false;
	      Result->push_back(theEmission->PerformEmission(aFragment));
	    }
	}
      // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
    } while (ThereIsTransition);   // end of do loop

    // stop if too many iterations
    if(count >= countmax) {
      G4ExceptionDescription ed;
      ed << "G4PreCompoundModel loop over " << countmax << " iterations; "
	 << "current G4Fragment: \n" << aFragment;
      G4Exception("G4PreCompoundModel::DeExcite()","had0034",JustWarning,
		  ed,"");
      PerformEquilibriumEmission(aFragment, Result);
      return Result;
    }
  } // end of for (;;) loop
  return Result;
}

////////////////////////////////////////////////////////////////////////////////
//       Initialisation
////////////////////////////////////////////////////////////////////////////////

void G4PreCompoundModel::UseHETCEmission() 
{
  PrintWarning("UseHETCEmission");
}

void G4PreCompoundModel::UseDefaultEmission() 
{ 
  PrintWarning("UseDefaultEmission");
}

void G4PreCompoundModel::UseGNASHTransition() 
{ 
  PrintWarning("UseGNASHTransition");
}

void G4PreCompoundModel::UseDefaultTransition() 
{ 
  PrintWarning("UseDefaultTransition");
}

void G4PreCompoundModel::SetOPTxs(G4int) 
{ 
  PrintWarning("UseOPTxs");
}

void G4PreCompoundModel::UseSICB() 
{ 
  PrintWarning("UseSICB");
}

void G4PreCompoundModel::UseNGB()  
{ 
  PrintWarning("UseNGB");
}

void G4PreCompoundModel::UseSCO()  
{ 
  PrintWarning("UseSCO");
}

void G4PreCompoundModel::UseCEMtr() 
{ 
  PrintWarning("UseCEMtr");
}

void G4PreCompoundModel::PrintWarning(const G4String& mname)
{
  G4ExceptionDescription ed;
  ed << "Obsolete method of the preCompound model is called: " 
     << mname << "() \n Instead a corresponding method of "
     << "G4DeexPrecoParameters class should be used";

  G4Exception("G4PreCompoundModel::ReadData()","had0803",JustWarning,ed);
}

////////////////////////////////////////////////////////////////////////////////
//       Documentation
////////////////////////////////////////////////////////////////////////////////

void G4PreCompoundModel::ModelDescription(std::ostream& outFile) const
{
  outFile 
    << "The GEANT4 precompound model is considered as an extension of the\n"
    <<	"hadron kinetic model. It gives a possibility to extend the low energy range\n"
    <<	"of the hadron kinetic model for nucleon-nucleus inelastic collision and it \n"
    <<	"provides a ”smooth” transition from kinetic stage of reaction described by the\n"
    <<	"hadron kinetic model to the equilibrium stage of reaction described by the\n"
    <<	"equilibrium deexcitation models.\n"
    <<	"The initial information for calculation of pre-compound nuclear stage\n"
    <<	"consists of the atomic mass number A, charge Z of residual nucleus, its\n"
    <<	"four momentum P0 , excitation energy U and number of excitons n, which equals\n"
    <<	"the sum of the number of particles p (from them p_Z are charged) and the number of\n"
    <<	"holes h.\n"
    <<	"At the preequilibrium stage of reaction, we follow the exciton model approach in ref. [1],\n"
    <<	"taking into account the competition among all possible nuclear transitions\n"
    <<	"with ∆n = +2, −2, 0 (which are deﬁned by their associated transition probabilities) and\n"
    <<	"the emission of neutrons, protons, deutrons, thritium and helium nuclei (also defined by\n"
    <<	"their associated emission  probabilities according to exciton model)\n"
    <<	"\n"
    <<	"[1] K.K. Gudima, S.G. Mashnik, V.D. Toneev, Nucl. Phys. A401 329 (1983)\n"
    << "\n";
}

void G4PreCompoundModel::DeExciteModelDescription(std::ostream& outFile) const
{
  outFile << "description of precompound model as used with DeExcite()" << "\n";
}

/////////////////////////////////////////////////////////////////////////////////////////
