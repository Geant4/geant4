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
// $Id: G4PreCompoundTransitions.cc 96603 2016-04-25 13:29:51Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PreCompoundTransitions
//
// Author:         V.Lara
//
// Modified:  
// 16.02.2008 J.M.Quesada fixed bugs 
// 06.09.2008 J.M.Quesada added external choices for:
//                      - "never go back"  hipothesis (useNGB=true) 
//                      -  CEM transition probabilities (useCEMtr=true)
// 30.10.2009 J.M.Quesada: CEM transition probabilities have been renormalized 
//                       (IAEA benchmark)
// 20.08.2010 V.Ivanchenko move constructor and destructor to the source and 
//                         optimise the code
// 30.08.2011 M.Kelsey - Skip CalculateProbability if no excitons

#include "G4PreCompoundTransitions.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4Fragment.hh"
#include "G4Proton.hh"
#include "G4Exp.hh"
#include "G4Log.hh"

G4PreCompoundTransitions::G4PreCompoundTransitions() 
{
  proton = G4Proton::Proton();
  G4DeexPrecoParameters* param = 
    G4NuclearLevelData::GetInstance()->GetParameters() ;
  FermiEnergy = param->GetFermiEnergy();
  r0 = param->GetTransitionsR0();
  aLDP = param->GetLevelDensity();
}

G4PreCompoundTransitions::~G4PreCompoundTransitions() 
{}

// Calculates transition probabilities with 
// DeltaN = +2 (Trans1) -2 (Trans2) and 0 (Trans3)
G4double G4PreCompoundTransitions::
CalculateProbability(const G4Fragment & aFragment)
{
  // Number of holes
  G4int H = aFragment.GetNumberOfHoles();
  // Number of Particles 
  G4int P = aFragment.GetNumberOfParticles();
  // Number of Excitons 
  G4int N = P+H;
  // Nucleus 
  G4int A = aFragment.GetA_asInt();
  G4int Z = aFragment.GetZ_asInt();
  G4double U = aFragment.GetExcitationEnergy();
  TransitionProb2 = 0.0;
  TransitionProb3 = 0.0;
  /*
  G4cout << "G4PreCompoundTransitions::CalculateProbability H/P/N/Z/A= " 
	 << H << " " << P << " " << N << " " << Z << " " << A <<G4endl;
  G4cout << aFragment << G4endl;
  */
  if(U < 10*eV || 0==N) { return 0.0; }
  
  //J. M. Quesada (Feb. 08) new physics
  // OPT=1 Transitions are calculated according to Gudima's paper 
  //       (original in G4PreCompound from VL) 
  // OPT=2 Transitions are calculated according to Gupta's formulae
  //
  static const G4double sixdpi2 = 6.0/CLHEP::pi2;
  if (useCEMtr) {
    // Relative Energy (T_{rel})
    G4double RelativeEnergy = 1.6*FermiEnergy + U/G4double(N);
    
    // Sample kind of nucleon-projectile 
    G4bool ChargedNucleon(false);
    if(G4int(P*G4UniformRand()) <= aFragment.GetNumberOfCharged()) {
      ChargedNucleon = true; 
    }
    
    // Relative Velocity: 
    // <V_{rel}>^2
    G4double RelativeVelocitySqr;
    if (ChargedNucleon) { 
      RelativeVelocitySqr = 2*RelativeEnergy/CLHEP::proton_mass_c2; 
    } else { 
      RelativeVelocitySqr = 2*RelativeEnergy/CLHEP::neutron_mass_c2; 
    }
    // <V_{rel}>
    G4double RelativeVelocity = std::sqrt(RelativeVelocitySqr);
    
    // Proton-Proton Cross Section
    G4double ppXSection = 
      (10.63/RelativeVelocitySqr - 29.92/RelativeVelocity + 42.9)
      * CLHEP::millibarn;
    // Proton-Neutron Cross Section
    G4double npXSection = 
      (34.10/RelativeVelocitySqr - 82.20/RelativeVelocity + 82.2)
      * CLHEP::millibarn;
    
    // Averaged Cross Section: \sigma(V_{rel})
    G4double AveragedXSection;
    if (ChargedNucleon)
      {
        //JMQ: small bug fixed
        AveragedXSection = ((Z-1)*ppXSection + (A-Z)*npXSection)/G4double(A-1);
      }
    else 
      {
        AveragedXSection = ((A-Z-1)*ppXSection + Z*npXSection)/G4double(A-1);
      }
    
    // Fermi relative energy ratio
    G4double FermiRelRatio = FermiEnergy/RelativeEnergy;
    
    // This factor is introduced to take into account the Pauli principle
    G4double PauliFactor = 1.0 - 1.4*FermiRelRatio;
    if (FermiRelRatio > 0.5) {
      G4double x = 2.0 - 1.0/FermiRelRatio;
      PauliFactor += 0.4*FermiRelRatio*x*x*std::sqrt(x);
    }
    // Interaction volume 
    G4double xx = 2*r0 + CLHEP::hbarc/(CLHEP::proton_mass_c2*RelativeVelocity);
    G4double Vint = CLHEP::pi*xx*xx*xx/0.75;
    
    // Transition probability for \Delta n = +2
    
    TransitionProb1 = std::max(0.0, AveragedXSection*PauliFactor
      *std::sqrt(2.0*RelativeEnergy/CLHEP::proton_mass_c2)/Vint);

    //JMQ 281009  phenomenological factor in order to increase 
    //   equilibrium contribution
    //   G4double factor=5.0;
    //   TransitionProb1 *= factor;
    
    // GE = g*E where E is Excitation Energy
    G4double GE = sixdpi2*aLDP*A*U;
    G4double Fph = G4double(P*P+H*H+P-3*H)*0.25;
    
    if(!useNGB) { 
        
      // F(p+1,h+1)
      G4double Fph1 = Fph + N*0.5;

      static const G4double plimit = 100;

      //JMQ/AH  bug fixed: if (U-Fph < 0.0) 
      if (GE-Fph1 > 0.0) { 
        G4double x0 = GE-Fph;
	G4double x1 = (N+1)*G4Log(x0/(GE-Fph1));
	if(x1 < plimit) {
	  x1 = G4Exp(x1)*TransitionProb1/x0;
    
	  // Transition probability for \Delta n = -2 (at F(p,h) = 0)
	  TransitionProb2 = std::max(0.0, (P*H*(N+1)*(N-2))*x1/x0);
        
	  // Transition probability for \Delta n = 0 (at F(p,h) = 0)
	  TransitionProb3 = std::max(0.0,((N+1)*(P*(P-1) + 4*P*H + H*(H-1)))*x1
				     /G4double(N));
	}
      }
    }

  } else {
    //JMQ: Transition probabilities from Gupta's work    
    // GE = g*E where E is Excitation Energy
    TransitionProb1 = std::max(0.0, U*(4.2e+12 - 3.6e+10*U/G4double(N+1)))
      /(16*CLHEP::c_light); 

    if (!useNGB && N > 1) {
      G4double GE = sixdpi2*aLDP*A*U; 
      TransitionProb2 = ((N-1)*(N-2)*P*H)*TransitionProb1/(GE*GE);  
    }
  }
  //  G4cout<<"U = "<<U<<G4endl;
  //  G4cout<<"N="<<N<<"  P="<<P<<"  H="<<H<<G4endl;
  //  G4cout<<"l+ ="<<TransitionProb1<<"  l- ="<< TransitionProb2
  //   <<"  l0 ="<< TransitionProb3<<G4endl; 
  return TransitionProb1 + TransitionProb2 + TransitionProb3;
}

void G4PreCompoundTransitions::PerformTransition(G4Fragment & result)
{
  G4double ChosenTransition = 
    G4UniformRand()*(TransitionProb1 + TransitionProb2 + TransitionProb3);
  G4int deltaN = 0;
  G4int Npart     = result.GetNumberOfParticles();
  G4int Ncharged  = result.GetNumberOfCharged();
  G4int Nholes    = result.GetNumberOfHoles();
  if (ChosenTransition <= TransitionProb1) 
    {
      // Number of excitons is increased on \Delta n = +2
      deltaN = 2;
    } 
  else if (ChosenTransition <= TransitionProb1+TransitionProb2) 
    {
      // Number of excitons is increased on \Delta n = -2
      deltaN = -2;
    }

  // AH/JMQ: Randomly decrease the number of charges if deltaN is -2 and  
  // in proportion to the number charges w.r.t. number of particles,  
  // PROVIDED that there are charged particles
  deltaN /= 2;

  //G4cout << "deltaN= " << deltaN << G4endl;

  // JMQ the following lines have to be before SetNumberOfCharged, 
  //     otherwise the check on number of charged vs. number of particles fails
  result.SetNumberOfParticles(Npart+deltaN);
  result.SetNumberOfHoles(Nholes+deltaN); 

  if(deltaN < 0) {
    if( (Ncharged == Npart) ||
	(Ncharged >= 1 && G4int(Npart*G4UniformRand()) <= Ncharged)) 
      { 
	result.SetNumberOfCharged(Ncharged+deltaN); // deltaN is negative!
      }

  } else if ( deltaN > 0 ) {
    // With weight Z/A, number of charged particles is increased with +1
    G4int A = result.GetA_asInt() - Npart;
    G4int Z = result.GetZ_asInt() - Ncharged;
    if((Z == A) ||  (Z > 0 && G4int(A*G4UniformRand()) <= Z)) 
      {
	result.SetNumberOfCharged(Ncharged+deltaN);
      }
  }
  
  // Number of charged can not be greater that number of particles
  if ( Npart < Ncharged ) 
    {
      result.SetNumberOfCharged(Npart);
    }
  //G4cout << "### After transition" << G4endl;
  //G4cout << result << G4endl;
}

