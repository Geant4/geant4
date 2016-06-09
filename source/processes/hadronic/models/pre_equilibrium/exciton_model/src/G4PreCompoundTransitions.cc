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
#include "G4Pow.hh"
#include "G4HadronicException.hh"
#include "G4PreCompoundParameters.hh"
#include "G4Proton.hh"

G4PreCompoundTransitions::G4PreCompoundTransitions() 
{
  proton = G4Proton::Proton();
  FermiEnergy = G4PreCompoundParameters::GetAddress()->GetFermiEnergy();
  r0 = G4PreCompoundParameters::GetAddress()->GetTransitionsr0();
  aLDP = G4PreCompoundParameters::GetAddress()->GetLevelDensity();
  g4pow = G4Pow::GetInstance();
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

  //G4cout << aFragment << G4endl;
  
  if(U < 10*eV || 0==N) { return 0.0; }
  
  //J. M. Quesada (Feb. 08) new physics
  // OPT=1 Transitions are calculated according to Gudima's paper 
  //       (original in G4PreCompound from VL) 
  // OPT=2 Transitions are calculated according to Gupta's formulae
  //
  if (useCEMtr){

    // Relative Energy (T_{rel})
    G4double RelativeEnergy = 1.6*FermiEnergy + U/G4double(N);
    
    // Sample kind of nucleon-projectile 
    G4bool ChargedNucleon(false);
    G4double chtest = 0.5;
    if (P > 0) { 
      chtest = G4double(aFragment.GetNumberOfCharged())/G4double(P); 
    }
    if (G4UniformRand() < chtest) { ChargedNucleon = true; }
    
    // Relative Velocity: 
    // <V_{rel}>^2
    G4double RelativeVelocitySqr(0.0);
    if (ChargedNucleon) { 
      RelativeVelocitySqr = 2.0*RelativeEnergy/CLHEP::proton_mass_c2; 
    } else { 
      RelativeVelocitySqr = 2.0*RelativeEnergy/CLHEP::neutron_mass_c2; 
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
    //  G4double AveragedXSection = (ppXSection+npXSection)/2.0;
    G4double AveragedXSection(0.0);
    if (ChargedNucleon)
      {
        //JMQ: small bug fixed
        //AveragedXSection=((Z-1.0) * ppXSection + (A-Z-1.0) * npXSection)/(A-1.0);
        AveragedXSection = ((Z-1)*ppXSection + (A-Z)*npXSection)/G4double(A-1);
      }
    else 
      {
        AveragedXSection = ((A-Z-1)*ppXSection + Z*npXSection)/G4double(A-1);
        //AveragedXSection = ((A-Z-1)*npXSection + Z*ppXSection)/G4double(A-1);
      }
    
    // Fermi relative energy ratio
    G4double FermiRelRatio = FermiEnergy/RelativeEnergy;
    
    // This factor is introduced to take into account the Pauli principle
    G4double PauliFactor = 1.0 - 1.4*FermiRelRatio;
    if (FermiRelRatio > 0.5) {
      G4double x = 2.0 - 1.0/FermiRelRatio;
      PauliFactor += 0.4*FermiRelRatio*x*x*std::sqrt(x);
      //PauliFactor += 
      //(2.0/5.0)*FermiRelRatio*std::pow(2.0 - (1.0/FermiRelRatio), 5.0/2.0);
    }
    // Interaction volume 
    //  G4double Vint = (4.0/3.0)
    //*pi*std::pow(2.0*r0 + hbarc/(proton_mass_c2*RelativeVelocity) , 3.0);
    G4double xx = 2.0*r0 + hbarc/(CLHEP::proton_mass_c2*RelativeVelocity);
    //    G4double Vint = (4.0/3.0)*CLHEP::pi*xx*xx*xx;
    G4double Vint = CLHEP::pi*xx*xx*xx/0.75;
    
    // Transition probability for \Delta n = +2
    
    TransitionProb1 = AveragedXSection*PauliFactor
      *std::sqrt(2.0*RelativeEnergy/CLHEP::proton_mass_c2)/Vint;

    //JMQ 281009  phenomenological factor in order to increase 
    //   equilibrium contribution
    //   G4double factor=5.0;
    //   TransitionProb1 *= factor;
    //
    if (TransitionProb1 < 0.0) { TransitionProb1 = 0.0; } 
    
    // GE = g*E where E is Excitation Energy
    G4double GE = (6.0/pi2)*aLDP*A*U;
    
    //G4double Fph = ((P*P+H*H+P-H)/4.0 - H/2.0);
    G4double Fph = G4double(P*P+H*H+P-3*H)/4.0;
    
    G4bool NeverGoBack(false);
    if(useNGB) { NeverGoBack=true; }
    
    //JMQ/AH  bug fixed: if (U-Fph < 0.0) NeverGoBack = true;
    if (GE-Fph < 0.0) { NeverGoBack = true; }
    
    // F(p+1,h+1)
    G4double Fph1 = Fph + N/2.0;
    
    G4double ProbFactor = g4pow->powN((GE-Fph)/(GE-Fph1),N+1);
    
    if (NeverGoBack)
      {
	TransitionProb2 = 0.0;
	TransitionProb3 = 0.0;
      }
    else 
      {
        // Transition probability for \Delta n = -2 (at F(p,h) = 0)
        TransitionProb2 = 
	  TransitionProb1 * ProbFactor * (P*H*(N+1)*(N-2))/((GE-Fph)*(GE-Fph));
        if (TransitionProb2 < 0.0) { TransitionProb2 = 0.0; } 
        
        // Transition probability for \Delta n = 0 (at F(p,h) = 0)
        TransitionProb3 = TransitionProb1*(N+1)* ProbFactor  
	  * (P*(P-1) + 4.0*P*H + H*(H-1))/(N*(GE-Fph));
        if (TransitionProb3 < 0.0) { TransitionProb3 = 0.0; }
      }
  } else {
    //JMQ: Transition probabilities from Gupta's work    
    // GE = g*E where E is Excitation Energy
    G4double GE = (6.0/pi2)*aLDP*A*U;
 
    G4double Kmfp=2.;
        
    //TransitionProb1=1./Kmfp*3./8.*1./c_light*1.0e-9*(1.4e+21*U-2./(N+1)*6.0e+18*U*U);
    TransitionProb1 = 3.0e-9*(1.4e+21*U - 1.2e+19*U*U/G4double(N+1))
      /(8*Kmfp*CLHEP::c_light);
    if (TransitionProb1 < 0.0) { TransitionProb1 = 0.0; }

    TransitionProb2=0.;
    TransitionProb3=0.;
    
    if (!useNGB && N > 1) {
      // TransitionProb2=1./Kmfp*3./8.*1./c_light*1.0e-9*(N-1.)*(N-2.)*P*H/(GE*GE)*(1.4e+21*U - 2./(N-1)*6.0e+18*U*U);      
      TransitionProb2 = 
	3.0e-9*(N-2)*P*H*(1.4e+21*U*(N-1) - 1.2e+19*U*U)/(8*Kmfp*c_light*GE*GE);      
      if (TransitionProb2 < 0.0) TransitionProb2 = 0.0; 
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
  //  G4int Nexcitons = result.GetNumberOfExcitons();
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

  // JMQ the following lines have to be before SetNumberOfCharged, otherwise the check on 
  // number of charged vs. number of particles fails
  result.SetNumberOfParticles(Npart+deltaN);
  result.SetNumberOfHoles(Nholes+deltaN); 

  if(deltaN < 0) {
    if( Ncharged >= 1 && G4int(Npart*G4UniformRand()) <= Ncharged) 
      { 
	result.SetNumberOfCharged(Ncharged+deltaN); // deltaN is negative!
      }

  } else if ( deltaN > 0 ) {
    // With weight Z/A, number of charged particles is increased with +1
    G4int A = result.GetA_asInt();
    G4int Z = result.GetZ_asInt();
    if( G4int(std::max(1, A - Npart)*G4UniformRand()) <= Z) 
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

