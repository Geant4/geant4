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
// $Id: G4PreCompoundTransitions.cc,v 1.22 2009/11/21 18:03:13 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PreCompoundIon
//
// Author:         V.Lara
//
// Modified:  
// 16.02.2008 J. M. Quesada fixed bugs 
// 06.09.2008 J. M. Quesada added external choices for:
//                      - "never go back"  hipothesis (useNGB=true) 
//                      -  CEM transition probabilities (useCEMtr=true)

// 30.10.09 J.M.Quesada: CEM transition probabilities have been renormalized 
//                       (IAEA benchmark)
//
#include "G4PreCompoundTransitions.hh"
#include "G4HadronicException.hh"

const G4PreCompoundTransitions & G4PreCompoundTransitions::
operator=(const G4PreCompoundTransitions &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4PreCompoundTransitions::operator= meant to not be accessable");
  return *this;
}


G4bool G4PreCompoundTransitions::operator==(const G4PreCompoundTransitions &) const
{
  return false;
}

G4bool G4PreCompoundTransitions::operator!=(const G4PreCompoundTransitions &) const
{
  return true;
}


G4double G4PreCompoundTransitions::
CalculateProbability(const G4Fragment & aFragment)
{
  //G4cout<<"In G4PreCompoundTransitions.cc  useNGB="<<useNGB<<G4endl;
  //G4cout<<"In G4PreCompoundTransitions.cc  useCEMtr="<<useCEMtr<<G4endl;

  // Fermi energy
  const G4double FermiEnergy = G4PreCompoundParameters::GetAddress()->GetFermiEnergy();
  
  // Nuclear radius
  const G4double r0 = G4PreCompoundParameters::GetAddress()->GetTransitionsr0();
   
  // In order to calculate the level density parameter
  // G4EvaporationLevelDensityParameter theLDP;

  // Number of holes
  G4double H = aFragment.GetNumberOfHoles();
  // Number of Particles 
  G4double P = aFragment.GetNumberOfParticles();
  // Number of Excitons 
  G4double N = P+H;
  // Nucleus 
  G4double A = aFragment.GetA();
  G4double Z = static_cast<G4double>(aFragment.GetZ());
  G4double U = aFragment.GetExcitationEnergy();
  
  if(U<10*eV) return 0.0;
  
  //J. M. Quesada (Feb. 08) new physics
  // OPT=1 Transitions are calculated according to Gudima's paper (original in G4PreCompound from VL) 
  // OPT=2 Transitions are calculated according to Gupta's formulae
  //
  
  
  
  if (useCEMtr){

    
    // Relative Energy (T_{rel})
    G4double RelativeEnergy = (8.0/5.0)*FermiEnergy + U/N;
    
    // Sample kind of nucleon-projectile 
    G4bool ChargedNucleon(false);
    G4double chtest = 0.5;
    if (P > 0) chtest = aFragment.GetNumberOfCharged()/P;
    if (G4UniformRand() < chtest) ChargedNucleon = true;
    
    // Relative Velocity: 
    // <V_{rel}>^2
    G4double RelativeVelocitySqr(0.0);
    if (ChargedNucleon) RelativeVelocitySqr = 2.0*RelativeEnergy/proton_mass_c2;
    else RelativeVelocitySqr = 2.0*RelativeEnergy/neutron_mass_c2;
    
    // <V_{rel}>
    G4double RelativeVelocity = std::sqrt(RelativeVelocitySqr);
    
    // Proton-Proton Cross Section
    G4double ppXSection = (10.63/RelativeVelocitySqr - 29.92/RelativeVelocity + 42.9)*millibarn;
    // Proton-Neutron Cross Section
    G4double npXSection = (34.10/RelativeVelocitySqr - 82.20/RelativeVelocity + 82.2)*millibarn;
    
    // Averaged Cross Section: \sigma(V_{rel})
    //  G4double AveragedXSection = (ppXSection+npXSection)/2.0;
    G4double AveragedXSection(0.0);
    if (ChargedNucleon)
      {
        //JMQ: small bug fixed
        //      AveragedXSection = ((Z-1.0) * ppXSection + (A-Z-1.0) * npXSection) / (A-1.0);
        AveragedXSection = ((Z-1.0) * ppXSection + (A-Z) * npXSection) / (A-1.0);
      }
    else 
      {
        AveragedXSection = ((A-Z-1.0) * ppXSection + Z * npXSection) / (A-1.0);
      }
    
    // Fermi relative energy ratio
    G4double FermiRelRatio = FermiEnergy/RelativeEnergy;
    
    // This factor is introduced to take into account the Pauli principle
    G4double PauliFactor = 1.0 - (7.0/5.0)*FermiRelRatio;
    if (FermiRelRatio > 0.5) PauliFactor += (2.0/5.0)*FermiRelRatio*std::pow(2.0 - (1.0/FermiRelRatio), 5.0/2.0);
    
    // Interaction volume 
    //  G4double Vint = (4.0/3.0)*pi*std::pow(2.0*r0 + hbarc/(proton_mass_c2*RelativeVelocity) , 3.0);
    G4double xx=2.0*r0 + hbarc/(proton_mass_c2*RelativeVelocity);
    G4double Vint = (4.0/3.0)*pi*xx*xx*xx;
    
    // Transition probability for \Delta n = +2
    
    TransitionProb1 = AveragedXSection*PauliFactor*std::sqrt(2.0*RelativeEnergy/proton_mass_c2)/Vint;

//JMQ 281009  phenomenological factor in order to increase equilibrium contribution
//   G4double factor=5.0;
//   TransitionProb1 *= factor;
//
    if (TransitionProb1 < 0.0) TransitionProb1 = 0.0; 
    
    G4double a = G4PreCompoundParameters::GetAddress()->GetLevelDensity();
    // GE = g*E where E is Excitation Energy
    G4double GE = (6.0/pi2)*a*A*U;
    
    G4double Fph = ((P*P+H*H+P-H)/4.0 - H/2.0);
    
    //G4bool NeverGoBack(false);
    G4bool NeverGoBack;
    if(useNGB)  NeverGoBack=true;
    else NeverGoBack=false;
    
    
    //JMQ/AH  bug fixed: if (U-Fph < 0.0) NeverGoBack = true;
    if (GE-Fph < 0.0) NeverGoBack = true;
    
    // F(p+1,h+1)
    G4double Fph1 = Fph + N/2.0;
    
    G4double ProbFactor = std::pow((GE-Fph)/(GE-Fph1),N+1.0);
    
    
    if (NeverGoBack)
      {
      TransitionProb2 = 0.0;
      TransitionProb3 = 0.0;
      }
    else 
      {
        // Transition probability for \Delta n = -2 (at F(p,h) = 0)
        TransitionProb2 = TransitionProb1 * ProbFactor * (P*H*(N+1.0)*(N-2.0))/((GE-Fph)*(GE-Fph));
        if (TransitionProb2 < 0.0) TransitionProb2 = 0.0; 
        
        // Transition probability for \Delta n = 0 (at F(p,h) = 0)
        TransitionProb3 = TransitionProb1* ((N+1.0)/N) * ProbFactor  * (P*(P-1.0) + 4.0*P*H + H*(H-1.0))/(GE-Fph);
        if (TransitionProb3 < 0.0) TransitionProb3 = 0.0; 
      }
    
    //  G4cout<<"U = "<<U<<G4endl;
    //  G4cout<<"N="<<N<<"  P="<<P<<"  H="<<H<<G4endl;
    //  G4cout<<"l+ ="<<TransitionProb1<<"  l- ="<< TransitionProb2<<"  l0 ="<< TransitionProb3<<G4endl; 
    return TransitionProb1 + TransitionProb2 + TransitionProb3;
  }
  
  else  {
    //JMQ: Transition probabilities from Gupta's work
    
    G4double a = G4PreCompoundParameters::GetAddress()->GetLevelDensity();
    // GE = g*E where E is Excitation Energy
    G4double GE = (6.0/pi2)*a*A*U;
 
    G4double Kmfp=2.;
        
    TransitionProb1=1./Kmfp*3./8.*1./c_light*1.0e-9*(1.4e+21*U-2./(N+1)*6.0e+18*U*U);
    if (TransitionProb1 < 0.0) TransitionProb1 = 0.0;
    
    if (useNGB){
      TransitionProb2=0.;
      TransitionProb3=0.;
    }
    else{        
      if (N<=1) TransitionProb2=0. ;    
      else  TransitionProb2=1./Kmfp*3./8.*1./c_light*1.0e-9*(N-1.)*(N-2.)*P*H/(GE*GE)*(1.4e+21*U - 2./(N-1)*6.0e+18*U*U);      
      if (TransitionProb2 < 0.0) TransitionProb2 = 0.0; 
      TransitionProb3=0.;
    }
    
      //  G4cout<<"U = "<<U<<G4endl;
    //  G4cout<<"N="<<N<<"  P="<<P<<"  H="<<H<<G4endl;
    //  G4cout<<"l+ ="<<TransitionProb1<<"  l- ="<< TransitionProb2<<"  l0 ="<< TransitionProb3<<G4endl; 
    return TransitionProb1 + TransitionProb2 + TransitionProb3;
  }
}


G4Fragment G4PreCompoundTransitions::PerformTransition(const G4Fragment & aFragment)
{
  G4Fragment result(aFragment);
  G4double ChosenTransition = G4UniformRand()*(TransitionProb1 + TransitionProb2 + TransitionProb3);
  G4int deltaN = 0;
  G4int Nexcitons = result.GetNumberOfExcitons();
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

  // AH/JMQ: Randomly decrease the number of charges if deltaN is -2 and in proportion 
  // to the number charges w.r.t. number of particles,  PROVIDED that there are charged particles
  if(deltaN < 0 && G4UniformRand() <= 
     static_cast<G4double>(result.GetNumberOfCharged())/static_cast<G4double>(result.GetNumberOfParticles()) 
     && (result.GetNumberOfCharged() >= 1)) {
    result.SetNumberOfCharged(result.GetNumberOfCharged()+deltaN/2); // deltaN is negative!
  }

  // JMQ the following lines have to be before SetNumberOfCharged, otherwise the check on 
  // number of charged vs. number of particles fails
  result.SetNumberOfParticles(result.GetNumberOfParticles()+deltaN/2);
  result.SetNumberOfHoles(result.GetNumberOfHoles()+deltaN/2); 

  // With weight Z/A, number of charged particles is increased with +1
  if ( ( deltaN > 0 ) &&
      (G4UniformRand() <= static_cast<G4double>(result.GetZ()-result.GetNumberOfCharged())/
       std::max(static_cast<G4double>(result.GetA()-Nexcitons),1.)))
    {
      result.SetNumberOfCharged(result.GetNumberOfCharged()+deltaN/2);
    }
  
  // Number of charged can not be greater that number of particles
  if ( result.GetNumberOfParticles() < result.GetNumberOfCharged() ) 
    {
      result.SetNumberOfCharged(result.GetNumberOfParticles());
    }
  
  return result;
}

