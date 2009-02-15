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
//J.M. Quesada (August2008). Based on:
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option
// JMQ (06 September 2008) Also external choices have been added for 
// superimposed Coulomb barrier (if useSICB is set true, by default is false) 
//
// JMQ (14 february 2009) bug fixed in emission width: hbarc instead of hbar_Planck in the denominator
//
#include <iostream>
using namespace std;

#include "G4EvaporationProbability.hh"
#include "G4PairingCorrection.hh"



G4EvaporationProbability::G4EvaporationProbability(const G4EvaporationProbability &) : G4VEmissionProbability()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4EvaporationProbability::copy_constructor meant to not be accessable");
}




const G4EvaporationProbability & G4EvaporationProbability::
operator=(const G4EvaporationProbability &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4EvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4EvaporationProbability::operator==(const G4EvaporationProbability &) const
{
    return false;
}

G4bool G4EvaporationProbability::operator!=(const G4EvaporationProbability &) const
{
    return true;
}
  
G4double G4EvaporationProbability::EmissionProbability(const G4Fragment & fragment, const G4double anEnergy)
{
    G4double EmissionProbability = 0.0;
    G4double MaximalKineticEnergy = anEnergy;

    if (MaximalKineticEnergy > 0.0 && fragment.GetExcitationEnergy() > 0.0) {
	EmissionProbability = CalculateProbability(fragment, MaximalKineticEnergy);

    }
    return EmissionProbability;
}

////////////////////////////////////

// Computes the integrated probability of evaporation channel
G4double G4EvaporationProbability::CalculateProbability(const G4Fragment & fragment, const G4double MaximalKineticEnergy)
{
    G4double ResidualA = fragment.GetA() - theA;
    G4double ResidualZ = fragment.GetZ() - theZ;
    G4double U = fragment.GetExcitationEnergy();
   
 if (OPTxs==0) {

	
    G4double NuclearMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetNucleusMass(theZ,theA);


    G4double delta0 = G4PairingCorrection::GetInstance()->GetPairingCorrection(static_cast<G4int>(fragment.GetA()),
									       static_cast<G4int>(fragment.GetZ()));

    G4double SystemEntropy = 2.0*std::sqrt(theEvapLDPptr->LevelDensityParameter(static_cast<G4int>(fragment.GetA()),
									   static_cast<G4int>(fragment.GetZ()),U)*
				      (U-delta0));
								  

    G4double RN = 1.5*fermi;

    G4double Alpha = CalcAlphaParam(fragment);
    G4double Beta = CalcBetaParam(fragment);
	
    G4double Rmax = MaximalKineticEnergy;
    G4double a = theEvapLDPptr->LevelDensityParameter(static_cast<G4int>(ResidualA),
						      static_cast<G4int>(ResidualZ),
						      Rmax);
    G4double GlobalFactor = static_cast<G4double>(Gamma) * (Alpha/(a*a)) *
	(NuclearMass*RN*RN*std::pow(ResidualA,2./3.))/
	(2.*pi* hbar_Planck*hbar_Planck);
    G4double Term1 = (2.0*Beta*a-3.0)/2.0 + Rmax*a;
    G4double Term2 = (2.0*Beta*a-3.0)*std::sqrt(Rmax*a) + 2.0*a*Rmax;
	
    G4double ExpTerm1 = 0.0;
    if (SystemEntropy <= 600.0) ExpTerm1 = std::exp(-SystemEntropy);
	
    G4double ExpTerm2 = 2.*std::sqrt(a*Rmax) - SystemEntropy;
    if (ExpTerm2 > 700.0) ExpTerm2 = 700.0;
    ExpTerm2 = std::exp(ExpTerm2);
	
    G4double Width = GlobalFactor*(Term1*ExpTerm1 + Term2*ExpTerm2);
	
    return Width;
             
 } else if (OPTxs==1 || OPTxs==2 ||OPTxs==3 || OPTxs==4) {

   G4double limit;
   if (useSICB) limit=theCoulombBarrierptr->GetCoulombBarrier(G4lrint(ResidualA),G4lrint(ResidualZ),U);
   else limit=0.;

   if (MaximalKineticEnergy <= limit)  return 0.0;


   // if Coulomb barrier cutoff is superimposed for all cross sections the limit is the Coulomb Barrier
   G4double LowerLimit= limit;

   //  MaximalKineticEnergy: asimptotic value (already accounted for in G4EvaporationChannel)     

   G4double UpperLimit = MaximalKineticEnergy;


   G4double Width = IntegrateEmissionProbability(fragment,LowerLimit,UpperLimit);

   return Width;
 } else{
   std::ostringstream errOs;
   errOs << "Bad option for cross sections at evaporation"  <<G4endl;
   throw G4HadronicException(__FILE__, __LINE__, errOs.str());
 }
  
}

/////////////////////////////////////////////////////////////////////

G4double G4EvaporationProbability::
IntegrateEmissionProbability(const G4Fragment & fragment, const G4double & Low, const G4double & Up )
{

  static const G4int N = 10;
  // 10-Points Gauss-Legendre abcisas and weights
  static const G4double w[N] = {
    0.0666713443086881,
    0.149451349150581,
    0.219086362515982,
    0.269266719309996,
    0.295524224714753,
    0.295524224714753,
    0.269266719309996,
    0.219086362515982,
    0.149451349150581,
    0.0666713443086881
  };
  static const G4double x[N] = {
    -0.973906528517172,
    -0.865063366688985,
    -0.679409568299024,
    -0.433395394129247,
    -0.148874338981631,
    0.148874338981631,
    0.433395394129247,
    0.679409568299024,
    0.865063366688985,
    0.973906528517172
  };

  G4double Total = 0.0;


  for (G4int i = 0; i < N; i++) 
    {

      G4double KineticE = ((Up-Low)*x[i]+(Up+Low))/2.0;

      Total += w[i]*ProbabilityDistributionFunction(fragment, KineticE);

    }
  Total *= (Up-Low)/2.0;
  return Total;
}


/////////////////////////////////////////////////////////
//New method (OPT=1,2,3,4)

G4double G4EvaporationProbability::ProbabilityDistributionFunction( const G4Fragment & fragment, const G4double K)
{ 

 


  G4double ResidualA = fragment.GetA() - theA;
  G4double ResidualZ = fragment.GetZ() - theZ;  
  G4double U = fragment.GetExcitationEnergy();

  //        if(K <= theCoulombBarrierptr->GetCoulombBarrier(G4lrint(ResidualA),G4lrint(ResidualZ),U)) return 0.0;   

  G4double delta1 = G4PairingCorrection::GetInstance()->GetPairingCorrection(static_cast<G4int>(ResidualA),static_cast<G4int>(ResidualZ));

 
  G4double delta0 = G4PairingCorrection::GetInstance()->GetPairingCorrection(static_cast<G4int>(fragment.GetA()),static_cast<G4int>(fragment.GetZ()));

  
  G4double ParticleMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetNucleusMass(theZ,theA);

  G4double theSeparationEnergy= G4NucleiProperties::GetMassExcess(static_cast<G4int>(theA),static_cast<G4int>(theZ)) +
    G4NucleiProperties::GetMassExcess(static_cast<G4int>(ResidualA),static_cast<G4int>(ResidualZ)) -
    G4NucleiProperties::GetMassExcess(static_cast<G4int>(fragment.GetA()),static_cast<G4int>(fragment.GetZ()));

  G4double a0 = theEvapLDPptr->LevelDensityParameter(static_cast<G4int>(fragment.GetA()),static_cast<G4int>(fragment.GetZ()),U - delta0);

  G4double a1 = theEvapLDPptr->LevelDensityParameter(static_cast<G4int>(ResidualA),static_cast<G4int>(ResidualZ),U - theSeparationEnergy - delta1);
  
  
  G4double E0=U-delta0;

  G4double E1=U-theSeparationEnergy-delta1-K;

  if (E1<0.) return 0.;

  //JMQ 14/02/09 BUG fixed: hbarc should be in the denominator instead of hbar_Planck 
  //Without 1/hbar_Panck remains as a width
  //  G4double  Prob=Gamma*ParticleMass/((pi*hbar_Planck)*(pi*hbar_Planck)*
  //std::exp(2*std::sqrt(a0*E0)))*K*CrossSection(fragment,K)*std::exp(2*std::sqrt(a1*E1))*millibarn;

  G4double Prob=Gamma*ParticleMass/((pi*hbarc)*(pi*hbarc)*std::exp(2*std::sqrt(a0*E0)))
    *K*CrossSection(fragment,K)*std::exp(2*std::sqrt(a1*E1))*millibarn;

  return Prob;
}


