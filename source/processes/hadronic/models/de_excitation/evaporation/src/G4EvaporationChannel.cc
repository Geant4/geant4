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
// $Id: G4EvaporationChannel.cc 85841 2014-11-05 15:35:06Z gcosmo $
//
//J.M. Quesada (August2008). Based on:
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modified:
// 03-09-2008 J.M. Quesada for external choice of inverse cross section option
// 06-09-2008 J.M. Quesada Also external choices have been added for superimposed 
//                 Coulomb barrier (if useSICB is set true, by default is false) 
// 17-11-2010 V.Ivanchenko in constructor replace G4VEmissionProbability by 
//            G4EvaporationProbability and do not new and delete probability
//            object at each call; use G4Pow

#include "G4EvaporationChannel.hh"
#include "G4PairingCorrection.hh"
#include "G4NucleiProperties.hh"
#include "G4Pow.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4EvaporationLevelDensityParameter.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Alpha.hh"

G4EvaporationChannel::G4EvaporationChannel(G4int anA, G4int aZ, 
					   const G4String & aName,
					   G4EvaporationProbability* aEmissionStrategy,
                                           G4VCoulombBarrier* aCoulombBarrier):
    G4VEvaporationChannel(aName),
    theA(anA),
    theZ(aZ),
    theEvaporationProbabilityPtr(aEmissionStrategy),
    theCoulombBarrierPtr(aCoulombBarrier),
    EmissionProbability(0.0),
    MaximalKineticEnergy(-1000.0)
{ 
  ResidualA = 0;
  ResidualZ = 0;
  ResidualMass = CoulombBarrier = 0.0;
  EvaporatedMass = G4NucleiProperties::GetNuclearMass(theA, theZ);
  theLevelDensityPtr = new G4EvaporationLevelDensityParameter;
  pairingCorrection = G4PairingCorrection::GetInstance();
}

G4EvaporationChannel::~G4EvaporationChannel()
{
  delete theLevelDensityPtr;
}

void G4EvaporationChannel::Initialise()
{
  //for inverse cross section choice
  theEvaporationProbabilityPtr->SetOPTxs(OPTxs);
  // for superimposed Coulomb Barrier for inverse cross sections
  theEvaporationProbabilityPtr->UseSICB(useSICB);

  G4VEvaporationChannel::Initialise();  
}

G4double G4EvaporationChannel::GetEmissionProbability(G4Fragment* fragment)
{
  G4int FragmentA = fragment->GetA_asInt();
  G4int FragmentZ = fragment->GetZ_asInt();
  ResidualA = FragmentA - theA;
  ResidualZ = FragmentZ - theZ;
  //G4cout << "G4EvaporationChannel::Initialize Z= " << theZ << " A= " << theA 
  //	 << " FragZ= " << FragmentZ << " FragA= " << FragmentA << G4endl;
  EmissionProbability = 0.0;

  // Only channels which are physically allowed are taken into account 
  if (ResidualA >= ResidualZ && ResidualZ > 0 && ResidualA >= theA) {
  
    //Effective excitation energy
    G4double ExEnergy = fragment->GetExcitationEnergy() - 
      pairingCorrection->GetPairingCorrection(FragmentA,FragmentZ);
    ResidualMass = G4NucleiProperties::GetNuclearMass(ResidualA, ResidualZ);
    G4double FragmentMass = fragment->GetGroundStateMass();
    G4double Etot = FragmentMass + ExEnergy;
  
    if(ExEnergy > 0.0 && Etot > ResidualMass + EvaporatedMass) {
  
      // Maximal Kinetic Energy
      MaximalKineticEnergy = ((Etot-ResidualMass)*(Etot+ResidualMass) 
	    + EvaporatedMass*EvaporatedMass)/(2.0*Etot) - EvaporatedMass;

      // Emission probability
      // Protection for the case Tmax<V. If not set in this way we could end up in an 
      // infinite loop in  the method GetKineticEnergy if OPTxs!=0 && useSICB=true. 
      // Of course for OPTxs=0 we have the Coulomb barrier 

      CoulombBarrier = 0.0;
      if (OPTxs==0 || (OPTxs!=0 && useSICB)) {
	CoulombBarrier = 
	  theCoulombBarrierPtr->GetCoulombBarrier(ResidualA,ResidualZ,ExEnergy);
      }
      // The threshold for charged particle emission must be  set to 0 if Coulomb 
      //cutoff  is included in the cross sections
      if (MaximalKineticEnergy > CoulombBarrier) {
	EmissionProbability = theEvaporationProbabilityPtr->
	  EmissionProbability(*fragment, MaximalKineticEnergy);
      }
    }
  }
  //G4cout << "G4EvaporationChannel:: probability= " << EmissionProbability << G4endl;   
  return EmissionProbability;
}

G4Fragment* G4EvaporationChannel::EmittedFragment(G4Fragment* theNucleus)
{
  G4Fragment* evFragment = 0;
  G4double evEnergy = SampleKineticEnergy(*theNucleus) + EvaporatedMass;

  G4ThreeVector momentum(IsotropicVector
    (std::sqrt((evEnergy - EvaporatedMass)*(evEnergy + EvaporatedMass))));
  
  G4LorentzVector EvaporatedMomentum(momentum, evEnergy);
  G4LorentzVector ResidualMomentum = theNucleus->GetMomentum();
  EvaporatedMomentum.boost(ResidualMomentum.boostVector());
  
  evFragment = new G4Fragment(theA,theZ,EvaporatedMomentum);
  ResidualMomentum -= EvaporatedMomentum;
  theNucleus->SetZandA_asInt(ResidualZ, ResidualA);
  theNucleus->SetMomentum(ResidualMomentum);

  return evFragment; 
} 

G4FragmentVector * G4EvaporationChannel::BreakUp(const G4Fragment & theNucleus)
{
  G4FragmentVector * theResult = new G4FragmentVector();
  G4Fragment* frag0 = new G4Fragment(theNucleus);
  G4Fragment* frag1 = EmittedFragment(frag0);
  if(frag1) { theResult->push_back(frag1); }
  theResult->push_back(frag0);
  return theResult;
} 

///////////////////////////////////////////
//JMQ: New method for MC sampling of kinetic energy. 
G4double G4EvaporationChannel::SampleKineticEnergy(const G4Fragment & aFragment)
{
  G4double T = 0.0;
  if (OPTxs==0) {
    // It uses Dostrovsky's approximation for the inverse reaction cross
    // in the probability for fragment emission
    // MaximalKineticEnergy energy in the original version (V.Lara) was calculated at 
    //the Coulomb barrier.
    
    G4double Rb = 4.0*theLevelDensityPtr->
      LevelDensityParameter(ResidualA+theA,ResidualZ+theZ,MaximalKineticEnergy)*
      MaximalKineticEnergy;
    G4double RbSqrt = std::sqrt(Rb);
    G4double PEX1 = 0.0;
    if (RbSqrt < 160.0) PEX1 = G4Exp(-RbSqrt);
    G4double Rk = 0.0;
    G4double FRk = 0.0;
    do {
      G4double RandNumber = G4UniformRand();
      Rk = 1.0 + (1./RbSqrt)*G4Log(RandNumber + (1.0-RandNumber)*PEX1);
      G4double Q1 = 1.0;
      G4double Q2 = 1.0;
      if (theZ == 0) { // for emitted neutron
        G4double Beta = (2.12/G4Pow::GetInstance()->Z23(ResidualA) - 0.05)*MeV/
          (0.76 + 2.2/G4Pow::GetInstance()->Z13(ResidualA));
        Q1 = 1.0 + Beta/(MaximalKineticEnergy);
        Q2 = Q1*std::sqrt(Q1);
      } 
      
      FRk = (3.0*std::sqrt(3.0)/2.0)/Q2 * Rk * (Q1 - Rk*Rk);
      
    } while (FRk < G4UniformRand());
    
    T =  MaximalKineticEnergy * (1.0-Rk*Rk) + CoulombBarrier;

  } else {    
    // Coulomb barrier is just included  in the cross sections
    G4double prob;
    do {  
      T=CoulombBarrier+G4UniformRand()*(MaximalKineticEnergy-CoulombBarrier);
      prob = theEvaporationProbabilityPtr->ProbabilityDistributionFunction(aFragment,T);
    } while (EmissionProbability*G4UniformRand() >= prob);
  }
  return T;
}

G4ThreeVector G4EvaporationChannel::IsotropicVector(G4double Magnitude)
    // Samples a isotropic random vectorwith a magnitud given by Magnitude.
    // By default Magnitude = 1.0
{
  G4double CosTheta = 1.0 - 2.0*G4UniformRand();
  G4double SinTheta = std::sqrt(1.0 - CosTheta*CosTheta);
  G4double Phi = twopi*G4UniformRand();
  G4ThreeVector Vector(Magnitude*std::cos(Phi)*SinTheta,
		       Magnitude*std::sin(Phi)*SinTheta,
		       Magnitude*CosTheta);
  return Vector;
}
