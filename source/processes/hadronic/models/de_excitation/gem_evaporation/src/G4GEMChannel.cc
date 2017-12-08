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
// $Id: G4GEMChannel.cc 107060 2017-11-01 15:00:04Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// J. M. Quesada (September 2009) bugs fixed in  probability distribution for kinetic 
//              energy sampling:
//  		-hbarc instead of hbar_Planck (BIG BUG)
//              -quantities for initial nucleus and residual are calculated separately
// V.Ivanchenko (September 2009) Added proper protection for unphysical final state 
// J. M. Quesada (October 2009) fixed bug in CoulombBarrier calculation 

#include "G4GEMChannel.hh"
#include "G4VCoulombBarrier.hh"
#include "G4GEMCoulombBarrier.hh"
#include "G4NuclearLevelData.hh"
#include "G4PairingCorrection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4RandomDirection.hh"

G4GEMChannel::G4GEMChannel(G4int theA, G4int theZ, const G4String & aName,
                           G4GEMProbability * aEmissionStrategy) :
  G4VEvaporationChannel(aName),
  A(theA),
  Z(theZ),
  theEvaporationProbabilityPtr(aEmissionStrategy),
  EmissionProbability(0.0),
  MaximalKineticEnergy(-CLHEP::GeV)
{ 
  theCoulombBarrierPtr = new G4GEMCoulombBarrier(theA, theZ);
  theEvaporationProbabilityPtr->SetCoulomBarrier(theCoulombBarrierPtr);
  theLevelDensityPtr = new G4EvaporationLevelDensityParameter;
  MyOwnLevelDensity = true;
  EvaporatedMass = G4NucleiProperties::GetNuclearMass(A, Z);
  ResidualMass = CoulombBarrier = 0.0;
  fG4pow = G4Pow::GetInstance(); 
  ResidualZ = ResidualA = 0;
  pairingCorrection = 
    G4NuclearLevelData::GetInstance()->GetPairingCorrection();
}

G4GEMChannel::~G4GEMChannel()
{
  if (MyOwnLevelDensity) { delete theLevelDensityPtr; }
  delete theCoulombBarrierPtr;
}

G4double G4GEMChannel::GetEmissionProbability(G4Fragment* fragment)
{
  G4int anA = fragment->GetA_asInt();
  G4int aZ  = fragment->GetZ_asInt();
  ResidualA = anA - A;
  ResidualZ = aZ - Z;
  /*
  G4cout << "G4GEMChannel: Z= " << Z << "  A= " << A 
	 << " FragmentZ= " << aZ << " FragmentA= " << anA
	 << " Zres= " << ResidualZ << " Ares= " << ResidualA 
	 << G4endl; 
  */
  // We only take into account channels which are physically allowed
  EmissionProbability = 0.0;

  // Only channels which are physically allowed are taken into account 
  if (ResidualA >= ResidualZ && ResidualZ > 0 && ResidualA >= A) {
  
    //Effective excitation energy
    G4double ExEnergy = fragment->GetExcitationEnergy()
      - pairingCorrection->GetPairingCorrection(anA, aZ);
    if(ExEnergy > 0.0) { 
      ResidualMass = G4NucleiProperties::GetNuclearMass(ResidualA, ResidualZ);
      G4double FragmentMass = fragment->GetGroundStateMass();
      G4double Etot = FragmentMass + ExEnergy;
      // Coulomb Barrier calculation
      CoulombBarrier = 
	theCoulombBarrierPtr->GetCoulombBarrier(ResidualA,ResidualZ,ExEnergy);
      /*  
      G4cout << "Eexc(MeV)= " << ExEnergy/MeV 
	     << " CoulBarrier(MeV)= " << CoulombBarrier/MeV << G4endl;
      */
      if(Etot > ResidualMass + EvaporatedMass + CoulombBarrier) {
  
	// Maximal Kinetic Energy
	MaximalKineticEnergy = ((Etot-ResidualMass)*(Etot+ResidualMass) 
	  + EvaporatedMass*EvaporatedMass)/(2.0*Etot) 
	  - EvaporatedMass - CoulombBarrier;

	//G4cout << "CBarrier(MeV)= " << CoulombBarrier/MeV << G4endl;

	if (MaximalKineticEnergy > 0.0) { 
	  // Total emission probability for this channel
	  EmissionProbability = theEvaporationProbabilityPtr->
	    EmissionProbability(*fragment, MaximalKineticEnergy);
	}
      }
    }
  }   
  //G4cout << "Prob= " << EmissionProbability << G4endl;
  return EmissionProbability;
}

G4Fragment* G4GEMChannel::EmittedFragment(G4Fragment* theNucleus)
{
  G4Fragment* evFragment = 0;
  G4double evEnergy = SampleKineticEnergy(*theNucleus) + EvaporatedMass;

  G4ThreeVector momentum = G4RandomDirection()*
    std::sqrt((evEnergy - EvaporatedMass)*(evEnergy + EvaporatedMass));
  
  G4LorentzVector EvaporatedMomentum(momentum, evEnergy);
  G4LorentzVector ResidualMomentum = theNucleus->GetMomentum();
  EvaporatedMomentum.boost(ResidualMomentum.boostVector());
  
  evFragment = new G4Fragment(A, Z, EvaporatedMomentum);
  ResidualMomentum -= EvaporatedMomentum;
  theNucleus->SetZandA_asInt(ResidualZ, ResidualA);
  theNucleus->SetMomentum(ResidualMomentum);

  return evFragment; 
} 

G4double G4GEMChannel::SampleKineticEnergy(const G4Fragment & fragment)
// Samples fragment kinetic energy.
{
  G4double U = fragment.GetExcitationEnergy();
  
  G4double Alpha = theEvaporationProbabilityPtr->CalcAlphaParam(fragment);
  G4double Beta = theEvaporationProbabilityPtr->CalcBetaParam(fragment);

  //                       ***RESIDUAL***
  //JMQ (September 2009) the following quantities  refer to the RESIDUAL:
  G4double delta0 = pairingCorrection->GetPairingCorrection(ResidualA,ResidualZ);
  G4double Ux = (2.5 + 150.0/ResidualA)*MeV;
  G4double Ex = Ux + delta0;
  G4double InitialLevelDensity;
  //                    ***end RESIDUAL ***
  
  //                       ***PARENT***
  //JMQ (September 2009) the following quantities   refer to the PARENT:
  
  G4double deltaCN = pairingCorrection->GetPairingCorrection(fragment.GetA_asInt(),
							     fragment.GetZ_asInt());
  G4double aCN = theLevelDensityPtr->LevelDensityParameter(fragment.GetA_asInt(),
							   fragment.GetZ_asInt(),
							   U-deltaCN);   
  G4double UxCN = (2.5 + 150.0/G4double(fragment.GetA_asInt()))*MeV;
  G4double ExCN = UxCN + deltaCN;
  G4double TCN = 1.0/(std::sqrt(aCN/UxCN) - 1.5/UxCN);
  //                       ***end PARENT***
  
  //JMQ quantities calculated for  CN in InitialLevelDensity
  if ( U < ExCN ) 
    {
      G4double E0CN = ExCN - TCN*(G4Log(TCN/MeV) - G4Log(aCN*MeV)/4.0 
				  - 1.25*G4Log(UxCN/MeV) + 2.0*std::sqrt(aCN*UxCN));
      InitialLevelDensity = (pi/12.0)*G4Exp((U-E0CN)/TCN)/TCN;
    } 
  else 
    {
      G4double x  = U-deltaCN;
      G4double x1 = std::sqrt(aCN*x);
      InitialLevelDensity = (pi/12.0)*G4Exp(2*x1)/(x*std::sqrt(x1));
    }
  
  G4double Spin = theEvaporationProbabilityPtr->GetSpin();
  //JMQ  BIG BUG fixed: hbarc instead of hbar_Planck !!!!
  //     it was fixed in total probability (for this channel) but remained still here!!
  //    G4double g = (2.0*Spin+1.0)*NuclearMass/(pi2* hbar_Planck*hbar_Planck);
  G4double gg = (2.0*Spin+1.0)*EvaporatedMass/(pi2* hbarc*hbarc);
  //
  //JMQ  fix on Rb and  geometrical cross sections according to Furihata's paper 
  //                      (JAERI-Data/Code 2001-105, p6)
  G4double Rb = 0.0; 
  if (A > 4) 
    {
      G4double Ad = fG4pow->Z13(ResidualA);
      G4double Aj = fG4pow->Z13(A);
      Rb = (1.12*(Aj + Ad) - 0.86*((Aj+Ad)/(Aj*Ad))+2.85)*fermi;
    }
  else if (A>1)
    {
      G4double Ad = fG4pow->Z13(ResidualA);
      G4double Aj = fG4pow->Z13(A);
      Rb=1.5*(Aj+Ad)*fermi;
    }
  else 
    {
      G4double Ad = fG4pow->Z13(ResidualA);
      Rb = 1.5*Ad*fermi;
    }
  G4double GeometricalXS = pi*Rb*Rb; 
    
  G4double ConstantFactor = gg*GeometricalXS*Alpha*pi/(InitialLevelDensity*12);
  //JMQ : this is the assimptotic maximal kinetic energy of the ejectile 
  //      (obtained by energy-momentum conservation). 
  //      In general smaller than U-theSeparationEnergy 
  G4double theEnergy = MaximalKineticEnergy + CoulombBarrier;
  G4double KineticEnergy;
  G4double Probability;

  for(G4int i=0; i<100; ++i) {
    KineticEnergy =  CoulombBarrier + G4UniformRand()*(MaximalKineticEnergy);
    G4double edelta = theEnergy-KineticEnergy-delta0;
    Probability = ConstantFactor*(KineticEnergy + Beta);
    G4double a = 
      theLevelDensityPtr->LevelDensityParameter(ResidualA,ResidualZ,edelta);
    G4double T = 1.0/(std::sqrt(a/Ux) - 1.5/Ux);
    //JMQ fix in units
	
    if (theEnergy - KineticEnergy < Ex) {
      G4double E0 = Ex - T*(G4Log(T) - G4Log(a)*0.25
			    - 1.25*G4Log(Ux) + 2.0*std::sqrt(a*Ux));
      Probability *= G4Exp((theEnergy-KineticEnergy-E0)/T)/T;
    } else {
      G4double e2 = edelta*edelta;
      Probability *= 
	G4Exp(2*std::sqrt(a*edelta) - 0.25*G4Log(a*edelta*e2*e2));
    }
    if(EmissionProbability*G4UniformRand() <= Probability) { break; }
  }
    
  return KineticEnergy;
} 

void G4GEMChannel::Dump() const
{
  theEvaporationProbabilityPtr->Dump();
}



