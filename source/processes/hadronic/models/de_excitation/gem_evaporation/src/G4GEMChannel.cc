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
// $Id: G4GEMChannel.cc,v 1.12 2010-11-18 16:21:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4PairingCorrection.hh"
#include "G4Pow.hh"

G4GEMChannel::G4GEMChannel(G4int theA, G4int theZ, const G4String & aName,
                           G4GEMProbability * aEmissionStrategy,
                           G4VCoulombBarrier * aCoulombBarrier) :
  G4VEvaporationChannel(aName),
  A(theA),
  Z(theZ),
  theEvaporationProbabilityPtr(aEmissionStrategy),
  theCoulombBarrierPtr(aCoulombBarrier),
  EmissionProbability(0.0),
  MaximalKineticEnergy(-1000.0)
{ 
  theLevelDensityPtr = new G4EvaporationLevelDensityParameter;
  MyOwnLevelDensity = true;
  EvaporatedMass = G4NucleiProperties::GetNuclearMass(A, Z);
  ResidualMass = CoulombBarrier = 0.0;
  fG4pow = G4Pow::GetInstance(); 
}

G4GEMChannel::G4GEMChannel() :
  G4VEvaporationChannel("dummy"),
  A(0),
  Z(0),
  theEvaporationProbabilityPtr(0),
  theCoulombBarrierPtr(0),
  EmissionProbability(0.0),
  MaximalKineticEnergy(-1000.0)
{ 
  theLevelDensityPtr = 0;
  MyOwnLevelDensity = true;
  EvaporatedMass = 0.0;
  ResidualMass = CoulombBarrier = 0.0;
  fG4pow = G4Pow::GetInstance(); 
}

G4GEMChannel::~G4GEMChannel()
{
  if (MyOwnLevelDensity) { delete theLevelDensityPtr; }
}

void G4GEMChannel::Initialize(const G4Fragment & fragment)
{
  G4int anA = fragment.GetA_asInt();
  G4int aZ  = fragment.GetZ_asInt();
  ResidualA = anA - A;
  ResidualZ = aZ - Z;
  //G4cout << "G4GEMChannel::Initialize: Z= " << aZ << " A= " << anA
  //	   << " Zres= " << ResidualZ << " Ares= " << ResidualA << G4endl; 

  // We only take into account channels which are physically allowed
  if (ResidualA <= 0 || ResidualZ <= 0 || ResidualA < ResidualZ ||
      (ResidualA == ResidualZ && ResidualA > 1)) 
    {
      CoulombBarrier = 0.0;
      MaximalKineticEnergy = -1000.0*MeV;
      EmissionProbability = 0.0;
    } 
  else 
    {
      // Effective excitation energy
      // JMQ 071009: pairing in ExEnergy should be the one of parent compound nucleus 
      // FIXED the bug causing reported crash by VI (negative Probabilities 
      // due to inconsistency in Coulomb barrier calculation (CoulombBarrier and -Beta 
      // param for protons must be the same)   
      //    G4double ExEnergy = fragment.GetExcitationEnergy() -
      //    G4PairingCorrection::GetInstance()->GetPairingCorrection(ResidualA,ResidualZ);
      G4double ExEnergy = fragment.GetExcitationEnergy() -
	G4PairingCorrection::GetInstance()->GetPairingCorrection(anA,aZ);

      //G4cout << "Eexc(MeV)= " << ExEnergy/MeV << G4endl;

      if( ExEnergy <= 0.0) {
	CoulombBarrier = 0.0;
	MaximalKineticEnergy = -1000.0*MeV;
	EmissionProbability = 0.0;

      } else {

	ResidualMass = G4NucleiProperties::GetNuclearMass(ResidualA, ResidualZ);

        // Coulomb Barrier calculation
	CoulombBarrier = theCoulombBarrierPtr->GetCoulombBarrier(ResidualA,ResidualZ,ExEnergy);
	//G4cout << "CBarrier(MeV)= " << CoulombBarrier/MeV << G4endl;

	//Maximal kinetic energy (JMQ : at the Coulomb barrier)
	MaximalKineticEnergy = 
	  CalcMaximalKineticEnergy(fragment.GetGroundStateMass()+ExEnergy);
	//G4cout << "MaxE(MeV)= " << MaximalKineticEnergy/MeV << G4endl;
		
	// Emission probability
	if (MaximalKineticEnergy <= 0.0) 
	  {
	    EmissionProbability = 0.0;
	  }
	else 
	  { 
	    // Total emission probability for this channel
	    EmissionProbability = 
	      theEvaporationProbabilityPtr->EmissionProbability(fragment,MaximalKineticEnergy);
	  }
      }
    }   
    //G4cout << "Prob= " << EmissionProbability << G4endl;
    return;
}

G4FragmentVector * G4GEMChannel::BreakUp(const G4Fragment & theNucleus)
{
  G4double EvaporatedKineticEnergy = CalcKineticEnergy(theNucleus);
  G4double EvaporatedEnergy = EvaporatedKineticEnergy + EvaporatedMass;
  
  G4ThreeVector momentum(IsotropicVector(std::sqrt(EvaporatedKineticEnergy*
						   (EvaporatedKineticEnergy+2.0*EvaporatedMass))));
    
  momentum.rotateUz(theNucleus.GetMomentum().vect().unit());

  G4LorentzVector EvaporatedMomentum(momentum,EvaporatedEnergy);
  EvaporatedMomentum.boost(theNucleus.GetMomentum().boostVector());
  G4Fragment * EvaporatedFragment = new G4Fragment(A,Z,EvaporatedMomentum);
  // ** And now the residual nucleus ** 
  G4double theExEnergy = theNucleus.GetExcitationEnergy();
  G4double theMass = theNucleus.GetGroundStateMass();
  G4double ResidualEnergy = theMass + (theExEnergy - EvaporatedKineticEnergy) - EvaporatedMass;
	
  G4LorentzVector ResidualMomentum(-momentum,ResidualEnergy);
  ResidualMomentum.boost(theNucleus.GetMomentum().boostVector());
	
  G4Fragment * ResidualFragment = new G4Fragment( ResidualA, ResidualZ, ResidualMomentum );
    
  G4FragmentVector * theResult = new G4FragmentVector;
    
#ifdef debug
  G4double Efinal = ResidualMomentum.e() + EvaporatedMomentum.e();
  G4ThreeVector Pfinal = ResidualMomentum.vect() + EvaporatedMomentum.vect();
  if (std::abs(Efinal-theNucleus.GetMomentum().e()) > 10.0*eV) {
    G4cout << "@@@@@@@@@@@@@@@@@@@@@ G4GEMChanel: ENERGY @@@@@@@@@@@@@@@@" << G4endl;
    G4cout << "Initial : " << theNucleus.GetMomentum().e()/MeV << " MeV    Final :" 
	   <<Efinal/MeV << " MeV    Delta: " <<  (Efinal-theNucleus.GetMomentum().e())/MeV 
	   << " MeV" << G4endl;
  }
  if (std::abs(Pfinal.x()-theNucleus.GetMomentum().x()) > 10.0*eV ||
      std::abs(Pfinal.y()-theNucleus.GetMomentum().y()) > 10.0*eV ||
      std::abs(Pfinal.z()-theNucleus.GetMomentum().z()) > 10.0*eV ) {
    G4cout << "@@@@@@@@@@@@@@@@@@@@@ G4GEMChanel: MOMENTUM @@@@@@@@@@@@@@@@" << G4endl;
    G4cout << "Initial : " << theNucleus.GetMomentum().vect() << " MeV    Final :" 
	   <<Pfinal/MeV << " MeV    Delta: " <<  Pfinal-theNucleus.GetMomentum().vect()
	   << " MeV" << G4endl;   
        
  }
#endif
  theResult->push_back(EvaporatedFragment);
  theResult->push_back(ResidualFragment);
  return theResult; 
} 

G4double G4GEMChannel::CalcMaximalKineticEnergy(const G4double NucleusTotalE)
// Calculate maximal kinetic energy that can be carried by fragment.
//JMQ this is not the assimptotic kinetic energy but the one at the Coulomb barrier
{
  G4double T = (NucleusTotalE*NucleusTotalE + EvaporatedMass*EvaporatedMass - ResidualMass*ResidualMass)/
    (2.0*NucleusTotalE) - EvaporatedMass - CoulombBarrier;
	
  return T;
}

G4double G4GEMChannel::CalcKineticEnergy(const G4Fragment & fragment)
// Samples fragment kinetic energy.
{
  G4double U = fragment.GetExcitationEnergy();
  
  G4double Alpha = theEvaporationProbabilityPtr->CalcAlphaParam(fragment);
  G4double Beta = theEvaporationProbabilityPtr->CalcBetaParam(fragment);

  G4double Normalization = theEvaporationProbabilityPtr->GetNormalization();

  //                       ***RESIDUAL***
  //JMQ (September 2009) the following quantities  refer to the RESIDUAL:
  G4double delta0 = G4PairingCorrection::GetInstance()->GetPairingCorrection(ResidualA,ResidualZ);
  G4double Ux = (2.5 + 150.0/ResidualA)*MeV;
  G4double Ex = Ux + delta0;
  G4double InitialLevelDensity;
  //                    ***end RESIDUAL ***
  
  //                       ***PARENT***
  //JMQ (September 2009) the following quantities   refer to the PARENT:
  
  G4double deltaCN = 
    G4PairingCorrection::GetInstance()->GetPairingCorrection(fragment.GetA_asInt(),
							     fragment.GetZ_asInt());
  G4double aCN = theLevelDensityPtr->LevelDensityParameter(fragment.GetA_asInt(),
							   fragment.GetZ_asInt(),U-deltaCN);   
  G4double UxCN = (2.5 + 150.0/fragment.GetA())*MeV;
  G4double ExCN = UxCN + deltaCN;
  G4double TCN = 1.0/(std::sqrt(aCN/UxCN) - 1.5/UxCN);
  //                       ***end PARENT***
  
  //JMQ quantities calculated for  CN in InitialLevelDensity
  if ( U < ExCN ) 
    {
      G4double E0CN = ExCN - TCN*(std::log(TCN/MeV) - std::log(aCN*MeV)/4.0 
				  - 1.25*std::log(UxCN/MeV) + 2.0*std::sqrt(aCN*UxCN));
      InitialLevelDensity = (pi/12.0)*std::exp((U-E0CN)/TCN)/TCN;
    } 
  else 
    {
      G4double x  = U-deltaCN;
      G4double x1 = std::sqrt(aCN*x);
      InitialLevelDensity = (pi/12.0)*std::exp(2*x1)/(x*std::sqrt(x1));
      //InitialLevelDensity = 
      //(pi/12.0)*std::exp(2*std::sqrt(aCN*(U-deltaCN)))/std::pow(aCN*std::pow(U-deltaCN,5.0),1.0/4.0);
    }
  
  const G4double Spin = theEvaporationProbabilityPtr->GetSpin();
  //JMQ  BIG BUG fixed: hbarc instead of hbar_Planck !!!!
  //     it was fixed in total probability (for this channel) but remained still here!!
  //    G4double g = (2.0*Spin+1.0)*NuclearMass/(pi2* hbar_Planck*hbar_Planck);
  G4double g = (2.0*Spin+1.0)*EvaporatedMass/(pi2* hbarc*hbarc);
  //
  //JMQ  fix on Rb and  geometrical cross sections according to Furihata's paper 
  //                      (JAERI-Data/Code 2001-105, p6)
  G4double Rb = 0.0; 
  if (A > 4) 
    {
      G4double Ad = fG4pow->Z13(ResidualA);
      G4double Aj = fG4pow->Z13(A);
      //        RN = 1.12*(R1 + R2) - 0.86*((R1+R2)/(R1*R2));
      Rb = 1.12*(Aj + Ad) - 0.86*((Aj+Ad)/(Aj*Ad))+2.85;
      Rb *= fermi;
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
  //   G4double GeometricalXS = pi*RN*RN*std::pow(G4double(fragment.GetA()-A),2./3.); 
  G4double GeometricalXS = pi*Rb*Rb; 
    
  G4double ConstantFactor = g*GeometricalXS*Alpha/InitialLevelDensity;
  ConstantFactor *= pi/12.0;
  //JMQ : this is the assimptotic maximal kinetic energy of the ejectile 
  //      (obtained by energy-momentum conservation). 
  //      In general smaller than U-theSeparationEnergy 
  G4double theEnergy = MaximalKineticEnergy + CoulombBarrier;
  G4double KineticEnergy;
  G4double Probability;

  do 
    {
      KineticEnergy =  CoulombBarrier + G4UniformRand()*MaximalKineticEnergy;
      Probability = ConstantFactor*(KineticEnergy + Beta);
      G4double a = 
	theLevelDensityPtr->LevelDensityParameter(ResidualA,ResidualZ,theEnergy-KineticEnergy-delta0);
      G4double T = 1.0/(std::sqrt(a/Ux) - 1.5/Ux);
      //JMQ fix in units
	
      if ( theEnergy-KineticEnergy < Ex) 
	{
	  G4double E0 = Ex - T*(std::log(T/MeV) - std::log(a*MeV)/4.0 
				- 1.25*std::log(Ux/MeV) + 2.0*std::sqrt(a*Ux));
	  Probability *= std::exp((theEnergy-KineticEnergy-E0)/T)/T;
	} 
      else 
	{
	  Probability *= std::exp(2*std::sqrt(a*(theEnergy-KineticEnergy-delta0)))/
	    std::pow(a*fG4pow->powN(theEnergy-KineticEnergy-delta0,5), 0.25);
	}
    }
  while (Normalization*G4UniformRand() > Probability);
    
  return KineticEnergy;
} 


G4ThreeVector G4GEMChannel::IsotropicVector(const G4double Magnitude)
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



