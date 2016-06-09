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
// $Id: G4EvaporationChannel.cc,v 1.19 2010-11-24 15:30:49 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4EvaporationLevelDensityParameter.hh"
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
  ResidualMass = CoulombBarrier=0.0;
  EvaporatedMass = G4NucleiProperties::GetNuclearMass(theA, theZ);
  theLevelDensityPtr = new G4EvaporationLevelDensityParameter;
}

G4EvaporationChannel::G4EvaporationChannel():
    G4VEvaporationChannel(""),
    theA(0),
    theZ(0),
    theEvaporationProbabilityPtr(0),
    theCoulombBarrierPtr(0),
    EmissionProbability(0.0),
    MaximalKineticEnergy(-1000.0)
{ 
  ResidualA = 0;
  ResidualZ = 0;
  EvaporatedMass = ResidualMass = CoulombBarrier = 0.0;
  theLevelDensityPtr = new G4EvaporationLevelDensityParameter;
}

G4EvaporationChannel::~G4EvaporationChannel()
{
  delete theLevelDensityPtr;
}

void G4EvaporationChannel::Initialize(const G4Fragment & fragment)
{
  //for inverse cross section choice
  theEvaporationProbabilityPtr->SetOPTxs(OPTxs);
  // for superimposed Coulomb Barrier for inverse cross sections
  theEvaporationProbabilityPtr->UseSICB(useSICB);
  
  G4int FragmentA = fragment.GetA_asInt();
  G4int FragmentZ = fragment.GetZ_asInt();
  ResidualA = FragmentA - theA;
  ResidualZ = FragmentZ - theZ;
  //G4cout << "G4EvaporationChannel::Initialize Z= " << theZ << " A= " << theA 
  //	 << " FragZ= " << FragmentZ << " FragA= " << FragmentA << G4endl;
  
  //Effective excitation energy
  G4double ExEnergy = fragment.GetExcitationEnergy() - 
    G4PairingCorrection::GetInstance()->GetPairingCorrection(FragmentA,FragmentZ);
  
  // Only channels which are physically allowed are taken into account 
  if (ResidualA <= 0 || ResidualZ <= 0 || ResidualA < ResidualZ ||
      (ResidualA == ResidualZ && ResidualA > 1) || ExEnergy <= 0.0) {
    CoulombBarrier = ResidualMass = 0.0;
    MaximalKineticEnergy = -1000.0*MeV;
    EmissionProbability = 0.0;
  } else {
    ResidualMass = G4NucleiProperties::GetNuclearMass(ResidualA, ResidualZ);
    G4double FragmentMass = fragment.GetGroundStateMass();
    CoulombBarrier = theCoulombBarrierPtr->GetCoulombBarrier(ResidualA,ResidualZ,ExEnergy);
    // Maximal Kinetic Energy
    MaximalKineticEnergy = CalcMaximalKineticEnergy(FragmentMass + ExEnergy);
    //MaximalKineticEnergy = ExEnergy + fragment.GetGroundStateMass() 
    //  - EvaporatedMass - ResidualMass;
    
    // Emission probability
    // Protection for the case Tmax<V. If not set in this way we could end up in an 
    // infinite loop in  the method GetKineticEnergy if OPTxs!=0 && useSICB=true. 
    // Of course for OPTxs=0 we have the Coulomb barrier 
    
    G4double limit;
    if (OPTxs==0 || (OPTxs!=0 && useSICB)) 
      limit= CoulombBarrier;
    else limit=0.;
    limit = std::max(limit, FragmentMass - ResidualMass - EvaporatedMass);
  
    // The threshold for charged particle emission must be  set to 0 if Coulomb 
    //cutoff  is included in the cross sections
    if (MaximalKineticEnergy <= limit) EmissionProbability = 0.0;
    else { 
      // Total emission probability for this channel
      EmissionProbability = theEvaporationProbabilityPtr->
        EmissionProbability(fragment, MaximalKineticEnergy);
    }
  }
  //G4cout << "G4EvaporationChannel:: probability= " << EmissionProbability << G4endl; 
  
  return;
}

G4FragmentVector * G4EvaporationChannel::BreakUp(const G4Fragment & theNucleus)
{
  /*
  G4double Ecm = GetKineticEnergy(theNucleus) + ResidualMass + EvaporatedMass;
  
  G4double EvaporatedEnergy = 
    ((Ecm-ResidualMass)*(Ecm+ResidualMass) + EvaporatedMass*EvaporatedMass)/(2*Ecm);
  */
  G4double EvaporatedEnergy = GetKineticEnergy(theNucleus) + EvaporatedMass;

  G4ThreeVector momentum(IsotropicVector
                         (std::sqrt((EvaporatedEnergy - EvaporatedMass)*
                                    (EvaporatedEnergy + EvaporatedMass))));
  
  G4LorentzVector EvaporatedMomentum(momentum,EvaporatedEnergy);
  G4LorentzVector ResidualMomentum = theNucleus.GetMomentum();
  EvaporatedMomentum.boost(ResidualMomentum.boostVector());
  
  G4Fragment * EvaporatedFragment = new G4Fragment(theA,theZ,EvaporatedMomentum);
  ResidualMomentum -= EvaporatedMomentum;

  G4Fragment * ResidualFragment = new G4Fragment(ResidualA, ResidualZ, ResidualMomentum);
 
  G4FragmentVector * theResult = new G4FragmentVector;
  
#ifdef debug
  G4double Efinal = ResidualMomentum.e() + EvaporatedMomentum.e();
  G4ThreeVector Pfinal = ResidualMomentum.vect() + EvaporatedMomentum.vect();
  if (std::abs(Efinal-theNucleus.GetMomentum().e()) > 1.0*keV) {
    G4cout << "@@@@@@@@@@@@@@@@@@@@@ G4Evaporation Chanel: ENERGY @@@@@@@@@@@@@@@@" << G4endl;
    G4cout << "Initial : " << theNucleus.GetMomentum().e()/MeV << " MeV    Final :" 
           <<Efinal/MeV << " MeV    Delta: " <<  (Efinal-theNucleus.GetMomentum().e())/MeV 
           << " MeV" << G4endl;
  }
  if (std::abs(Pfinal.x()-theNucleus.GetMomentum().x()) > 1.0*keV ||
      std::abs(Pfinal.y()-theNucleus.GetMomentum().y()) > 1.0*keV ||
      std::abs(Pfinal.z()-theNucleus.GetMomentum().z()) > 1.0*keV ) {
    G4cout << "@@@@@@@@@@@@@@@@@@@@@ G4Evaporation Chanel: MOMENTUM @@@@@@@@@@@@@@@@" << G4endl;
    G4cout << "Initial : " << theNucleus.GetMomentum().vect() << " MeV    Final :" 
           <<Pfinal/MeV << " MeV    Delta: " <<  Pfinal-theNucleus.GetMomentum().vect()
           << " MeV" << G4endl;   
    
  }
#endif
  theResult->push_back(EvaporatedFragment);
  theResult->push_back(ResidualFragment);
  return theResult; 
} 

/////////////////////////////////////////
// Calculates the maximal kinetic energy that can be carried by fragment.
G4double G4EvaporationChannel::CalcMaximalKineticEnergy(G4double NucleusTotalE)
{
  // This is the "true" assimptotic kinetic energy (from energy conservation)	
  G4double Tmax = 
    ((NucleusTotalE-ResidualMass)*(NucleusTotalE+ResidualMass) + EvaporatedMass*EvaporatedMass)
    /(2.0*NucleusTotalE) - EvaporatedMass;
  
  //JMQ (13-09-08) bug fixed: in the original version the Tmax is calculated
  //at the Coulomb barrier
  //IMPORTANT: meaning of Tmax differs in OPTxs=0 and OPTxs!=0
  //When OPTxs!=0 Tmax is the TRUE (assimptotic) maximal kinetic energy
  
  if(OPTxs==0) 
    Tmax=Tmax- CoulombBarrier;
  
  return Tmax;
}

///////////////////////////////////////////
//JMQ: New method for MC sampling of kinetic energy. Substitutes old CalcKineticEnergy
G4double G4EvaporationChannel::GetKineticEnergy(const G4Fragment & aFragment)
{
  
  if (OPTxs==0) {
    // It uses Dostrovsky's approximation for the inverse reaction cross
    // in the probability for fragment emission
    // MaximalKineticEnergy energy in the original version (V.Lara) was calculated at 
    //the Coulomb barrier.
    
    
    if (MaximalKineticEnergy < 0.0) {
      throw G4HadronicException(__FILE__, __LINE__, 
                                "G4EvaporationChannel::CalcKineticEnergy: maximal kinetic at the Coulomb barrier is less than 0");
    }
    G4double Rb = 4.0*theLevelDensityPtr->
      LevelDensityParameter(ResidualA+theA,ResidualZ+theZ,MaximalKineticEnergy)*
      MaximalKineticEnergy;
    G4double RbSqrt = std::sqrt(Rb);
    G4double PEX1 = 0.0;
    if (RbSqrt < 160.0) PEX1 = std::exp(-RbSqrt);
    G4double Rk = 0.0;
    G4double FRk = 0.0;
    do {
      G4double RandNumber = G4UniformRand();
      Rk = 1.0 + (1./RbSqrt)*std::log(RandNumber + (1.0-RandNumber)*PEX1);
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
    
    G4double result =  MaximalKineticEnergy * (1.0-Rk*Rk) + CoulombBarrier;
    return result;
  } else if (OPTxs==1 || OPTxs==2 || OPTxs==3 || OPTxs==4) {
    
    // Coulomb barrier is just included  in the cross sections
    G4double V = 0;
    if(useSICB) { V= CoulombBarrier; }

    V = std::max(V, aFragment.GetGroundStateMass()-EvaporatedMass-ResidualMass);

    G4double Tmax=MaximalKineticEnergy;
    G4double T(0.0);
    G4double NormalizedProbability(1.0);

    // VI: This is very ineffective - create new objects at each call to the method    
    /*
    // A pointer is created in order to access the distribution function.
    G4EvaporationProbability * G4EPtemp = 0;
    
    if (theA==1 && theZ==0) G4EPtemp=new G4NeutronEvaporationProbability();
    else if (theA==1 && theZ==1) G4EPtemp=new G4ProtonEvaporationProbability();
    else if (theA==2 && theZ==1 ) G4EPtemp=new G4DeuteronEvaporationProbability();
    else if (theA==3 && theZ==1 ) G4EPtemp=new G4TritonEvaporationProbability();
    else if (theA==3 && theZ==2 ) G4EPtemp=new G4He3EvaporationProbability();
    else if (theA==4 && theZ==2) G4EPtemp=new G4AlphaEvaporationProbability(); 
    else {
      std::ostringstream errOs;
      errOs << "ejected particle out of range in G4EvaporationChannel"  << G4endl;
      throw G4HadronicException(__FILE__, __LINE__, errOs.str());
    }

      //for cross section selection and superimposed Coulom Barrier for xs
      G4EPtemp->SetOPTxs(OPTxs);
      G4EPtemp->UseSICB(useSICB);
    */

    // use local pointer and not create a new one
    do
      {  
        T=V+G4UniformRand()*(Tmax-V);
        NormalizedProbability = 
	  theEvaporationProbabilityPtr->ProbabilityDistributionFunction(aFragment,T)/
          GetEmissionProbability();
        
      }
    while (G4UniformRand() > NormalizedProbability);
    //   delete G4EPtemp;
    return T;
  } else{
    std::ostringstream errOs;
    errOs << "Bad option for energy sampling in evaporation"  <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, errOs.str());
  }
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





   


  

