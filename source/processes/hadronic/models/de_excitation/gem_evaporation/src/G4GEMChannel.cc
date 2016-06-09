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
// $Id: G4GEMChannel.cc,v 1.5 2006/06/29 20:22:07 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//

#include "G4GEMChannel.hh"
#include "G4PairingCorrection.hh"


G4GEMChannel::G4GEMChannel(const G4int theA, const G4int theZ,
                           G4GEMProbability * aEmissionStrategy,
                           G4VCoulombBarrier * aCoulombBarrier):
    A(theA),
    Z(theZ),
    theEvaporationProbabilityPtr(aEmissionStrategy),
    theCoulombBarrierPtr(aCoulombBarrier),
    EmissionProbability(0.0),
    MaximalKineticEnergy(-1000.0)
{ 
    theLevelDensityPtr = new G4EvaporationLevelDensityParameter;
    MyOwnLevelDensity = true;
}

G4GEMChannel::G4GEMChannel(const G4int theA, const G4int theZ, const G4String & aName,
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
}

G4GEMChannel::G4GEMChannel(const G4int theA, const G4int theZ, const G4String * aName,
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
}


G4GEMChannel::~G4GEMChannel()
{
    if (MyOwnLevelDensity) delete theLevelDensityPtr;
}




G4GEMChannel::G4GEMChannel(const G4GEMChannel & ) : G4VEvaporationChannel()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4GEMChannel::copy_costructor meant to not be accessable");
}

const G4GEMChannel & G4GEMChannel::operator=(const G4GEMChannel & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4GEMChannel::operator= meant to not be accessable");
    return *this;
}

G4bool G4GEMChannel::operator==(const G4GEMChannel & right) const 
{
    return (this == (G4GEMChannel *) &right);
    //  return false;
}

G4bool G4GEMChannel::operator!=(const G4GEMChannel & right) const 
{
    return (this != (G4GEMChannel *) &right);
    //  return true;
}



void G4GEMChannel::Initialize(const G4Fragment & fragment)
{

    G4int anA = static_cast<G4int>(fragment.GetA());
    G4int aZ = static_cast<G4int>(fragment.GetZ());
    AResidual = anA - A;
    ZResidual = aZ - Z;
    
    // Effective excitation energy
    G4double ExEnergy = fragment.GetExcitationEnergy() -
      G4PairingCorrection::GetInstance()->GetPairingCorrection(anA,aZ);
    
    // We only take into account channels which are physically allowed
    if (AResidual <= 0 || ZResidual <= 0 || AResidual < ZResidual ||
        (AResidual == ZResidual && AResidual > 1) || ExEnergy <= 0.0) 
      {
        CoulombBarrier = 0.0;
        MaximalKineticEnergy = -1000.0*MeV;
        EmissionProbability = 0.0;
      } 
    else 
      {
        // Coulomb Barrier calculation
        CoulombBarrier = theCoulombBarrierPtr->GetCoulombBarrier(AResidual,ZResidual,ExEnergy);

//  	std::cout << "\tfragment (" << A << ',' << Z << ") residual (" << AResidual << ',' << ZResidual << ')'; 
//  	std::cout << " U = " << fragment.GetExcitationEnergy();
//  	std::cout << " delta = " << G4PairingCorrection::GetInstance()->GetPairingCorrection(anA,aZ);
//  	std::cout << " U-delta = " << ExEnergy;
//  	std::cout << " V = " << CoulombBarrier;
        
        // Maximal Kinetic Energy
        MaximalKineticEnergy = CalcMaximalKineticEnergy(G4ParticleTable::GetParticleTable()->
                                                        GetIonTable()->GetNucleusMass(aZ,anA)+ExEnergy);

//  	std::cout << " Tmax-V = " << MaximalKineticEnergy << '\n';
		
        // Emission probability
        if (MaximalKineticEnergy <= 0.0) 
	  {
	    EmissionProbability = 0.0;
	  }
	else 
	  { 
	    // Total emission probability for this channel
	    EmissionProbability = theEvaporationProbabilityPtr->EmissionProbability(fragment,MaximalKineticEnergy);
	  }
      }
    
    return;
}


G4FragmentVector * G4GEMChannel::BreakUp(const G4Fragment & theNucleus)
{
    
    G4double EvaporatedKineticEnergy = CalcKineticEnergy(theNucleus);
    G4double EvaporatedMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z,A);
    G4double EvaporatedEnergy = EvaporatedKineticEnergy + EvaporatedMass;
    
    G4ThreeVector momentum(IsotropicVector(std::sqrt(EvaporatedKineticEnergy*
                                                (EvaporatedKineticEnergy+2.0*EvaporatedMass))));
    
    momentum.rotateUz(theNucleus.GetMomentum().vect().unit());

    G4LorentzVector EvaporatedMomentum(momentum,EvaporatedEnergy);
    EvaporatedMomentum.boost(theNucleus.GetMomentum().boostVector());
    
    G4Fragment * EvaporatedFragment = new G4Fragment(A,Z,EvaporatedMomentum);
#ifdef PRECOMPOUND_TEST
    EvaporatedFragment->SetCreatorModel(G4String("G4Evaporation"));
#endif
    // ** And now the residual nucleus ** 
    G4double theExEnergy = theNucleus.GetExcitationEnergy();
    G4double theMass = G4ParticleTable::GetParticleTable()->
        GetIonTable()->GetNucleusMass(static_cast<G4int>(theNucleus.GetZ()),
				      static_cast<G4int>(theNucleus.GetA()));
    G4double ResidualEnergy = theMass + (theExEnergy - EvaporatedKineticEnergy) - EvaporatedMass;
	
    G4LorentzVector ResidualMomentum(-momentum,ResidualEnergy);
    ResidualMomentum.boost(theNucleus.GetMomentum().boostVector());
	
    G4Fragment * ResidualFragment = new G4Fragment( AResidual, ZResidual, ResidualMomentum );

#ifdef PRECOMPOUND_TEST
    ResidualFragment->SetCreatorModel(G4String("ResidualNucleus"));
#endif
    
    
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
{
    G4double ResidualMass = G4ParticleTable::GetParticleTable()->
        GetIonTable()->GetNucleusMass( ZResidual, AResidual );
    G4double EvaporatedMass = G4ParticleTable::GetParticleTable()->
        GetIonTable()->GetNucleusMass( Z, A );
	
    G4double T = (NucleusTotalE*NucleusTotalE + EvaporatedMass*EvaporatedMass - ResidualMass*ResidualMass)/
        (2.0*NucleusTotalE) -
        EvaporatedMass - CoulombBarrier;
	
    return T;
}




G4double G4GEMChannel::CalcKineticEnergy(const G4Fragment & fragment)
    // Samples fragment kinetic energy.
{
    G4double U = fragment.GetExcitationEnergy();
	
    G4double NuclearMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetNucleusMass(Z,A);

    G4double a = theLevelDensityPtr->LevelDensityParameter(static_cast<G4int>(fragment.GetA()),
							   static_cast<G4int>(fragment.GetZ()),U);
    G4double delta0 = G4PairingCorrection::GetInstance()->GetPairingCorrection(static_cast<G4int>(fragment.GetA()),
									       static_cast<G4int>(fragment.GetZ()));
    
    G4double Alpha = theEvaporationProbabilityPtr->CalcAlphaParam(fragment);
    G4double Beta = theEvaporationProbabilityPtr->CalcBetaParam(fragment);
    G4double Normalization = theEvaporationProbabilityPtr->GetNormalization();
    
    G4double Ux = (2.5 + 150.0/AResidual)*MeV;
    G4double Ex = Ux + delta0;
    G4double T = 1.0/(std::sqrt(a/Ux) - 1.5/Ux);
    G4double E0 = Ex - T*(std::log(T*MeV) - std::log(a/MeV)/4.0 - 1.25*std::log(Ux*MeV) + 2.0*std::sqrt(a*Ux));

    G4double InitialLevelDensity;
    if ( U < Ex ) 
      {
        InitialLevelDensity = (pi/12.0)*std::exp((U-E0)/T)/T;
      } 
    else 
      {
        InitialLevelDensity = (pi/12.0)*std::exp(2*std::sqrt(a*(U-delta0)))/std::pow(a*std::pow(U-delta0,5.0),1.0/4.0);
      }

    const G4double Spin = theEvaporationProbabilityPtr->GetSpin();
    G4double g = (2.0*Spin+1.0)*NuclearMass/(pi2* hbar_Planck*hbar_Planck);
    G4double RN = 0.0;
    if (A > 4) 
    {
	G4double R1 = std::pow(G4double(fragment.GetA()-A),1.0/3.0);
	G4double R2 = std::pow(G4double(A),1.0/3.0);
	RN = 1.12*(R1 + R2) - 0.86*((R1+R2)/(R1*R2));
	RN *= fermi;
    }
    else 
    {
	RN = 1.5*fermi;
    }
    G4double GeometricalXS = pi*RN*RN*std::pow(G4double(fragment.GetA()-A),2./3.); 
    

    G4double ConstantFactor = g*GeometricalXS*Alpha/InitialLevelDensity;
    ConstantFactor *= pi/12.0;
    G4double theEnergy = MaximalKineticEnergy + CoulombBarrier;
    G4double KineticEnergy;
    G4double Probability;

//      std::cout << "\t\tEjectile (" << A << ',' << Z << ") V = " << CoulombBarrier 
//  	      << " Beta = " << Beta << " V+Beta = " << CoulombBarrier+Beta << '\n';

    do 
      {
        KineticEnergy =  CoulombBarrier + G4UniformRand()*MaximalKineticEnergy;
        Probability = ConstantFactor*(KineticEnergy + Beta);
        if ( theEnergy-KineticEnergy < Ex) 
	  {
            Probability *= std::exp((theEnergy-KineticEnergy-E0)/T)/T;
	  } 
	else 
	  {
            Probability *= std::exp(2*std::sqrt(a*(theEnergy-KineticEnergy-delta0)))/
	      std::pow(a*std::pow(theEnergy-KineticEnergy-delta0,5.0),1.0/4.0);
	  }
        Probability /= Normalization;
    }
    while (G4UniformRand() > Probability);


    //  ---------------------------------------------------------------------------------------------------
    
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



