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
//JMQ & MAC 7/12/07 : new MonteCarlo sampling  of the kinetic energy
//
// $Id: G4EvaporationChannel.cc,v 1.9 2008-05-24 16:34:33 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
//J.M. Quesada (Apr.2008). Rebuilt class. Mayor changes. Unused items have been removed (constructors..). New methods

#include "G4EvaporationChannel.hh"
#include "G4PairingCorrection.hh"



G4EvaporationChannel::G4EvaporationChannel(const G4int theA, const G4int theZ, const G4String & aName,
					   G4VEmissionProbability * aEmissionStrategy,
                                           G4VCoulombBarrier * aCoulombBarrier):
    G4VEvaporationChannel(aName),
    A(theA),
    Z(theZ),
    theEvaporationProbabilityPtr(aEmissionStrategy),
    theCoulombBarrierPtr(aCoulombBarrier),
    EmissionProbability(0.0),
    MaximalKineticEnergy(-1000.0)
{ 
    
}

G4EvaporationChannel::~G4EvaporationChannel()
{
  
}

G4EvaporationChannel::G4EvaporationChannel(const G4EvaporationChannel & ) : G4VEvaporationChannel()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4EvaporationChannel::copy_costructor meant to not be accessable");
}

const G4EvaporationChannel & G4EvaporationChannel::operator=(const G4EvaporationChannel & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4EvaporationChannel::operator= meant to not be accessable");
    return *this;
}

G4bool G4EvaporationChannel::operator==(const G4EvaporationChannel & right) const 
{
    return (this == (G4EvaporationChannel *) &right);
    //  return false;
}

G4bool G4EvaporationChannel::operator!=(const G4EvaporationChannel & right) const 
{
    return (this != (G4EvaporationChannel *) &right);
    //  return true;
}



void G4EvaporationChannel::Initialize(const G4Fragment & fragment)
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
	(AResidual == ZResidual && AResidual > 1) || ExEnergy <= 0.0) {
	MaximalKineticEnergy = -1000.0*MeV;
	EmissionProbability = 0.0;
    } else {
	// Coulomb Barrier calculation
	CoulombBarrier = theCoulombBarrierPtr->GetCoulombBarrier(AResidual,ZResidual,ExEnergy);
	// Maximal Kinetic Energy
	MaximalKineticEnergy = CalcMaximalKineticEnergy(G4ParticleTable::GetParticleTable()->
							GetIonTable()->GetNucleusMass(aZ,anA)+ExEnergy);		
	// Emission probability
	if (MaximalKineticEnergy <= 0.0) EmissionProbability = 0.0;
	else { 
	    // Total emission probability for this channel
	    EmissionProbability = theEvaporationProbabilityPtr->EmissionProbability(fragment, MaximalKineticEnergy);
	}
    }

//JMQ & MAC (14/12/07): protection against non-physical situation
if(EmissionProbability>0 && MaximalKineticEnergy<0)
{
	std::ostringstream errOs;
  	errOs << "Non-physical situation: EmissionProbability >0 & Tmax <0" <<G4endl;
        throw G4HadronicException(__FILE__, __LINE__, errOs.str());
}
    return;
}

G4FragmentVector * G4EvaporationChannel::BreakUp(const G4Fragment & theNucleus)
{
// JMQ: 7/12/07  new MonteCarlo sampling method of kinetic energy
      G4double EvaporatedKineticEnergy=GetKineticEnergy(theNucleus);

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
    G4double theMass = G4ParticleTable::GetParticleTable()->GetIonTable()->
	GetNucleusMass(static_cast<G4int>(theNucleus.GetZ()),
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



G4double G4EvaporationChannel::CalcMaximalKineticEnergy(const G4double NucleusTotalE)
    // Calculate maximal kinetic energy that can be carried by fragment.
{
    G4double ResidualMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetNucleusMass( ZResidual, AResidual );
    G4double EvaporatedMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetNucleusMass( Z, A );
	
    G4double Tmax = (NucleusTotalE*NucleusTotalE + EvaporatedMass*EvaporatedMass - 		
	ResidualMass*ResidualMass)/(2.0*NucleusTotalE) - EvaporatedMass;

    return Tmax;
}



G4double G4EvaporationChannel::GetKineticEnergy(const G4Fragment & aFragment)
//JMQ & MAC : 7/12/07 new method for MC sampling. Substitutes old CalcKineticEnergy
{
 G4double V=CoulombBarrier;
 G4double Tmax=MaximalKineticEnergy;
  G4double T(0.0);
  G4double NormalizedProbability(1.0);
  G4double Gamma;
  if (A==1) {Gamma=2;}
  else if (A==2 && Z==1) {Gamma=1;}
  else if (A==3 ) {Gamma=2;}
  else if (A==4 && Z==2){Gamma=1;}
  else {
 	std::ostringstream errOs;
      errOs << "ejected particle out of range"  << G4endl;
    throw G4HadronicException(__FILE__, __LINE__, errOs.str());
    }

// JMQ & MAC. new pointer to object is created in order to access the distribution function.

G4EvaporationProbability * G4EPtemp=new G4EvaporationProbability(A,Z,Gamma,theCoulombBarrierPtr);

  do
     {  
      T=V+G4UniformRand()*(Tmax-V);
       NormalizedProbability=G4EPtemp->ProbabilityDistributionFunction(T,aFragment)/
(this->GetEmissionProbability());

}
   while (G4UniformRand() > NormalizedProbability);
   delete G4EPtemp;
   return T;
}

 

G4ThreeVector G4EvaporationChannel::IsotropicVector(const G4double Magnitude)
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



   


  

