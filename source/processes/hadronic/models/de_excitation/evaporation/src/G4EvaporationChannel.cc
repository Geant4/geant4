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
// $Id: G4EvaporationChannel.cc 90273 2015-05-22 10:20:32Z gcosmo $
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
					   G4EvaporationProbability* aprob,
                                           G4VCoulombBarrier* barrier):
    G4VEvaporationChannel(aName),
    theA(anA),
    theZ(aZ),
    theProbability(aprob),
    theCoulombBarrier(barrier),
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
  theProbability->SetOPTxs(OPTxs);
  // for superimposed Coulomb Barrier for inverse cross sections
  theProbability->UseSICB(useSICB);

  G4VEvaporationChannel::Initialise();  
}

G4double G4EvaporationChannel::GetEmissionProbability(G4Fragment* fragment)
{
  G4int FragA = fragment->GetA_asInt();
  G4int FragZ = fragment->GetZ_asInt();
  ResidualA = FragA - theA;
  ResidualZ = FragZ - theZ;
  //G4cout << "G4EvaporationChannel::Initialize Z= " << theZ << " A= " << theA 
  //	 << " FragZ= " << FragZ << " FragA= " << FragA << G4endl;
  EmissionProbability = 0.0;

  // Only channels which are physically allowed are taken into account 
  if (ResidualA >= ResidualZ && ResidualZ > 0 && ResidualA >= theA) {
  
    //Effective excitation energy
    G4double ExEnergy = fragment->GetExcitationEnergy();
    G4double delta0 = 
      std::max(0.0,pairingCorrection->GetPairingCorrection(FragA,FragZ));
    G4double delta1 = 
      std::max(0.0,pairingCorrection->GetPairingCorrection(ResidualA,ResidualZ));
    ResidualMass = G4NucleiProperties::GetNuclearMass(ResidualA, ResidualZ);
    G4double FragmentMass = fragment->GetGroundStateMass();
    G4double Etot = FragmentMass + ExEnergy;
    G4double ResMass = ResidualMass + delta1;  

    if(ExEnergy >= delta0 && Etot > ResMass + EvaporatedMass) {
  
      // Maximal Kinetic Energy
      MaximalKineticEnergy = ((Etot-ResMass)*(Etot+ResMass) 
	    + EvaporatedMass*EvaporatedMass)/(2.0*Etot) - EvaporatedMass;

      // The threshold for charged particle emission must be set to 
      // 0 if Coulomb cutoff  is included in the cross sections
      // Of course for OPTxs=0 we have the Coulomb barrier 

      CoulombBarrier = 0.0;
      if (OPTxs==0 || useSICB) {
	CoulombBarrier = 
	  theCoulombBarrier->GetCoulombBarrier(ResidualA,ResidualZ,ExEnergy);
      }
      if (MaximalKineticEnergy > CoulombBarrier) {
	EmissionProbability = theProbability->
	  TotalProbability(*fragment, CoulombBarrier, MaximalKineticEnergy);
      }
    }
  }
  //G4cout << "G4EvaporationChannel:: probability= " 
  // << EmissionProbability << G4endl;   
  return EmissionProbability;
}

G4Fragment* G4EvaporationChannel::EmittedFragment(G4Fragment* theNucleus)
{
  G4Fragment* evFragment = 0;
  G4double evEnergy = EvaporatedMass +
    theProbability->SampleKineticEnergy(CoulombBarrier,
					MaximalKineticEnergy);

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
