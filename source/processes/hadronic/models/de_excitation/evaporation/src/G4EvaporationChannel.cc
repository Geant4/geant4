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
// $Id: G4EvaporationChannel.cc 103162 2017-03-20 09:40:58Z gcosmo $
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
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4Alpha.hh"

G4EvaporationChannel::G4EvaporationChannel(G4int anA, G4int aZ, 
					   const G4String & aName,
					   G4EvaporationProbability* aprob,
                                           G4VCoulombBarrier* barrier):
    G4VEvaporationChannel(aName),
    theA(anA),
    theZ(aZ),
    theProbability(aprob),
    theCoulombBarrier(barrier)
{ 
  ResA = ResZ = 0;
  Mass = CoulombBarrier = MinKinEnergy = MaxKinEnergy = EmissionProbability = 0.0; 
  EvapMass = G4NucleiProperties::GetNuclearMass(theA, theZ);
  pairingCorrection = G4PairingCorrection::GetInstance();
}

G4EvaporationChannel::~G4EvaporationChannel()
{}

void G4EvaporationChannel::Initialise()
{
  theProbability->Initialise();
  G4VEvaporationChannel::Initialise();  
}

G4double G4EvaporationChannel::GetEmissionProbability(G4Fragment* fragment)
{
  G4int FragA = fragment->GetA_asInt();
  G4int FragZ = fragment->GetZ_asInt();
  ResA = FragA - theA;
  ResZ = FragZ - theZ;

  G4double FragmentMass = fragment->GetGroundStateMass();
  G4double ExEnergy = fragment->GetExcitationEnergy();
  Mass = FragmentMass + ExEnergy;
  //G4cout << "G4EvaporationChannel::Initialize Z= " << theZ << " A= " << theA 
  //	 << " FragZ= " << FragZ << " FragA= " << FragA << G4endl;
  EmissionProbability = 0.0;

  // Only channels which are physically allowed are taken into account 
  if (ResA >= ResZ && ResZ > 0 && ResA >= theA) {
  
    //Effective excitation energy
    G4double ResMass = G4NucleiProperties::GetNuclearMass(ResA, ResZ);

    CoulombBarrier = (0 == theZ) ? 0.0 : 
      theCoulombBarrier->GetCoulombBarrier(ResA,ResZ,ExEnergy);

    G4double delta0 = 
      std::max(0.0,pairingCorrection->GetPairingCorrection(FragA,FragZ));
    G4double delta1 = 
      std::max(0.0,pairingCorrection->GetPairingCorrection(ResA,ResZ));
    ResMass += delta1;  
    /*
    G4cout << "ExEnergy= " << ExEnergy << " Ec= " << CoulombBarrier
	   << " delta0= " << delta0 << " delta1= " << delta1
	   << " Free= " << Mass - ResMass - EvapMass 
	   << G4endl;
    */
    // for OPTxs >0 penetration under the barrier is taken into account
    // G4double elim = (0 == OPTxs) ? CoulombBarrier : CoulombBarrier*0.5;
    static const G4double dCB = 3.5*CLHEP::MeV;
    G4double elim = (0 == OPTxs) ? CoulombBarrier : CoulombBarrier - dCB*theZ;
    if(ExEnergy >= delta0 && Mass >= ResMass + EvapMass + elim) {
      G4double xm2 = (Mass - EvapMass)*(Mass - EvapMass);
      G4double xm  = Mass - EvapMass - elim;
      MinKinEnergy = (0.0 >= elim) ? 0.0 : std::max(0.5*(xm2 - xm*xm)/Mass, 0.0);
      MaxKinEnergy = std::max(0.5*(xm2 - ResMass*ResMass)/Mass, 0.0);
      //G4cout << "Emin= " << MinKinEnergy << " Emax= " << MaxKinEnergy 
      //     << "  xm= " << xm  << G4endl;
      EmissionProbability = theProbability->
	TotalProbability(*fragment, MinKinEnergy, MaxKinEnergy, CoulombBarrier);
    }
  }
  //G4cout << "G4EvaporationChannel:: probability= " 
  //    << EmissionProbability << G4endl;   
  return EmissionProbability;
}

G4Fragment* G4EvaporationChannel::EmittedFragment(G4Fragment* theNucleus)
{
  G4Fragment* evFragment = nullptr;
  G4double ekin = 0.0;
  if(ResA <= 4 && 
    ((ResA == 4 && ResZ == 2) || (ResA == 3 && ResZ == 2) ||
     (ResA == 3 && ResZ == 1) || (ResA == 2 && ResZ == 1) ||
     (ResA == 1 && ResZ == 1) || (ResA == 1 && ResZ == 0) )) {
    G4double mres = G4NucleiProperties::GetNuclearMass(ResA, ResZ);
    ekin = 0.5*(Mass*Mass - mres*mres + EvapMass*EvapMass)/Mass - EvapMass;
  } else {
    ekin = theProbability->SampleKineticEnergy(MinKinEnergy, MaxKinEnergy,
					       CoulombBarrier);
  }
  G4LorentzVector lv0 = theNucleus->GetMomentum();
  G4LorentzVector lv(std::sqrt(ekin*(ekin + 2.0*EvapMass))*G4RandomDirection(), 
		     ekin + EvapMass);
  lv.boost(lv0.boostVector());

  evFragment = new G4Fragment(theA, theZ, lv);
  lv0 -= lv;
  theNucleus->SetZandA_asInt(ResZ, ResA);
  theNucleus->SetMomentum(lv0);

  return evFragment; 
} 
