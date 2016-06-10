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
// $Id: G4MuMinusCapturePrecompound.cc 91836 2015-08-07 07:25:54Z gcosmo $
//
//-----------------------------------------------------------------------------
//
// GEANT4 Class file 
//
// File name:  G4MuMinusCapturePrecompound
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 22 April 2012 on base of G4MuMinusCaptureCascade
//
//
//-----------------------------------------------------------------------------
//
// Modifications: 
//
//-----------------------------------------------------------------------------

#include "G4MuMinusCapturePrecompound.hh"
#include "Randomize.hh" 
#include "G4RandomDirection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4MuonMinus.hh"
#include "G4NeutrinoMu.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Triton.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4NucleiProperties.hh"
#include "G4VPreCompoundModel.hh"
#include "G4PreCompoundModel.hh"
#include "G4HadronicInteractionRegistry.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuMinusCapturePrecompound::G4MuMinusCapturePrecompound(
    G4VPreCompoundModel* ptr)
  : G4HadronicInteraction("muMinusNuclearCapture")
{ 
  fMuMass = G4MuonMinus::MuonMinus()->GetPDGMass(); 
  fProton = G4Proton::Proton();
  fNeutron = G4Neutron::Neutron();
  fThreshold = 10*MeV;
  fTime = 0.0;
  fPreCompound = ptr;
  if(!ptr) { 
    G4HadronicInteraction* p =
      G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
    ptr = static_cast<G4VPreCompoundModel*>(p); 
    fPreCompound = ptr;
    if(!ptr) { fPreCompound = new G4PreCompoundModel(); }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuMinusCapturePrecompound::~G4MuMinusCapturePrecompound()
{
  result.Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4HadFinalState* 
G4MuMinusCapturePrecompound::ApplyYourself(const G4HadProjectile& projectile, 
					   G4Nucleus& targetNucleus)
{
  result.Clear();
  result.SetStatusChange(stopAndKill);
  fTime = projectile.GetGlobalTime();
  G4double time0 = fTime;

  G4double muBindingEnergy = projectile.GetBoundEnergy();

  G4int Z = targetNucleus.GetZ_asInt(); 
  G4int A = targetNucleus.GetA_asInt();
  G4double massA = G4NucleiProperties::GetNuclearMass(A, Z);

  /*
  G4cout << "G4MuMinusCapturePrecompound::ApplyYourself: Emu= "
	 << muBindingEnergy << G4endl;
  */
  // Energy on K-shell
  G4double muEnergy = fMuMass + muBindingEnergy;
  G4double muMom =std::sqrt(muBindingEnergy*(muBindingEnergy + 2.0*fMuMass));
  G4double availableEnergy = massA + fMuMass - muBindingEnergy;
  G4double residualMass = G4NucleiProperties::GetNuclearMass(A, Z - 1);

  G4ThreeVector vmu = muMom*G4RandomDirection();
  G4LorentzVector aMuMom(vmu, muEnergy);

  const G4double nenergy = keV;    

  // p or 3He as a target 
  // two body reaction mu- + A(Z,A) -> nuMu + A(Z-1,A)
  if((1 == Z && 1 == A) || (2 == Z && 3 == A)) {

    const G4ParticleDefinition* pd = 0;
    if(1 == Z) { pd = fNeutron; }
    else { pd = G4Triton::Triton(); }

    //
    //  Computation in assumption of CM reaction
    //  
    G4double e = 0.5*(availableEnergy - 
		      residualMass*residualMass/availableEnergy);

    G4ThreeVector nudir = G4RandomDirection();
    AddNewParticle(G4NeutrinoMu::NeutrinoMu(), nudir, e);
    nudir *= -1.0;
    AddNewParticle(pd, nudir, availableEnergy - e - residualMass);

  // d or 4He as a target 
  // three body reaction mu- + A(Z,A) -> nuMu + n + A(Z-1,A)
  // extra neutron produced at rest
  } else if((1 == Z && 2 == A) || (2 == Z && 4 == A)) {

    const G4ParticleDefinition* pd = 0;
    if(1 == Z) { pd = fNeutron; }
    else { pd = G4Triton::Triton(); }

    availableEnergy -= neutron_mass_c2 - nenergy;
    residualMass = pd->GetPDGMass();

    //
    //  Computation in assumption of CM reaction
    //  
    G4double e = 0.5*(availableEnergy - 
		      residualMass*residualMass/availableEnergy);

    G4ThreeVector nudir = G4RandomDirection();
    AddNewParticle(G4NeutrinoMu::NeutrinoMu(), nudir, e);
    nudir *= -1.0;
    AddNewParticle(pd, nudir, availableEnergy - e - residualMass);

    // extra low-energy neutron
    nudir = G4RandomDirection();
    AddNewParticle(fNeutron, nudir, nenergy);

  } else {
    // sample mu- + p -> nuMu + n reaction in CM of muonic atom

    // nucleus
    G4LorentzVector momInitial(0.0,0.0,0.0,availableEnergy);
    G4LorentzVector momResidual, momNu;

    // pick random proton inside nucleus 
    G4double eEx;
    fNucleus.Init(A, Z);
    const std::vector<G4Nucleon>& nucleons= fNucleus.GetNucleons();
    const G4ParticleDefinition* pDef;

    G4int reentryCount = 0;
  
    do {
      ++reentryCount;
      G4int index = 0;
      do {
	index=G4int(A*G4UniformRand());
	pDef = nucleons[index].GetDefinition();
      } while(pDef != fProton);
      G4LorentzVector momP = nucleons[index].Get4Momentum();

      // Get CMS kinematics
      G4LorentzVector theCMS = momP + aMuMom;
      G4ThreeVector bst = theCMS.boostVector();

      G4double Ecms = theCMS.mag();
      G4double Enu  = 0.5*(Ecms - neutron_mass_c2*neutron_mass_c2/Ecms);
      eEx = 0.0;

      if(Enu > 0.0) {
	// make the nu, and transform to lab;
	momNu.set(Enu*G4RandomDirection(), Enu);

	// nu in lab.
	momNu.boost(bst);
	momResidual = momInitial - momNu;
	eEx = momResidual.mag() - residualMass;
        if(eEx < 0.0 && eEx + nenergy >= 0.0) {
          momResidual.set(0.0, 0.0, 0.0, residualMass);
          eEx = 0.0;
	}
      }
      // in the case of many iterations stop the loop
      // with zero excitation energy
      if(reentryCount > 100 && eEx < 0.0) {
	G4ExceptionDescription ed;
	ed << "Call for " << GetModelName() << G4endl;
	ed << "Target  Z= " << Z  
	   << "  A= " << A << "  Eex(MeV)= " << eEx/MeV << G4endl;
	ed << " ApplyYourself does not completed after 100 attempts -"
	   << " excitation energy is set to zero";
	G4Exception("G4MuMinusCapturePrecompound::ApplyYourself", "had006", 
		    JustWarning, ed);
	momResidual.set(0.0, 0.0, 0.0, residualMass);
	eEx = 0.0;
      }
      // Loop checking, 06-Aug-2015, Vladimir Ivanchenko
    } while(eEx <= 0.0);

    G4ThreeVector dir = momNu.vect().unit();
    AddNewParticle(G4NeutrinoMu::NeutrinoMu(), dir, momNu.e());

    G4Fragment initialState(A, Z-1, momResidual);
    initialState.SetNumberOfExcitedParticle(2,0);
    initialState.SetNumberOfHoles(1,1);

    // decay time for pre-compound/de-excitation starts from zero
    G4ReactionProductVector* rpv = fPreCompound->DeExcite(initialState);
    size_t n = rpv->size();
    for(size_t i=0; i<n; ++i) {
      G4ReactionProduct* rp = (*rpv)[i];

      // reaction time
      fTime = time0 + rp->GetTOF();
      G4ThreeVector direction = rp->GetMomentum().unit();
      AddNewParticle(rp->GetDefinition(), direction, rp->GetKineticEnergy());
      delete rp;
    }
    delete rpv;
  } 
  if(verboseLevel > 1)
    G4cout << "G4MuMinusCapturePrecompound::ApplyYourself:  Nsec= " 
	   << result.GetNumberOfSecondaries() 
	   <<" E0(MeV)= " <<availableEnergy/MeV
	   <<" Mres(GeV)= " <<residualMass/GeV
	   <<G4endl;

  return &result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuMinusCapturePrecompound::ModelDescription(std::ostream& outFile) const
{
  outFile << "Sampling of mu- capture by atomic nucleus from K-shell"
	  << " mesoatom orbit.\n"
	  << "Primary reaction mu- + p -> n + neutrino, neutron providing\n"
	  << "  initial excitation of the target nucleus and PreCompound"
	  << " model samples final state\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
