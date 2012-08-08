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
// $Id: G4HadronElasticProcess.cc,v 1.5 2010-11-19 18:50:03 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Geant4 Hadron Elastic Scattering Process 
// 
// Created 26 July 2012 V.Ivanchenko from G4WHadronElasticProcess
//  
// Modified:
//

#include "G4HadronElasticProcess.hh"
#include "G4Nucleus.hh"
#include "G4ProcessManager.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4HadronElasticDataSet.hh"
#include "G4ProductionCutsTable.hh"
#include "G4HadronicException.hh"
#include <iostream>
#include <typeinfo>

G4HadronElasticProcess::G4HadronElasticProcess(const G4String& pName)
  : G4HadronicProcess(pName) 
{
  SetProcessSubType(fHadronElastic);
  AddDataSet(new G4HadronElasticDataSet);
  theNeutron  = G4Neutron::Neutron();
  lowestEnergy = 1.*keV;
  lowestEnergyNeutron = 1.e-6*eV;
}

G4HadronElasticProcess::~G4HadronElasticProcess()
{}

void G4HadronElasticProcess::Description() const
{
  char* dirName = getenv("G4PhysListDocDir");
  if (dirName) {
    std::ofstream outFile;
    G4String outFileName = GetProcessName() + ".html";
    G4String pathName = G4String(dirName) + "/" + outFileName;
    outFile.open(pathName);
    outFile << "<html>\n";
    outFile << "<head>\n";

    outFile << "<title>Description of G4HadronElasticProcess</title>\n";
    outFile << "</head>\n";
    outFile << "<body>\n";

    outFile << "G4HadronElasticProcess handles the elastic scattering of\n"
            << "hadrons by invoking one or more hadronic models and one or\n"
            << "more hadronic cross sections.\n";

    outFile << "</body>\n";
    outFile << "</html>\n";
    outFile.close();
  }
}

G4VParticleChange* G4HadronElasticProcess::PostStepDoIt(
				  const G4Track& track, 
				  const G4Step& step)
{
  aParticleChange.Initialize(track);
  G4double weight = track.GetWeight();
  aParticleChange.ProposeWeight(weight);

  G4double kineticEnergy = track.GetKineticEnergy();
  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();
  const G4ParticleDefinition* part = dynParticle->GetDefinition();

  // protection against numerical problems
  // primary energy cannot be very low
  if(part == theNeutron) {
    if(kineticEnergy <= lowestEnergyNeutron) 
      { return G4VDiscreteProcess::PostStepDoIt(track,step); }
  } else if(kineticEnergy <= lowestEnergy)
      { return G4VDiscreteProcess::PostStepDoIt(track,step); }

  G4Material* material = track.GetMaterial();
  G4Nucleus* targNucleus = GetTargetNucleusPointer();

  // Select element
  G4Element* elm = 0;
  try
    {
      elm = GetCrossSectionDataStore()->SampleZandA(dynParticle, material, 
						    *targNucleus);
    }
  catch(G4HadronicException & aR)
    {
      G4ExceptionDescription ed;
      DumpState(track,"SampleZandA",ed); 
      ed << " PostStepDoIt failed on element selection" << G4endl;
      G4Exception("G4HadronElasticProcess::PostStepDoIt", "had003", 
		  FatalException, ed);
    }
  G4HadronicInteraction* hadi = 0;
  try
    {
      hadi = ChooseHadronicInteraction( kineticEnergy, material, elm );
    }
  catch(G4HadronicException & aE)
    {
      G4ExceptionDescription ed;
      ed << "Target element "<< elm->GetName()<<"  Z= " 
	 << targNucleus->GetZ_asInt() << "  A= " 
	 << targNucleus->GetA_asInt() << G4endl;
      DumpState(track,"ChooseHadronicInteraction",ed);
      ed << " No HadronicInteraction found out" << G4endl;
      G4Exception("G4HadronElasticProcess::PostStepDoIt", "had005", 
		  FatalException, ed);
    }

  size_t idx = track.GetMaterialCutsCouple()->GetIndex();
  G4double tcut = (*(G4ProductionCutsTable::GetProductionCutsTable()
		     ->GetEnergyCutsVector(3)))[idx];
  hadi->SetRecoilEnergyThreshold(tcut);

  // Initialize the hadronic projectile from the track
  //  G4cout << "track " << track.GetDynamicParticle()->Get4Momentum()<<G4endl;
  G4HadProjectile theProj(track);
  if(verboseLevel>1) {
    G4cout << "G4HadronElasticProcess::PostStepDoIt for " 
	   << part->GetParticleName()
	   << " in " << material->GetName() 
	   << " Target Z= " << targNucleus->GetZ_asInt() 
	   << " A= " << targNucleus->GetA_asInt() << G4endl; 
  }

  G4HadFinalState* result = 0;
  try
    {
      result = hadi->ApplyYourself( theProj, *targNucleus);
    }
  catch(G4HadronicException aR)
    {
      G4ExceptionDescription ed;
      ed << "Call for " << hadi->GetModelName() << G4endl;
      ed << "Target element "<< elm->GetName()<<"  Z= " 
	 << targNucleus->GetZ_asInt() 
	 << "  A= " << targNucleus->GetA_asInt() << G4endl;
      DumpState(track,"ApplyYourself",ed);
      ed << " ApplyYourself failed" << G4endl;
      G4Exception("G4HadronElasticProcess::PostStepDoIt", "had006", 
		  FatalException, ed);
    }

  // Check the result for catastrophic energy non-conservation
  // cannot be applied because is not guranteed that recoil 
  // nucleus is created
  // result = CheckResult(theProj, targNucleus, result);

  // directions
  G4ThreeVector indir = track.GetMomentumDirection();
  G4double phi = CLHEP::twopi*G4UniformRand();
  G4ThreeVector it(0., 0., 1.);
  G4ThreeVector outdir = result->GetMomentumChange();

  if(verboseLevel>1) {
    G4cout << "Efin= " << result->GetEnergyChange()
	   << " de= " << result->GetLocalEnergyDeposit()
	   << " nsec= " << result->GetNumberOfSecondaries()
	   << " dir= " << outdir
	   << G4endl;
  }

  // energies  
  G4double edep = result->GetLocalEnergyDeposit();
  G4double efinal = result->GetEnergyChange();
  if(efinal < 0.0) { efinal = 0.0; }
  if(edep < 0.0)   { edep = 0.0; }

  // protection against numerical problems
  if((part == theNeutron && efinal <= lowestEnergyNeutron) || 
     (part != theNeutron && efinal <= lowestEnergy)) {
    edep += efinal;
    efinal = 0.0;
  }

  // primary change
  aParticleChange.ProposeEnergy(efinal);

  G4TrackStatus status = track.GetTrackStatus();
  if(efinal > 0.0) {
    outdir.rotate(phi, it);
    outdir.rotateUz(indir);
    aParticleChange.ProposeMomentumDirection(outdir);
  } else {
    if(part->GetProcessManager()->GetAtRestProcessVector()->size() > 0)
         { status = fStopButAlive; }
    else { status = fStopAndKill; }
    aParticleChange.ProposeTrackStatus(status);
  }

  //G4cout << "Efinal= " << efinal << "  TrackStatus= " << status << G4endl;

  aParticleChange.SetNumberOfSecondaries(0);

  // recoil
  if(result->GetNumberOfSecondaries() > 0) {
    G4DynamicParticle* p = result->GetSecondary(0)->GetParticle();

    if(p->GetKineticEnergy() > lowestEnergy) {
      aParticleChange.SetNumberOfSecondaries(1);
      G4ThreeVector pdir = p->GetMomentumDirection();
      // G4cout << "recoil " << pdir << G4endl;
      //!! is not needed for models inheriting G4HadronElastic
      pdir.rotate(phi, it);
      pdir.rotateUz(indir);
      // G4cout << "recoil rotated " << pdir << G4endl;
      p->SetMomentumDirection(pdir);

      // in elastic scattering time and weight are not changed
      G4Track* t = new G4Track(p, track.GetGlobalTime(), 
			       track.GetPosition());
      t->SetWeight(weight);
      t->SetTouchableHandle(track.GetTouchableHandle());
      aParticleChange.AddSecondary(t);

    } else {
      edep += p->GetKineticEnergy();
      delete p;
    }
  }
  aParticleChange.ProposeLocalEnergyDeposit(edep);
  aParticleChange.ProposeNonIonizingEnergyDeposit(edep);
  result->Clear();

  return G4VDiscreteProcess::PostStepDoIt(track,step);
}

