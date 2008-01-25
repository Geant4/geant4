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
// Created:
// 16.01.08 V.Ivanchenko use initialization similar to other CHIPS models
//

#include "G4ChiralInvariantPhaseSpace.hh"
#include "G4ParticleTable.hh"
#include "G4LorentzVector.hh"
#include "G4DynamicParticle.hh"
#include "G4IonTable.hh"
#include "G4Neutron.hh"

G4ChiralInvariantPhaseSpace::G4ChiralInvariantPhaseSpace()
{}

G4ChiralInvariantPhaseSpace::~G4ChiralInvariantPhaseSpace()
{}

G4HadFinalState * G4ChiralInvariantPhaseSpace::ApplyYourself(
	   const G4HadProjectile& aTrack, 
	   G4Nucleus& aTargetNucleus, 
	   G4HadFinalState * aChange)
{
  G4HadFinalState * aResult;
  if(aChange != 0)
  {
    aResult = aChange;
  }
  else
  {
    aResult = & theResult;
    aResult->Clear();
    aResult->SetStatusChange(stopAndKill);
  }
  //projectile properties needed in constructor of quasmon
  G4LorentzVector proj4Mom;
  proj4Mom = aTrack.Get4Momentum();
  G4int projectilePDGCode = aTrack.GetDefinition()
				  ->GetPDGEncoding();
  
  //target properties needed in constructor of quasmon
  G4int targetZ = G4int(aTargetNucleus.GetZ()+0.5);
  G4int targetA = G4int(aTargetNucleus.GetN()+0.5);
  G4int targetPDGCode = 90000000 + 1000*targetZ + (targetA-targetZ);
  // NOT NECESSARY ______________
  G4double targetMass = G4ParticleTable::GetParticleTable()->GetIonTable()
                                                           ->GetIonMass(targetZ, targetA);
  G4LorentzVector targ4Mom(0.,0.,0.,targetMass);  
  // END OF NOT NECESSARY^^^^^^^^

  // V.Ivanchenko set the same value as for other CHIPS models
  G4int nop = 152; 
  //G4int nop = 164; // nuclear clusters up to A=21
  G4double fractionOfSingleQuasiFreeNucleons = 0.4;
  G4double fractionOfPairedQuasiFreeNucleons = 0.0;
  if(targetA>27) fractionOfPairedQuasiFreeNucleons = 0.04;
  G4double clusteringCoefficient = 4.;
  G4double temperature = 180.;
  G4double halfTheStrangenessOfSee = 0.1; // = s/d = s/u
  G4double etaToEtaPrime = 0.3;
  
  // construct and fragment the quasmon
  //G4QCHIPSWorld aWorld(nop);              // Create CHIPS World of nop particles
  G4QCHIPSWorld::Get()->GetParticles(nop);  // Create CHIPS World of nop particles
  G4QNucleus::SetParameters(fractionOfSingleQuasiFreeNucleons,
                            fractionOfPairedQuasiFreeNucleons,
			    clusteringCoefficient);
  G4Quasmon::SetParameters(temperature,
                           halfTheStrangenessOfSee,
			   etaToEtaPrime);
  //  G4QEnvironment::SetParameters(solidAngle);
//  G4cout << "Input info "<< projectilePDGCode << " " 
//         << targetPDGCode <<" "
//	 << 1./MeV*proj4Mom<<" "
//	 << 1./MeV*targ4Mom << " "
//	 << nop << G4endl;
  G4QHadronVector projHV;
  G4QHadron* iH = new G4QHadron(projectilePDGCode, 1./MeV*proj4Mom);
  projHV.push_back(iH);
  G4QEnvironment* pan= new G4QEnvironment(projHV, targetPDGCode);
  //G4Quasmon* pan= new G4Quasmon(projectilePDGCode, targetPDGCode, 1./MeV*proj4Mom, 1./MeV*targ4Mom, nop);
  G4QHadronVector* output=0;
  try
  {
    output = pan->Fragment();
  }
  catch(G4HadronicException & aR)
  {
    G4cerr << "Exception thrown passing through G4ChiralInvariantPhaseSpace "<<G4endl;
    G4cerr << " targetPDGCode = "<< targetPDGCode <<G4endl;
    G4cerr << " Dumping the information in the pojectile list"<<G4endl;
    for(size_t i=0; i< projHV.size(); i++)
    {
      G4cerr <<"  Incoming 4-momentum and PDG code of "<<i<<"'th hadron: "
             <<" "<< projHV[i]->Get4Momentum()<<" "<<projHV[i]->GetPDGCode()<<G4endl;
    }
    throw;
  }
  std::for_each(projHV.begin(), projHV.end(), DeleteQHadron());
  projHV.clear();
  delete pan;
  
  // Fill the particle change.
  G4DynamicParticle * theSec;
#ifdef CHIPSdebug
  G4cout << "G4ChiralInvariantPhaseSpace: NEW EVENT #ofHadrons="
	 <<output->size()<<G4endl;
#endif
  unsigned int particle;
  for( particle = 0; particle < output->size(); particle++)
  {
    if(output->operator[](particle)->GetNFragments() != 0) 
    {
      delete output->operator[](particle);
      continue;
    }
    theSec = new G4DynamicParticle;  
    G4int pdgCode = output->operator[](particle)->GetPDGCode();
#ifdef CHIPSdebug
    G4cout << "G4ChiralInvariantPhaseSpace: h#"<<particle
	   <<", PDG="<<pdgCode<<G4endl;
#endif
    G4ParticleDefinition * theDefinition;
    // Note that I still have to take care of strange nuclei
    // For this I need the mass calculation, and a changed interface
    // for ion-tablel ==> work for Hisaya @@@@@@@
    // Then I can sort out the pdgCode. I also need a decau process 
    // for strange nuclei; may be another chips interface
    if(pdgCode>90000000) 
    {
      G4int aZ = (pdgCode-90000000)/1000;
      G4int anN = pdgCode-90000000-1000*aZ;
      theDefinition = G4ParticleTable::GetParticleTable()->FindIon(aZ,anN+aZ,0,aZ);
      if(aZ == 0 && anN == 1) theDefinition = G4Neutron::Neutron();
    }
    else
    {
      theDefinition = G4ParticleTable::GetParticleTable()
	->FindParticle(output->operator[](particle)->GetPDGCode());
    }
    theSec->SetDefinition(theDefinition);
    theSec->SetMomentum(output->operator[](particle)->Get4Momentum().vect());
    aResult->AddSecondary(theSec); 
    delete output->operator[](particle);
  }
  delete output;
  return aResult;
}

