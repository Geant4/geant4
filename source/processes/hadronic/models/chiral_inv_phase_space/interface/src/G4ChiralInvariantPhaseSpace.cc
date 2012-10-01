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

// Modified:
// 24.08.10 A. Dotti (andrea.dotti@cern.ch) handle exceptions
//          thrown by Q4QEnvironment::Fragment retying interaction
// 17.06.10 A. Dotti (andrea.dotti@cern.ch) handle case in which 
//          Q4QEnvironment returns a 90000000 fragment (see code comments)

// Created:
// 16.01.08 V.Ivanchenko use initialization similar to other CHIPS models
//

#include "G4ChiralInvariantPhaseSpace.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4LorentzVector.hh"
#include "G4DynamicParticle.hh"
#include "G4IonTable.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"

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
  G4int targetZ = aTargetNucleus.GetZ_asInt();
  G4int targetA = aTargetNucleus.GetA_asInt();
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

  //AND
  //A. Dotti 24 Aug. : Trying to handle situation when G4QEnvironment::Fragment throws an exception
  //                   Seen by ATLAS for gamma+Nuclear interaction for Gamma@2.4GeV on Al
  //                   The poor-man solution here is to re-try interaction if a G4QException is catched
  //                   Warning: G4QExcpetion does NOT inherit from base class G4HadException
  G4QHadronVector* output=0;

  bool retry=true;
  int retryes=0;
  int maxretries=10;

  while ( retry && retryes < maxretries ) 
    {
      G4QEnvironment* pan= new G4QEnvironment(projHV, targetPDGCode);

      try
	{
	  ++retryes;
	  output = pan->Fragment();
	  retry=false;//If here, Fragment did not throw! (AND)
	}
      catch( ... )
	{
	  G4cerr << "***WARNING*** Exception thrown passing through G4ChiralInvariantPhaseSpace "<<G4endl;
	  G4cerr << " targetPDGCode = "<< targetPDGCode <<G4endl;
	  G4cerr << " Dumping the information in the pojectile list"<<G4endl;
	  for(size_t i=0; i< projHV.size(); i++)
	    {
	      G4cerr <<"  Incoming 4-momentum and PDG code of "<<i<<"'th hadron: "
		     <<" "<< projHV[i]->Get4Momentum()<<" "<<projHV[i]->GetPDGCode()<<G4endl;
	    }
	  G4cerr << "Retrying interaction "<<G4endl; //AND
	  //throw; //AND
	}
      delete pan;
    } //AND
  if ( retryes >= maxretries ) //AND
    {
      G4cerr << "***ERROR*** Maximum number of retries ("<<maxretries<<") reached for G4QEnvironment::Fragment(), exception is being thrown" << G4endl;
      throw G4HadronicException(__FILE__,__LINE__,"G4ChiralInvariantPhaseSpace::ApplyYourself(...) - Maximum number of re-tries reached, abandon interaction");
    }
  
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
    //AND
    else if ( pdgCode == 90000000 && output->operator[](particle)->Get4Momentum().m()<1*MeV )
      {
	//A. Dotti: 
	//We cure the case the model returns a (A,Z)=0,0 G4QHadron with a very small mass
	//We convert this to a gamma. According to the author of the model this is the 
	//correct thing to do and it is done also in other parts of the CHIPS package
	theDefinition = G4Gamma::Gamma();
      }
    //AND
    else
    {
      theDefinition = G4ParticleTable::GetParticleTable()
	->FindParticle(output->operator[](particle)->GetPDGCode());
    }

    //AND
    //A. Dotti: Handle problematic cases in which one of the products has not been recognized
    //          This should never happen but we want to add an extra-protection
    if ( theDefinition == NULL )
      {
	//If we arrived here something bad really happened. We do not know how to handle the products of the interaction, we give up, resetting the 
	//result and keeping the primary alive.
	G4cerr<<"**WARNING*** G4ChiralInvariantPhaseSpace::ApplyYourself(...) : G4QEnvironment::Fragment() returns an invalid fragment\n with fourMom(MeV)=";
	G4cerr<<output->operator[](particle)->Get4Momentum()<<" and mass(MeV)="<<output->operator[](particle)->Get4Momentum().m();
	G4cerr<<". Offending PDG is:"<<pdgCode<<" abandon interaction. \n Taget PDG was:"<<targetPDGCode<<" \n Dumping the information in the projectile list:\n";
	for(size_t i=0; i< projHV.size(); i++)
	  {
	    G4cerr <<" Incoming 4-momentum and PDG code of "<<i<<"'th hadron: "
		 <<" "<< projHV[i]->Get4Momentum()<<" "<<projHV[i]->GetPDGCode()<<"\n";
	  }
	G4cerr<<"\n Please report as bug \n***END OF MESSAGE***"<<G4endl;

	for ( unsigned int cparticle=0 ; cparticle<output->size();++cparticle)
	  delete output->operator[](cparticle);
	delete output;
	std::for_each(projHV.begin(), projHV.end(), DeleteQHadron());
	projHV.clear();

	aResult->Clear();
	aResult->SetStatusChange(isAlive);
	aResult->SetEnergyChange(aTrack.GetKineticEnergy());
        aResult->SetMomentumChange(aTrack.Get4Momentum().vect().unit());
	return aResult;
      }
    //AND


    theSec->SetDefinition(theDefinition);
    theSec->SetMomentum(output->operator[](particle)->Get4Momentum().vect());
    aResult->AddSecondary(theSec); 
    delete output->operator[](particle);
  }
  delete output;
  std::for_each(projHV.begin(), projHV.end(), DeleteQHadron());
  projHV.clear();
  return aResult;
}

