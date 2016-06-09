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
// Implementation of the HETC88 code into Geant4.
// Evaporation and De-excitation parts
// T. Lampen, Helsinki Institute of Physics, May-2000

#include "globals.hh"
#include "G4BertiniEvaporation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4BENeutronChannel.hh"
#include "G4BEProtonChannel.hh"
#include "G4BEDeuteronChannel.hh"
#include "G4BETritonChannel.hh"
#include "G4BEHe3Channel.hh"
#include "G4BEHe4Channel.hh"
#include "G4BEGammaDeexcitation.hh"

// verboseLevels : 4 inform about basic probabilities and branchings
//                 6 some details about calculations
//                 8 inform about final values
//                10 inform everything

G4BertiniEvaporation::G4BertiniEvaporation()
{
    G4cout << "Info G4BertiniEvaporation: This is still very fresh code."<< G4endl;
    G4cout << "     G4BertiniEvaporation: feed-back for improvement is very wellcome."<< G4endl;
    verboseLevel = 0;
 
    channelVector.push_back( new G4BENeutronChannel );
    channelVector.push_back( new G4BEProtonChannel  );
    channelVector.push_back( new G4BEDeuteronChannel);
    channelVector.push_back( new G4BETritonChannel  );
    channelVector.push_back( new G4BEHe3Channel     );
    channelVector.push_back( new G4BEHe4Channel     );
}


G4BertiniEvaporation::~G4BertiniEvaporation()
{
  while ( !channelVector.empty() )
    {
      delete channelVector.back();
      channelVector.pop_back();
    }
}


void G4BertiniEvaporation::setVerboseLevel( const G4int verbose )
{
  verboseLevel = verbose;

  // Update verbose level to all evaporation channels.
  std::vector< G4BertiniEvaporationChannel * >::iterator iChannel = channelVector.begin();  
  for ( ; iChannel != channelVector.end() ; *iChannel++ )                  
    ( *iChannel )->setVerboseLevel( verboseLevel );
}


G4FragmentVector * G4BertiniEvaporation::BreakItUp( G4LayeredNucleus & nucleus )
{
  G4int nucleusA;
  G4int nucleusZ;
  G4int i;
  G4double totalProbability;
  G4double excE;
  G4double newExcitation;
  G4double nucleusTotalMomentum;
  G4double nucleusKineticEnergy;
  G4double mRes; // Mass of residual nucleus.
  G4ThreeVector nucleusMomentumVector;
  G4DynamicParticle *pEmittedParticle;
  std::vector< G4DynamicParticle * > secondaryParticleVector;
  G4FragmentVector * result = new G4FragmentVector;
  
  // Read properties of the nucleus.
  nucleusA = nucleus.GetA_asInt();
  nucleusZ = nucleus.GetZ_asInt();
  excE                  = nucleus.GetEnergyDeposit();
  nucleusMomentumVector = nucleus.GetMomentum();

  // Move to CMS frame, save initial velocity of the nucleus to boostToLab vector.
  G4ThreeVector boostToLab( ( 1/G4NucleiProperties::GetNuclearMass( nucleusA, nucleusZ ) ) 
			    * nucleusMomentumVector ); // xx mass ok?

  if ( verboseLevel >= 10 )
    G4cout << " G4BertiniEvaporation : initial kinematics : boostToLab vector = " << boostToLab << G4endl
	   << "                     excitation energy  : " << excE<< G4endl;

  if ( nucleusA == 8 && nucleusZ == 4 )
    {
      splitBe8( excE, boostToLab, secondaryParticleVector );
      fillResult( secondaryParticleVector, result);
      return result;
    }
  
  // Initialize evaporation channels and calculate sum of emission
  // probabilities.
  std::vector< G4BertiniEvaporationChannel * >::iterator iChannel = channelVector.begin();  
  totalProbability = 0;
  for ( ; iChannel != channelVector.end() ; *iChannel++ )                  
     {
       ( *iChannel )->setNucleusA( nucleusA );
       ( *iChannel )->setNucleusZ( nucleusZ );
       ( *iChannel )->setExcitationEnergy( excE );
       ( *iChannel )->setPairingCorrection( 1 );
       ( *iChannel )->calculateProbability();
       totalProbability += ( *iChannel )->getProbability();
     }

  // If no evaporation process is possible, try without pairing energy.
  if ( totalProbability == 0 )
    {
      if ( verboseLevel >= 4 )
	G4cout << "G4BertiniEvaporation : no emission possible with pairing correction, trying without it" << G4endl;
      for ( iChannel = channelVector.begin() ; iChannel != channelVector.end() ; *iChannel++ )
	{
	  ( *iChannel )->setPairingCorrection(0);
	  ( *iChannel )->calculateProbability();
	  totalProbability += ( *iChannel )->getProbability();
	}
      if ( verboseLevel >= 4 )
	G4cout << "                       probability without correction " << totalProbability << G4endl;

    }

  // Normal evaporation cycle follows.
  // Particles are evaporated till all probabilities are zero.
  while ( totalProbability > 0 )
    {
      G4BertiniEvaporationChannel *pSelectedChannel;

      // Sample active emission channel.
      G4double sampledProbability = G4UniformRand() * totalProbability;
      if ( verboseLevel >= 10 ) G4cout << "          RandomProb : " << sampledProbability << G4endl;
      for ( iChannel = channelVector.begin() ; iChannel != channelVector.end() ; *iChannel++ )
	{
	  sampledProbability -= ( *iChannel )->getProbability();
	  if ( sampledProbability < 0 ) break;
	}
      pSelectedChannel = ( *iChannel );
      if ( iChannel == channelVector.end() )
	throw G4HadronicException(__FILE__, __LINE__,  "G4BertiniEvaporation : Error while sampling evaporation particle" );

      if ( verboseLevel >= 4 ) 
	G4cout << "G4BertiniEvaporation : particle " << pSelectedChannel->getName() << " selected " << G4endl;

      // Max 10 tries to get a physically acceptable particle energy
      // in this emission channel.
      i = 0;

      do
	{
	  pEmittedParticle = pSelectedChannel->emit();
	  // This loop checks that particle with too large energy is not emitted.
	  // CMS frame is considered in this loop. Nonrelativistic treatment. xxx


	  const G4int zRes = nucleusZ - pSelectedChannel->getParticleZ(); 
	  const G4int aRes  = nucleusA - pSelectedChannel->getParticleA(); 
	  // const G4double eBind = G4NucleiProperties::GetBindingEnergy( aRes, zRes );  // Binding energy of the nucleus.
	  mRes  = G4NucleiProperties::GetNuclearMass( aRes, zRes ); // Mass of the target nucleus
	  //      In HETC88:
	  //   	  eBind = Z * (-0.78244) + A * 8.36755 - cameron ( A , Z );
	  //  	  mRes =  zRes * 938.79304 + ( aRes - zRes ) * 939.57548 - eBind; 
	  //      where cameron ( A, Z ) is the mass excess by Cameron, see Canadian Journal of Physics,
	  //      vol. 35, 1957, p.1022

	  nucleusTotalMomentum = pEmittedParticle->GetTotalMomentum(); // CMS frame
	  nucleusKineticEnergy = std::pow( nucleusTotalMomentum, 2 ) / ( 2 * mRes );
	  newExcitation = excE - pEmittedParticle->GetKineticEnergy() - nucleusKineticEnergy - pSelectedChannel->getQ();

 	  if ( verboseLevel >= 10)
	    G4cout << "G4BertiniEvaporation : Kinematics " << G4endl
		   << "                    part kinetic E in CMS = " 
		   << pEmittedParticle->GetKineticEnergy() << G4endl
		   << "                    new excitation E      =  " 
		   << newExcitation << G4endl;

	  i++;
	  if ( !( newExcitation > 0 ) && i < 10) delete pEmittedParticle;
	} while ( !( newExcitation > 0 )  &&  i < 10 );

      if ( i >= 10 ) 
	{
	  // No appropriate particle energy found.
	  // Set probability of this channel to zero 
	  // and try to sample another particle on the 
	  // next round.
	  delete pEmittedParticle;

	  if ( verboseLevel >= 4 )
	    G4cout << "G4BertiniEvaporation : No appropriate energy for particle " 
		   << pSelectedChannel->getName() << " found." << G4endl;

	  pSelectedChannel->setProbability( 0 );

	  totalProbability = 0;
	  for ( ; iChannel != channelVector.end() ; *iChannel++ )                  
	    totalProbability += ( *iChannel )->getProbability();
	} // Treatment of physically unacceptable particle ends.
      else
	{
	  // Good particle found. 

	  // Transform particle properties to lab frame and save it.
	  G4LorentzVector particleFourVector = pEmittedParticle->Get4Momentum();
	  G4LorentzVector nucleusFourVector(  - pEmittedParticle->GetMomentum(), // CMS Frame.
					      nucleusKineticEnergy + mRes ); // Total energy.
	  particleFourVector.boost( boostToLab );
	  nucleusFourVector.boost(  boostToLab );
	  G4DynamicParticle *pPartLab  = new G4DynamicParticle( pEmittedParticle->GetDefinition(),
								particleFourVector.e(), // Total energy
								particleFourVector.vect() ); // momentum vector
	  secondaryParticleVector.push_back( pPartLab );

	  // Update residual nucleus and boostToLab vector.
	  nucleusA = nucleusA - pSelectedChannel->getParticleA();
	  nucleusZ = nucleusZ - pSelectedChannel->getParticleZ();
	  excE = newExcitation;
	  boostToLab = G4ThreeVector ( ( 1/mRes ) * nucleusFourVector.vect() );

	  if ( verboseLevel >= 10 )
	    G4cout << "  particle mom in cms " << pEmittedParticle->GetMomentum() << G4endl 
		   << "  particle mom in cms " << pEmittedParticle->GetTotalMomentum() << G4endl 
		   << "          Etot in cms " << pEmittedParticle->GetTotalEnergy() << G4endl
		   << "          Ekin in cms " << pEmittedParticle->GetKineticEnergy() << G4endl
		   << "        particle mass " << pEmittedParticle->GetMass() << G4endl
		   << "          boost  vect " << boostToLab << G4endl
		   << "          boosted 4v  " << particleFourVector << G4endl
		   << "          nucleus 4v  " << nucleusFourVector << G4endl
		   << "          nucl cm mom " << nucleusTotalMomentum << G4endl
		   << " particle k.e. in lab " << pPartLab->GetKineticEnergy() << G4endl
		   << "     new boost vector " << boostToLab << G4endl;

	  delete pEmittedParticle;

	  // If the residual nucleus is Be8, split it.
	  if ( nucleusA == 8 && nucleusZ == 4 )
	    {
	      splitBe8( excE, boostToLab, secondaryParticleVector );
              fillResult( secondaryParticleVector, result);
	      return result;
	    }

	  // Update evaporation channels and 
	  // total emission probability.
	  totalProbability = 0;
	  for ( iChannel = channelVector.begin() ; iChannel != channelVector.end() ; *iChannel++ )  
	    {
	      ( *iChannel )->setNucleusA( nucleusA ); 
	      ( *iChannel )->setNucleusZ( nucleusZ );
	      ( *iChannel )->setExcitationEnergy( excE ); 
	      ( *iChannel )->calculateProbability();
	      totalProbability += ( *iChannel )->getProbability();
	    }
	} // Treatment of physically acceptable particle ends.
    } // End of evaporation cycle.

  // Gamma de-excitation.
  G4BEGammaDeexcitation *pGammaChannel = new G4BEGammaDeexcitation;
  pGammaChannel->setNucleusA( nucleusA );
  pGammaChannel->setNucleusZ( nucleusZ );
  pGammaChannel->setVerboseLevel( verboseLevel );
  pGammaChannel->setExcitationEnergy( excE );

  // Emit equally sampled photons until all the excitation energy is
  // used; the last photon gets the remaining de-excitation energy.
  // Change of momentum of the nucleus is neglected.
  G4double gammaEnergy;
  G4double totalDeExcEnergy = 0;

  while ( excE > 0 )
    {
      pEmittedParticle = pGammaChannel->emit();
      gammaEnergy = pEmittedParticle->GetKineticEnergy();
      
      if ( ( totalDeExcEnergy + gammaEnergy ) > excE ) 
	{
	  // All the remaining energy is used here.
	  gammaEnergy = excE - totalDeExcEnergy;
	  excE = 0;
	}

      totalDeExcEnergy += gammaEnergy; 

      if ( verboseLevel >= 10 )
	G4cout << " G4BertiniEvaporation : gamma de-excitation, g of " << gammaEnergy << " MeV " << G4endl;
      
      G4LorentzVector gammaFourVector( pEmittedParticle->GetMomentum(),
				       pEmittedParticle->GetKineticEnergy() );
      gammaFourVector.boost( boostToLab );
      pEmittedParticle->SetKineticEnergy( gammaFourVector.e() );
      pEmittedParticle->SetMomentumDirection( gammaFourVector.vect().unit() );

      secondaryParticleVector.push_back( pEmittedParticle );

      } 

  delete pGammaChannel;

  // Residual nucleus is not returned.

  fillResult( secondaryParticleVector, result);
  
  return result;
}


void G4BertiniEvaporation::splitBe8( const G4double E,
				     const G4ThreeVector boostToLab,
				     std::vector< G4DynamicParticle * > & secondaryParticleVector )
{
  G4double kineticEnergy; 
  G4double u;
  G4double v;
  G4double w;
  const G4double Be8DecayEnergy = 0.093 * MeV; 

  if ( E <= 0 ) throw G4HadronicException(__FILE__, __LINE__,  "G4BertiniEvaporation : excitation energy < 0 " );
  if ( verboseLevel >= 4 ) G4cout << "     Be8 split to 2 x He4" << G4endl;

  kineticEnergy = 0.5 * ( E + Be8DecayEnergy ); 
  
  // Create two alpha particles in CMS frame.
  isotropicCosines( u , v , w );
  G4ThreeVector momentumDirectionCMS( u, v, w ); 

  G4DynamicParticle *pP1Cms = new G4DynamicParticle( G4Alpha::Alpha(),
						     momentumDirectionCMS,
						     kineticEnergy );
  G4DynamicParticle *pP2Cms = new G4DynamicParticle( G4Alpha::Alpha(),
						     -momentumDirectionCMS,
						     kineticEnergy );
  G4LorentzVector fourVector1( pP1Cms->GetMomentum(), 
			       pP1Cms->GetTotalEnergy() );
  G4LorentzVector fourVector2( pP2Cms->GetMomentum(), 
			       pP2Cms->GetTotalEnergy() );

  // Transform into lab frame by transforming the four vectors. Then
  // add to the vector of secondary particles.
  fourVector1.boost( boostToLab );
  fourVector2.boost( boostToLab );
  G4DynamicParticle * pP1Lab  = new G4DynamicParticle( G4Alpha::Alpha(),
						      fourVector1.vect(),
						      fourVector1.e() );
  G4DynamicParticle * pP2Lab  = new G4DynamicParticle( G4Alpha::Alpha(),
						      fourVector2.vect(),
						      fourVector2.e() );
  secondaryParticleVector.push_back( pP1Lab );
  secondaryParticleVector.push_back( pP2Lab );

  delete pP1Cms;
  delete pP2Cms;

  return;
}


void G4BertiniEvaporation::fillResult( std::vector<G4DynamicParticle *> secondaryParticleVector,
				      G4FragmentVector * aResult )
{
  // Fill the vector pParticleChange with secondary particles stored in vector.
  for ( size_t i = 0 ; i < secondaryParticleVector.size() ; i++ )
  {
    G4int aZ = static_cast<G4int> (secondaryParticleVector[i]->GetDefinition()->GetPDGCharge() );
    G4int aA = static_cast<G4int> (secondaryParticleVector[i]->GetDefinition()->GetBaryonNumber());
    G4LorentzVector aMomentum = secondaryParticleVector[i]->Get4Momentum();
    if(aA>0) 
    {
      aResult->push_back( new G4Fragment(aA, aZ, aMomentum) ); 
    }
    else 
    {
      aResult->push_back( new G4Fragment(aMomentum, secondaryParticleVector[i]->GetDefinition()) ); 
    }
  }
  return;
}


void G4BertiniEvaporation::isotropicCosines( G4double & u, G4double & v, G4double & w )
{
  // Samples isotropic random direction cosines.
  G4double CosTheta = 1.0 - 2.0 * G4UniformRand();
  G4double SinTheta = std::sqrt( 1.0 - CosTheta * CosTheta );
  G4double Phi = twopi * G4UniformRand();

  u = std::cos( Phi ) * SinTheta;
  v = std::cos( Phi ) * CosTheta,
  w = std::sin( Phi );

  return;
}
