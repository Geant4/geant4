// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Implementation of the HETC88 code into Geant4.
// Evaporation and De-excitation parts
// T. Lampen, Helsinki Institute of Physics, May-2000

#include "globals.hh"
#include "G4BertiniEvaporation.hh"
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
  vector< G4BertiniEvaporationChannel * >::iterator iChannel = channelVector.begin();  
  for ( ; iChannel != channelVector.end() ; *iChannel++ )                  
    ( *iChannel )->setVerboseLevel( verboseLevel );
}


G4VParticleChange * G4BertiniEvaporation::BreakItUp( G4LayeredNucleus & nucleus )
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
  vector< G4DynamicParticle * > secondaryParticleVector;
  G4ParticleChange *pParticleChange;

  pParticleChange = new G4ParticleChange(); 

  //************
  //    fillParticleChange( secondaryParticleVector, pParticleChange );
  //    return pParticleChange;
 //***********




  // Read properties of the nucleus.
  nucleusA = ( G4int ) nucleus.GetN(); // GetN should in fact get GetA
  nucleusZ = ( G4int ) nucleus.GetZ();
  excE                  = nucleus.GetEnergyDeposit();
  nucleusMomentumVector = nucleus.GetMomentum();

  // Move to CMS frame, save initial velocity of the nucleus to boostToLab vector.
  G4ThreeVector boostToLab( ( 1/G4NucleiProperties::GetAtomicMass( nucleusA, nucleusZ ) ) 
			    * nucleusMomentumVector ); // xx mass ok?

  if ( verboseLevel >= 10 )
    G4cout << " G4BertiniEvaporation : initial kinematics : boostToLab vector = " << boostToLab << endl
	   << "                     excitation energy  : " << excE<< endl;

  if ( nucleusA == 8 && nucleusZ == 4 )
    {
      splitBe8( excE, boostToLab, secondaryParticleVector );
      fillParticleChange( secondaryParticleVector, pParticleChange );
      return pParticleChange;
    }
  
  // Initialize evaporation channels and calculate sum of emission
  // probabilities.
  vector< G4BertiniEvaporationChannel * >::iterator iChannel = channelVector.begin();  
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
	G4cout << "G4BertiniEvaporation : no emission possible with pairing correction, trying without it" << endl;
      for ( iChannel = channelVector.begin() ; iChannel != channelVector.end() ; *iChannel++ )
	{
	  ( *iChannel )->setPairingCorrection(0);
	  ( *iChannel )->calculateProbability();
	  totalProbability += ( *iChannel )->getProbability();
	}
      if ( verboseLevel >= 4 )
	G4cout << "                       probability without correction " << totalProbability << endl;

    }

  // Normal evaporation cycle follows.
  // Particles are evaporated till all probabilities are zero.
  while ( totalProbability > 0 )
    {
      G4BertiniEvaporationChannel *pSelectedChannel;

      // Sample active emission channel.
      G4double sampledProbability = G4UniformRand() * totalProbability;
      if ( verboseLevel >= 10 ) G4cout << "          RandomProb : " << sampledProbability << endl;
      for ( iChannel = channelVector.begin() ; iChannel != channelVector.end() ; *iChannel++ )
	{
	  sampledProbability -= ( *iChannel )->getProbability();
	  if ( sampledProbability < 0 ) break;
	}
      pSelectedChannel = ( *iChannel );
      if ( iChannel == channelVector.end() )
	G4Exception( "G4BertiniEvaporation : Error while sampling evaporation particle" );

      if ( verboseLevel >= 4 ) 
	G4cout << "G4BertiniEvaporation : particle " << pSelectedChannel->getName() << " selected " << endl;

      // Max 10 tries to get a physically acceptable particle energy
      // in this emission channel.
      i = 0;

      do
	{
	  // This loop checks that particle with too large energy is not emitted.
	  // CMS frame is considered in this loop. Nonrelativistic treatment. xxx
	  pEmittedParticle = pSelectedChannel->emit();

	  const G4int zRes = nucleusZ - pSelectedChannel->getParticleZ(); 
	  const G4int aRes  = nucleusA - pSelectedChannel->getParticleA(); 
	  const G4double eBind = G4NucleiProperties::GetBindingEnergy( aRes, zRes );  // Binding energy of the nucleus.
	  mRes  = G4NucleiProperties::GetAtomicMass( aRes, zRes ); // Mass of the target nucleus
	  //      In HETC88:
	  //   	  eBind = Z * (-0.78244) + A * 8.36755 - cameron ( A , Z );
	  //  	  mRes =  zRes * 938.79304 + ( aRes - zRes ) * 939.57548 - eBind; 
	  //      where cameron ( A, Z ) is the mass excess by Cameron, see Canadian Journal of Physics,
	  //      vol. 35, 1957, p.1022

	  nucleusTotalMomentum = pEmittedParticle->GetTotalMomentum(); // CMS frame
	  nucleusKineticEnergy = pow( nucleusTotalMomentum, 2 ) / ( 2 * mRes );
	  newExcitation = excE - pEmittedParticle->GetKineticEnergy() - nucleusKineticEnergy - pSelectedChannel->getQ();

 	  if ( verboseLevel >= 10)
	    G4cout << "G4BertiniEvaporation : Kinematics " << endl
		   << "                    part kinetic E in CMS = " 
		   << pEmittedParticle->GetKineticEnergy() << endl
		   << "                    new excitation E      =  " 
		   << newExcitation << endl;

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
		   << pSelectedChannel->getName() << " found." << endl;

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
	    G4cout << "  particle mom in cms " << pEmittedParticle->GetMomentum() << endl 
		   << "  particle mom in cms " << pEmittedParticle->GetTotalMomentum() << endl 
		   << "          Etot in cms " << pEmittedParticle->GetTotalEnergy() << endl
		   << "          Ekin in cms " << pEmittedParticle->GetKineticEnergy() << endl
		   << "        particle mass " << pEmittedParticle->GetMass() << endl
		   << "          boost  vect " << boostToLab << endl
		   << "          boosted 4v  " << particleFourVector << endl
		   << "          nucleus 4v  " << nucleusFourVector << endl
		   << "          nucl cm mom " << nucleusTotalMomentum << endl
		   << " particle k.e. in lab " << pPartLab->GetKineticEnergy() << endl
		   << "     new boost vector " << boostToLab << endl;

	  delete pEmittedParticle;

	  // If the residual nucleus is Be8, split it.
	  if ( nucleusA == 8 && nucleusZ == 4 )
	    {
	      splitBe8( excE, boostToLab, secondaryParticleVector );
	      fillParticleChange( secondaryParticleVector, pParticleChange );
	      return pParticleChange;
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
	G4cout << " G4BertiniEvaporation : gamma de-excitation, g of " << gammaEnergy << " MeV " << endl;
      
      G4LorentzVector gammaFourVector( pEmittedParticle->GetMomentum(),
				       pEmittedParticle->GetKineticEnergy() );
      gammaFourVector.boost( boostToLab );
      pEmittedParticle->SetKineticEnergy( gammaFourVector.e() );
      pEmittedParticle->SetMomentumDirection( gammaFourVector.vect().unit() );

      secondaryParticleVector.push_back( pEmittedParticle );

      } 

  delete pGammaChannel;

  // Residual nucleus is not returned.

  fillParticleChange( secondaryParticleVector, pParticleChange );
  
  return pParticleChange;
}


void G4BertiniEvaporation::splitBe8( const G4double E,
				     const G4ThreeVector boostToLab,
				     vector< G4DynamicParticle * > & secondaryParticleVector )
{
  G4double kineticEnergy; 
  G4double u;
  G4double v;
  G4double w;
  const G4double Be8DecayEnergy = 0.093 * MeV; 

  if ( E <= 0 ) G4Exception( "G4BertiniEvaporation : excitation energy < 0 " );
  if ( verboseLevel >= 4 ) G4cout << "     Be8 split to 2 x He4" << endl;

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


void G4BertiniEvaporation::fillParticleChange( vector<G4DynamicParticle *> secondaryParticleVector,
					    G4ParticleChange * pParticleChange )
{
  // Fill the vector pParticleChange with secondary particles stored in vector.
  pParticleChange->SetNumberOfSecondaries( secondaryParticleVector.size() );
  for ( G4int i = 0 ; i < secondaryParticleVector.size() ; i++ )
    pParticleChange->AddSecondary( secondaryParticleVector[i] ); 
  return;
}


void G4BertiniEvaporation::isotropicCosines( G4double & u, G4double & v, G4double & w )
{
  // Samples isotropic random direction cosines.
  G4double CosTheta = 1.0 - 2.0 * G4UniformRand();
  G4double SinTheta = sqrt( 1.0 - CosTheta * CosTheta );
  G4double Phi = twopi * G4UniformRand();

  u = cos( Phi ) * SinTheta;
  v = cos( Phi ) * CosTheta,
  w = sin( Phi );

  return;
}
