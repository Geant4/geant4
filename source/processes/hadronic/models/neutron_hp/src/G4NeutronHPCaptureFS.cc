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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 12-April-06 Enable IC electron emissions T. Koi 
// 26-January-07 Add G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION flag
// 081024 G4NucleiPropertiesTable:: to G4NucleiProperties::
// 101203 Bugzilla/Geant4 Problem 1155 Lack of residual in some case 
//
#include "G4NeutronHPCaptureFS.hh"
#include "G4Gamma.hh"
#include "G4ReactionProduct.hh"
#include "G4Nucleus.hh"
#include "G4PhotonEvaporation.hh"
#include "G4Fragment.hh"
#include "G4ParticleTable.hh" 
#include "G4NeutronHPDataUsed.hh"

  G4HadFinalState * G4NeutronHPCaptureFS::ApplyYourself(const G4HadProjectile & theTrack)
  {

    G4int i;
    theResult.Clear();
// prepare neutron
    G4double eKinetic = theTrack.GetKineticEnergy();
    const G4HadProjectile *incidentParticle = &theTrack;
    G4ReactionProduct theNeutron( const_cast<G4ParticleDefinition *>(incidentParticle->GetDefinition()) );
    theNeutron.SetMomentum( incidentParticle->Get4Momentum().vect() );
    theNeutron.SetKineticEnergy( eKinetic );

// prepare target
    G4ReactionProduct theTarget; 
    G4Nucleus aNucleus;
    G4double eps = 0.0001;
    if(targetMass<500*MeV)
      targetMass = ( G4NucleiProperties::GetNuclearMass( static_cast<G4int>(theBaseA+eps) , static_cast<G4int>(theBaseZ+eps) )) /
                     G4Neutron::Neutron()->GetPDGMass();
    G4ThreeVector neutronVelocity = 1./G4Neutron::Neutron()->GetPDGMass()*theNeutron.GetMomentum();
    G4double temperature = theTrack.GetMaterial()->GetTemperature();
    theTarget = aNucleus.GetBiasedThermalNucleus(targetMass, neutronVelocity, temperature);

// go to nucleus rest system
    theNeutron.Lorentz(theNeutron, -1*theTarget);
    eKinetic = theNeutron.GetKineticEnergy();

// dice the photons

    G4ReactionProductVector * thePhotons = 0;
    if ( HasFSData() && !getenv ( "G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION" ) ) 
    { 
      thePhotons = theFinalStatePhotons.GetPhotons(eKinetic);
    }
    else
    {
      G4ThreeVector aCMSMomentum = theNeutron.GetMomentum()+theTarget.GetMomentum();
      G4LorentzVector p4(aCMSMomentum, theTarget.GetTotalEnergy() + theNeutron.GetTotalEnergy());
      G4Fragment nucleus(static_cast<G4int>(theBaseA+1), static_cast<G4int>(theBaseZ) ,p4);
      G4PhotonEvaporation photonEvaporation;
      // T. K. add
      photonEvaporation.SetICM( TRUE );
      G4FragmentVector* products = photonEvaporation.BreakItUp(nucleus);
      G4FragmentVector::iterator i;
      thePhotons = new G4ReactionProductVector;
      for(i=products->begin(); i!=products->end(); i++)
      {
        G4ReactionProduct * theOne = new G4ReactionProduct;
        // T. K. add 
        if ( (*i)->GetParticleDefinition() != 0 ) 
           theOne->SetDefinition( (*i)->GetParticleDefinition() );
        else
           theOne->SetDefinition( G4Gamma::Gamma() ); // this definiion will be over writen
        
        // T. K. comment out below line
        //theOne->SetDefinition( G4Gamma::Gamma() );
        G4ParticleTable* theTable = G4ParticleTable::GetParticleTable();
        if((*i)->GetMomentum().mag() > 10*MeV) 
                 theOne->SetDefinition( 
                 theTable->FindIon(static_cast<G4int>(theBaseZ), static_cast<G4int>(theBaseA+1), 0, static_cast<G4int>(theBaseZ)) );
        theOne->SetMomentum( (*i)->GetMomentum().vect() ) ;
        theOne->SetTotalEnergy( (*i)->GetMomentum().t() );
        thePhotons->push_back(theOne);
        delete *i;
      } 
      delete products;
    }



// add them to the final state

    G4int nPhotons = 0;
    if(thePhotons!=0) nPhotons=thePhotons->size();
    G4int nParticles = nPhotons;
    if(1==nPhotons) nParticles = 2;


//Make at least one photon  
//101203 TK
    if ( nPhotons == 0 )
    {
       G4ReactionProduct * theOne = new G4ReactionProduct;
       theOne->SetDefinition( G4Gamma::Gamma() ); 
       G4double theta = pi*G4UniformRand();
       G4double phi = twopi*G4UniformRand();
       G4double sinth = std::sin(theta);
       G4ThreeVector direction( sinth*std::cos(phi), sinth*std::sin(phi), std::cos(theta) );
       theOne->SetMomentum( direction ) ;
       thePhotons->push_back(theOne);
       nPhotons++; // 0 -> 1
    }
//One photon case: energy set to Q-value 
//101203 TK
    if ( nPhotons == 1 )
    {
       G4ThreeVector direction = thePhotons->operator[](0)->GetMomentum().unit();
       G4double Q = G4ParticleTable::GetParticleTable()->FindIon(static_cast<G4int>(theBaseZ), static_cast<G4int>(theBaseA), 0, static_cast<G4int>(theBaseZ))->GetPDGMass() + G4Neutron::Neutron()->GetPDGMass()
         - G4ParticleTable::GetParticleTable()->FindIon(static_cast<G4int>(theBaseZ), static_cast<G4int>(theBaseA+1), 0, static_cast<G4int>(theBaseZ))->GetPDGMass();
       thePhotons->operator[](0)->SetMomentum( Q*direction );
    } 
//

    // back to lab system
    for(i=0; i<nPhotons; i++)
    {
      thePhotons->operator[](i)->Lorentz(*(thePhotons->operator[](i)), theTarget);
    }
    
    // Recoil, if only one gamma
    if (1==nPhotons)
    {
       G4DynamicParticle * theOne = new G4DynamicParticle;
       G4ParticleDefinition * aRecoil = G4ParticleTable::GetParticleTable()
                                        ->FindIon(static_cast<G4int>(theBaseZ), static_cast<G4int>(theBaseA+1), 0, static_cast<G4int>(theBaseZ));
       theOne->SetDefinition(aRecoil);
       // Now energy; 
       // Can be done slightly better @
       G4ThreeVector aMomentum =  theTrack.Get4Momentum().vect()
                                 +theTarget.GetMomentum()
				 -thePhotons->operator[](0)->GetMomentum();

       G4ThreeVector theMomUnit = aMomentum.unit();
       G4double aKinEnergy =  theTrack.GetKineticEnergy()
                             +theTarget.GetKineticEnergy(); // gammas come from Q-value
       G4double theResMass = aRecoil->GetPDGMass();
       G4double theResE = aRecoil->GetPDGMass()+aKinEnergy;
       G4double theAbsMom = std::sqrt(theResE*theResE - theResMass*theResMass);
       G4ThreeVector theMomentum = theAbsMom*theMomUnit;
       theOne->SetMomentum(theMomentum);
       theResult.AddSecondary(theOne);
    }

    // Now fill in the gammas.
    for(i=0; i<nPhotons; i++)
    {
      // back to lab system
      G4DynamicParticle * theOne = new G4DynamicParticle;
      theOne->SetDefinition(thePhotons->operator[](i)->GetDefinition());
      theOne->SetMomentum(thePhotons->operator[](i)->GetMomentum());
      theResult.AddSecondary(theOne);
      delete thePhotons->operator[](i);
    }
    delete thePhotons; 

//101203TK
    G4bool residual = false;
    G4ParticleDefinition * aRecoil = G4ParticleTable::GetParticleTable()
                                   ->FindIon(static_cast<G4int>(theBaseZ), static_cast<G4int>(theBaseA+1), 0, static_cast<G4int>(theBaseZ));
    for ( G4int i = 0 ; i != theResult.GetNumberOfSecondaries() ; i++ )
    {
       if ( theResult.GetSecondary(i)->GetParticle()->GetDefinition() == aRecoil ) residual = true;
    }

    if ( residual == false )
    {
       G4ParticleDefinition * aRecoil = G4ParticleTable::GetParticleTable()
                                        ->FindIon(static_cast<G4int>(theBaseZ), static_cast<G4int>(theBaseA+1), 0, static_cast<G4int>(theBaseZ));
       G4int nNonZero = 0;
       G4LorentzVector p_photons(0,0,0,0);
       for ( G4int i = 0 ; i != theResult.GetNumberOfSecondaries() ; i++ )
       {
          p_photons += theResult.GetSecondary(i)->GetParticle()->Get4Momentum();
          // To many 0 momentum photons -> Check PhotonDist 
          if ( theResult.GetSecondary(i)->GetParticle()->Get4Momentum() > 0 ) nNonZero++;
       }

       // Can we include kinetic energy here?
       G4double deltaE = ( theTrack.Get4Momentum().e() + theTarget.GetTotalEnergy() )
                       - ( p_photons.e() + aRecoil->GetPDGMass() );

//Add photons
       if ( nPhotons - nNonZero > 0 ) 
       {
              //G4cout << "TKDB G4NeutronHPCaptureFS::ApplyYourself we will create additional " << nPhotons - nNonZero << " photons" << G4endl;
          std::vector<G4double> vRand;
          vRand.push_back( 0.0 );
          for ( G4int i = 0 ; i != nPhotons - nNonZero - 1 ; i++ )
          { 
             vRand.push_back( G4UniformRand() );
          }
          vRand.push_back( 1.0 );
          std::sort( vRand.begin(), vRand.end() );

          std::vector<G4double> vEPhoton;
          for ( G4int i = 0 ; i < (G4int)vRand.size() - 1 ; i++ )
          {
             vEPhoton.push_back( deltaE * ( vRand[i+1] - vRand[i] ) );
          }
          std::sort( vEPhoton.begin(), vEPhoton.end() );

          for ( G4int i = 0 ; i < nPhotons - nNonZero - 1 ; i++ )
          {
             //Isotopic in LAB OK?
             G4double theta = pi*G4UniformRand();
             G4double phi = twopi*G4UniformRand();
             G4double sinth = std::sin(theta);
             G4double en = vEPhoton[i];
             G4ThreeVector tempVector(en*sinth*std::cos(phi), en*sinth*std::sin(phi), en*std::cos(theta) );
              
             p_photons += G4LorentzVector ( tempVector, tempVector.mag() );
             G4DynamicParticle * theOne = new G4DynamicParticle;
             theOne->SetDefinition( G4Gamma::Gamma() );
             theOne->SetMomentum( tempVector );
             theResult.AddSecondary(theOne);
          }

//        Add last photon 
          G4DynamicParticle * theOne = new G4DynamicParticle;
          theOne->SetDefinition( G4Gamma::Gamma() );
//        For better momentum conservation 
          G4ThreeVector lastPhoton = -p_photons.vect().unit()*vEPhoton.back();
          p_photons += G4LorentzVector( lastPhoton , lastPhoton.mag() );
          theOne->SetMomentum( lastPhoton );
          theResult.AddSecondary(theOne);
       }

//Add residual 
       G4DynamicParticle * theOne = new G4DynamicParticle;
       G4ThreeVector aMomentum = theTrack.Get4Momentum().vect() + theTarget.GetMomentum()
			       - p_photons.vect();
       theOne->SetDefinition(aRecoil);
       theOne->SetMomentum( aMomentum );
       theResult.AddSecondary(theOne);

    }
//101203TK END

// clean up the primary neutron
    theResult.SetStatusChange(stopAndKill);
    return &theResult;
  }

  void G4NeutronHPCaptureFS::Init (G4double A, G4double Z, G4String & dirName, G4String & )
  {
    G4String tString = "/FS/";
    G4bool dbool;
    G4NeutronHPDataUsed aFile = theNames.GetName(static_cast<G4int>(A), static_cast<G4int>(Z), dirName, tString, dbool);
    G4String filename = aFile.GetName();
    theBaseA = A;
    theBaseZ = G4int(Z+.5);
    if(!dbool || ( Z<2.5 && ( std::abs(theBaseZ - Z)>0.0001 || std::abs(theBaseA - A)>0.0001)))
    {
      hasAnyData = false;
      hasFSData = false; 
      hasXsec = false;
      return;
    }
    std::ifstream theData(filename, std::ios::in);
    
    hasFSData = theFinalStatePhotons.InitMean(theData); 
    if(hasFSData)
    {
      targetMass = theFinalStatePhotons.GetTargetMass();
      theFinalStatePhotons.InitAngular(theData); 
      theFinalStatePhotons.InitEnergies(theData); 
    }
    theData.close();
  }
