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
// particle_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 100413 Fix bug in incidence energy by T. Koi  
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPEnAngCorrelation.hh"
#include "G4LorentzRotation.hh"
#include "G4LorentzVector.hh"
#include "G4RotationMatrix.hh"
#include "G4IonTable.hh"

G4ReactionProduct * G4ParticleHPEnAngCorrelation::SampleOne(G4double anEnergy)
{ 
  G4ReactionProduct * result = new G4ReactionProduct;
  
  // do we have an appropriate distribution
   if(nProducts!=1) throw G4HadronicException(__FILE__, __LINE__, "More than one product in SampleOne");
  
  // get the result
  G4ReactionProductVector * temp=0;
  G4int i=0;

  G4int icounter=0;
  G4int icounter_max=1024;
  while(temp == 0) {
      icounter++;
      if ( icounter > icounter_max ) {
	 G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
         break;
      }
     temp = theProducts[i++].Sample(anEnergy,1); 
  }
  
  // is the multiplicity correct
  if(temp->size()!=1) throw G4HadronicException(__FILE__, __LINE__, "SampleOne: Yield not correct");
  
  // fill result
  result = temp->operator[](0);
  
  // some garbage collection
  delete temp;
  
  // return result
  return result;
}

G4ReactionProductVector * G4ParticleHPEnAngCorrelation::Sample(G4double anEnergy)
{
  G4ReactionProductVector * result = new G4ReactionProductVector;
  G4int i;
  G4ReactionProductVector * it;
  G4ReactionProduct theCMS;
  G4LorentzRotation toZ;

  if(frameFlag==2 
     || frameFlag==3) // Added for particle HP
  {
    // simplify and double check @
    G4ThreeVector the3IncidentPart = fCache.Get().theProjectileRP->GetMomentum(); //theProjectileRP has value in LAB
    G4double nEnergy = fCache.Get().theProjectileRP->GetTotalEnergy();
    G4ThreeVector the3Target = fCache.Get().theTarget->GetMomentum();  //theTarget has value in LAB
    G4double tEnergy = fCache.Get().theTarget->GetTotalEnergy();
    G4double totE = nEnergy+tEnergy;
    G4ThreeVector the3CMS = the3Target+the3IncidentPart;
    theCMS.SetMomentum(the3CMS);
    G4double cmsMom = std::sqrt(the3CMS*the3CMS);
    G4double sqrts = std::sqrt((totE-cmsMom)*(totE+cmsMom));
    theCMS.SetMass(sqrts);
    theCMS.SetTotalEnergy(totE);
    G4ReactionProduct aIncidentPart;
    aIncidentPart.Lorentz(*fCache.Get().theProjectileRP, theCMS);
    //TKDB 100413 
    //ENDF-6 Formats Manual ENDF-102
    //CHAPTER 6. FILE 6: PRODUCT ENERGY-ANGLE DISTRIBUTIONS
    //LCT Reference system for secondary energy and angle (incident energy is always given in the LAB system)
    //anEnergy = aIncidentPart.GetKineticEnergy();
    anEnergy = fCache.Get().theProjectileRP->GetKineticEnergy(); //should be same argumment of "anEnergy"

    G4LorentzVector Ptmp (aIncidentPart.GetMomentum(), aIncidentPart.GetTotalEnergy());

    toZ.rotateZ(-1*Ptmp.phi());
    toZ.rotateY(-1*Ptmp.theta());
  }
  fCache.Get().theTotalMeanEnergy=0;
  G4LorentzRotation toLab(toZ.inverse()); //toLab only change axis NOT to LAB system
  //- get first number of particles, to check if sum of Z and N is not bigger than target values
  std::vector<int> nParticles;
  bool bNPOK = true;
//TKDB_PHP_150507
#ifdef PHP_AS_HP
#endif
//TKDB_PHP_161107
  G4int iTry(0);
//TKDB_PHP_161107
//TKDB_PHP_150507
  do {
    G4int sumZ = 0;
    G4int sumA = 0;
    nParticles.clear();
    for(i=0; i<nProducts; i++) 
      {
	G4int massCode = G4int(theProducts[i].GetMassCode());
	G4int nPart;
	nPart = theProducts[i].GetMultiplicity(anEnergy);
	sumZ += massCode/1000 * nPart;
	sumA += massCode % 1000 * nPart;
#ifdef G4PHPDEBUG
	if( getenv("G4ParticleHPDebug") ) G4cout << i << " G4ParticleHPEnAngCorrelation::MULTIPLICITY " << massCode << " sumZ " << sumZ << " sumA " << sumA << " NPART " << nPart << G4endl;
#endif
	nParticles.push_back( nPart );
      }
    bNPOK = true;
    double targetZ = fCache.Get().theTarget->GetDefinition()->GetAtomicNumber();
    double targetA = fCache.Get().theTarget->GetDefinition()->GetAtomicMass();
    targetZ += fCache.Get().theProjectileRP->GetDefinition()->GetAtomicNumber();
    targetA += fCache.Get().theProjectileRP->GetDefinition()->GetAtomicMass();
    if ( bAdjustFinalState ) {
/*
G4cout << "TKDB G4ParticleHPEnAngCorrelation::Sample 1" << G4endl;
G4cout << "TKDB "
<< "targetZ = " << targetZ
<< ", targetA = " << targetA
<< ", sumZ = " << sumZ
<< ", sumA = " << sumA
<< ", int( targetZ-sumZ ) = " << int( targetZ-sumZ )
<< ", int( targetA-sumA ) = " << int( targetA-sumA )
//<< ", G4IonTable::GetIonTable()->GetIon ( int(targetZ - sumZ), (int)(targetA - sumA), 0.0 ) = " << G4IonTable::GetIonTable()->GetIon ( int(targetZ - sumZ), (int)(targetA - sumA), 0.0 )
<< G4endl;
*/
      //if ( (sumZ != targetZ || sumA != targetA ) && 
      //   (sumZ > targetZ || sumA > targetA  
      //    || ! G4IonTable::GetIonTable()->GetIon ( int(targetZ - sumZ), (int)(targetA - sumA), 0.0 ) ) ){  // e.g. Z=3, A=2
      if ( ( sumZ != targetZ || sumA != targetA ) 
        && ( sumZ  > targetZ || sumA  > targetA || (targetZ-sumZ) >= (targetA-sumA) ) ) {  
                                                        // e.g. Z=3, A=2
	 bNPOK = false;
	 //nParticles.clear();
#ifdef G4PHPDEBUG
	  if ( getenv("G4ParticleHPDebug") ) 
	     G4cerr << " WRONG MULTIPLICITY Z= " << sumZ 
		 << " > " << targetZ
		 << " A= " <<  sumA 
		 << " > " << targetA << G4endl;
#endif
      }
    }
//TKDB_PHP_150507
#ifdef PHP_AS_HP
#endif
//TKDB_PHP_161107
   iTry++;
   if ( iTry > 1024 ) {
      G4Exception("G4ParticleHPEnAngCorrelation::Sample",
                  "Warning",
                  JustWarning,
                  "Too many trials were done. Exiting current loop by force. You may have Probably, the result violating (baryon number) conservation law will be obtained.");
      bNPOK=true;
   }
//TKDB_PHP_161107
//TKDB_PHP_150507

  }while(!bNPOK); // Loop checking, 11.05.2015, T. Koi

  for(i=0; i<nProducts; i++)
  {
    //-    if( nParticles[i] == 0 ) continue;
    it = theProducts[i].Sample(anEnergy,nParticles[i]); 
    G4double aMeanEnergy = theProducts[i].MeanEnergyOfThisInteraction();
    //    if( getenv("G4PHPTEST") ) G4cout << " EnAnG energy sampled " << it->operator[](0)->GetKineticEnergy() << " aMeanEnergy " << aMeanEnergy << G4endl; // GDEB
    //if(aMeanEnergy>0)
    //151120 TK Modified for solving reproducibility problem 
    //This change may have side effect. 
    if(aMeanEnergy>=0)
    {
      fCache.Get().theTotalMeanEnergy += aMeanEnergy;
    }
    else
    {
      fCache.Get().theTotalMeanEnergy = anEnergy/nProducts+theProducts[i].GetQValue();
    }
    if(it!=0)
    {
      for(unsigned int ii=0; ii<it->size(); ii++)
      {
	//if(!getenv("G4PHP_NO_LORENTZ_BOOST")) {
	G4LorentzVector pTmp1 (it->operator[](ii)->GetMomentum(),
			       it->operator[](ii)->GetTotalEnergy());
	pTmp1 = toLab*pTmp1;
    if( getenv("G4PHPTEST") )	G4cout << " G4particleHPEnAngCorrelation COS THETA " <<  std::cos(it->operator[](ii)->GetMomentum().theta()) << G4endl;
	it->operator[](ii)->SetMomentum(pTmp1.vect());
	it->operator[](ii)->SetTotalEnergy(pTmp1.e());
	if( getenv("G4PHPTEST") ) G4cout << " G4particleHPEnAngCorrelation COS THETA after toLab " <<  std::cos(it->operator[](ii)->GetMomentum().theta()) << G4endl;

	if(frameFlag==1) // target rest //TK 100413 should be LAB?
	{
	  it->operator[](ii)->Lorentz(*(it->operator[](ii)), -1.*(*fCache.Get().theTarget)); //TK 100413 Is this really need?
	}
	else if(frameFlag==2 ) // CMS
        {
#ifdef G4PHPDEBUG
	  if( getenv("G4ParticleHPDebug") ) 
	    G4cout <<"G4ParticleHPEnAngCorrelation: before Lorentz boost "<<
	      it->at(ii)->GetKineticEnergy()<<" "<<
	      it->at(ii)->GetMomentum()<<G4endl;
#endif
	  it->operator[](ii)->Lorentz(*(it->operator[](ii)), -1.*theCMS);
#ifdef G4PHPDEBUG
	  if( getenv("G4ParticleHPDebug") ) 
	    G4cout <<"G4ParticleHPEnAngCorrelation: after Lorentz boost "<<
		it->at(ii)->GetKineticEnergy()<<" "<<
	      it->at(ii)->GetMomentum()<<G4endl;
#endif
	}
        //TK120515 migrate frameFlag (MF6 LCT) = 3 
	else if(frameFlag==3) // CMS A<=4 other LAB
        {
           if ( theProducts[i].GetMassCode() > 4 ) //Alpha AWP 3.96713
           {
              //LAB
              it->operator[](ii)->Lorentz(*(it->operator[](ii)), -1.*(*fCache.Get().theTarget)); //TK 100413 Is this really need?
#ifdef G4PHPDEBUG
	  if( getenv("G4ParticleHPDebug") ) 
	    G4cout <<"G4ParticleHPEnAngCorrelation: after Lorentz boost "<<
		it->at(ii)->GetKineticEnergy()<<" "<<
	      it->at(ii)->GetMomentum()<<G4endl;
#endif
           }
           else
           {
              //CMS
              it->operator[](ii)->Lorentz(*(it->operator[](ii)), -1.*theCMS);
#ifdef G4PHPDEBUG
	  if( getenv("G4ParticleHPDebug") ) 
	    G4cout <<"G4ParticleHPEnAngCorrelation: after Lorentz boost "<<
		it->at(ii)->GetKineticEnergy()<<" "<<
	      it->at(ii)->GetMomentum()<<G4endl;
#endif
           }
        }
	else
	{
	  throw G4HadronicException(__FILE__, __LINE__, "G4ParticleHPEnAngCorrelation::Sample: The frame of the finalstate is not specified");
	}
	      if( getenv("G4PHPTEST") ) G4cout << frameFlag << " G4particleHPEnAngCorrelation COS THETA after Lorentz " <<  std::cos(it->operator[](ii)->GetMomentum().theta()) << G4endl;

	// }//getenv("G4PHP_NO_LORENTZ_BOOST")) 
	//	G4cout <<  ii << " EnAnG energy after boost " << it->operator[](ii)->GetKineticEnergy() << G4endl; //GDEB
	result->push_back(it->operator[](ii));
      }
      delete it;
    }
  }   

  return result;
}

