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
// Historic fragment from M.Komogorov; clean-up still necessary @@@
#include "G4ExcitedStringDecay.hh"
#include "G4KineticTrack.hh"

G4ExcitedStringDecay::G4ExcitedStringDecay()
{
	theStringDecay=NULL;
}

G4ExcitedStringDecay::G4ExcitedStringDecay(G4VLongitudinalStringDecay * aStringDecay)
: theStringDecay(aStringDecay)
{}

G4ExcitedStringDecay::G4ExcitedStringDecay(const G4ExcitedStringDecay &) : G4VStringFragmentation()
{
} 

G4ExcitedStringDecay::~G4ExcitedStringDecay()
{
}

const G4ExcitedStringDecay & G4ExcitedStringDecay::operator=(const G4ExcitedStringDecay &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4ExcitedStringDecay::operator= meant to not be accessable");
  return *this;
}

int G4ExcitedStringDecay::operator==(const G4ExcitedStringDecay &) const
{
  return 0;
}

int G4ExcitedStringDecay::operator!=(const G4ExcitedStringDecay &) const
{
  return 1;
}

G4KineticTrackVector *G4ExcitedStringDecay::FragmentString
				(const G4ExcitedString &theString)
{
	if ( theStringDecay == NULL ) 

	    theStringDecay=new G4LundStringFragmentation();
    
	return theStringDecay->FragmentString(theString);
}
	
G4KineticTrackVector *G4ExcitedStringDecay::FragmentStrings
				(const G4ExcitedStringVector * theStrings)
{
  G4LorentzVector KTsum(0.,0.,0.,0.);
  for ( unsigned int astring=0; astring < theStrings->size(); astring++)
  {
	KTsum+= theStrings->operator[](astring)->Get4Momentum();

	if( !(KTsum.e()<1) && !(KTsum.e()>-1) )
	{
          throw G4HadronicException(__FILE__, __LINE__, 
	                           "G4ExcitedStringDecay::FragmentStrings received nan string...");
	}
  }

  G4KineticTrackVector * theResult = new G4KineticTrackVector;
  G4int attempts(0);
  G4bool success=false;
  do {

	std::for_each(theResult->begin() , theResult->end() , DeleteKineticTrack());
	theResult->clear();

	attempts++;
	G4LorentzVector KTsecondaries(0.,0.,0.,0.);
	G4bool NeedEnergyCorrector=false;

	for ( unsigned int astring=0; astring < theStrings->size(); astring++)
	{
//G4cout<<G4endl<<"String No "<<astring+1<<" "<<theStrings->operator[](astring)->Get4Momentum().mag()<<G4endl;

          G4KineticTrackVector * generatedKineticTracks = NULL;
  
	  if ( theStrings->operator[](astring)->IsExcited() )
	  {
  	     generatedKineticTracks=FragmentString(*theStrings->operator[](astring));
	  } else {
	     generatedKineticTracks = new G4KineticTrackVector;
	     generatedKineticTracks->push_back(theStrings->operator[](astring)->GetKineticTrack());
	  }    

	  if (generatedKineticTracks == NULL) 
	  {
	     G4cerr << "G4VPartonStringModel:No KineticTracks produced" << G4endl;
	     continue;
	  }

          G4LorentzVector KTsum1(0.,0.,0.,0.);
          for ( unsigned int aTrack=0; aTrack<generatedKineticTracks->size();aTrack++)
	  {
//G4cout<<" Prod part "<<(*generatedKineticTracks)[aTrack]->GetDefinition()->GetParticleName()<<" "<<(*generatedKineticTracks)[aTrack]->Get4Momentum()<<G4endl;
             theResult->push_back(generatedKineticTracks->operator[](aTrack));
             KTsum1+= (*generatedKineticTracks)[aTrack]->Get4Momentum();
	  }
	  KTsecondaries+=KTsum1;
	
	  if  ( KTsum1.e() > 0 && std::abs((KTsum1.e()-theStrings->operator[](astring)->Get4Momentum().e()) / KTsum1.e()) > perMillion ) 
	  {
//--debug--           G4cout << "String secondaries(" <<generatedKineticTracks->size()<< ")  momentum: " 
//--debug--	          << theStrings->operator[](astring)->Get4Momentum() << " " << KTsum1 << G4endl;
	    NeedEnergyCorrector=true;
 	  }

//        clean up
	  delete generatedKineticTracks;
	}
//--DEBUG  G4cout << "Strings/secs total  4 momentum " << KTsum << " " <<KTsecondaries << G4endl;

        success=true;
	if ( NeedEnergyCorrector ) success=EnergyAndMomentumCorrector(theResult, KTsum);

  } while(!success && (attempts < 100));

#ifdef debug_ExcitedStringDecay
  G4LorentzVector  KTsum1=0;
  for ( unsigned int aTrack=0; aTrack<theResult->size();aTrack++)
  {
      G4cout << " corrected tracks .. " << (*theResult)[aTrack]->GetDefinition()->GetParticleName()
      <<"  " << (*theResult)[aTrack]->Get4Momentum() << G4endl;
      KTsum1+= (*theResult)[aTrack]->Get4Momentum();
  }
  G4cout << "Needcorrector/success " << NeedEnergyCorrector << "/" << success << ", Corrected total  4 momentum " << KTsum1  << G4endl;
  if ( ! success ) G4cout << "failed to correct E/p" << G4endl;  
#endif

  return theResult;
}

G4bool G4ExcitedStringDecay::EnergyAndMomentumCorrector
		(G4KineticTrackVector* Output, G4LorentzVector& TotalCollisionMom)   
  {
    const int    nAttemptScale = 500;
    const double ErrLimit = 1.E-5;
    if (Output->empty())
       return TRUE;
    G4LorentzVector SumMom;
    G4double        SumMass = 0;     
    G4double        TotalCollisionMass = TotalCollisionMom.m();

    if( !(TotalCollisionMass<1) && !(TotalCollisionMass>-1) )
    {
      std::cout << "TotalCollisionMomentum = "<<TotalCollisionMom<<G4endl;
      throw G4HadronicException(__FILE__, __LINE__, "G4ExcitedStringDecay received nan mass...");
    }

//G4cout<<G4endl<<"EnergyAndMomentumCorrector "<<Output->size()<<G4endl;
    // Calculate sum hadron 4-momenta and summing hadron mass
    unsigned int cHadron;
    for(cHadron = 0; cHadron < Output->size(); cHadron++)
    {
        SumMom  += Output->operator[](cHadron)->Get4Momentum();
	if( !(SumMom<1) && !(SumMom>-1) )
	{
          throw G4HadronicException(__FILE__, __LINE__, "G4ExcitedStringDecay::EnergyAndMomentumCorrector() received nan momentum...");
	}
        SumMass += Output->operator[](cHadron)->GetDefinition()->GetPDGMass();
	if( !(SumMass<1) && !(SumMass>-1) )
	{
          throw G4HadronicException(__FILE__, __LINE__, "G4ExcitedStringDecay::EnergyAndMomentumCorrector() received nan mass...");
	}
    }

//G4cout<<"SumMass TotalCollisionMass "<<SumMass<<" "<<TotalCollisionMass<<G4endl;

    // Cannot correct a single particle
    if (Output->size() < 2) return FALSE;

    if (SumMass > TotalCollisionMass) return FALSE;
    SumMass = SumMom.m2();
    if (SumMass < 0) return FALSE;
    SumMass = std::sqrt(SumMass);

     // Compute c.m.s. hadron velocity and boost KTV to hadron c.m.s.
    G4ThreeVector Beta = -SumMom.boostVector();
    Output->Boost(Beta);

    // Scale total c.m.s. hadron energy (hadron system mass).
    // It should be equal interaction mass
    G4double Scale = 1;
    G4int cAttempt = 0;
    G4double Sum = 0;
    G4bool success = false;
    for(cAttempt = 0; cAttempt < nAttemptScale; cAttempt++)
    {
      Sum = 0;
      for(cHadron = 0; cHadron < Output->size(); cHadron++)
      {
        G4LorentzVector HadronMom = Output->operator[](cHadron)->Get4Momentum();
        HadronMom.setVect(Scale*HadronMom.vect());
        G4double E = std::sqrt(HadronMom.vect().mag2() + sqr(Output->operator[](cHadron)->GetDefinition()->GetPDGMass()));
        HadronMom.setE(E);
        Output->operator[](cHadron)->Set4Momentum(HadronMom);
        Sum += E;
      } 
      Scale = TotalCollisionMass/Sum;    
#ifdef debug_G4ExcitedStringDecay 
      G4cout << "Scale-1=" << Scale -1 
                << ",  TotalCollisionMass=" << TotalCollisionMass
		<< ",  Sum=" << Sum
		<< G4endl;
#endif     
      if (std::fabs(Scale - 1) <= ErrLimit) 
      {
        success = true;
	break;
      }
    }
#ifdef debug_G4ExcitedStringDecay     
    if(!success)
    {
      G4cout << "G4ExcitedStringDecay::EnergyAndMomentumCorrector - Warning"<<G4endl;
      G4cout << "   Scale not unity at end of iteration loop: "<<TotalCollisionMass<<" "<<Sum<<" "<<Scale<<G4endl;
      G4cout << "   Number of secondaries: " << Output->size() << G4endl;
      G4cout << "   Wanted total energy: " <<  TotalCollisionMom.e() << G4endl; 
      G4cout << "   Increase number of attempts or increase ERRLIMIT"<<G4endl;
//       throw G4HadronicException(__FILE__, __LINE__, "G4ExcitedStringDecay failed to correct...");
    }
#endif     
    // Compute c.m.s. interaction velocity and KTV back boost   
    Beta = TotalCollisionMom.boostVector();
    Output->Boost(Beta);

    return success;
  }
