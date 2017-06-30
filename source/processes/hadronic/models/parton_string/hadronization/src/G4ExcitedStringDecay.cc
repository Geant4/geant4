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
#include "G4SystemOfUnits.hh"
#include "G4KineticTrack.hh"

//#define debug_G4ExcitedStringDecay
//#define debug_G4ExcitedStringCorr

G4ExcitedStringDecay::G4ExcitedStringDecay() : G4VStringFragmentation(),theStringDecay(0)
{}

G4ExcitedStringDecay::G4ExcitedStringDecay(G4VLongitudinalStringDecay * aStringDecay)
: G4VStringFragmentation(),
  theStringDecay(aStringDecay)
{}

G4ExcitedStringDecay::G4ExcitedStringDecay(const G4ExcitedStringDecay &)
: G4VStringFragmentation(),
  theStringDecay(0)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4ExcitedStringDecay::copy ctor not accessible");
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

G4KineticTrackVector *G4ExcitedStringDecay::FragmentString(const G4ExcitedString &theString)
{
 if ( theStringDecay == NULL ) theStringDecay=new G4LundStringFragmentation();
 return theStringDecay->FragmentString(theString);
}
	
G4KineticTrackVector *G4ExcitedStringDecay::FragmentStrings(const G4ExcitedStringVector * theStrings)
{
  G4LorentzVector KTsum(0.,0.,0.,0.);

#ifdef debug_G4ExcitedStringDecay
  G4cout<<G4endl;
  G4cout<<"--------------------------- G4ExcitedStringDecay ----------------------"<<G4endl;
  G4cout<<"Hadronization of Excited Strings: theStrings->size() "<<theStrings->size()<<G4endl;
#endif

  for ( unsigned int astring=0; astring < theStrings->size(); astring++)
//  for ( unsigned int astring=0; astring < 1; astring++)
  {
   if ( theStrings->operator[](astring)->IsExcited() )
        {KTsum+= theStrings->operator[](astring)->Get4Momentum();}
   else {KTsum+=theStrings->operator[](astring)->GetKineticTrack()->Get4Momentum();}
  }

  G4LorentzRotation toCms( -1 * KTsum.boostVector() );
  G4LorentzRotation toLab(toCms.inverse());
  G4LorentzVector Ptmp;
  KTsum=G4LorentzVector(0.,0.,0.,0.);

  for ( unsigned int astring=0; astring < theStrings->size(); astring++)
//  for ( unsigned int astring=0; astring < 1; astring++)
  {
   if ( theStrings->operator[](astring)->IsExcited() )
   {
    Ptmp=toCms * theStrings->operator[](astring)->GetLeftParton()->Get4Momentum();
    theStrings->operator[](astring)->GetLeftParton()->Set4Momentum(Ptmp);

    Ptmp=toCms * theStrings->operator[](astring)->GetRightParton()->Get4Momentum();
    theStrings->operator[](astring)->GetRightParton()->Set4Momentum(Ptmp);

    KTsum+= theStrings->operator[](astring)->Get4Momentum();
   }
   else
   {
    Ptmp=toCms * theStrings->operator[](astring)->GetKineticTrack()->Get4Momentum();
    theStrings->operator[](astring)->GetKineticTrack()->Set4Momentum(Ptmp);
    KTsum+= theStrings->operator[](astring)->GetKineticTrack()->Get4Momentum();
   }
  }

  G4KineticTrackVector * theResult = new G4KineticTrackVector;
  G4int attempts(0);
  G4bool success=false;
  G4bool NeedEnergyCorrector=false;
  do {
#ifdef debug_G4ExcitedStringDecay  
        G4cout<<"New try No "<<attempts<<" to hadronize strings"<<G4endl;
#endif

	std::for_each(theResult->begin() , theResult->end() , DeleteKineticTrack());
	theResult->clear();

	attempts++;

	G4LorentzVector KTsecondaries(0.,0.,0.,0.);
	NeedEnergyCorrector=false;

	for ( unsigned int astring=0; astring < theStrings->size(); astring++)
//	for ( unsigned int astring=0; astring < 1; astring++)                       // Uzhi
	{
#ifdef debug_G4ExcitedStringDecay  
          G4cout<<"String No "<<astring+1<<" Excited?  "<<theStrings->operator[](astring)->IsExcited()<<G4endl;

          G4cout<<"String No "<<astring+1<<" 4Momentum "<<theStrings->operator[](astring)->Get4Momentum()
          <<" "<<theStrings->operator[](astring)->Get4Momentum().mag()<<G4endl;
#endif

          G4KineticTrackVector * generatedKineticTracks = NULL;
	  if ( theStrings->operator[](astring)->IsExcited() )
	  {
#ifdef debug_G4ExcitedStringDecay  
             G4cout<<"Fragment String with partons: "
                   <<theStrings->operator[](astring)->GetLeftParton()->GetPDGcode() <<" "
                   <<theStrings->operator[](astring)->GetRightParton()->GetPDGcode()<<" "
                   <<"Direction "<<theStrings->operator[](astring)->GetDirection()<<G4endl;
#endif
  	     generatedKineticTracks=FragmentString(*theStrings->operator[](astring));
#ifdef debug_G4ExcitedStringDecay  
            G4cout<<"(G4ExcitedStringDecay) Number of produced hadrons = "
                  <<generatedKineticTracks->size()<<G4endl;
#endif
	  } else {
#ifdef debug_G4ExcitedStringDecay  
              G4cout<<"   GetTrack from the String"<<G4endl;
#endif
             G4LorentzVector Mom=theStrings->operator[](astring)->GetKineticTrack()->Get4Momentum();
             G4KineticTrack * aTrack= new G4KineticTrack(
                            theStrings->operator[](astring)->GetKineticTrack()->GetDefinition(),
                            theStrings->operator[](astring)->GetKineticTrack()->GetFormationTime(),
                            G4ThreeVector(0), Mom);

             aTrack->SetPosition(theStrings->operator[](astring)->GetKineticTrack()->GetPosition());

#ifdef debug_G4ExcitedStringDecay  
             G4cout<<"   A particle stored in the track is "<<aTrack->GetDefinition()->GetParticleName()<<G4endl;
#endif

	     generatedKineticTracks = new G4KineticTrackVector;
	     generatedKineticTracks->push_back(aTrack);
	  }    

	  if (generatedKineticTracks->size() == 0)
	  {
//	     G4cerr << "G4VPartonStringModel:No KineticTracks produced" << G4endl; 
//	     continue;                                                             
             success=false; NeedEnergyCorrector=false; break;                      
	  }

          G4LorentzVector KTsum1(0.,0.,0.,0.);
          for ( unsigned int aTrack=0; aTrack<generatedKineticTracks->size();aTrack++)
	  {
#ifdef debug_G4ExcitedStringDecay  
             G4cout<<"Prod part No. "<<aTrack+1<<" "
                   <<(*generatedKineticTracks)[aTrack]->GetDefinition()->GetParticleName()<<" "
                   <<(*generatedKineticTracks)[aTrack]->Get4Momentum()<<G4endl;
#endif
             theResult->push_back(generatedKineticTracks->operator[](aTrack));
             KTsum1+= (*generatedKineticTracks)[aTrack]->Get4Momentum();
	  }
	  KTsecondaries+=KTsum1;
	
#ifdef debug_G4ExcitedStringDecay  
          G4cout << "String secondaries(" <<generatedKineticTracks->size()<< ")"<<G4endl
                 <<"Init  string  momentum: "<< theStrings->operator[](astring)->Get4Momentum()<<G4endl
                 <<"Final hadrons momentum: "<< KTsum1 << G4endl;
#endif

	  if  ( KTsum1.e() > 0 && std::abs((KTsum1.e()-theStrings->operator[](astring)->Get4Momentum().e()) / KTsum1.e()) > perMillion )
	  {
	    NeedEnergyCorrector=true;
 	  }

#ifdef debug_G4ExcitedStringDecay  
          G4cout<<"NeedEnergyCorrection yes/no "<<NeedEnergyCorrector<<G4endl;
#endif

//        clean up
	  delete generatedKineticTracks;
	  success=true;                          
	}

	if ( NeedEnergyCorrector ) success=EnergyAndMomentumCorrector(theResult, KTsum);
  } while(!success && (attempts < 100));  /* Loop checking, 07.08.2015, A.Ribon */

  for ( unsigned int aTrack=0; aTrack<theResult->size();aTrack++)
  {
   Ptmp=(*theResult)[aTrack]->Get4Momentum();
   Ptmp.transform( toLab);
   (*theResult)[aTrack]->Set4Momentum(Ptmp);
  }

#ifdef debug_G4ExcitedStringDecay  
  G4cout<<"End of the strings fragmentation (G4ExcitedStringDecay)"<<G4endl;

  G4LorentzVector  KTsum1(0.,0.,0.,0.); 

  for ( unsigned int aTrack=0; aTrack<theResult->size();aTrack++)
  {
      G4cout << " corrected tracks .. " << (*theResult)[aTrack]->GetDefinition()->GetParticleName()
      <<"  " << (*theResult)[aTrack]->Get4Momentum() 
      <<"  " << (*theResult)[aTrack]->Get4Momentum().mag()<< G4endl;
      KTsum1+= (*theResult)[aTrack]->Get4Momentum();
  }

  G4cout << "Needcorrector/success " << NeedEnergyCorrector << "/" << success 
         << ", Corrected total  4 momentum " << KTsum1  << G4endl;
  if ( ! success ) G4cout << "failed to correct E/p" << G4endl;  

  G4cout<<"End of the Hadronization (G4ExcitedStringDecay)"<<G4endl;
#endif

if(!success)                                                                   // Uzhi 28 June 2016
{                                                                              // Uzhi 28 June 2016
 if(theResult->size() != 0)                                                    // Uzhi 7 Sept. 2016 it was theResult != 0
 {std::for_each(theResult->begin() , theResult->end() , DeleteKineticTrack()); // Uzhi 28 June 2016
  theResult->clear();                                                          // Uzhi 28 June 2016
  delete theResult; theResult=0;                                               // Uzhi 28 June 2016
 }
  for ( unsigned int astring=0; astring < theStrings->size(); astring++)       // Uzhi 28 June 2016
//  for ( unsigned int astring=0; astring < 1; astring++) // Uzhi 24 Oct. 2014                    Uzhi 2016 Need more correct 
  {
   if ( theStrings->operator[](astring)->IsExcited() )
   {
    Ptmp=theStrings->operator[](astring)->GetLeftParton()->Get4Momentum();
    Ptmp.transform( toLab);
    theStrings->operator[](astring)->GetLeftParton()->Set4Momentum(Ptmp);

    Ptmp=theStrings->operator[](astring)->GetRightParton()->Get4Momentum();
    Ptmp.transform( toLab);
    theStrings->operator[](astring)->GetRightParton()->Set4Momentum(Ptmp);
   }
   else
   {
    Ptmp=theStrings->operator[](astring)->GetKineticTrack()->Get4Momentum();
    Ptmp.transform( toLab);
    theStrings->operator[](astring)->GetKineticTrack()->Set4Momentum(Ptmp);
   }
  }                                                                      // Uzhi 24 Oct. 2014
}                                                                              // Uzhi 28 June 2016
  return theResult;
}

G4bool G4ExcitedStringDecay::EnergyAndMomentumCorrector
		(G4KineticTrackVector* Output, G4LorentzVector& TotalCollisionMom)   
  {
    const int    nAttemptScale = 500;
    const double ErrLimit = 1.E-5;
    if (Output->empty()) return TRUE;
    G4LorentzVector SumMom;
    G4double        SumMass = 0;
    G4double        TotalCollisionMass = TotalCollisionMom.m();

    std::vector<G4double> HadronMass; G4double HadronM(0.);                                       // Uzhi 8 Sept. 2016

#ifdef debug_G4ExcitedStringCorr
    G4cout<<G4endl<<"EnergyAndMomentumCorrector. Number of particles: "<<Output->size()<<G4endl;
#endif
    // Calculate sum hadron 4-momenta and summing hadron mass
    unsigned int cHadron;
    for(cHadron = 0; cHadron < Output->size(); cHadron++)
    {
        SumMom  += Output->operator[](cHadron)->Get4Momentum();
        HadronM=Output->operator[](cHadron)->Get4Momentum().mag(); HadronMass.push_back(HadronM); // Uzhi 8 Sept. 2016
        SumMass += Output->operator[](cHadron)->Get4Momentum().mag();  //GetDefinition()->GetPDGMass();
    }

#ifdef debug_G4ExcitedStringCorr
   G4cout<<"Sum part mom "<<SumMom<<" "<<SumMom.mag()<<G4endl
         <<"Sum str  mom "<<TotalCollisionMom<<" "<<TotalCollisionMom.mag()<<G4endl;
   G4cout<<"SumMass TotalCollisionMass "<<SumMass<<" "<<TotalCollisionMass<<G4endl;
#endif

    // Cannot correct a single particle
    if (Output->size() < 2) return FALSE;

    if (SumMass > TotalCollisionMass) return FALSE;
    SumMass = SumMom.m2();
    if (SumMass < 0) return FALSE;
    SumMass = std::sqrt(SumMass);

     // Compute c.m.s. hadron velocity and boost KTV to hadron c.m.s.
//    G4ThreeVector Beta = -SumMom.boostVector();          
    G4ThreeVector Beta = -TotalCollisionMom.boostVector(); 
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
        HadronM = HadronMass.at(cHadron);                                                         // Uzhi 8 Sept. 2016
        G4LorentzVector HadronMom = Output->operator[](cHadron)->Get4Momentum();
        HadronMom.setVect(Scale*HadronMom.vect());
        G4double E = std::sqrt(HadronMom.vect().mag2() + sqr(HadronM));                           // Uzhi 8 Sept. 2016
                                                       //sqr(Output->operator[](cHadron)->GetDefinition()->GetPDGMass()));
        HadronMom.setE(E);
        Output->operator[](cHadron)->Set4Momentum(HadronMom);
        Sum += E;
      } 
      Scale = TotalCollisionMass/Sum;    
#ifdef debug_G4ExcitedStringCorr
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

#ifdef debug_G4ExcitedStringCorr
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
