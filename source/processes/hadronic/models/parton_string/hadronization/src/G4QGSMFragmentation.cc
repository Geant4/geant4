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
// $Id$
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Maxim Komogorov, 10-Jul-1998
// -----------------------------------------------------------------------------
#include "G4QGSMFragmentation.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4FragmentingString.hh"
#include "G4DiQuarks.hh"
#include "G4Quarks.hh"

// Class G4QGSMFragmentation 
//****************************************************************************************
 
G4QGSMFragmentation::G4QGSMFragmentation() :
arho(0.5), aphi(0.), an(-0.5), ala(-0.75), aksi(-1.), alft(0.5)
   {
   }

G4QGSMFragmentation::~G4QGSMFragmentation()
   {
   }

//----------------------------------------------------------------------------------------------------------

G4KineticTrackVector* G4QGSMFragmentation::FragmentString(const G4ExcitedString& theString)
{
//    Can no longer modify Parameters for Fragmentation.
	PastInitPhase=true;
	
// 	check if string has enough mass to fragment...
	G4KineticTrackVector * LeftVector=LightFragmentationTest(&theString);
	if ( LeftVector != 0 ) return LeftVector;
	
	LeftVector = new G4KineticTrackVector;
	G4KineticTrackVector * RightVector=new G4KineticTrackVector;

// this should work but its only a semi deep copy. %GF	G4ExcitedString theStringInCMS(theString);
        G4ExcitedString *theStringInCMS=CPExcited(theString);
	G4LorentzRotation toCms=theStringInCMS->TransformToAlignedCms();

	G4bool success=false, inner_sucess=true;
	G4int attempt=0;
	while ( !success && attempt++ < StringLoopInterrupt )
	{
		G4FragmentingString *currentString=new G4FragmentingString(*theStringInCMS);

		std::for_each(LeftVector->begin(), LeftVector->end(), DeleteKineticTrack());
		LeftVector->clear();
		std::for_each(RightVector->begin(), RightVector->end(), DeleteKineticTrack());
		RightVector->clear();
		
		inner_sucess=true;  // set false on failure..
		while (! StopFragmenting(currentString) )
		{  // Split current string into hadron + new string
			G4FragmentingString *newString=0;  // used as output from SplitUp...
			G4KineticTrack * Hadron=Splitup(currentString,newString);
			if ( Hadron != 0 && IsFragmentable(newString)) 
			{
			   if ( currentString->GetDecayDirection() > 0 )
				   LeftVector->push_back(Hadron);
       			   else
	  			   RightVector->push_back(Hadron);
			   delete currentString;
			   currentString=newString;
			} else {
			 // abandon ... start from the beginning
			   if (newString) delete newString;          // Uzhi restore 20.06.08
			   if (Hadron)    delete Hadron;
			   inner_sucess=false;
			   break;
			}
		} 
		// Split current string into 2 final Hadrons
		if ( inner_sucess && 
		     SplitLast(currentString,LeftVector, RightVector) ) 
		{
			success=true;
		}
		delete currentString;
	}
	
	delete theStringInCMS;
	
	if ( ! success )
	{
		std::for_each(LeftVector->begin(), LeftVector->end(), DeleteKineticTrack());
		LeftVector->clear();
		std::for_each(RightVector->begin(), RightVector->end(), DeleteKineticTrack());
		delete RightVector;
		return LeftVector;
	}
		
	// Join Left- and RightVector into LeftVector in correct order.
	while(!RightVector->empty())
	{
	    LeftVector->push_back(RightVector->back());
	    RightVector->erase(RightVector->end()-1);
	}
	delete RightVector;

	CalculateHadronTimePosition(theString.Get4Momentum().mag(), LeftVector);

	G4LorentzRotation toObserverFrame(toCms.inverse());

	for(size_t C1 = 0; C1 < LeftVector->size(); C1++)
	{
	   G4KineticTrack* Hadron = LeftVector->operator[](C1);
	   G4LorentzVector Momentum = Hadron->Get4Momentum();
	   Momentum = toObserverFrame*Momentum;
	   Hadron->Set4Momentum(Momentum);
	   G4LorentzVector Coordinate(Hadron->GetPosition(), Hadron->GetFormationTime());
	   Momentum = toObserverFrame*Coordinate;
	   Hadron->SetFormationTime(Momentum.e());
	   G4ThreeVector aPosition(Momentum.vect());
	   Hadron->SetPosition(theString.GetPosition()+aPosition);
	}
	return LeftVector;
		


}

//----------------------------------------------------------------------------------------------------------

G4double G4QGSMFragmentation::GetLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,  G4ParticleDefinition* , G4double , G4double )
{    
  G4double z;    
  G4double theA(0), d1, d2, yf;
  G4int absCode = std::abs( PartonEncoding );
  if (absCode < 10)
  { 
    if(absCode == 1 || absCode == 2) theA = arho;
    else if(absCode == 3) theA = aphi;
    else throw G4HadronicException(__FILE__, __LINE__, "Unknown PDGencoding in G4QGSMFragmentation::G4LightConeZ");

    do 	
    {
      z  = zmin + G4UniformRand() * (zmax - zmin);
      d1 =  (1. - z);
      d2 =  (alft - theA);
      yf = std::pow(d1, d2);
    } 
    while (G4UniformRand() > yf);
  }
  else
  {       
    if(absCode == 1103 || absCode == 2101 || 
       absCode == 2203 || absCode == 2103)
    {
      d2 =  (alft - (2.*an - arho));
    }
    else if(absCode == 3101 || absCode == 3103 ||
            absCode == 3201 || absCode == 3203)
    {
      d2 =  (alft - (2.*ala - arho));
    }
    else
    {
      d2 =  (alft - (2.*aksi - arho));
    }

    do  
    {
      z = zmin + G4UniformRand() * (zmax - zmin);
      d1 =  (1. - z);
      yf = std::pow(d1, d2);
    } 
    while (G4UniformRand() > yf); 
  }
  return z;
}
//-----------------------------------------------------------------------------------------

G4LorentzVector * G4QGSMFragmentation::SplitEandP(G4ParticleDefinition * pHadron,
                                                  G4FragmentingString * string,    // Uzhi
                                                  G4FragmentingString * ) // Uzhi
{
       G4double HadronMass = pHadron->GetPDGMass();

       // calculate and assign hadron transverse momentum component HadronPx andHadronPy
       G4ThreeVector thePt;
       thePt=SampleQuarkPt();
       G4ThreeVector HadronPt = thePt +string->DecayPt();
       HadronPt.setZ(0);
       //...  sample z to define hadron longitudinal momentum and energy
       //... but first check the available phase space
       G4double DecayQuarkMass2  = sqr(string->GetDecayParton()->GetPDGMass());
       G4double HadronMass2T = sqr(HadronMass) + HadronPt.mag2();
       if (DecayQuarkMass2 + HadronMass2T >= SmoothParam*(string->Mass2()) ) 
          return 0;		// have to start all over!

       //... then compute allowed z region  z_min <= z <= z_max 
 
       G4double zMin = HadronMass2T/(string->Mass2());
       G4double zMax = 1. - DecayQuarkMass2/(string->Mass2());
       if (zMin >= zMax) return 0;		// have to start all over!
	
       G4double z = GetLightConeZ(zMin, zMax,
		       string->GetDecayParton()->GetPDGEncoding(), pHadron,
		       HadronPt.x(), HadronPt.y());      
       
       //... now compute hadron longitudinal momentum and energy
       // longitudinal hadron momentum component HadronPz

        HadronPt.setZ(0.5* string->GetDecayDirection() *
			(z * string->LightConeDecay() - 
			 HadronMass2T/(z * string->LightConeDecay())));
        G4double HadronE  = 0.5* (z * string->LightConeDecay() + 
				  HadronMass2T/(z * string->LightConeDecay()));

       G4LorentzVector * a4Momentum= new G4LorentzVector(HadronPt,HadronE);

       return a4Momentum;
}


//-----------------------------------------------------------------------------------------

G4bool G4QGSMFragmentation::SplitLast(G4FragmentingString * string,
					     G4KineticTrackVector * LeftVector,
    					     G4KineticTrackVector * RightVector)
{
    //... perform last cluster decay
    G4ThreeVector ClusterVel =string->Get4Momentum().boostVector();
    G4double ResidualMass    =string->Mass(); 
    G4double ClusterMassCut = ClusterMass;
    G4int cClusterInterrupt = 0;
    G4ParticleDefinition * LeftHadron, * RightHadron;
    do
    {
        if (cClusterInterrupt++ >= ClusterLoopInterrupt)
        {
          return false;
        }
	G4ParticleDefinition * quark = NULL;
	string->SetLeftPartonStable(); // to query quark contents..
	if (string->DecayIsQuark() && string->StableIsQuark() ) 
	{
	   //... there are quarks on cluster ends
		LeftHadron= QuarkSplitup(string->GetLeftParton(), quark);
	} else {
	   //... there is a Diquark on cluster ends
		G4int IsParticle;
		if ( string->StableIsQuark() ) {
		  IsParticle=(string->GetLeftParton()->GetPDGEncoding()>0) ? -1 : +1; 
		} else {
		  IsParticle=(string->GetLeftParton()->GetPDGEncoding()>0) ? +1 : -1;
		}
      		pDefPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
      		quark = QuarkPair.second;
      		LeftHadron=hadronizer->Build(QuarkPair.first, string->GetLeftParton());
	}
        RightHadron = hadronizer->Build(string->GetRightParton(), quark);

       //... repeat procedure, if mass of cluster is too low to produce hadrons
       //... ClusterMassCut = 0.15*GeV model parameter
	if ( quark->GetParticleSubType()== "quark" ) {ClusterMassCut = 0.;}
	else {ClusterMassCut = ClusterMass;}
    } 
    while (ResidualMass <= LeftHadron->GetPDGMass() + RightHadron->GetPDGMass()  + ClusterMassCut);

    //... compute hadron momenta and energies   
    G4LorentzVector  LeftMom, RightMom;
    G4ThreeVector    Pos;
    Sample4Momentum(&LeftMom, LeftHadron->GetPDGMass(), &RightMom, RightHadron->GetPDGMass(), ResidualMass);
    LeftMom.boost(ClusterVel);
    RightMom.boost(ClusterVel);
    LeftVector->push_back(new G4KineticTrack(LeftHadron, 0, Pos, LeftMom));
    RightVector->push_back(new G4KineticTrack(RightHadron, 0, Pos, RightMom));

    return true;

}

//----------------------------------------------------------------------------------------------------------

G4bool G4QGSMFragmentation::IsFragmentable(const G4FragmentingString * const string)
{
	return sqr(FragmentationMass(string)+MassCut) <
			string->Mass2();
}

//----------------------------------------------------------------------------------------------------------

G4bool G4QGSMFragmentation::StopFragmenting(const G4FragmentingString * const string)
{
	return
         sqr(FragmentationMass(string,&G4HadronBuilder::BuildHighSpin)+MassCut) >
         string->Get4Momentum().mag2();
}

//----------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------

void G4QGSMFragmentation::Sample4Momentum(G4LorentzVector* Mom, G4double Mass, G4LorentzVector* AntiMom, G4double AntiMass, G4double InitialMass) 
    {
    G4double r_val = sqr(InitialMass*InitialMass - Mass*Mass - AntiMass*AntiMass) - sqr(2.*Mass*AntiMass);
    G4double Pabs = (r_val > 0.)? std::sqrt(r_val)/(2.*InitialMass) : 0;

    //... sample unit vector       
    G4double pz = 1. - 2.*G4UniformRand();  
    G4double st     = std::sqrt(1. - pz * pz)*Pabs;
    G4double phi    = 2.*pi*G4UniformRand();
    G4double px = st*std::cos(phi);
    G4double py = st*std::sin(phi);
    pz *= Pabs;
    
    Mom->setPx(px); Mom->setPy(py); Mom->setPz(pz);
    Mom->setE(std::sqrt(Pabs*Pabs + Mass*Mass));

    AntiMom->setPx(-px); AntiMom->setPy(-py); AntiMom->setPz(-pz);
    AntiMom->setE (std::sqrt(Pabs*Pabs + AntiMass*AntiMass));
    }

    
//*********************************************************************************************
