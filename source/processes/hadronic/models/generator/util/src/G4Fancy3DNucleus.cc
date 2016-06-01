// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Fancy3DNucleus.cc,v 1.5 1998/11/19 17:10:15 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4Fancy3DNucleus ----------------
//             by Gunter Folger, May 1998.
//       class for a 3D nucleus, arranging nucleons in space and momentum.
// ------------------------------------------------------------

#include "G4Fancy3DNucleus.hh"
#include "G4NuclearFermiDensity.hh"
#include "G4NuclearShellModelDensity.hh"
#include "G4NucleiPropertiesTable.hh"
#include "Randomize.hh"
#include "G4ios.hh"

G4Fancy3DNucleus::G4Fancy3DNucleus()
 : nucleondistance(0.8*fermi)
{
	theDensity=NULL;
	theNucleons=NULL;
	currentNucleon=-1;
	myA=0;
	myZ=0;
//G4cout <<"G4Fancy3DNucleus::G4Fancy3DNucleus()"<<endl;
}

/* No use for these
 *
 *G4Fancy3DNucleus::G4Fancy3DNucleus(const G4Fancy3DNucleus &right)
 * : nucleondistance(0.8*fermi)  {}
 *const G4Fancy3DNucleus & G4Fancy3DNucleus::operator=(const G4Fancy3DNucleus &right)
 *{
 *}
 *
 *
 *int G4Fancy3DNucleus::operator==(const G4Fancy3DNucleus &right) const
 *{
 *}
 *
 *int G4Fancy3DNucleus::operator!=(const G4Fancy3DNucleus &right) const
 *{
 *}
 *
 */

G4Fancy3DNucleus::~G4Fancy3DNucleus()
{
  if(theNucleons!=NULL) delete [] theNucleons;
  if(theDensity!=NULL) delete theDensity;
}


void G4Fancy3DNucleus::Init(G4double theA, G4double theZ)
{
  G4int i;
//  G4cout << "G4Fancy3DNucleus::Init(theA, theZ) called"<<endl;
  currentNucleon=-1;
  if(theNucleons!=NULL) delete [] theNucleons;

  theRWNucleons.clear();

  myZ = theZ;
  myA= ( G4UniformRand()>theA-G4int(theA) ) ? G4int(theA) : G4int(theA)+1;

  theNucleons = new G4Nucleon[myA];
  
//  G4cout << "myA, myZ" << myA << ", " << myZ << endl;

  if(theDensity!=NULL) delete theDensity;
  if ( myA < 17 ) {
     theDensity = new G4NuclearShellModelDensity(myA, myZ);
  } else {
     theDensity = new G4NuclearFermiDensity(myA, myZ);
  }

  theFermi.Init(myA, myZ);
  
  ChooseNucleons();
  
  ChoosePositions();
  
//  CenterNucleons();

  ChooseFermiMomenta();
  
  return;
}

G4bool G4Fancy3DNucleus::StartLoop()
{
	currentNucleon=0;
	return theNucleons != NULL;
}	

G4Nucleon * G4Fancy3DNucleus::GetNextNucleon()
{
  return ( currentNucleon>=0 && currentNucleon<myA ) ? 
			theNucleons+currentNucleon++  : NULL;
}


const RWTPtrOrderedVector<G4Nucleon> & G4Fancy3DNucleus::GetNucleons()
{
	if ( theRWNucleons.isEmpty() )
	{
	    for (G4int i=0; i< myA; i++)
	    {
	        theRWNucleons.append(theNucleons+i);
	    }
	}
	return theRWNucleons;
}

G4double G4Fancy3DNucleus::BindingEnergy()
{
	return G4NucleiPropertiesTable::GetBindingEnergy(myZ,myA);
}

G4double G4Fancy3DNucleus::GetNuclearRadius()
{
	return GetNuclearRadius(0.5);
}

G4double G4Fancy3DNucleus::GetNuclearRadius(const G4double maxRelativeDensity)
{	
	return theDensity->GetRadius(maxRelativeDensity);
}

G4double G4Fancy3DNucleus::GetOuterRadius()
{
	G4double maxradius2=0;
	
	for (int i=0; i<myA; i++)
	{
	   if ( theNucleons[i].GetPosition().mag2() > maxradius2 )
	   {
	   	maxradius2=theNucleons[i].GetPosition().mag2();
	   }
	}
	return sqrt(maxradius2)+nucleondistance;
}

G4double G4Fancy3DNucleus::GetMass()
{
        return   myZ*G4Proton::Proton()->GetPDGMass() + 
                 (myA-myZ)*G4Neutron::Neutron()->GetPDGMass() -
                 BindingEnergy();
}



void G4Fancy3DNucleus::DoLorentzBoost(const G4LorentzVector & theBoost)
{
	for (G4int i=0; i<myA; i++){
	    theNucleons[i].Boost(theBoost);
	}
}

void G4Fancy3DNucleus::DoLorentzBoost(const G4ThreeVector & theBeta)
{
	for (G4int i=0; i<myA; i++){
	    theNucleons[i].Boost(theBeta);
	}
}

void G4Fancy3DNucleus::DoLorentzContraction(const G4ThreeVector & theBeta)
{
	G4double factor=(1-sqrt(1-theBeta.mag2()))/theBeta.mag2(); // (gamma-1)/gamma/beta**2
	for (G4int i=0; i< myA; i++)
	{
	     G4ThreeVector rprime=theNucleons[i].GetPosition() - 
	         factor * (theBeta*theNucleons[i].GetPosition()) * 
	         theNucleons[i].GetPosition();
	     theNucleons[i].SetPosition(rprime);
	}    
}

void G4Fancy3DNucleus::DoLorentzContraction(const G4LorentzVector & theBoost)
{
	G4ThreeVector beta= 1/theBoost.e() * theBoost.vect();
	DoLorentzBoost(beta);
}



void G4Fancy3DNucleus::CenterNucleons()
{
	G4ThreeVector	center;
	
	for (G4int i=0; i<myA; i++ )
	{
	   center+=theNucleons[i].GetPosition();
	}   
	center *= -1./myA;
	DoTranslation(center);
}

void G4Fancy3DNucleus::DoTranslation(const G4ThreeVector & theShift)
{
	for (G4int i=0; i<myA; i++ )
	{
	   G4ThreeVector tempV = theNucleons[i].GetPosition() + theShift;
	   theNucleons[i].SetPosition(tempV);
	}   
}

//----------------------- private Implementation Methods-------------

void G4Fancy3DNucleus::ChooseNucleons()
{
	G4int protons=0,nucleons=0;
	
	while (nucleons < myA )
	{
	  if ( protons < myZ && G4UniformRand() < (G4double)(myZ-protons)/(G4double)(myA-nucleons) )
	  {
	     protons++;
	     theNucleons[nucleons++].SetParticleType(G4Proton::Proton());
	  }
	  else if ( (nucleons-protons) < (myA-myZ) )
	  {
	       theNucleons[nucleons++].SetParticleType(G4Neutron::Neutron());
	  }
	  else G4cout << "G4Fancy3DNucleus::ChooseNucleons not efficient" << endl;
	}
	return;
}

void G4Fancy3DNucleus::ChoosePositions()
{
	G4int i=0;
	G4ThreeVector	aPos,center;
	G4bool		freeplace;
	G4double maxR=GetNuclearRadius(0.01);   //  there are no nucleons at a 
	                                        //  relative Density of 0.01
	
	while ( i < myA )
	{
	   do 
	   {   aPos=G4ThreeVector( (2*G4UniformRand()-1.),
	   				   (2*G4UniformRand()-1.),
					   (2*G4UniformRand()-1.));
	   } while (aPos.mag2() > 1. );
	   aPos *=maxR;
	   if (theDensity->GetRelativeDensity(aPos) >  G4UniformRand() )
	   {
	      freeplace= true;
	      for( int j=0; j<i && freeplace; j++) 
	      {
	        freeplace= freeplace && 
	              (theNucleons[j].GetPosition()-aPos).mag() > nucleondistance;
	      }
	      if ( freeplace ) theNucleons[i++].SetPosition(aPos);
	   }
	}
	 
}

void G4Fancy3DNucleus::ChooseFermiMomenta()
{
	G4int i;
	G4ThreeVector * momentum=new G4ThreeVector[myA];
	G4ThreeVector sum,sum2;

   for (G4int hardtry=0; hardtry<50; hardtry++)
   {
   for (G4int ntry=0; ntry<50 ; ntry ++ )
      {
      if ( ntry > 49 ) G4cout << "G4Nucleus: Difficulties finding nucleon momenta, number of try " << ntry<< endl;
	
	for (i=0; i < myA-1; i++ )    // momenta for all but the last nucleon
	{
	   momentum[i]= theFermi.GetMomentum(theDensity->GetDensity(theNucleons[i].GetPosition()));
	   sum+=momentum[i];
	}
	G4int best;
	G4double testsum;
	do 
	{
	   sum=0;
	   for (i=0; i < myA-1 ; i++ )
	   	{ sum+=momentum[i]; }
//	   G4cout << "Momenum Sum " << sum.mag() << endl;
	   G4ThreeVector test;
	   testsum=sum.mag() * (1-perCent); // forces improvement!
	   if ( testsum > theFermi.GetFermiMomentum(theDensity->GetDensity(theNucleons[myA-1].GetPosition())) ) 
	   {
	      best=-1;
	      for (i=0; i < myA-1 ;i++) 
	      {   
	         test=-1. * sum.unit() * momentum[i].mag();
	         if ( (sum+test-momentum[i]).mag() < testsum )
	         {
		   best=i;
		   testsum=(sum+test-momentum[i]).mag();	      
		 }
		 if ( testsum < theFermi.GetFermiMomentum(theDensity->GetDensity(theNucleons[myA-1].GetPosition())))
		 {  break; }
	      }
	      if (best != -1 )
	      {
		  momentum[best]=-1.*sum.unit() * momentum[best].mag();
	      } else
	      {
	   	 G4Nucleon swap= theNucleons[ntry%(myA-1)];
	   	 theNucleons[ntry%(myA-1)]=theNucleons[myA-1];
	   	 theNucleons[myA-1]=swap;
	      }
	   }   
	} while ( best != -1 && 
	          testsum > theFermi.GetFermiMomentum(theDensity->GetDensity(theNucleons[myA-1].GetPosition())));
	
	if ( best != -1 )  // Success 
	{
	   for (sum=0, i=0; i < myA-1 ; i++ ) sum+=momentum[i]; 
	   ;
	   momentum[myA-1]= -1 * sum;
	   break;
	}
     }
     ChoosePositions();
     }

	G4double energy;
	for ( i=0; i< myA ; i++ )
	{
	   energy=theNucleons[i].GetParticleType()->GetPDGMass() 
	          + BindingEnergy()/myA;
	   G4LorentzVector tempV(momentum[i],energy);
	   theNucleons[i].SetMomentum(tempV);
	}
	   
	delete [] momentum;
}
