// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Fancy3DNucleus.cc,v 1.7 2000/01/27 12:42:13 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
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
#include "g4rw/tvvector.h"

G4Fancy3DNucleus::G4Fancy3DNucleus()
 : nucleondistance(0.8*fermi)
{
	theDensity=NULL;
	theNucleons=NULL;
	currentNucleon=-1;
	myA=0;
	myZ=0;
//G4cout <<"G4Fancy3DNucleus::G4Fancy3DNucleus()"<<G4endl;
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
//  G4cout << "G4Fancy3DNucleus::Init(theA, theZ) called"<<G4endl;
  currentNucleon=-1;
  if(theNucleons!=NULL) delete [] theNucleons;

  theRWNucleons.clear();

  myZ = G4int(theZ);
  myA= ( G4UniformRand()>theA-G4int(theA) ) ? G4int(theA) : G4int(theA)+1;

  theNucleons = new G4Nucleon[myA];
  
//  G4cout << "myA, myZ" << myA << ", " << myZ << G4endl;

  if(theDensity!=NULL) delete theDensity;
  if ( myA < 17 ) {
     theDensity = new G4NuclearShellModelDensity(myA, myZ);
  } else {
     theDensity = new G4NuclearFermiDensity(myA, myZ);
  }

  theFermi.Init(myA, myZ);
  
  ChooseNucleons();
  
  ChoosePositions();
  
//  CenterNucleons();   // This would introduce a bias 

  ChooseFermiMomenta();
  
  G4double Ebinding= BindingEnergy()/myA;
  
  for (G4int aNucleon=0; aNucleon < myA; aNucleon++)
  {
	theNucleons[aNucleon].SetBindingEnergy(Ebinding);
  }
  
  
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


const G4RWTPtrOrderedVector<G4Nucleon> & G4Fancy3DNucleus::GetNucleons()
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
	  else G4cout << "G4Fancy3DNucleus::ChooseNucleons not efficient" << G4endl;
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
    
    G4double * fermiM=new G4double[myA];
	
    for (G4int ntry=0; ntry<1 ; ntry ++ )
    {	
	for (i=0; i < myA; i++ )    // momenta for all, including last, in case we swa
	{
	   fermiM[i]=theFermi.GetFermiMomentum(theDensity->GetDensity(theNucleons[i].GetPosition()));
	   momentum[i]= theFermi.GetMomentum(theDensity->GetDensity(theNucleons[i].GetPosition()));
	}
	
	if (ReduceSum(momentum,fermiM) ) 
	  break;
//       G4cout <<" G4FancyNucleus: iterating to find momenta: "<< ntry<< G4endl;
    }

//     G4ThreeVector sum;
//     for (G4int index=0; index<myA;sum+=momentum[index++])
//     ;
//     cout << "final sum / mag() " << sum << " / " << sum.mag() << G4endl;
    
    
    G4double energy;
    for ( i=0; i< myA ; i++ )
    {
       energy=theNucleons[i].GetParticleType()->GetPDGMass() 
	      + BindingEnergy()/myA;
       G4LorentzVector tempV(momentum[i],energy);
       theNucleons[i].SetMomentum(tempV);
    }

    delete [] momentum;
    delete [] fermiM;
}

  class G4Fancy3DNucleusHelper // Helper class 
  {
    public:
    	G4Fancy3DNucleusHelper(const G4ThreeVector &vec,const G4double size,const G4int index)
    	: Vector(vec), Size(size), anInt(index) {}
    	int operator ==(const G4Fancy3DNucleusHelper &right) const
    	{
    		return this==&right;
    	}
    	int operator < (const G4Fancy3DNucleusHelper &right) const
    	{
    		return this->Size<right.Size;
    	}
    	const G4ThreeVector& vector() const
    	{
    		return Vector;
    	}
    	const G4double size() const
    	{
    		return Size;
    	}
    	const G4int index() const
    	{
    		return anInt;
    	}
    private:
    	G4Fancy3DNucleusHelper operator =(const G4Fancy3DNucleusHelper &right) const 
    	{
    		G4cout <<" G4Fancy3DNucleus::G4Fancy3DNucleusHelper op = called------------------" << G4endl;
		return G4Fancy3DNucleusHelper();
    	}
    	G4Fancy3DNucleusHelper(): Vector(0), Size(0), anInt(0) {G4cout << "def ctor for MixMasch" << G4endl;}
	const G4ThreeVector Vector;
	const G4double Size;
	const G4int anInt;
  };

G4bool G4Fancy3DNucleus::ReduceSum(G4ThreeVector * momentum, G4double *pFermiM)
{
	G4ThreeVector sum;
	G4double PFermi=theFermi.GetFermiMomentum(theDensity->GetDensity(theNucleons[myA-1].GetPosition()));
	
	for (G4int i=0; i < myA-1 ; i++ )
	     { sum+=momentum[i]; }

// check if have to do anything at all..
	if ( sum.mag() <= PFermi ) 
	{
 		momentum[myA-1]=-sum;
		return true;
	}
	
// find all possible changes in momentum, changing only the component parallel to sum
	G4ThreeVector testDir=sum.unit();
	G4RWTPtrSortedVector<G4Fancy3DNucleusHelper> testSums(myA-1);		// Sorted on delta.mag()
	for ( G4int aNucleon=0; aNucleon < myA-1; aNucleon++){
		G4ThreeVector delta=2*((momentum[aNucleon]*testDir)* testDir);
		testSums.insert(new G4Fancy3DNucleusHelper(delta,delta.mag(),aNucleon));		
	}
	
//    reduce Momentum Sum until the next would be allowed.
	G4int index=testSums.entries();
	while ( (sum-testSums[--index]->vector()).mag()>PFermi && index>0) 
	{
		// Only take one which improve, ie. don't change sign and overshoot...
		if ( sum.mag() > (sum-testSums[index]->vector()).mag() ) {
// 		   momentum[testSums[index]->index()]-=testSums[index]->vector();
// 		   sum-=testSums[index]->vector();
		}
	}

	if ( (sum-testSums[index]->vector()).mag() <= PFermi ) 
	{
		G4int best=-1;
		G4double pBest=2*PFermi; // anything larger than PFermi
		for ( G4int aNucleon=0; aNucleon<=index; aNucleon++)
		{
			// find the momentum closest to choosen momentum for last Nucleon.
			G4double pTry=(testSums[aNucleon]->vector()-sum).mag();
			if ( pTry < PFermi 
			 &&  abs(momentum[myA-1].mag() - pTry ) < pBest )
		        {
			   pBest=abs(momentum[myA-1].mag() - pTry );	 
			   best=aNucleon;
			}
		}
		if ( best < 0 )  G4Exception( " Logic error in Fancy3DNucleus");
// 		momentum[testSums[best]->index()]-=testSums[best]->vector();
// 		momentum[myA-1]=testSums[best]->vector()-sum;

		testSums.clearAndDestroy();
		return true;
		
	} 
	testSums.clearAndDestroy();

	// try to compensate momentum using another Nucleon....
	
	G4int swapit=-1;
	while (swapit< myA-1
	   && theFermi.GetFermiMomentum(theDensity->GetDensity(theNucleons[++swapit].GetPosition())) < PFermi )
	   ; 
	
	if (swapit == myA-1 ) return false;
	
	// Now we have a nucleon with a bigger Fermi Momentum.
	// Exchange with last nucleon.. and iterate.
// 	G4cout << " Nucleon to swap with : " << swapit << G4endl;
// 	G4cout << " Fermi momentum test, and better.. " << PFermi << " / "
// 	    << theFermi.GetFermiMomentum(theDensity->GetDensity(theNucleons[swapit].GetPosition())) << G4endl;
//	cout << theNucleons[swapit]<< G4endl << theNucleons[myA-1] << G4endl;
//	cout << momentum[swapit] << G4endl << momentum[myA-1] << G4endl;
	G4Nucleon swap= theNucleons[swapit];
	G4ThreeVector mom_swap=momentum[swapit];
	G4double pf=pFermiM[swapit];
	theNucleons[swapit]=theNucleons[myA-1];
	momentum[swapit]=momentum[myA-1];
	pFermiM[swapit]=pFermiM[myA-1];
	theNucleons[myA-1]=swap;
	momentum[myA-1]=mom_swap;
	pFermiM[myA-1]=pf;
//	cout << "after swap" <<G4endl<< theNucleons[swapit] << G4endl << theNucleons[myA-1] << G4endl;
//	cout << momentum[swapit] << G4endl << momentum[myA-1] << G4endl;
	return ReduceSum(momentum,pFermiM);	
}
