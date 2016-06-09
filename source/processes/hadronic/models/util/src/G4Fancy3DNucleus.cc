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
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4Fancy3DNucleus ----------------
//             by Gunter Folger, May 1998.
//       class for a 3D nucleus, arranging nucleons in space and momentum.
// ------------------------------------------------------------

#include "G4Fancy3DNucleus.hh"
#include "G4NuclearFermiDensity.hh"
#include "G4NuclearShellModelDensity.hh"
#include "G4NucleiProperties.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <algorithm>
#include "G4HadronicException.hh"


G4Fancy3DNucleus::G4Fancy3DNucleus()
 : nucleondistance(0.8*fermi)
{
	theDensity=0;
	theNucleons=0;
	currentNucleon=-1;
	myA=0;
	myZ=0;
//G4cout <<"G4Fancy3DNucleus::G4Fancy3DNucleus()"<<G4endl;
}

G4Fancy3DNucleus::~G4Fancy3DNucleus()
{
  if(theNucleons) delete [] theNucleons;
  if(theDensity) delete theDensity;
}

#if defined(NON_INTEGER_A_Z)
void G4Fancy3DNucleus::Init(G4double theA, G4double theZ)
{
  G4int intZ = G4int(theZ);
  G4int intA= ( G4UniformRand()>theA-G4int(theA) ) ? G4int(theA) : G4int(theA)+1;
   // forward to integer Init()
  Init(intA, intZ);

}
#endif

void G4Fancy3DNucleus::Init(G4int theA, G4int theZ)
{
//  G4cout << "G4Fancy3DNucleus::Init(theA, theZ) called"<<G4endl;
  currentNucleon=-1;
  if(theNucleons) delete [] theNucleons;

  theRWNucleons.clear();

  myZ = theZ;
  myA= theA;

  theNucleons = new G4Nucleon[myA];
  
//  G4cout << "myA, myZ" << myA << ", " << myZ << G4endl;

  if(theDensity) delete theDensity;
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
	return theNucleons;
}	

G4Nucleon * G4Fancy3DNucleus::GetNextNucleon()
{
  return ( currentNucleon>=0 && currentNucleon<myA ) ? 
			theNucleons+currentNucleon++  : 0;
}

const std::vector<G4Nucleon *> & G4Fancy3DNucleus::GetNucleons()
{
	if ( theRWNucleons.size()==0 )
	{
	    for (G4int i=0; i< myA; i++)
	    {
	        theRWNucleons.push_back(theNucleons+i);
	    }
	}
	return theRWNucleons;
}

//void G4Fancy3DNucleus::SortNucleonsIncZ() // on increased Z-coordinates Uzhi 29.08.08

   bool G4Fancy3DNucleusHelperForSortInZ(const G4Nucleon* nuc1, const G4Nucleon* nuc2)
{
	return nuc1->GetPosition().z() < nuc2->GetPosition().z();
}    

//void G4Fancy3DNucleus::SortNucleonsInZ()
void G4Fancy3DNucleus::SortNucleonsIncZ() // on increased Z-coordinates Uzhi 29.08.08
{
	
	GetNucleons();   // make sure theRWNucleons is initialised

	if (theRWNucleons.size() < 2 ) return; 

	std::sort( theRWNucleons.begin(),theRWNucleons.end(),G4Fancy3DNucleusHelperForSortInZ); 

// now copy sorted nucleons to theNucleons array. TheRWNucleons are pointers in theNucleons
//  so we need to copy to new, and then swap. 
        G4Nucleon * sortedNucleons = new G4Nucleon[myA];
	for ( unsigned int i=0; i<theRWNucleons.size(); i++ )
	{
	   sortedNucleons[i]= *(theRWNucleons[i]);
	}

	theRWNucleons.clear();   // about to delete array these point to....
	delete [] theNucleons;
	
	theNucleons=sortedNucleons;

	return;
}

void G4Fancy3DNucleus::SortNucleonsDecZ() // on decreased Z-coordinates Uzhi 29.08.08
{
        G4Nucleon * sortedNucleons = new G4Nucleon[myA];
	
	GetNucleons();   // make sure theRWNucleons is initialised

	if (theRWNucleons.size() < 2 ) return; 
	std::sort( theRWNucleons.begin(),theRWNucleons.end(),G4Fancy3DNucleusHelperForSortInZ); 

// now copy sorted nucleons to theNucleons array. TheRWNucleons are pointers in theNucleons
//  so we need to copy to new, and then swap. 
	for ( unsigned int i=0; i<theRWNucleons.size(); i++ )
	{
	   sortedNucleons[i]= *(theRWNucleons[myA-1-i]);  // Uzhi 29.08.08
	}
	theRWNucleons.clear();   // about to delete array elements these point to....
	delete [] theNucleons;
	theNucleons=sortedNucleons;

	return;
}

G4double G4Fancy3DNucleus::BindingEnergy()
{
  return G4NucleiProperties::GetBindingEnergy(myA,myZ);
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
	return std::sqrt(maxradius2)+nucleondistance;
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
	G4double factor=(1-std::sqrt(1-theBeta.mag2()))/theBeta.mag2(); // (gamma-1)/gamma/beta**2
	for (G4int i=0; i< myA; i++)
	{
	     G4ThreeVector rprime=theNucleons[i].GetPosition() - 
	         factor * (theBeta*theNucleons[i].GetPosition()) * 
	       // theNucleons[i].GetPosition();
	       theBeta;  
	     theNucleons[i].SetPosition(rprime);
	}    
}

void G4Fancy3DNucleus::DoLorentzContraction(const G4LorentzVector & theBoost)
{
	G4ThreeVector beta= 1/theBoost.e() * theBoost.vect();
	// DoLorentzBoost(beta); 
	DoLorentzContraction(beta);
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

const G4VNuclearDensity * G4Fancy3DNucleus::GetNuclearDensity() const
{
	return theDensity;
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
	G4ThreeVector	aPos, delta;
        std::vector<G4ThreeVector> places;
	places.reserve(myA);
	G4bool		freeplace;
	static G4double nd2 = sqr(nucleondistance);
	G4double maxR=GetNuclearRadius(0.001);   //  there are no nucleons at a
	                                        //  relative Density of 0.01
	G4int jr=0;
	G4int jx,jy;
	G4double arand[600];
	G4double *prand=arand;
	while ( i < myA )
	{
	   do
	   {   	
		if ( jr < 3 ) 
		{
		    jr=std::min(600,9*(myA - i));
		    CLHEP::RandFlat::shootArray(jr, prand );
		}
		jx=--jr;
		jy=--jr;
		aPos=G4ThreeVector( (2*arand[jx]-1.),
	   				   (2*arand[jy]-1.),
					   (2*arand[--jr]-1.));
	   } while (aPos.mag2() > 1. );
	   aPos *=maxR;
	   G4double density=theDensity->GetRelativeDensity(aPos);
	   if (G4UniformRand() < density)
	   {
	      freeplace= true;
	      std::vector<G4ThreeVector>::iterator iplace;
	      for( iplace=places.begin(); iplace!=places.end() && freeplace;++iplace)
	      {
	        delta = *iplace - aPos;
		freeplace= delta.mag2() > nd2;
	      }
	      
	      if ( freeplace )
	      {
		  G4double pFermi=theFermi.GetFermiMomentum(theDensity->GetDensity(aPos));
		    // protons must at least have binding energy of CoulombBarrier, so
		    //  assuming the Fermi energy corresponds to a potential, we must place these such
		    //  that the Fermi Energy > CoulombBarrier
		  if (theNucleons[i].GetDefinition() == G4Proton::Proton())
		  {
	             G4double eFermi= std::sqrt( sqr(pFermi) + sqr(theNucleons[i].GetDefinition()->GetPDGMass()) )
		                      - theNucleons[i].GetDefinition()->GetPDGMass();
	             if (eFermi <= CoulombBarrier() ) freeplace=false;
		  }
	      }
	      if ( freeplace )
	      {
		  theNucleons[i].SetPosition(aPos);
		  places.push_back(aPos);
		  ++i;
	      }
	   }
	}

}

void G4Fancy3DNucleus::ChooseFermiMomenta()
{
    G4int i;
    G4double density;
    G4ThreeVector * momentum=new G4ThreeVector[myA];

    G4double * fermiM=new G4double[myA];

    for (G4int ntry=0; ntry<1 ; ntry ++ )
    {
	for (i=0; i < myA; i++ )    // momenta for all, including last, in case we swap nucleons
	{
	   density = theDensity->GetDensity(theNucleons[i].GetPosition());
	   fermiM[i] = theFermi.GetFermiMomentum(density);
	   G4ThreeVector mom=theFermi.GetMomentum(density);
	   if (theNucleons[i].GetDefinition() == G4Proton::Proton())
	   {
	      G4double eMax = std::sqrt(sqr(fermiM[i]) +sqr(theNucleons[i].GetDefinition()->GetPDGMass()) )
	                      - CoulombBarrier();
	      if ( eMax > theNucleons[i].GetDefinition()->GetPDGMass() )
	      {
	          G4double pmax2= sqr(eMax) - sqr(theNucleons[i].GetDefinition()->GetPDGMass());
		  fermiM[i] = std::sqrt(pmax2);
		  while ( mom.mag2() > pmax2 )
		  {
		      mom=theFermi.GetMomentum(density, fermiM[i]);
		  }
	      }  else
	      {
	          G4cerr << "G4Fancy3DNucleus: difficulty finding proton momentum" << G4endl;
		  mom=G4ThreeVector(0,0,0);
	      }

	   }
	   momentum[i]= mom;
	}

	if (ReduceSum(momentum,fermiM) )
	  break;
//       G4cout <<" G4FancyNucleus: iterating to find momenta: "<< ntry<< G4endl;
    }

//     G4ThreeVector sum;
//     for (G4int index=0; index<myA;sum+=momentum[index++])
//     ;
//     G4cout << "final sum / mag() " << sum << " / " << sum.mag() << G4endl;

    G4double energy;
    for ( i=0; i< myA ; i++ )
    {
       energy = theNucleons[i].GetParticleType()->GetPDGMass()
	        - BindingEnergy()/myA;
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
    		return size()<right.size();
    	}
    	const G4ThreeVector& vector() const
    	{
    		return Vector;
    	}
    	G4double size() const
    	{
    		return Size;
    	}
    	G4int index() const
    	{
    		return anInt;
    	}
    	G4Fancy3DNucleusHelper operator =(const G4Fancy3DNucleusHelper &right)
    	{
	  Vector = right.Vector;
	  Size = right.Size;
	  anInt = right.anInt;
	  return *this;
    	}

    private:
    	G4Fancy3DNucleusHelper(): Vector(0), Size(0), anInt(0) {G4cout << "def ctor for MixMasch" << G4endl;}
	G4ThreeVector Vector;
	G4double Size;
	G4int anInt;
  };

G4bool G4Fancy3DNucleus::ReduceSum(G4ThreeVector * momentum, G4double *pFermiM)
{
	G4ThreeVector sum;
	G4double PFermi=pFermiM[myA-1];

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
	std::vector<G4Fancy3DNucleusHelper> testSums;		// Sorted on delta.mag()

	for ( G4int aNucleon=0; aNucleon < myA-1; aNucleon++){
		G4ThreeVector delta=2*((momentum[aNucleon]*testDir)* testDir);
		testSums.push_back(G4Fancy3DNucleusHelper(delta,delta.mag(),aNucleon));
	}
	std::sort(testSums.begin(), testSums.end());

//    reduce Momentum Sum until the next would be allowed.
	G4int index=testSums.size();
	while ( (sum-testSums[--index].vector()).mag()>PFermi && index>0)
	{
		// Only take one which improve, ie. don't change sign and overshoot...
		if ( sum.mag() > (sum-testSums[index].vector()).mag() ) {
		   momentum[testSums[index].index()]-=testSums[index].vector();
		   sum-=testSums[index].vector();
		}
	}

	if ( (sum-testSums[index].vector()).mag() <= PFermi )
	{
		G4int best=-1;
		G4double pBest=2*PFermi; // anything larger than PFermi
		for ( G4int aNucleon=0; aNucleon<=index; aNucleon++)
		{
			// find the momentum closest to choosen momentum for last Nucleon.
			G4double pTry=(testSums[aNucleon].vector()-sum).mag();
			if ( pTry < PFermi
			 &&  std::abs(momentum[myA-1].mag() - pTry ) < pBest )
		        {
			   pBest=std::abs(momentum[myA-1].mag() - pTry );
			   best=aNucleon;
			}
		}
		if ( best < 0 )  
		{
		  G4String text = "G4Fancy3DNucleus.cc: Logic error in ReduceSum()";
  	          throw G4HadronicException(__FILE__, __LINE__, text);
		}
		momentum[testSums[best].index()]-=testSums[best].vector();
		momentum[myA-1]=testSums[best].vector()-sum;

		testSums.clear();
		return true;

	}
	testSums.clear();

	// try to compensate momentum using another Nucleon....

	G4int swapit=-1;
	while (swapit<myA-1)
	{
	  if ( pFermiM[++swapit] > PFermi ) break;
	}
	if (swapit == myA-1 ) return false;

	// Now we have a nucleon with a bigger Fermi Momentum.
	// Exchange with last nucleon.. and iterate.
	G4Nucleon swap= theNucleons[swapit];
	G4ThreeVector mom_swap=momentum[swapit];
	G4double pf=pFermiM[swapit];
	theNucleons[swapit]=theNucleons[myA-1];
	momentum[swapit]=momentum[myA-1];
	pFermiM[swapit]=pFermiM[myA-1];
	theNucleons[myA-1]=swap;
	momentum[myA-1]=mom_swap;
	pFermiM[myA-1]=pf;
	return ReduceSum(momentum,pFermiM);
}

G4double G4Fancy3DNucleus::CoulombBarrier()
{
  G4double coulombBarrier = (1.44/1.14) * MeV * myZ / (1.0 + std::pow(G4double(myA),1./3.));
  return coulombBarrier;
}
