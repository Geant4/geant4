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
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4Fancy3DNucleus ----------------
//             by Gunter Folger, May 1998.
//       class for a 3D nucleus, arranging nucleons in space and momentum.
// ------------------------------------------------------------
// 20110805  M. Kelsey -- Remove C-style array (pointer) of G4Nucleons,
//		make vector a container of objects.  Move Helper class
//		to .hh.  Move testSums, places, momentum and fermiM to
//		class data members for reuse.

#include <algorithm>

#include "G4Fancy3DNucleus.hh"
#include "G4Fancy3DNucleusHelper.hh"
#include "G4NuclearFermiDensity.hh"
#include "G4NuclearShellModelDensity.hh"
#include "G4NucleiProperties.hh"
#include "G4HyperNucleiProperties.hh"
#include "G4Nucleon.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4Pow.hh"
#include "G4HadronicException.hh"

#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4RandomDirection.hh"
#include "G4LorentzRotation.hh"
#include "G4RotationMatrix.hh"             
#include "G4PhysicalConstants.hh"

G4Fancy3DNucleus::G4Fancy3DNucleus()
  : myA(0), myZ(0), myL(0), theNucleons(250), currentNucleon(-1), theDensity(0), 
    nucleondistance(0.8*fermi),excitationEnergy(0.),
    places(250), momentum(250), fermiM(250), testSums(250)
{
}

G4Fancy3DNucleus::~G4Fancy3DNucleus()
{
  if(theDensity) delete theDensity;
}

#if defined(NON_INTEGER_A_Z)
void G4Fancy3DNucleus::Init(G4double theA, G4double theZ, G4int numberOfLambdas)
{
  G4int intZ = G4int(theZ);
  G4int intA= ( G4UniformRand()>theA-G4int(theA) ) ? G4int(theA) : G4int(theA)+1;
   // forward to integer Init()
  Init(intA, intZ, std::max(numberOfLambdas, 0));

}
#endif

void G4Fancy3DNucleus::Init(G4int theA, G4int theZ, G4int numberOfLambdas)
{
  currentNucleon=-1;
  theNucleons.clear();
  nucleondistance = 0.8*fermi;
  places.clear();
  momentum.clear();
  fermiM.clear();
  testSums.clear();

  myZ = theZ;
  myA = theA;
  myL = std::max(numberOfLambdas, 0);  // Cannot be negative
  excitationEnergy=0;

  theNucleons.resize(myA);	// Pre-loads vector with empty elements

  // For simplicity, we neglect eventual Lambdas in the nucleus as far as the
  // density of nucler levels and the Fermi level are concerned.
  
  if(theDensity) delete theDensity;
  if ( myA < 17 ) {
     theDensity = new G4NuclearShellModelDensity(myA, myZ);
     if( myA == 12 ) nucleondistance=0.9*fermi; 
  } else {
     theDensity = new G4NuclearFermiDensity(myA, myZ);
  }

  theFermi.Init(myA, myZ);
  
  ChooseNucleons();
  
  ChoosePositions();
  
  if( myA == 12 ) CenterNucleons();   // This would introduce a bias

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
	return (theNucleons.size()>0);
}	

// Returns by pointer; null pointer indicates end of loop
G4Nucleon * G4Fancy3DNucleus::GetNextNucleon()
{
  return ( (currentNucleon>=0 && currentNucleon<myA) ? 
	   &theNucleons[currentNucleon++] : 0 );
}

const std::vector<G4Nucleon> & G4Fancy3DNucleus::GetNucleons()
{
  return theNucleons;
}


// Class-scope function to sort nucleons by Z coordinate
bool G4Fancy3DNucleusHelperForSortInZ(const G4Nucleon& nuc1, const G4Nucleon& nuc2)
{
	return nuc1.GetPosition().z() < nuc2.GetPosition().z();
}

void G4Fancy3DNucleus::SortNucleonsIncZ()
{
  if (theNucleons.size() < 2 ) return;	 	// Avoid unnecesary work

  std::sort(theNucleons.begin(), theNucleons.end(),
	    G4Fancy3DNucleusHelperForSortInZ); 
}

void G4Fancy3DNucleus::SortNucleonsDecZ()
{
  if (theNucleons.size() < 2 ) return;		// Avoid unnecessary work
  SortNucleonsIncZ();

  std::reverse(theNucleons.begin(), theNucleons.end());
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
  if ( myL <= 0 ) return myZ*G4Proton::Proton()->GetPDGMass() + 
                         (myA-myZ)*G4Neutron::Neutron()->GetPDGMass() -
                         BindingEnergy();
  else            return G4HyperNucleiProperties::GetNuclearMass(myA, myZ, myL); 
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
	G4double beta2=theBeta.mag2();
	if (beta2 > 0) {
	   G4double factor=(1-std::sqrt(1-beta2))/beta2; // (gamma-1)/gamma/beta**2
	   G4ThreeVector rprime;
	   for (G4int i=0; i< myA; i++) {
	     rprime = theNucleons[i].GetPosition() - 
	       factor * (theBeta*theNucleons[i].GetPosition()) * theBeta;  
	     theNucleons[i].SetPosition(rprime);
	   }
	}    
}

void G4Fancy3DNucleus::DoLorentzContraction(const G4LorentzVector & theBoost)
{
	if (theBoost.e() !=0 ) {
	   G4ThreeVector beta = theBoost.vect()/theBoost.e();
	   DoLorentzContraction(beta);
	}
}



void G4Fancy3DNucleus::CenterNucleons()
{
	G4ThreeVector	center;
	
	for (G4int i=0; i<myA; i++ )
	{
	   center+=theNucleons[i].GetPosition();
	}   
	center /= -myA;
	DoTranslation(center);
}

void G4Fancy3DNucleus::DoTranslation(const G4ThreeVector & theShift)
{
  G4ThreeVector tempV;
  for (G4int i=0; i<myA; i++ )
    {
      tempV = theNucleons[i].GetPosition() + theShift;
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
  G4int protons=0, nucleons=0, lambdas=0;
  G4double probProton = ( G4double(myZ) )/( G4double(myA) );
  G4double probLambda = myL > 0 ? ( G4double(myL) )/( G4double(myA) ) : 0.0;
  while ( nucleons < myA ) {  /* Loop checking, 30-Oct-2015, G.Folger */
    G4double rnd = G4UniformRand();
    if ( rnd < probProton ) {
      if ( protons < myZ ) {
	protons++;
	theNucleons[nucleons++].SetParticleType(G4Proton::Proton());
      }
    } else if ( rnd < probProton + probLambda ) {
      if ( lambdas < myL ) {
	lambdas++;
	theNucleons[nucleons++].SetParticleType(G4Lambda::Lambda());
      }
    } else {
      if ( (nucleons - protons - lambdas) < (myA - myZ - myL) ) {
	theNucleons[nucleons++].SetParticleType(G4Neutron::Neutron());
      }
    }
  }
  return;
}

void G4Fancy3DNucleus::ChoosePositions()
{
  if( myA != 12) {

    G4int i=0;
    G4ThreeVector	aPos, delta;
    G4bool		freeplace;
    const G4double nd2=sqr(nucleondistance);
    G4double maxR=GetNuclearRadius(0.001);   //  there are no nucleons at a
	                                     //  relative Density of 0.01
    G4int jr=0;
    G4int jx,jy;
    G4double arand[600];
    G4double *prand=arand;
    places.clear();				// Reset data buffer
    G4int interationsLeft=1000*myA;
    while ( (i < myA) && (--interationsLeft>0))  /* Loop checking, 30-Oct-2015, G.Folger */
    {
      do
      {   	
        if ( jr < 3 ) 
	{
          jr=std::min(600,9*(myA - i));
          G4RandFlat::shootArray(jr,prand);
	  //CLHEP::RandFlat::shootArray(jr, prand );
	}
	jx=--jr;
	jy=--jr;
	aPos.set((2*arand[jx]-1.), (2*arand[jy]-1.), (2*arand[--jr]-1.));
      } while (aPos.mag2() > 1. );  /* Loop checking, 30-Oct-2015, G.Folger */
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
	if ( freeplace ) {
          G4double pFermi=theFermi.GetFermiMomentum(theDensity->GetDensity(aPos));
          // protons must at least have binding energy of CoulombBarrier, so
          //  assuming the Fermi energy corresponds to a potential, we must place these such
          //  that the Fermi Energy > CoulombBarrier
	  if (theNucleons[i].GetDefinition() == G4Proton::Proton())
          {
            G4double nucMass = theNucleons[i].GetDefinition()->GetPDGMass();
            G4double eFermi= std::sqrt( sqr(pFermi) + sqr(nucMass) ) - nucMass;
	    if (eFermi <= CoulombBarrier() ) freeplace=false;
	  }
	} 
	if ( freeplace ) {
          theNucleons[i].SetPosition(aPos);
	  places.push_back(aPos);
	  ++i;
	}
      }
    }
    if (interationsLeft<=0) {
      G4Exception("model/util/G4Fancy3DNucleus.cc", "mod_util001", FatalException,
                  "Problem to place nucleons");
    }

  } else {
    // Start insertion
    // Alpha cluster structure of carbon nuclei, C-12, is implemented according to
    //       P. Bozek, W. Broniowski, E.R. Arriola and M. Rybczynski
    //                Phys. Rev. C90, 064902 (2014) 
    const G4double Lbase=3.05*fermi;
    const G4double Disp=0.552;        // 0.91^2*2/3 fermi^2
    const G4double nd2=sqr(nucleondistance);
    const G4ThreeVector Corner1=G4ThreeVector( Lbase/2.,         0., 0.);
    const G4ThreeVector Corner2=G4ThreeVector(-Lbase/2.,         0., 0.);
    const G4ThreeVector Corner3=G4ThreeVector(       0.,Lbase*0.866, 0.);  // 0.866=sqrt(3)/2
    G4ThreeVector R1;
    R1=G4ThreeVector(G4RandGauss::shoot(0.,Disp), G4RandGauss::shoot(0.,Disp), G4RandGauss::shoot(0.,Disp))*fermi + Corner1;
    theNucleons[0].SetPosition(R1); // First nucleon of the first He-4
    G4int loopCounterLeft = 10000;
    for(G4int ii=1; ii<4; ii++)     // 2 - 4 nucleons of the first He-4
    {
      G4bool Continue;
      do
      {
        R1=G4ThreeVector(G4RandGauss::shoot(0.,Disp), G4RandGauss::shoot(0.,Disp), G4RandGauss::shoot(0.,Disp))*fermi + Corner1;
        theNucleons[ii].SetPosition(R1);
        Continue=false;
        for(G4int jj=0; jj < ii; jj++)
        {
          if( (theNucleons[ii].GetPosition() - theNucleons[jj].GetPosition()).mag2() <= nd2 ) {Continue = true; break;}
        }
      } while( Continue && --loopCounterLeft > 0 );  /* Loop checking, 12-Dec-2017, A.Ribon */
    }
    if ( loopCounterLeft <= 0 ) {
      G4Exception("model/util/G4Fancy3DNucleus.cc", "mod_util002", FatalException,
                  "Unable to find a good position for the first alpha cluster");
    }
    loopCounterLeft = 10000;
    for(G4int ii=4; ii<8; ii++)     // 5 - 8 nucleons of the second He-4
    {
      G4bool Continue;
      do
      {
        R1=G4ThreeVector(G4RandGauss::shoot(0.,Disp), G4RandGauss::shoot(0.,Disp), G4RandGauss::shoot(0.,Disp))*fermi + Corner2;
        theNucleons[ii].SetPosition(R1);
        Continue=false;
        for(G4int jj=0; jj < ii; jj++)
        {
          if( (theNucleons[ii].GetPosition() - theNucleons[jj].GetPosition()).mag2() <= nd2 ) {Continue = true; break;}
        }
      } while( Continue && --loopCounterLeft > 0 );  /* Loop checking, 12-Dec-2017, A.Ribon */
    }
    if ( loopCounterLeft <= 0 ) {
      G4Exception("model/util/G4Fancy3DNucleus.cc", "mod_util003", FatalException,
                  "Unable to find a good position for the second alpha cluster");
    }
    loopCounterLeft = 10000;
    for(G4int ii=8; ii<12; ii++)    // 9 - 12 nucleons of the third He-4
    {
      G4bool Continue;
      do
      {
        R1=G4ThreeVector(G4RandGauss::shoot(0.,Disp), G4RandGauss::shoot(0.,Disp), G4RandGauss::shoot(0.,Disp))*fermi + Corner3;
        theNucleons[ii].SetPosition(R1);
        Continue=false;
        for(G4int jj=0; jj < ii; jj++)
        {
          if( (theNucleons[ii].GetPosition() - theNucleons[jj].GetPosition()).mag2() <= nd2 ) {Continue = true; break;}
        }
      } while( Continue && --loopCounterLeft > 0 );  /* Loop checking, 12-Dec-2017, A.Ribon */
    }
    if ( loopCounterLeft <= 0 ) {
      G4Exception("model/util/G4Fancy3DNucleus.cc", "mod_util004", FatalException,
                  "Unable to find a good position for the third alpha cluster");
    }
    G4LorentzRotation RandomRotation;
    RandomRotation.rotateZ(2.*pi*G4UniformRand());
    RandomRotation.rotateY(std::acos(2.*G4UniformRand()-1.));
    // Randomly rotation of the created nucleus
    G4LorentzVector Pos;
    for(G4int ii=0; ii<myA; ii++ )
    {
      Pos=G4LorentzVector(theNucleons[ii].GetPosition(),0.); Pos *=RandomRotation;
      G4ThreeVector NewPos = Pos.vect();
      theNucleons[ii].SetPosition(NewPos);
    }

  }
}

void G4Fancy3DNucleus::ChooseFermiMomenta()
{
    G4int i;
    G4double density;

    // Pre-allocate buffers for filling by index
    momentum.resize(myA, G4ThreeVector(0.,0.,0.));
    fermiM.resize(myA, 0.*GeV);

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
		  while ( mom.mag2() > pmax2 )  /* Loop checking, 30-Oct-2015, G.Folger */
		  {
		      mom=theFermi.GetMomentum(density, fermiM[i]);
		  }
	      }  else
	      {
                  //AR-21Dec2017 : emit a "JustWarning" exception instead of writing on the error stream.
	          //G4cerr << "G4Fancy3DNucleus: difficulty finding proton momentum" << G4endl;
                  G4ExceptionDescription ed;
                  ed << "Nucleus Z A " << myZ << " " << myA << G4endl;
                  ed << "proton with eMax=" << eMax << G4endl;
                  G4Exception( "G4Fancy3DNucleus::ChooseFermiMomenta(): difficulty finding proton momentum, set it to (0,0,0)",
                               "HAD_FANCY3DNUCLEUS_001", JustWarning, ed );
		  mom=G4ThreeVector(0,0,0);
	      }

	   }
	   momentum[i]= mom;
	}

	if ( ReduceSum() ) break;
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
       // GF 11-05-2011: set BindingEnergy to be T of Nucleon with p , ~ p**2/2m
       //theNucleons[i].SetBindingEnergy(
       //     0.5*sqr(fermiM[i])/theNucleons[i].GetParticleType()->GetPDGMass());
    }
}


G4bool G4Fancy3DNucleus::ReduceSum()
{
	G4ThreeVector sum;
	G4double PFermi=fermiM[myA-1];

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
	testSums.clear();
	testSums.resize(myA-1);		// Allocate block for filling below

	G4ThreeVector delta;
	for (G4int aNucleon=0; aNucleon < myA-1; ++aNucleon) {
	  delta = 2.*((momentum[aNucleon]*testDir)*testDir);

	  testSums[aNucleon].Fill(delta, delta.mag(), aNucleon);
	}

	std::sort(testSums.begin(), testSums.end());

//    reduce Momentum Sum until the next would be allowed.
	G4int index=(G4int)testSums.size();
	while ( (sum-testSums[--index].Vector).mag()>PFermi && index>0)  /* Loop checking, 30-Oct-2015, G.Folger */
	{
	  // Only take one which improve, ie. don't change sign and overshoot...
	  if ( sum.mag() > (sum-testSums[index].Vector).mag() ) {
	    momentum[testSums[index].Index]-=testSums[index].Vector;
	    sum-=testSums[index].Vector;
	  }
	}

	if ( (sum-testSums[index].Vector).mag() <= PFermi )
	{
		G4int best=-1;
		G4double pBest=2*PFermi; // anything larger than PFermi
		for ( G4int aNucleon=0; aNucleon<=index; aNucleon++)
		{
			// find the momentum closest to choosen momentum for last Nucleon.
			G4double pTry=(testSums[aNucleon].Vector-sum).mag();
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
		momentum[testSums[best].Index]-=testSums[best].Vector;
		momentum[myA-1]=testSums[best].Vector-sum;

		return true;

	}

	// try to compensate momentum using another Nucleon....
	G4int swapit=-1;
	while (swapit<myA-1)  /* Loop checking, 30-Oct-2015, G.Folger */
	{
	  if ( fermiM[++swapit] > PFermi ) break;
	}
	if (swapit == myA-1 ) return false;

	// Now we have a nucleon with a bigger Fermi Momentum.
	// Exchange with last nucleon.. and iterate.
	std::swap(theNucleons[swapit], theNucleons[myA-1]);
	std::swap(momentum[swapit], momentum[myA-1]);
	std::swap(fermiM[swapit], fermiM[myA-1]);
	return ReduceSum();
}

G4double G4Fancy3DNucleus::CoulombBarrier()
{
  static const G4double cfactor = (1.44/1.14) * MeV;
  return cfactor*myZ/(1.0 + G4Pow::GetInstance()->Z13(myA));
}
