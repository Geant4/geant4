//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4AtomicTransitionManager.hh,v 1.2 ????
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Authors: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//          Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
//  
//  16 Sept 2001  First committed to cvs
//
// -------------------------------------------------------------------

#include "G4AtomicDeexcitation.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"

G4AtomicDeexcitation::G4AtomicDeexcitation()
{ }

G4AtomicDeexcitation::~G4AtomicDeexcitation()

{ }

G4std::vector<G4DynamicParticle*>* G4AtomicDeexcitation::GenerateParticles(G4int Z,G4int shellId)
{ 
  G4std::vector<G4DynamicParticle*>* vectorOfParticles = new G4std::vector<G4DynamicParticle*>;
  G4DynamicParticle* aParticle;
  G4int provShellId = 0;
  G4int counter = 0;
  
  // The aim of this loop is to generate more than one fluorecence photon 
  // from the same ionizing event 
  while (provShellId >= 0)
    {
      if (counter == 0) 
	// First call to GenerateParticles(...):
	// shellId is given by the process
	{
	  provShellId = SelectTypeOfTransition(Z, shellId);
	  
	  if  ( provShellId >0)
	    {
	      aParticle = GenerateFluorescence(Z,shellId,provShellId);  
	    }
	  else if ( provShellId ==-1)
	    {
	      aParticle = GenerateAuger(Z, shellId);
	    }
	  else
	    {
	      G4Exception("G4AtomicDeexcitation: starting shell uncorrect: check it");
	    }
	}
      else 
	// Following calls to GenerateParticles(...):
	// newShellId is given by GenerateFluorescence(...)
	{
	  provShellId = SelectTypeOfTransition(Z,newShellId);
	  if  ( provShellId >0)
	    {
	      aParticle = GenerateFluorescence(Z,newShellId,provShellId);
	    }
	  else if ( provShellId ==-1)
	    {
	      aParticle = GenerateAuger(Z, newShellId);
	    }
	  else
	    {
	      G4Exception("G4AtomicDeexcitation: starting shell uncorrect: check it");
	    }
	}
      counter++;
      vectorOfParticles->push_back(aParticle);
    } 
 
  return vectorOfParticles;
}

const G4int G4AtomicDeexcitation::SelectTypeOfTransition(G4int Z, G4int shellId)
{
  if (shellId <=0 ) 
    {G4Exception("G4AtomicDeexcitation: zero or negative shellId");}
  
  G4AtomicTransitionManager*  transitionManager = G4AtomicTransitionManager::Instance();
  G4int provShellId = -1;
  G4int shellNum = 0;
  G4int maxNumOfShells = transitionManager->NumberOfReachableShells(Z);  
  
  const G4AtomicTransition* refShell = transitionManager->ReachableShell(Z,maxNumOfShells-1);

  // This loop gives shellNum the value of the index of shellId
  // in the vector storing the list of the shells reachable through
  // a radiative transition
  if ( shellId <= refShell->FinalShellId())
    {
      while (shellId != transitionManager->ReachableShell(Z,shellNum)->FinalShellId())
	{
	  if(shellNum ==maxNumOfShells-1)
	    {
	      break;
	    }
	  shellNum++;
	}
      G4int transProb = 1;
   
      G4double partialProb = G4UniformRand();      
      G4double partSum = 0;
      const G4AtomicTransition* aShell = transitionManager->ReachableShell(Z,shellNum);      
      G4int trSize =  (aShell->TransitionProbabilities()).size();
    
      // Loop over the shells wich can provide an electron for a 
      // radiative transition towards shellId:
      // in every loop the partial sum of the first transProb shells
      // is calculated and compared with a random number [0,1].
      // If the partial sum is greater the shell whose index transProb
      // is chosen as the starting shell for a radiative transition
      // and its identity is returned
      // Else, terminateded the loop, -1 is returned
      while(transProb < trSize){
	
	 partSum += aShell->TransitionProbability(transProb);

	 if(partialProb <= partSum)
	   {
	     provShellId = aShell->OriginatingShellId(transProb);
	     break;
	   }
	 transProb++;
      }
    }
  else 
    {
      provShellId = -1;

    }
  return provShellId;
}

G4DynamicParticle* G4AtomicDeexcitation::GenerateFluorescence(G4int Z, 
							      G4int shellId,
							      G4int provShellId )
{ 
  G4AtomicTransitionManager*  transitionManager = G4AtomicTransitionManager::Instance();
  //  G4int provenienceShell = provShellId;

  //isotropic angular distribution for the outcoming photon
  G4double newcosTh = 1.-2.*G4UniformRand();
  G4double  newsinTh = sqrt(1.-newcosTh*newcosTh);
  G4double newPhi = twopi*G4UniformRand();
  
  G4double xDir =  newsinTh*sin(newPhi);
  G4double yDir = newsinTh*cos(newPhi);
  G4double zDir = newcosTh;
  
  G4ThreeVector newGammaDirection(xDir,yDir,zDir);
  
  G4int shellNum = 0;
  G4int maxNumOfShells = transitionManager->NumberOfReachableShells(Z);
  
  // find the index of the shell named shellId
  while (shellId != transitionManager->
	 ReachableShell(Z,shellNum)->FinalShellId())
    {
      if(shellNum == maxNumOfShells-1)
	{
	  break;
	}
      shellNum++;
    }
  // number of shell from wich an electron can reach shellId
  size_t transitionSize = transitionManager->
    ReachableShell(Z,shellNum)->OriginatingShellIds().size();
  
  size_t index = 0;
  
  // find the index of the shell named provShellId in the vector
  // storing the shells from which shellId can be reached 
  while (provShellId != transitionManager->
	 ReachableShell(Z,shellNum)->OriginatingShellId(index))
    {
      if(index ==  transitionSize-1)
	{
	  break;
	}
      index++;
    }
  // energy of the gamma leaving provShellId for shellId
  G4double transitionEnergy = transitionManager->
    ReachableShell(Z,shellNum)->TransitionEnergy(index);
  
  // This is the shell where the new vacancy is: it is the same
  // shell where the electron came from
  newShellId = transitionManager->
    ReachableShell(Z,shellNum)->OriginatingShellId(index);
  
  
  G4DynamicParticle* newPart = new G4DynamicParticle(G4Gamma::Gamma(), 
						     newGammaDirection,
						     transitionEnergy);
  return newPart;
}

G4DynamicParticle* G4AtomicDeexcitation::GenerateAuger(G4int Z, G4int shellId)
{
  return 0;
}







