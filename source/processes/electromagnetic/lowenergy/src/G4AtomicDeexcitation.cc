#include "g4std/fstream"
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
  while (provShellId >= 0)
    {
      if (counter == 0)
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

const G4int G4AtomicDeexcitation::SelectTypeOfTransition(G4int Z, G4int shellId){
   
  G4AtomicTransitionManager*  transitionManager = G4AtomicTransitionManager::Instance();
  G4int provShellId = -1;
  G4int shellNum = 0;
  G4int maxNumOfShells = transitionManager->NumberOfReachableShells(Z);
  
  if (shellId <=0 ) 
    {G4Exception("G4AtomicDeexcitation: zero or negative shellId");}
  
  const G4AtomicTransition* refShell = 
    transitionManager->ReachableShell(Z,maxNumOfShells-1);
 
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
      
      const G4AtomicTransition* aShell = 
	transitionManager->ReachableShell(Z,shellNum);
      
      G4int trSize =  (aShell->TransitionProbabilities()).size();
    
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
//questo significa che la vacanza e' su una shell piu' esterna della piu'
      //esterna tra le raggiungibili
    }
  return provShellId;
}

G4DynamicParticle* G4AtomicDeexcitation::GenerateFluorescence(G4int Z, 
							      G4int shellId,
							      G4int provShellId )
{ 
  G4AtomicTransitionManager*  transitionManager = 
    G4AtomicTransitionManager::Instance();
  G4int provenienceShell=0;
  provenienceShell = provShellId;
  //isotropic angular distribution
  
  G4double newcosTh = 1-2*G4UniformRand();
  G4double  newsinTh = sqrt(1-newcosTh*newcosTh);
  G4double newPhi = twopi*G4UniformRand();
  G4double xDir =  newsinTh*sin(newPhi);
  
  G4double yDir = newsinTh*cos(newPhi);
  
  G4double zDir = newcosTh;
  
    G4ThreeVector newGammaDirection(xDir,yDir,zDir);
    
      G4int shellNum = 0;
      G4int maxNumOfShells = transitionManager->NumberOfReachableShells(Z);
      
      //     const G4AtomicTransition* refShell = 
      //transitionManager->ReachableShell(Z,maxNumOfShells-1);
      
      while (shellId != transitionManager->ReachableShell(Z,shellNum)->FinalShellId())
	{
	  if(shellNum == maxNumOfShells-1)
	    {
	      break;
	    }
	  shellNum++;
	}
      size_t transitionSize = transitionManager->ReachableShell(Z,shellNum)->OriginatingShellIds().size();
      
      size_t index = 0;

      while (provShellId !=transitionManager->ReachableShell(Z,shellNum)->OriginatingShellId(index))
	{
	  if(index ==  transitionSize-1)
	    {
	      break;
	    }
	  index++;
	}
      G4double transitionEnergy = transitionManager->ReachableShell(Z,shellNum)->TransitionEnergy(index)*MeV;
     
      newShellId = transitionManager->ReachableShell(Z,shellNum)->OriginatingShellId(index);
    
      G4DynamicParticle* newPart = new G4DynamicParticle(G4Gamma::Gamma(), 
							 newGammaDirection,
							 transitionEnergy);
      return newPart;
      /* }
    else{
      return 0;
      G4cout<<" Ho restituito uno zero"<<G4endl;
      //return 0 if the atom is in the fundamental state
    }
      */
}

G4DynamicParticle* G4AtomicDeexcitation::GenerateAuger(G4int Z, G4int shellId)
{
  return 0;
}







