#ifndef GhadNuclearGeometry_h
#define GhadNuclearGeometry_h

#include "G4Sphere.hh"
#include "G4NucleusShell.hh"
#include "G4V3DNucleus.hh"
#include "GhadParticles.hh"

class GhadNuclearGeometry
{
  public:
    void Init(G4V3DNucleus * aNuc)
    {
      theNuc = aNuc;
      theShells.clear();
      G4double d=0.99;
      while(d>0.0)
      {
        G4NucleusShell aShell(theNuc, d);
	theShells.push_back(aShell);
	d-=0.05;
      }
    }
    
    G4double DistanceToBoundary(G4double aStep, GhadParticles::iterator aPro)
    {
      if (IsOnBoundary(aPro)) return aStep;
      
      G4double result = aStep;
      G4ThreeVector pos = aPro->GetPosition();
      G4ThreeVector dir = aPro->GetMom().vect().unit();
      std::vector<G4NucleusShell>::iterator i=theShells.begin();
      std::vector<G4NucleusShell>::iterator innterShell=theShells.end();
      G4double dMin= DBL_MAX;

      while(i!=theShells.end() && !i->Contains(pos.mag())) {i++;};
      innterShell = i;

      G4double aStepCandidate = DBL_MAX;
      if(innterShell==theShells.begin())
      {
        // outside
	G4Sphere theShell("", 0, innterShell->GetOuterRadius(), 0, twopi, 0, pi);
	aStepCandidate = theShell.DistanceToIn(pos, dir);
      }
      else if(innterShell==theShells.end())
      {
        // center
	G4Sphere theShell("", 0, (innterShell-1)->GetOuterRadius(), 0, twopi, 0, pi);
	aStepCandidate = theShell.DistanceToOut(pos, dir);
      }
      else
      {
        // between
	G4Sphere theShell("", innterShell->GetOuterRadius(), (innterShell-1)->GetOuterRadius(), 0, twopi, 0, pi);
	aStepCandidate = theShell.DistanceToOut(pos, dir);
      }
      result = std::min(result, aStepCandidate);
      return result;
    }
  
    G4bool IsOnBoundary(GhadParticles::iterator aPro)
    {
      G4bool result = false;
      G4ThreeVector pos = aPro->GetPosition();
      std::vector<G4NucleusShell>::iterator i;
      G4int counter = 0;
      for(i=theShells.begin(); i!= theShells.end(); i++)
      {
        counter++;
	if(i->IsOnShell(pos.mag())) 
	{
	  result = true;
	  break;
	}
      }
      return false;
    }
    
    G4double GetPotentialStep(GhadParticles::iterator thePro)
    {
      G4double result = 0; 
      G4ThreeVector pos = thePro->GetPosition();
      std::vector<G4NucleusShell>::iterator i=theShells.begin();
      while(i!=theShells.end() && !i->IsOnShell(pos.mag())) {i++;};
      G4double outer = 0;
      G4double inner = 0;
      if(i!=theShells.begin()) outer = (i-1)->GetPotential(thePro->GetDefinition());
      if(i!=theShells.end())   inner = i->GetPotential(thePro->GetDefinition());
      result = outer-inner;
      G4cout << "Potential step is "<<result<<" "<<i<<G4endl;
      return result/MeV;
    }
    
  private:
    std::vector<G4NucleusShell> theShells;
    G4V3DNucleus * theNuc;
};

#endif
