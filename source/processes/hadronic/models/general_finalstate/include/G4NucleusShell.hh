#ifndef G4NucleusShell_h
#define G4NucleusShell_h

#include "globals.hh"
#include "GhadPotential.hh"
#include "G4V3DNucleus.hh"

class G4NucleusShell
{
  public:
    G4NucleusShell(G4V3DNucleus * theNuc, G4double d)
    {
      outerRadius = theNuc->GetNuclearRadius(d);
      thePot = new GhadPotential(theNuc->GetMassNumber(), theNuc->GetCharge());
      Count();
    }
    G4NucleusShell(const G4NucleusShell & aS)
    {
      outerRadius = aS.outerRadius;
      thePot = aS.thePot;
      Count();
    }
    ~G4NucleusShell() {;}
    G4bool Contains(G4double aRad) {return aRad>outerRadius;}
    G4bool IsOnShell(G4double aRad) {return abs(aRad-outerRadius)<thickness;}
    G4double GetPotential(G4ParticleDefinition * aP) 
    {
      return thePot->GetPotential(aP, outerRadius);
    }
    G4double GetOuterRadius() {return outerRadius;}
  
  private:
  
  G4double outerRadius;
  const static G4double thickness = 0.01*fermi;
  GhadPotential * thePot;
  
  void Count();
  static G4int total;
  G4int my;
};

#endif
