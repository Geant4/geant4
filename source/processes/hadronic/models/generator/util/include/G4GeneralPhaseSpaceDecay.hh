// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GeneralPhaseSpaceDecay.hh,v 1.3 1999/12/15 14:52:50 gunter Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
//$Id: G4GeneralPhaseSpaceDecay.hh,v 1.1 1997/05/21
// ----------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, A. Feliciello, 20th May 1998
//
//      Note: this class is a generalization of the 
//            G4PhaseSpaceDecayChannel one
// ----------------------------------------------------------------
#ifndef G4GeneralPhaseSpaceDecay_h
#define G4GeneralPhaseSpaceDecay_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDecayChannel.hh"

class G4GeneralPhaseSpaceDecay : public G4VDecayChannel
{
  public:
    //Constructors 
      G4GeneralPhaseSpaceDecay(G4int Verbose = 1);

      G4GeneralPhaseSpaceDecay(const G4String& theParentName,
			       G4double        theBR,
			       G4int           theNumberOfDaughters,
			       const G4String& theDaughterName1,
			       const G4String& theDaughterName2 = "",
			       const G4String& theDaughterName3 = "");

      G4GeneralPhaseSpaceDecay(const G4String& theParentName,
                               G4double        theParentMass,
			       G4double        theBR,
			       G4int           theNumberOfDaughters,
			       const G4String& theDaughterName1,
			       const G4String& theDaughterName2 = "",
			       const G4String& theDaughterName3 = "");

    //  Destructor
      virtual ~G4GeneralPhaseSpaceDecay();

  public:
     G4double GetParentMass() const;
     void SetParentMass(const G4double aParentMass);
     virtual G4DecayProducts* DecayIt(G4double mass=0.0);   
     static G4double Pmx(G4double e, G4double p1, G4double p2);

  protected:
     G4DecayProducts* OneBodyDecayIt();
     G4DecayProducts* TwoBodyDecayIt();
     G4DecayProducts* ThreeBodyDecayIt();
     G4DecayProducts* ManyBodyDecayIt();
     
  private:
     G4double parentmass;
    
};  



inline G4double G4GeneralPhaseSpaceDecay::GetParentMass() const
{
  return parentmass;
}

inline void G4GeneralPhaseSpaceDecay::SetParentMass(const G4double aParentMass)
{
  parentmass = aParentMass;
}



inline
 G4double G4GeneralPhaseSpaceDecay::Pmx(G4double e, G4double p1, G4double p2)
{
   // calculate momentum of daughter particles in two-body decay
   G4double ppp = (e+p1+p2)*(e+p1-p2)*(e-p1+p2)*(e-p1-p2)/(4.0*e*e);
   if (ppp>0) return sqrt(ppp);
   else       return -1.;
}

#endif
