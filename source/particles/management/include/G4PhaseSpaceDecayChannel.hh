// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhaseSpaceDecayChannel.hh,v 1.2 1999-10-28 23:24:13 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      27 July 1996 H.Kurashige
//      30 May 1997 H.Kurashige
// ------------------------------------------------------------
#ifndef G4PhaseSpaceDecayChannel_h
#define G4PhaseSpaceDecayChannel_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDecayChannel.hh"

class G4PhaseSpaceDecayChannel :public G4VDecayChannel
{
  public:  // With Description
    //Constructors 
      G4PhaseSpaceDecayChannel(G4int Verbose = 1);
      G4PhaseSpaceDecayChannel(const G4String& theParentName,
			       G4double        theBR,
			       G4int           theNumberOfDaughters,
			       const G4String& theDaughterName1,
			       const G4String& theDaughterName2 = "",
			       const G4String& theDaughterName3 = "",
			       const G4String& theDaughterName4 = ""   );

  public: 
   //  Destructor
      virtual ~G4PhaseSpaceDecayChannel();

  public:  // With Description
     virtual G4DecayProducts *DecayIt(G4double);   

  public: 
     static G4double Pmx(G4double e, G4double p1, G4double p2);

  private: 
     G4DecayProducts *OneBodyDecayIt();
     G4DecayProducts *TwoBodyDecayIt();
     G4DecayProducts *ThreeBodyDecayIt();
     G4DecayProducts *ManyBodyDecayIt();
};  

inline
 G4double G4PhaseSpaceDecayChannel::Pmx(G4double e, G4double p1, G4double p2)
{
   // calcurate momentum of daughter particles in two-body decay
   G4double ppp = (e+p1+p2)*(e+p1-p2)*(e-p1+p2)*(e-p1-p2)/(4.0*e*e);
   if (ppp>0) return sqrt(ppp);
   else       return -1.;
}

#endif
