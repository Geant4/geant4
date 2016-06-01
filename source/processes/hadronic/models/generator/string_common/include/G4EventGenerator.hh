// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EventGenerator.hh,v 1.1 1998/08/22 08:57:14 hpw Exp $
// GEANT4 tag $Name: geant4-00 $
//
#ifndef G4EventGenerator_h
#define G4EventGenerator_h 1

#include "G4HadronicInteraction.hh"
#include "G4VertexCode.hh"
#include "G4InteractionCode.hh"

class G4EventGenerator : public G4HadronicInteraction
{
  public:
      G4EventGenerator();
      ~G4EventGenerator();

  private:
      G4EventGenerator(const G4EventGenerator &right);
      const G4EventGenerator & operator=(const G4EventGenerator &right);
      int operator==(const G4EventGenerator &right) const;
      int operator!=(const G4EventGenerator &right) const;

  public:
      virtual G4double GetWidth(G4VertexCode &theCode) = 0;
      virtual G4double GetCrossSection(G4InteractionCode &theCode) = 0;
      virtual G4ReactionProduct * GetFinalState(G4InteractionCode &theCode) = 0;
     // please also define ApplyYourself for scattering off Hydrogen
};


#endif


