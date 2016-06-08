// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VHighEnergyGenerator.hh,v 1.3 1999/12/15 14:52:37 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
#ifndef G4VHighEnergyGenerator_h
#define G4VHighEnergyGenerator_h 1

#include "G4Nucleus.hh"
#include "G4Track.hh"
class G4KineticTrackVector;
#include "G4ReactionProduct.hh"
#include "G4V3DNucleus.hh"

class G4VHighEnergyGenerator 
{
  public:
      G4VHighEnergyGenerator();
      virtual ~G4VHighEnergyGenerator();

  private:
      G4VHighEnergyGenerator(const G4VHighEnergyGenerator &right);
      const G4VHighEnergyGenerator & operator=(const G4VHighEnergyGenerator &right);
      int operator==(const G4VHighEnergyGenerator &right) const;
      int operator!=(const G4VHighEnergyGenerator &right) const;

 public:
      virtual G4V3DNucleus * GetWoundedNucleus() const = 0;
      virtual G4KineticTrackVector * Scatter(const G4Nucleus &theNucleus, 
                                             const G4DynamicParticle &thePrimary) = 0;

  private:

};
#endif


