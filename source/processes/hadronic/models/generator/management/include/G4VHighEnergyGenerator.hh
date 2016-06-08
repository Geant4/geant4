// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VHighEnergyGenerator.hh,v 1.4 2000/12/14 09:36:03 hpw Exp $
// GEANT4 tag $Name: geant4-03-01 $
//
#ifndef G4VHighEnergyGenerator_h
#define G4VHighEnergyGenerator_h 1

// Class Description
// Base class for high energy interaction models in geant4. By merit of inheriting
// from this class a high energy interaction model can be used in conjunction with
// any cascade, precompound model and evaporation phase in the
// generation of complete final states for inelastic scattering.
// Class Description - End

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


