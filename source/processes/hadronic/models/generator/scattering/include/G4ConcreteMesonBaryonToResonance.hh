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
// $Id: G4ConcreteMesonBaryonToResonance.hh,v 1.3 2002/12/12 19:17:39 gunter Exp $ //

#ifndef G4ConcreteMesonBaryonToResonance_h
#define G4ConcreteMesonBaryonToResonance_h

#include "globals.hh"
#include "G4VAnnihilationCollision.hh"
#include "G4VCrossSectionSource.hh"
#include "G4VAngularDistribution.hh"
#include "G4KineticTrackVector.hh"
#include "g4std/vector"
#include "G4XNDeltaTable.hh"
#include "G4ParticleTypeConverter.hh"
#include "G4BaryonWidth.hh"
#include "G4BaryonPartialWidth.hh"

//class G4KineticTrack;

class G4ConcreteMesonBaryonToResonance : public G4VAnnihilationCollision
{

public:

  G4ConcreteMesonBaryonToResonance(const G4ParticleDefinition* aPrimary,
				   const G4ParticleDefinition* bPriamry,
				   const G4ParticleDefinition* aSecondary,
				   const G4String& partWidthLabel);

  virtual ~G4ConcreteMesonBaryonToResonance();

  virtual G4bool IsInCharge(const G4KineticTrack& trk1, 
			    const G4KineticTrack& trk2) const;

  virtual const G4std::vector<G4String>& GetListOfColliders(G4int whichOne) const
  {
    G4Exception("Tried to call G4ConcreteNNToNDelta::GetListOfColliders. Please find out why!");
    G4std::vector<G4String> * aList = new G4std::vector<G4String>;
    return *aList;
  } 
  
  virtual G4String GetName() const
  {
    return "ConcreteMesonBaryonToResonance";
  }

  G4bool operator==(const G4ConcreteMesonBaryonToResonance &right) const;
  G4bool operator!=(const G4ConcreteMesonBaryonToResonance &right) const;


protected:

  virtual const G4VCrossSectionSource* GetCrossSectionSource() const 
    { return crossSectionSource; }

  virtual const G4ParticleDefinition* GetOutgoingParticle(const G4KineticTrack& trk1, 
							  const G4KineticTrack& trk2) const;

private:  

  G4VCrossSectionSource* crossSectionSource;
  const G4ParticleDefinition* thePrimary1;
  const G4ParticleDefinition* thePrimary2;
  const G4ParticleDefinition* theSecondary;
  static G4BaryonWidth theBaryonWidth;
  static G4BaryonPartialWidth theBaryonPartialWidth;
  static G4ParticleTypeConverter myConv;
};

#endif
