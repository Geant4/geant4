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

#ifndef G4ConcreteNNTwoBodyResonance_HH
#define G4ConcreteNNTwoBodyResonance_HH

#include "globals.hh"
#include "G4VScatteringCollision.hh"
#include "G4VCrossSectionSource.hh"
#include "G4VAngularDistribution.hh"
#include "G4KineticTrackVector.hh"
#include <vector>
#include "G4XDeltaDeltaTable.hh"

class G4ConcreteNNTwoBodyResonance : public G4VScatteringCollision
{

public:

  G4ConcreteNNTwoBodyResonance(const G4ParticleDefinition* aPrimary,
		       const G4ParticleDefinition* bPriamry,
		       const G4ParticleDefinition* aSecondary,
		       const G4ParticleDefinition* bSecondary,
		       const G4VXResonanceTable& sigmaTable);

  virtual ~G4ConcreteNNTwoBodyResonance();

  virtual G4bool IsInCharge(const G4KineticTrack& trk1, 
			    const G4KineticTrack& trk2) const;

  virtual const std::vector<G4String>& GetListOfColliders(G4int ) const
  {
    G4Exception("Tried to call G4ConcreteNNTwoBodyResonance::GetListOfColliders. Please find out why!");
    std::vector<G4String> * aList = new std::vector<G4String>;
    return *aList;
  } 
  
  virtual G4String GetName() const
  {
    return "G4ConcreteNNTwoBodyResonance";
  }

  G4bool operator==(const G4ConcreteNNTwoBodyResonance &right) const;
  G4bool operator!=(const G4ConcreteNNTwoBodyResonance &right) const;


protected:

  virtual const G4VCrossSectionSource* GetCrossSectionSource() const { return crossSectionSource; }

  virtual const std::vector<const G4ParticleDefinition*> & GetOutgoingParticles() const
  {
    return theOutGoing;
  }

private:  

  G4VCrossSectionSource* crossSectionSource;
  const G4ParticleDefinition* thePrimary1;
  const G4ParticleDefinition* thePrimary2;
  std::vector<const G4ParticleDefinition*> theOutGoing;

};

#endif
