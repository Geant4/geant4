//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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
void establish_G4MT_TLS_G4ConcreteNNTwoBodyResonance(const G4ParticleDefinition* aPrimary,
		       const G4ParticleDefinition* bPriamry,
		       const G4ParticleDefinition* aSecondary,
		       const G4ParticleDefinition* bSecondary,
		       const G4VXResonanceTable& sigmaTable);

  G4ConcreteNNTwoBodyResonance(void *s1, void *s2, void *s3, void *s4, void *s5, void *s6, void *s7);

  virtual ~G4ConcreteNNTwoBodyResonance();

  virtual G4bool IsInCharge(const G4KineticTrack& trk1, 
			    const G4KineticTrack& trk2) const;

  virtual const std::vector<G4String>& GetListOfColliders(G4int ) const
  {
    throw G4HadronicException(__FILE__, __LINE__, "Tried to call G4ConcreteNNTwoBodyResonance::GetListOfColliders. Please find out why!");
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
  G4ConcreteNNTwoBodyResonance(const G4ConcreteNNTwoBodyResonance &);
  G4ConcreteNNTwoBodyResonance & operator= (const G4ConcreteNNTwoBodyResonance &);

private:  

  G4VCrossSectionSource* crossSectionSource;
  const G4ParticleDefinition* thePrimary1;
  const G4ParticleDefinition* thePrimary2;
  std::vector<const G4ParticleDefinition*> theOutGoing;

};

#endif
