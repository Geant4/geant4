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
// $Id: G4ConcreteMesonBaryonToResonance.hh,v 1.3 2006-06-29 20:34:47 gunter Exp $ //

#ifndef G4ConcreteMesonBaryonToResonance_h
#define G4ConcreteMesonBaryonToResonance_h

#include "globals.hh"
#include "G4VAnnihilationCollision.hh"
#include "G4VCrossSectionSource.hh"
#include "G4VAngularDistribution.hh"
#include "G4KineticTrackVector.hh"
#include <vector>
#include "G4XNDeltaTable.hh"
#include "G4ParticleTypeConverter.hh"
#include "G4BaryonWidth.hh"
#include "G4BaryonPartialWidth.hh"
#include "G4Threading.hh"

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

  virtual const std::vector<G4String>& GetListOfColliders(G4int ) const
  {
    throw G4HadronicException(__FILE__, __LINE__, "Tried to call G4ConcreteMesonBaryonToResonance::GetListOfColliders. Please find out why!");
    std::vector<G4String> * aList = new std::vector<G4String>;
    return *aList;
  } 
  
  virtual G4String GetName() const
  {
    return "ConcreteMesonBaryonToResonance";
  }

  G4bool operator==(const G4ConcreteMesonBaryonToResonance &right) const;
  G4bool operator!=(const G4ConcreteMesonBaryonToResonance &right) const;

private:
  G4ConcreteMesonBaryonToResonance(const G4ConcreteMesonBaryonToResonance &);
  G4ConcreteMesonBaryonToResonance & operator= (const G4ConcreteMesonBaryonToResonance &);

protected:

  virtual const G4VCrossSectionSource* GetCrossSectionSource() const 
    { return crossSectionSource; }

  virtual const G4ParticleDefinition* GetOutgoingParticle(const G4KineticTrack& trk1, 
							  const G4KineticTrack& trk2) const;

private:  

  static void InitialisePointers();

  G4VCrossSectionSource* crossSectionSource;
  const G4ParticleDefinition* thePrimary1;
  const G4ParticleDefinition* thePrimary2;
  const G4ParticleDefinition* theSecondary;

  static G4BaryonWidth & theBaryonWidth();
  static G4BaryonPartialWidth & theBaryonPartialWidth();
  static G4ParticleTypeConverter & myConv();

  static G4BaryonWidth*           baryonWidth;
  static G4BaryonPartialWidth*    baryonPartialWidth;
  static G4ParticleTypeConverter* particleTypeConverter;

#ifdef G4MULTITHREADED
  static G4Mutex concreteMesonBaryonToResonanceMutex;
#endif
};

#endif
