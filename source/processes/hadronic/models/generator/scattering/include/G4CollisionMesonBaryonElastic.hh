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
// $Id: G4CollisionMesonBaryonElastic.hh,v 1.8 2003/06/16 17:07:55 gunter Exp $ //

#ifndef G4CollisionMesonBaryonElastic_h
#define G4CollisionMesonBaryonElastic_h

#include "globals.hh"
#include "G4VCollision.hh"
#include "G4VCrossSectionSource.hh"
#include "G4VAngularDistribution.hh"
#include "G4KineticTrackVector.hh"
#include <vector>
#include "G4VElasticCollision.hh"

class G4KineticTrack;

class G4CollisionMesonBaryonElastic : public G4VElasticCollision
{

public:

  G4CollisionMesonBaryonElastic();
  virtual ~G4CollisionMesonBaryonElastic();

  G4bool operator==(const G4CollisionMesonBaryonElastic &right) const;
  G4bool operator!=(const G4CollisionMesonBaryonElastic &right) const;
  
  virtual G4String GetName() const;
  virtual G4bool IsInCharge(const G4KineticTrack& trk1, 
			    const G4KineticTrack& trk2) const;
protected:

  virtual const G4VCrossSectionSource* GetCrossSectionSource() const 
  { return crossSectionSource; }
  virtual const G4VAngularDistribution* GetAngularDistribution() const 
  { return angularDistribution; }
  virtual const std::vector<G4String>& GetListOfColliders(G4int ) const;

private:  

  G4VCrossSectionSource* crossSectionSource;
  G4VAngularDistribution* angularDistribution;

  std::vector<G4String> dummy;


};

#endif
