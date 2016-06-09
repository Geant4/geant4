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
// $Id: G4CollisionNNElastic.hh,v 1.4 2003/06/16 17:07:58 gunter Exp $ //

#ifndef G4CollisionNNElastic_h
#define G4CollisionNNElastic_h

#include "globals.hh"
#include "G4VCollision.hh"
#include "G4VCrossSectionSource.hh"
#include "G4VAngularDistribution.hh"
#include "G4KineticTrackVector.hh"
#include "G4VElasticCollision.hh"
#include <vector>

class G4KineticTrack;

class G4CollisionNNElastic : public G4VElasticCollision
{

public:

  G4CollisionNNElastic();

  virtual ~G4CollisionNNElastic();

  G4bool operator==(const G4CollisionNNElastic &right) const;
  G4bool operator!=(const G4CollisionNNElastic &right) const;

  virtual G4bool IsInCharge(const G4KineticTrack& trk1, 
			    const G4KineticTrack& trk2) const;

  virtual G4String GetName() const { return "NN Elastic Collision"; }
  

protected:

  virtual const G4VCrossSectionSource* GetCrossSectionSource() const { return crossSectionSource; }
  virtual const G4VAngularDistribution* GetAngularDistribution() const { return angularDistribution; }

  virtual const std::vector<G4String>& GetListOfColliders(G4int whichOne) const;  


private:  

  G4VCrossSectionSource* crossSectionSource;
  G4VAngularDistribution* angularDistribution;

  std::vector<G4String> colliders1;
  std::vector<G4String> colliders2;

};

#endif
