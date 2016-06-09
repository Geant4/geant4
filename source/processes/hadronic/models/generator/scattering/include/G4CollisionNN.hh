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
// $Id: G4CollisionNN.hh,v 1.5 2003/06/16 17:07:57 gunter Exp $ //

#ifndef G4CollisionNN_h
#define G4CollisionNN_h

#include "globals.hh"
#include "G4GeneralNNCollision.hh"
#include "G4CollisionVector.hh"
#include "G4VCrossSectionSource.hh"
#include <vector>

class G4KineticTrack;

class G4CollisionNN : public G4GeneralNNCollision
{

public:

  G4CollisionNN();

  virtual ~G4CollisionNN();

  G4bool operator==(const G4CollisionNN &right) const;
  G4bool operator!=(const G4CollisionNN &right) const;

  virtual G4String GetName() const { return "NN CollisionComposite"; }

  virtual G4double CrossSection(const G4KineticTrack& trk1, 
				const G4KineticTrack& trk2) const;

protected:

  virtual const G4VCrossSectionSource* GetCrossSectionSource() const { return crossSectionSource; }
  virtual const G4VAngularDistribution* GetAngularDistribution() const { return 0; }

  virtual const G4CollisionVector* GetComponents() const { return components; } 

  virtual const std::vector<G4String>& GetListOfColliders(G4int whichOne) const;  

private:  

  G4CollisionVector* components;

  G4VCrossSectionSource* crossSectionSource;

  std::vector<G4String> colliders1;
  std::vector<G4String> colliders2;
};

#endif
