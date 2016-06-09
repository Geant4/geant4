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
// $Id: G4CollisionNN.hh,v 1.2 2006-06-29 20:32:55 gunter Exp $ //

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

private:
  G4CollisionNN(const G4CollisionNN &);
  G4CollisionNN & operator= (const G4CollisionNN &);


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
