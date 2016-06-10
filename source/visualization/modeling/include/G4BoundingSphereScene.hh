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
// $Id: G4BoundingSphereScene.hh 81056 2014-05-20 09:02:16Z gcosmo $
//
// 
// John Allison  7th June 1997
// An artificial scene to reuse G4VScene code to calculate a bounding sphere.

#ifndef G4BOUNDINGSPHERESCENE_HH
#define G4BOUNDINGSPHERESCENE_HH

#include "G4PseudoScene.hh"
#include "G4VisExtent.hh"

class G4VModel;

class G4BoundingSphereScene: public G4PseudoScene {

public:

  G4BoundingSphereScene (G4VModel* pModel = 0);

  virtual ~G4BoundingSphereScene ();

  G4VisExtent GetBoundingSphereExtent ();

  const G4Point3D& GetCentre() const {return fCentre;}

  G4double GetRadius() const {return fRadius;}

  void SetCentre(const G4Point3D& centre) {fCentre = centre;}

  ////////////////////////////////////////////////////////////////
  // The following 2 functions can be used by any code which wishes to
  // accrue a bounding sphere.  Just instantiate a
  // G4BoundingSphereScene and use AccrueBoundingSphere.

  void ResetBoundingSphere ();
  void AccrueBoundingSphere (const G4Point3D& centre,
			     G4double radius);

private:

  void ProcessVolume (const G4VSolid& solid);

  G4VModel* fpModel;  // Instantiating code may optionally set this.
  G4Point3D fCentre;
  G4double fRadius;
};

#endif
