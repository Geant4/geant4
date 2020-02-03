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
//
// 
// John Allison  15th February 2019
// An artificial scene to reuse G4VScene code to calculate a bounding extent.

#ifndef G4BOUNDINGEXTENTSCENE_HH
#define G4BOUNDINGEXTENTSCENE_HH

#include "G4PseudoScene.hh"
#include "G4VisExtent.hh"

class G4VModel;

class G4BoundingExtentScene: public G4PseudoScene {

public:

  G4BoundingExtentScene (G4VModel* pModel = 0);

  virtual ~G4BoundingExtentScene ();

  const G4VisExtent& GetExtent () const
  {return fExtent;}
  const G4VisExtent& GetBoundingExtent () const
  {return fExtent;}

  ////////////////////////////////////////////////////////////////
  // The following 2 functions can be used by any code which wishes to
  // accrue a bounding sphere.  Just instantiate a
  // G4BoundingExtentScene and use AccrueBoundingExtent.

  void ResetBoundingExtent ();
  void AccrueBoundingExtent (const G4VisExtent&);

private:

  void ProcessVolume (const G4VSolid& solid);

  G4VModel* fpModel;
  G4VisExtent fExtent;
};

#endif
