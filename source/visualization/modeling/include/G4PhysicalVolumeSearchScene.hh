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
// $Id: G4PhysicalVolumeSearchScene.hh 81056 2014-05-20 09:02:16Z gcosmo $
//
// 
// John Allison  10th August 1998.
// An artificial scene to find physical volumes.

#ifndef G4PHYSICALVOLUMESEARCHSCENE_HH
#define G4PHYSICALVOLUMESEARCHSCENE_HH

#include "G4PseudoScene.hh"
#include "G4VisExtent.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeModel.hh"

class G4PhysicalVolumeSearchScene: public G4PseudoScene {

public:

  G4PhysicalVolumeSearchScene
  (G4PhysicalVolumeModel*,
   const G4String& requiredPhysicalVolumeName,
   G4int requiredCopyNo,
   G4int verbosity = 99);

  virtual ~G4PhysicalVolumeSearchScene ();

  const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>&
  GetFoundFullPVPath() const
  {return fFoundFullPVPath;}

  G4int GetFoundDepth() const {return fFoundDepth;}

  G4VPhysicalVolume* GetFoundVolume() const {return fpFoundPV;}

  const G4Transform3D&  GetFoundTransformation () const
  {return fFoundObjectTransformation;}

private:

  void ProcessVolume(const G4VSolid&);
  
  const G4PhysicalVolumeModel* fpPVModel;
  G4String             fRequiredPhysicalVolumeName;
  G4int                fRequiredCopyNo;
  std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>
                       fFoundFullPVPath;            // Found full path.
  G4int                fFoundDepth;                 // Found depth.
  G4VPhysicalVolume*   fpFoundPV;                   // Found physical volume.
  G4LogicalVolume*     fpFoundLV;                   // Found logical volume.
  G4Transform3D        fFoundObjectTransformation;  // Found transformation.
  G4int                fVerbosity;
  G4bool               fMultipleOccurrence;
};

#endif
