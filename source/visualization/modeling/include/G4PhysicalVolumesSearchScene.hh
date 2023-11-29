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
// John Allison  5th September 2018
// Originally G4PhysicalVolumeSearchScene, 10th August 1998, which only
// found the first occurrence of a volume.
// This class (note the extra 's' in the name) produces a vector of all
// occurrences. It can match a physical volume name with the required
// match. The latter can be of the form "/regexp/", where regexp is a
// regular expression (see C++ regex), or a plain string, in which case
// there must be an exact match.

#ifndef G4PHYSICALVOLUMESSEARCHSCENE_HH
#define G4PHYSICALVOLUMESSEARCHSCENE_HH

#include "G4PseudoScene.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeModel.hh"

class G4PhysicalVolumesSearchScene: public G4PseudoScene
{
public:

  G4PhysicalVolumesSearchScene
  (G4PhysicalVolumeModel* pSearchVolumeModel,    // usually a world
   const G4String&        requiredPhysicalVolumeName,
   G4int                  requiredCopyNo = -1); // -1 means any copy no

  virtual ~G4PhysicalVolumesSearchScene () {}

  struct Findings
  {
    Findings
    (G4VPhysicalVolume* pSearchPV,
     G4VPhysicalVolume* pFoundPV,
     G4int foundPVCopyNo = 0,
     G4int foundDepth = 0,
     std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>
     foundBasePVPath =
     std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>(),
     std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>
     foundFullPVPath =
     std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>(),
     G4Transform3D foundObjectTransformation = G4Transform3D())
    : fpSearchPV(pSearchPV)
    , fpFoundPV(pFoundPV)
    , fFoundPVCopyNo(foundPVCopyNo)
    , fFoundDepth(foundDepth)
    , fFoundBasePVPath(foundBasePVPath)
    , fFoundFullPVPath(foundFullPVPath)
    , fFoundObjectTransformation(foundObjectTransformation) {}
    Findings(const G4PhysicalVolumeModel::TouchableProperties& tp)
    : fpSearchPV(nullptr)
    , fpFoundPV(tp.fpTouchablePV)
    , fFoundPVCopyNo(tp.fCopyNo)
    , fFoundDepth(0)
    , fFoundBasePVPath(tp.fTouchableBaseFullPVPath)
    , fFoundFullPVPath(tp.fTouchableFullPVPath)
    , fFoundObjectTransformation(tp.fTouchableGlobalTransform) {}
    G4VPhysicalVolume*   fpSearchPV;   // Searched physical volume.
    G4VPhysicalVolume*   fpFoundPV;    // Found physical volume.
    G4int                fFoundPVCopyNo;  // Found Copy number.
    G4int                fFoundDepth;  // Found depth.
    std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>
    fFoundBasePVPath;    // Base path (e.g., empty for world volume)
    std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>
    fFoundFullPVPath;    // Full path of found volume
    G4Transform3D        fFoundObjectTransformation;  // Found transformation.
  };

  const std::vector<Findings>& GetFindings() const
  {return fFindings;}

private:

  class Matcher {
  public:
    Matcher(const G4String& requiredMatch);
    G4bool Match(const G4String&);
    // Match the string with the required match. The latter can be of the form
    // "/regexp/", where regexp is a regular expression (see C++ regex),
    // or a plain string, in which case there must be an exact match.
  private:
    G4bool   fRegexFlag;  // True if fRequiredMatch is of the form "/.../".
    G4String fRequiredMatch;
  };

  void ProcessVolume(const G4VSolid&);
  
  const G4PhysicalVolumeModel* fpSearchVolumesModel;
  Matcher                      fMatcher;
  G4int                        fRequiredCopyNo;
  std::vector<Findings>        fFindings;
};

#endif
