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
// $Id: G4ModelingParameters.hh 103627 2017-04-19 13:30:26Z gcosmo $
//
// 
// John Allison  31st December 1997.
//
// Class Description:
//
// Parameters associated with the modeling of GEANT4 objects.

#ifndef G4MODELINGPARAMETERS_HH
#define G4MODELINGPARAMETERS_HH

#include "globals.hh"
#include "G4VisExtent.hh"
#include "G4VisAttributes.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeModel.hh"

#include <vector>
#include <utility>

class G4LogicalVolume;
class G4VisAttributes;
class G4VSolid;
class G4Event;

class G4ModelingParameters {

public: // With description

  // Currently requested drawing style.
  enum DrawingStyle {
    wf,         // Draw edges    - no hidden line removal (wireframe).
    hlr,        // Draw edges    - hidden lines removed.
    hsr,        // Draw surfaces - hidden surfaces removed.
    hlhsr       // Draw surfaces and edges - hidden removed.
  };

  // enums and nested class for communicating a modification to the vis
  // attributes for a specfic touchable defined by PVNameCopyNoPath.
  enum VisAttributesSignifier {
    VASVisibility,
    VASDaughtersInvisible,
    VASColour,
    VASLineStyle,
    VASLineWidth,
    VASForceWireframe,
    VASForceSolid,
    VASForceAuxEdgeVisible,
    VASForceLineSegmentsPerCircle
  };

  class PVNameCopyNo {
  public:
    // Normal constructor
    PVNameCopyNo(G4String name, G4int copyNo)
    : fName(name), fCopyNo(copyNo) {}
    // Constructor from G4PhysicalVolumeModel::G4PhysicalVolumeNodeID
    PVNameCopyNo
    (const G4PhysicalVolumeModel::G4PhysicalVolumeNodeID& nodeID)
    : fName(nodeID.GetPhysicalVolume()->GetName())
    , fCopyNo(nodeID.GetCopyNo()) {}
    const G4String& GetName() const {return fName;}
    G4int GetCopyNo() const {return fCopyNo;}
    G4bool operator!=(const PVNameCopyNo&) const;
    G4bool operator==(const PVNameCopyNo& rhs) const {return !operator!=(rhs);}
  private:
    G4String fName;
    G4int fCopyNo;
  };
  typedef std::vector<PVNameCopyNo> PVNameCopyNoPath;
  typedef PVNameCopyNoPath::const_iterator PVNameCopyNoPathConstIterator;

  class PVPointerCopyNo {
  public:
    // Normal constructor
    PVPointerCopyNo(G4VPhysicalVolume* pPV, G4int copyNo)
    : fpPV(pPV), fCopyNo(copyNo) {}
    // Constructor from G4PhysicalVolumeModel::G4PhysicalVolumeNodeID
    PVPointerCopyNo
    (const G4PhysicalVolumeModel::G4PhysicalVolumeNodeID& nodeID)
    : fpPV(nodeID.GetPhysicalVolume())
    , fCopyNo(nodeID.GetCopyNo()) {}
    const G4String& GetName() const;
    const G4VPhysicalVolume* GetPVPointer() const {return fpPV;}
    G4int GetCopyNo() const {return fCopyNo;}
    G4bool operator!=(const PVPointerCopyNo&) const;
    G4bool operator==(const PVPointerCopyNo& rhs) const {return !operator!=(rhs);}
  private:
    G4VPhysicalVolume* fpPV;
    G4int fCopyNo;
  };
  typedef std::vector<PVPointerCopyNo> PVPointerCopyNoPath;
  typedef PVPointerCopyNoPath::const_iterator PVPointerCopyNoPathConstIterator;

  class VisAttributesModifier {
  public:
    VisAttributesModifier
    (const G4VisAttributes& visAtts,
     VisAttributesSignifier signifier,
     const PVNameCopyNoPath& path):
    fVisAtts(visAtts), fSignifier(signifier), fPVNameCopyNoPath(path) {}
    VisAttributesModifier
    (const G4VisAttributes& visAtts,
     VisAttributesSignifier signifier,
     const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>& path);
    const G4VisAttributes& GetVisAttributes() const
    {return fVisAtts;}
    VisAttributesSignifier GetVisAttributesSignifier() const
    {return fSignifier;}
    const PVNameCopyNoPath& GetPVNameCopyNoPath() const
    {return fPVNameCopyNoPath;}
    void SetVisAttributes(const G4VisAttributes& visAtts)
    {fVisAtts = visAtts;}
    void SetVisAttributesSignifier(VisAttributesSignifier signifier)
    {fSignifier = signifier;}
    void SetPVNameCopyNoPath(const PVNameCopyNoPath& PVNameCopyNoPath)
    {fPVNameCopyNoPath = PVNameCopyNoPath;}
    G4bool operator!=(const VisAttributesModifier&) const;
    G4bool operator==(const VisAttributesModifier& rhs) const
    {return !operator!=(rhs);}
  private:
    G4VisAttributes fVisAtts;
    VisAttributesSignifier fSignifier;
    PVNameCopyNoPath fPVNameCopyNoPath;
  };

  G4ModelingParameters ();

  G4ModelingParameters (const G4VisAttributes* pDefaultVisAttributes,
                        DrawingStyle drawingStyle,
			G4bool isCulling,
			G4bool isCullingInvisible,
			G4bool isDensityCulling,
			G4double visibleDensity,
			G4bool isCullingCovered,
			G4int noOfSides);
  // Culling and clipping policy for G4PhysicalVolumeModel.

  ~G4ModelingParameters ();

  // Note: uses default assignment operator and copy constructor.

  G4bool operator != (const G4ModelingParameters&) const;

  // Get and Is functions...
  G4bool           IsWarning                     () const;
  const G4VisAttributes* GetDefaultVisAttributes () const;
  DrawingStyle     GetDrawingStyle               () const;
  G4bool           IsCulling                     () const;
  G4bool           IsCullingInvisible            () const;
  G4bool           IsDensityCulling              () const;
  G4double         GetVisibleDensity             () const;
  G4bool           IsCullingCovered              () const;
  G4bool           IsExplode                     () const;
  G4double         GetExplodeFactor              () const;
  const G4Point3D& GetExplodeCentre              () const;
  G4int            GetNoOfSides                  () const;
  G4VSolid*        GetSectionSolid               () const;
  G4VSolid*        GetCutawaySolid               () const;
  const G4Event*   GetEvent                      () const;
  const std::vector<VisAttributesModifier>& GetVisAttributesModifiers() const;

  // Set functions...
  void SetWarning              (G4bool);
  void SetDefaultVisAttributes (const G4VisAttributes* pDefaultVisAttributes);
  void SetDrawingStyle         (DrawingStyle);
  void SetCulling              (G4bool);
  void SetCullingInvisible     (G4bool);
  void SetDensityCulling       (G4bool);
  void SetVisibleDensity       (G4double);
  void SetCullingCovered       (G4bool);
  void SetExplodeFactor        (G4double explodeFactor);
  void SetExplodeCentre        (const G4Point3D& explodeCentre);
  G4int SetNoOfSides           (G4int);  // Returns actual number set.
  void SetSectionSolid         (G4VSolid* pSectionSolid);
  void SetCutawaySolid         (G4VSolid* pCutawaySolid);
  void SetEvent                (const G4Event* pEvent);
  void SetVisAttributesModifiers(const std::vector<VisAttributesModifier>&);

  friend std::ostream& operator <<
  (std::ostream& os, const G4ModelingParameters&);
  
  friend std::ostream& operator <<
  (std::ostream& os, const PVNameCopyNoPath&);

  friend std::ostream& operator <<
  (std::ostream& os, const PVPointerCopyNoPath&);

  friend std::ostream& operator <<
  (std::ostream& os,
   const std::vector<VisAttributesModifier>&);

private:

  // Data members...
  G4bool       fWarning;         // Print warnings if true.
  const G4VisAttributes* fpDefaultVisAttributes;
  DrawingStyle fDrawingStyle;    // Drawing style.
  G4bool       fCulling;         // Culling requested.
  G4bool       fCullInvisible;   // Cull (don't Draw) invisible objects.
  G4bool       fDensityCulling;  // Density culling requested.  If so...
  G4double     fVisibleDensity;  // ...density lower than this not drawn.
  G4bool       fCullCovered;     // Cull daughters covered by opaque mothers.
  G4double     fExplodeFactor;   // Explode along radius by this factor...
  G4Point3D    fExplodeCentre;   // ...about this centre.
  G4int        fNoOfSides;       // ...if polygon approximates circle.
  G4VSolid*    fpSectionSolid;   // For generic section (DCUT).
  G4VSolid*    fpCutawaySolid;   // For generic cutaways.
  const G4Event* fpEvent;        // Event being processed.
  std::vector<VisAttributesModifier> fVisAttributesModifiers;
};

std::ostream& operator <<
(std::ostream& os, const G4ModelingParameters&);

std::ostream& operator <<
(std::ostream& os, const G4ModelingParameters::PVNameCopyNoPath&);

std::ostream& operator <<
(std::ostream& os, const G4ModelingParameters::PVPointerCopyNoPath&);

std::ostream& operator <<
(std::ostream& os,
 const std::vector<G4ModelingParameters::VisAttributesModifier>&);

#include "G4ModelingParameters.icc"

#endif
