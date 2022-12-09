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
// class G4PVPlacement Implementation
//
// ----------------------------------------------------------------------

#include "G4PVPlacement.hh"
#include "G4AffineTransform.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"

// ----------------------------------------------------------------------
// Constructor
//
G4PVPlacement::G4PVPlacement( G4RotationMatrix* pRot,
                        const G4ThreeVector& tlate,
                        const G4String& pName,
                              G4LogicalVolume* pLogical,
                              G4VPhysicalVolume* pMother,
                              G4bool pMany,
                              G4int pCopyNo,
                              G4bool pSurfChk )
  : G4VPhysicalVolume(pRot, tlate, pName, pLogical, pMother),
    fmany(pMany), fcopyNo(pCopyNo)
{
  if (pMother)
  {
    G4LogicalVolume* motherLogical = pMother->GetLogicalVolume();
    if (pLogical == motherLogical)
    {
      G4Exception("G4PVPlacement::G4PVPlacement()", "GeomVol0002",
                  FatalException, "Cannot place a volume inside itself!");
    }
    SetMotherLogical(motherLogical);
    motherLogical->AddDaughter(this);
    if (pSurfChk) { CheckOverlaps(); }
  }
}

// ----------------------------------------------------------------------
// Constructor
//
G4PVPlacement::G4PVPlacement( const G4Transform3D& Transform3D,
                              const G4String& pName,
                                    G4LogicalVolume* pLogical,
                                    G4VPhysicalVolume* pMother,
                                    G4bool pMany,
                                    G4int pCopyNo,
                                    G4bool pSurfChk )
  : G4VPhysicalVolume(NewPtrRotMatrix(Transform3D.getRotation().inverse()),
                      Transform3D.getTranslation(), pName, pLogical, pMother),
    fmany(pMany), fcopyNo(pCopyNo)
{
  fallocatedRotM = (GetRotation() != 0);
  if (pMother)
  {
    G4LogicalVolume* motherLogical = pMother->GetLogicalVolume();
    if (pLogical == motherLogical)
      G4Exception("G4PVPlacement::G4PVPlacement()", "GeomVol0002",
                  FatalException, "Cannot place a volume inside itself!");
    SetMotherLogical(motherLogical);
    motherLogical->AddDaughter(this);
    if (pSurfChk) { CheckOverlaps(); }
  }
}

// ----------------------------------------------------------------------
// Constructor
//
// The logical volume of the mother is utilised (not the physical)
//
G4PVPlacement::G4PVPlacement( G4RotationMatrix* pRot,
                        const G4ThreeVector& tlate,
                              G4LogicalVolume* pCurrentLogical,
                        const G4String& pName,
                              G4LogicalVolume* pMotherLogical,
                              G4bool pMany,
                              G4int pCopyNo,
                              G4bool pSurfChk )
  : G4VPhysicalVolume(pRot, tlate, pName, pCurrentLogical, nullptr),
    fmany(pMany), fcopyNo(pCopyNo)
{
  if (pCurrentLogical == pMotherLogical)
  {
    G4Exception("G4PVPlacement::G4PVPlacement()", "GeomVol0002",
                FatalException, "Cannot place a volume inside itself!");
  }
  SetMotherLogical(pMotherLogical);
  if (pMotherLogical) { pMotherLogical->AddDaughter(this); }
  if ((pSurfChk) && (pMotherLogical)) { CheckOverlaps(); }
}


// ----------------------------------------------------------------------
// Constructor
//
G4PVPlacement::G4PVPlacement( const G4Transform3D& Transform3D,
                                    G4LogicalVolume* pCurrentLogical,
                              const G4String& pName,
                                    G4LogicalVolume* pMotherLogical,
                                    G4bool pMany,
                                    G4int pCopyNo,
                                    G4bool pSurfChk )
  : G4VPhysicalVolume(nullptr, Transform3D.getTranslation(),
                      pName, pCurrentLogical, nullptr),
    fmany(pMany), fcopyNo(pCopyNo)
{
  if (pCurrentLogical == pMotherLogical)
  {
    G4Exception("G4PVPlacement::G4PVPlacement()", "GeomVol0002",
                FatalException, "Cannot place a volume inside itself!");
  }
  SetRotation( NewPtrRotMatrix(Transform3D.getRotation().inverse()) );
  fallocatedRotM = (GetRotation() != nullptr);
  SetMotherLogical(pMotherLogical);
  if (pMotherLogical) { pMotherLogical->AddDaughter(this); }
  if ((pSurfChk) && (pMotherLogical)) { CheckOverlaps(); }
}

// ----------------------------------------------------------------------
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4PVPlacement::G4PVPlacement( __void__& a )
  : G4VPhysicalVolume(a)
{
}

// ----------------------------------------------------------------------
// Destructor
//
G4PVPlacement::~G4PVPlacement()
{
  if( fallocatedRotM ){ delete this->GetRotation() ; }
}

// ----------------------------------------------------------------------
// IsMany
//
G4bool G4PVPlacement::IsMany() const
{
  return fmany;
}

// ----------------------------------------------------------------------
// SetCopyNo
//
void G4PVPlacement::SetCopyNo(G4int newCopyNo)
{
  fcopyNo = newCopyNo;
}

// ----------------------------------------------------------------------
// IsReplicated
//
G4bool G4PVPlacement::IsReplicated() const
{
  return false;
}

// ----------------------------------------------------------------------
// IsParameterised
//
G4bool G4PVPlacement::IsParameterised() const
{
  return false;
}

// ----------------------------------------------------------------------
// GetParameterisation
//
G4VPVParameterisation* G4PVPlacement::GetParameterisation() const
{
  return nullptr;
}

// ----------------------------------------------------------------------
// GetReplicationData
//
void G4PVPlacement::
GetReplicationData( EAxis&, G4int&, G4double&, G4double&, G4bool& ) const
{
  // No-operations
}

// ----------------------------------------------------------------------
// IsRegularRepeatedStructure
//
// This is for specialised repeated volumes (replicas, parameterised vol.)
//
G4bool G4PVPlacement::IsRegularStructure() const
{
  return false;
}

// ----------------------------------------------------------------------
// IsRegularRepeatedStructure
//
// This is for specialised repeated volumes (replicas, parameterised vol.)
//
G4int G4PVPlacement::GetRegularStructureId() const
{
  return 0;
}

// ----------------------------------------------------------------------
// VolumeType
//
// Information to help identify sub-navigator which will be used
//
EVolume G4PVPlacement::VolumeType() const
{
  return kNormal;
}

// ----------------------------------------------------------------------
// CheckOverlaps
//
G4bool G4PVPlacement::CheckOverlaps(G4int res, G4double tol,
                                    G4bool verbose, G4int maxErr)
{
  if (res <= 0) { return false; }

  G4VSolid* solid = GetLogicalVolume()->GetSolid();
  G4LogicalVolume* motherLog = GetMotherLogical();
  if (motherLog == nullptr) { return false; }

  G4int trials = 0;
  G4bool retval = false;

  if (verbose)
  {
    G4cout << "Checking overlaps for volume "
           << GetName() << ':' << GetCopyNo()
           << " (" << solid->GetEntityType() << ") ... ";
  }

  // Check that random points are gererated correctly
  //
  G4ThreeVector ptmp = solid->GetPointOnSurface();
  if (solid->Inside(ptmp) != kSurface)
  {
    G4String position[3] = { "outside", "surface", "inside" };
    std::ostringstream message;
    message << "Sample point is not on the surface !" << G4endl
            << "          The issue is detected for volume "
            << GetName() << ':' << GetCopyNo()
            << " (" << solid->GetEntityType() << ")" << G4endl
            << "          generated point " << ptmp
            << " is " << position[solid->Inside(ptmp)];
    G4Exception("G4PVPlacement::CheckOverlaps()",
                "GeomVol1002", JustWarning, message);
    return false;
  }

  // Generate random points on the surface of the solid,
  // transform them into the mother volume coordinate system
  // and find the bonding box
  //
  std::vector<G4ThreeVector> points(res);
  G4double xmin =  kInfinity, ymin =  kInfinity, zmin =  kInfinity;
  G4double xmax = -kInfinity, ymax = -kInfinity, zmax = -kInfinity;
  G4AffineTransform Tm(GetRotation(), GetTranslation());
  for (G4int i = 0; i < res; ++i)
  {
    points[i] = Tm.TransformPoint(solid->GetPointOnSurface());
    xmin = std::min(xmin, points[i].x());
    ymin = std::min(ymin, points[i].y());
    zmin = std::min(zmin, points[i].z());
    xmax = std::max(xmax, points[i].x());
    ymax = std::max(ymax, points[i].y());
    zmax = std::max(zmax, points[i].z());
  }
  G4ThreeVector scenter(0.5*(xmax+xmin), 0.5*(ymax+ymin), 0.5*(zmax+zmin));
  G4double sradius = 0.5*G4ThreeVector(xmax-xmin, ymax-ymin, zmax-zmin).mag();

  // Check overlap with the mother volume
  //
  G4int overlapCount = 0;
  G4double overlapSize = -kInfinity;
  G4ThreeVector overlapPoint;
  G4VSolid* motherSolid = motherLog->GetSolid();
  for (G4int i = 0; i < res; ++i)
  {
    G4ThreeVector mp = points[i];
    if (motherSolid->Inside(mp) != kOutside) continue;
    G4double distin = motherSolid->DistanceToIn(mp);
    if (distin < tol) continue; // too small overlap
    ++overlapCount;
    if (distin <= overlapSize) continue;
    overlapSize = distin;
    overlapPoint = mp;
  }

  // Print information on overlap
  //
  if (overlapCount > 0)
  {
    ++trials;
    retval = true;
    std::ostringstream message;
    message << "Overlap with mother volume !" << G4endl
            << "          Overlap is detected for volume "
            << GetName() << ':' << GetCopyNo()
            << " (" << solid->GetEntityType() << ")"
            << " with its mother volume " << motherLog->GetName()
            << " (" << motherSolid->GetEntityType() << ")" << G4endl
            << "          protrusion at mother local point " << overlapPoint
            << " by " << G4BestUnit(overlapSize, "Length")
            << " (max of " << overlapCount << " cases)";
    if (trials >= maxErr)
    {
      message << G4endl
              << "NOTE: Reached maximum fixed number -" << maxErr
              << "- of overlaps reports for this volume !";
    }
    G4Exception("G4PVPlacement::CheckOverlaps()",
                "GeomVol1002", JustWarning, message);
    if (trials >= maxErr)  { return true; }
  }

  // Checking overlaps with each 'sister' volumes
  //
  G4VSolid* previous = nullptr;
  G4ThreeVector pmin_local(0.,0.,0.), pmax_local(0.,0.,0.);

  for (std::size_t k = 0; k < motherLog->GetNoDaughters(); ++k)
  {
    G4VPhysicalVolume* daughter = motherLog->GetDaughter((G4int)k);
    if (daughter == this) continue;
    G4bool check_encapsulation = true;

    G4AffineTransform Td(daughter->GetRotation(), daughter->GetTranslation());
    G4VSolid* daughterSolid = daughter->GetLogicalVolume()->GetSolid();
    if (previous != daughterSolid)
    {
      daughterSolid->BoundingLimits(pmin_local, pmax_local);
      previous = daughterSolid;
    }
    overlapCount = 0;
    overlapSize = -kInfinity;
    if (!Td.IsRotated()) { // no rotation, only translation
      G4ThreeVector offset = Td.NetTranslation();
      G4ThreeVector pmin(pmin_local + offset);
      G4ThreeVector pmax(pmax_local + offset);
      if (pmin.x() >= xmax) continue;
      if (pmin.y() >= ymax) continue;
      if (pmin.z() >= zmax) continue;
      if (pmax.x() <= xmin) continue;
      if (pmax.y() <= ymin) continue;
      if (pmax.z() <= zmin) continue;
      for (G4int i = 0; i < res; ++i)
      {
        G4ThreeVector p = points[i];
        if (p.x() <= pmin.x()) continue;
        if (p.x() >= pmax.x()) continue;
        if (p.y() <= pmin.y()) continue;
        if (p.y() >= pmax.y()) continue;
        if (p.z() <= pmin.z()) continue;
        if (p.z() >= pmax.z()) continue;
        G4ThreeVector md = p - offset;
        if (daughterSolid->Inside(md) == kInside)
        {
          check_encapsulation = false;
          G4double distout = daughterSolid->DistanceToOut(md);
          if (distout < tol) continue; // too small overlap
          ++overlapCount;
          if (distout <= overlapSize) continue;
          overlapSize = distout;
          overlapPoint = md;
        }
      }
    }
    else // transformation with rotation
    {
      G4ThreeVector pmin(pmin_local), pmax(pmax_local);
      G4ThreeVector dcenter = Td.TransformPoint(0.5*(pmin + pmax));
      G4double dradius = 0.5*((pmax - pmin).mag());
      if ((scenter - dcenter).mag2() >= (sradius + dradius)*(sradius + dradius)) continue;
      if (dcenter.x() - dradius >= xmax) continue;
      if (dcenter.y() - dradius >= ymax) continue;
      if (dcenter.z() - dradius >= zmax) continue;
      if (dcenter.x() + dradius <= xmin) continue;
      if (dcenter.y() + dradius <= ymin) continue;
      if (dcenter.z() + dradius <= zmin) continue;

      G4ThreeVector pbox[8] = {
        G4ThreeVector(pmin.x(), pmin.y(), pmin.z()),
        G4ThreeVector(pmax.x(), pmin.y(), pmin.z()),
        G4ThreeVector(pmin.x(), pmax.y(), pmin.z()),
        G4ThreeVector(pmax.x(), pmax.y(), pmin.z()),
        G4ThreeVector(pmin.x(), pmin.y(), pmax.z()),
        G4ThreeVector(pmax.x(), pmin.y(), pmax.z()),
        G4ThreeVector(pmin.x(), pmax.y(), pmax.z()),
        G4ThreeVector(pmax.x(), pmax.y(), pmax.z())
      };
      G4double dxmin =  kInfinity, dymin =  kInfinity, dzmin =  kInfinity;
      G4double dxmax = -kInfinity, dymax = -kInfinity, dzmax = -kInfinity;
      for (G4int i = 0; i < 8; ++i)
      {
        G4ThreeVector p = Td.TransformPoint(pbox[i]);
        dxmin = std::min(dxmin, p.x());
        dymin = std::min(dymin, p.y());
        dzmin = std::min(dzmin, p.z());
        dxmax = std::max(dxmax, p.x());
        dymax = std::max(dymax, p.y());
        dzmax = std::max(dzmax, p.z());
      }
      if (dxmin >= xmax) continue;
      if (dymin >= ymax) continue;
      if (dzmin >= zmax) continue;
      if (dxmax <= xmin) continue;
      if (dymax <= ymin) continue;
      if (dzmax <= zmin) continue;
      for (G4int i = 0; i < res; ++i)
      {
        G4ThreeVector p = points[i];
        if (p.x() >= dxmax) continue;
        if (p.x() <= dxmin) continue;
        if (p.y() >= dymax) continue;
        if (p.y() <= dymin) continue;
        if (p.z() >= dzmax) continue;
        if (p.z() <= dzmin) continue;
        G4ThreeVector md = Td.InverseTransformPoint(p);
        if (daughterSolid->Inside(md) == kInside)
        {
          check_encapsulation = false;
          G4double distout = daughterSolid->DistanceToOut(md);
          if (distout < tol) continue; // too small overlap
          ++overlapCount;
          if (distout <= overlapSize) continue;
          overlapSize = distout;
          overlapPoint = md;
        }
      }
    }

    // Print information on overlap
    //
    if (overlapCount > 0)
    {
      ++trials;
      retval = true;
      std::ostringstream message;
      message << "Overlap with volume already placed !" << G4endl
              << "          Overlap is detected for volume "
              << GetName() << ':' << GetCopyNo()
              << " (" << solid->GetEntityType() << ") with "
              << daughter->GetName() << ':' << daughter->GetCopyNo()
              << " (" << daughterSolid->GetEntityType() << ")" << G4endl
              << "          overlap at local point " << overlapPoint
              << " by " << G4BestUnit(overlapSize, "Length")
              << " (max of " << overlapCount << " cases)";
      if (trials >= maxErr)
      {
        message << G4endl
                << "NOTE: Reached maximum fixed number -" << maxErr
                << "- of overlaps reports for this volume !";
      }
      G4Exception("G4PVPlacement::CheckOverlaps()",
                  "GeomVol1002", JustWarning, message);
      if (trials >= maxErr)  { return true; }
    }
    else if (check_encapsulation)
    {
      // Now checking that 'sister' volume is not totally included
      // and overlapping. Generate a single point inside of
      // the 'sister' volume and verify that the point in NOT inside
      // the current volume
      //
      G4ThreeVector pSurface = daughterSolid->GetPointOnSurface();
      G4ThreeVector normal = daughterSolid->SurfaceNormal(pSurface);
      G4ThreeVector pInside = pSurface - normal*1.e-4; // move point to inside
      G4ThreeVector dPoint = (daughterSolid->Inside(pInside) == kInside) ?
        pInside : pSurface;

      // Transform the generated point to the mother's coordinate system
      // and then to current volume's coordinate system
      //
      G4ThreeVector mp2 = Td.TransformPoint(dPoint);
      G4ThreeVector msi = Tm.InverseTransformPoint(mp2);

      if (solid->Inside(msi) == kInside)
      {
        ++trials;
        retval = true;
        std::ostringstream message;
        message << "Overlap with volume already placed !" << G4endl
                << "          Overlap is detected for volume "
                << GetName() << ':' << GetCopyNo()
                << " (" << solid->GetEntityType() << ")" << G4endl
                << "          apparently fully encapsulating volume "
                << daughter->GetName() << ':' << daughter->GetCopyNo()
                << " (" << daughterSolid->GetEntityType() << ")"
                << " at the same level!";
        if (trials >= maxErr)
        {
          message << G4endl
                  << "NOTE: Reached maximum fixed number -" << maxErr
                  << "- of overlaps reports for this volume !";
        }
        G4Exception("G4PVPlacement::CheckOverlaps()",
                    "GeomVol1002", JustWarning, message);
        if (trials >= maxErr)  { return true; }
      }
    }
  }

  if (verbose && trials == 0) { G4cout << "OK! " << G4endl; }
  return retval;
}

// ----------------------------------------------------------------------
// NewPtrRotMatrix
//
// Auxiliary function for 2nd & 4th constructors (those with G4Transform3D)
// Creates a new rotation matrix on the heap (using "new") and copies its
// argument into it.
//
// NOTE: Ownership of the returned pointer is left to the caller !
//       No entity is currently responsible to delete this memory.
//
G4RotationMatrix*
G4PVPlacement::NewPtrRotMatrix(const G4RotationMatrix &RotMat)
{
  G4RotationMatrix* pRotMatrix;
  if ( RotMat.isIdentity() )
  {
     pRotMatrix = nullptr;
  }
  else
  {
     pRotMatrix = new G4RotationMatrix(RotMat);
  }
  return pRotMatrix;
}
