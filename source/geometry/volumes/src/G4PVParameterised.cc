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
// $Id: G4PVParameterised.cc 73250 2013-08-22 13:22:23Z gcosmo $
//
// 
// class G4PVParameterised
//
// Implementation
//
// ----------------------------------------------------------------------

#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G4AffineTransform.hh"
#include "G4UnitsTable.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"

// ----------------------------------------------------------------------
// Constructor
//
G4PVParameterised::G4PVParameterised( const G4String& pName,
                                            G4LogicalVolume* pLogical,
                                            G4VPhysicalVolume* pMother,
                                      const EAxis pAxis,
                                      const G4int nReplicas,
                                            G4VPVParameterisation *pParam,
                                            G4bool pSurfChk )
  : G4PVReplica(pName, pLogical, pMother, pAxis, nReplicas, 0, 0),
    fparam(pParam)
{
#ifdef G4VERBOSE  
  if ((pMother) && (pMother->IsParameterised()))
  {
    std::ostringstream message, hint;
    message << "A parameterised volume is being placed" << G4endl
            << "inside another parameterised volume !";
    hint << "To make sure that no overlaps are generated," << G4endl
         << "you should verify the mother replicated shapes" << G4endl
         << "are of the same type and dimensions." << G4endl
         << "   Mother physical volume: " << pMother->GetName() << G4endl
         << "   Parameterised volume: " << pName << G4endl
         << "  (To switch this warning off, compile with G4_NO_VERBOSE)";
    G4Exception("G4PVParameterised::G4PVParameterised()", "GeomVol1002",
                JustWarning, message, G4String(hint.str()));
  }
#endif
  if (pSurfChk) { CheckOverlaps(); }
}

// ----------------------------------------------------------------------
// Constructor
//
G4PVParameterised::G4PVParameterised( const G4String& pName,
                                            G4LogicalVolume* pLogical,
                                            G4LogicalVolume* pMotherLogical,
                                      const EAxis pAxis,
                                      const G4int nReplicas,
                                            G4VPVParameterisation *pParam,
                                            G4bool pSurfChk )
  : G4PVReplica(pName, pLogical, pMotherLogical, pAxis, nReplicas, 0, 0),
    fparam(pParam)
{
  if (pSurfChk) { CheckOverlaps(); }
}

// ----------------------------------------------------------------------
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4PVParameterised::G4PVParameterised( __void__& a )
  : G4PVReplica(a), fparam(0)
{
}

// ----------------------------------------------------------------------
// Destructor
//
G4PVParameterised::~G4PVParameterised()
{
}

// ----------------------------------------------------------------------
// GetParameterisation
//
G4VPVParameterisation* G4PVParameterised::GetParameterisation() const
{
  return fparam;
}

// ----------------------------------------------------------------------
// IsParameterised
//
G4bool G4PVParameterised::IsParameterised() const
{
  return true;
}

// ----------------------------------------------------------------------
// GetReplicationData
//
void G4PVParameterised::GetReplicationData( EAxis& axis,
                                            G4int& nReplicas,
                                            G4double& width,
                                            G4double& offset,
                                            G4bool& consuming) const
{
  axis = faxis;
  nReplicas = fnReplicas;
  width = fwidth;
  offset = foffset;
  consuming = false;
}

// ----------------------------------------------------------------------
// SetRegularStructureId
//
void  G4PVParameterised::SetRegularStructureId( G4int Code )
{
  G4PVReplica::SetRegularStructureId( Code );
  // To undertake additional preparation, a derived volume must
  //   redefine this method, while calling also the above method.
}


// ----------------------------------------------------------------------
// CheckOverlaps
//
G4bool
G4PVParameterised::CheckOverlaps(G4int res, G4double tol,
                                 G4bool verbose, G4int maxErr)
{
  if (res<=0) { return false; }

  G4int trials = 0;
  G4bool retval = false;
  G4VSolid *solidA = 0, *solidB = 0;
  G4LogicalVolume *motherLog = GetMotherLogical();
  G4VSolid *motherSolid = motherLog->GetSolid();
  std::vector<G4ThreeVector> points;

  if (verbose)
  {
    G4cout << "Checking overlaps for parameterised volume "
           << GetName() << " ... ";
  }

  for (G4int i=0; i<GetMultiplicity(); i++)
  {
    solidA = fparam->ComputeSolid(i, this);
    solidA->ComputeDimensions(fparam, i, this);
    fparam->ComputeTransformation(i, this);

    // Create the transformation from daughter to mother
    //
    G4AffineTransform Tm( GetRotation(), GetTranslation() );

    // Generate random points on surface according to the given resolution,
    // transform them to the mother's coordinate system and if no overlaps
    // with the mother volume, cache them in a vector for later use with
    // the daughters
    //
    for (G4int n=0; n<res; n++)
    {
      G4ThreeVector mp = Tm.TransformPoint(solidA->GetPointOnSurface());

      // Checking overlaps with the mother volume
      //
      if (motherSolid->Inside(mp)==kOutside)
      {
        G4double distin = motherSolid->DistanceToIn(mp);
        if (distin > tol)
        {
          trials++; retval = true;
          std::ostringstream message;
          message << "Overlap with mother volume !" << G4endl
                  << "         Overlap is detected for volume "
                  << GetName() << ", parameterised instance: " << i << G4endl
                  << "          with its mother volume "
                  << motherLog->GetName() << G4endl
                  << "          at mother local point " << mp << ", "
                  << "overlapping by at least: "
                  << G4BestUnit(distin, "Length");
          if (trials>=maxErr)
          {
            message << G4endl
                    << "NOTE: Reached maximum fixed number -" << maxErr
                    << "- of overlaps reports for this volume !";
          }
          G4Exception("G4PVParameterised::CheckOverlaps()",
                      "GeomVol1002", JustWarning, message);
          if (trials>=maxErr)  { return true; }
        }
      }
      points.push_back(mp);
    }

    // Checking overlaps with each other parameterised instance
    //
    std::vector<G4ThreeVector>::iterator pos;
    for (G4int j=i+1; j<GetMultiplicity(); j++)
    {
      solidB = fparam->ComputeSolid(j,this);
      solidB->ComputeDimensions(fparam, j, this);
      fparam->ComputeTransformation(j, this);

      // Create the transformation for daughter volume
      //
      G4AffineTransform Td( GetRotation(), GetTranslation() );

      for (pos=points.begin(); pos!=points.end(); pos++)
      {
        // Transform each point according to daughter's frame
        //
        G4ThreeVector md = Td.Inverse().TransformPoint(*pos);

        if (solidB->Inside(md)==kInside)
        {
          G4double distout = solidB->DistanceToOut(md);
          if (distout > tol)
          {
            trials++; retval = true;
            std::ostringstream message;
            message << "Overlap within parameterised volumes !" << G4endl
                    << "          Overlap is detected for volume "
                    << GetName() << ", parameterised instance: " << i << G4endl
                    << "          with parameterised volume instance: " << j
                    << G4endl
                    << "          at local point " << md << ", "
                    << "overlapping by at least: "
                    << G4BestUnit(distout, "Length")
                    << ", related to volume instance: " << j << ".";
            if (trials>=maxErr)
            {
              message << G4endl
                      << "NOTE: Reached maximum fixed number -" << maxErr
                      << "- of overlaps reports for this volume !";
            }
            G4Exception("G4PVParameterised::CheckOverlaps()",
                        "GeomVol1002", JustWarning, message);
            if (trials>=maxErr)  { return true; }
          }
        }
      }
    }
  }
  if (verbose)
  {
    G4cout << "OK! " << G4endl;
  }

  return retval;
}
