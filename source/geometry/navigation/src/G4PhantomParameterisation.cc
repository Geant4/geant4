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
// class G4PhantomParameterisation implementation
//
// May 2007 Pedro Arce, first version
//
// --------------------------------------------------------------------

#include "G4PhantomParameterisation.hh"

#include "globals.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VVolumeMaterialScanner.hh"
#include "G4GeometryTolerance.hh"

//------------------------------------------------------------------
G4PhantomParameterisation::G4PhantomParameterisation()
{
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}

 
//------------------------------------------------------------------
G4PhantomParameterisation::~G4PhantomParameterisation()
{
}


//------------------------------------------------------------------
void G4PhantomParameterisation::
BuildContainerSolid( G4VPhysicalVolume* pMotherPhysical )
{
  fContainerSolid = pMotherPhysical->GetLogicalVolume()->GetSolid();
  fContainerWallX = fNoVoxelsX * fVoxelHalfX;
  fContainerWallY = fNoVoxelsY * fVoxelHalfY;
  fContainerWallZ = fNoVoxelsZ * fVoxelHalfZ;

  // CheckVoxelsFillContainer();
}

//------------------------------------------------------------------
void G4PhantomParameterisation::
BuildContainerSolid( G4VSolid* pMotherSolid )
{
  fContainerSolid = pMotherSolid;
  fContainerWallX = fNoVoxelsX * fVoxelHalfX;
  fContainerWallY = fNoVoxelsY * fVoxelHalfY;
  fContainerWallZ = fNoVoxelsZ * fVoxelHalfZ;

  // CheckVoxelsFillContainer();
}


//------------------------------------------------------------------
void G4PhantomParameterisation::
ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  // Voxels cannot be rotated, return translation
  //
  G4ThreeVector trans = GetTranslation( copyNo );

  physVol->SetTranslation( trans );
}


//------------------------------------------------------------------
G4ThreeVector G4PhantomParameterisation::
GetTranslation(const G4int copyNo ) const
{
  CheckCopyNo( copyNo );

  std::size_t nx;
  std::size_t ny;
  std::size_t nz;

  ComputeVoxelIndices( copyNo, nx, ny, nz );

  G4ThreeVector trans( (2*nx+1)*fVoxelHalfX - fContainerWallX,
                       (2*ny+1)*fVoxelHalfY - fContainerWallY,
                       (2*nz+1)*fVoxelHalfZ - fContainerWallZ);
  return trans;
}


//------------------------------------------------------------------
G4VSolid* G4PhantomParameterisation::
ComputeSolid(const G4int, G4VPhysicalVolume* pPhysicalVol) 
{
  return pPhysicalVol->GetLogicalVolume()->GetSolid();
}
       

//------------------------------------------------------------------
G4Material* G4PhantomParameterisation::
ComputeMaterial(const G4int copyNo, G4VPhysicalVolume *, const G4VTouchable *) 
{ 
  CheckCopyNo( copyNo );
  std::size_t matIndex = GetMaterialIndex(copyNo);

  return fMaterials[ matIndex ];
}


//------------------------------------------------------------------
std::size_t G4PhantomParameterisation::
GetMaterialIndex( std::size_t copyNo ) const
{
  CheckCopyNo( copyNo );

  if( fMaterialIndices == nullptr ) { return 0; }
  return *(fMaterialIndices+copyNo);
}


//------------------------------------------------------------------
std::size_t G4PhantomParameterisation::
GetMaterialIndex( std::size_t nx, std::size_t ny, std::size_t nz ) const
{
  std::size_t copyNo = nx + fNoVoxelsX*ny + fNoVoxelsXY*nz;
  return GetMaterialIndex( copyNo );
}


//------------------------------------------------------------------
G4Material*
G4PhantomParameterisation::GetMaterial( std::size_t nx, std::size_t ny, std::size_t nz) const
{
  return fMaterials[GetMaterialIndex(nx,ny,nz)];
}


//------------------------------------------------------------------
G4Material* G4PhantomParameterisation::GetMaterial( std::size_t copyNo ) const
{
  return fMaterials[GetMaterialIndex(copyNo)];
}


//------------------------------------------------------------------
void G4PhantomParameterisation::
ComputeVoxelIndices(const G4int copyNo, std::size_t& nx,
                          std::size_t& ny, std::size_t& nz ) const
{
  CheckCopyNo( copyNo );
  nx = std::size_t(copyNo%fNoVoxelsX);
  ny = std::size_t( (copyNo/fNoVoxelsX)%fNoVoxelsY );
  nz = std::size_t(copyNo/fNoVoxelsXY);
}


//------------------------------------------------------------------
void G4PhantomParameterisation::
CheckVoxelsFillContainer( G4double contX, G4double contY, G4double contZ ) const
{
  G4double toleranceForWarning = 0.25*kCarTolerance;

  // Any bigger value than 0.25*kCarTolerance will give a warning in
  // G4NormalNavigation::ComputeStep(), because the Inverse of a container
  // translation that is Z+epsilon gives -Z+epsilon (and the maximum tolerance
  // in G4Box::Inside is 0.5*kCarTolerance
  //
  G4double toleranceForError = 1.*kCarTolerance;  

  // Any bigger value than kCarTolerance will give an error in GetReplicaNo()
  //
  if( std::fabs(contX-fNoVoxelsX*fVoxelHalfX) >= toleranceForError
   || std::fabs(contY-fNoVoxelsY*fVoxelHalfY) >= toleranceForError  
   || std::fabs(contZ-fNoVoxelsZ*fVoxelHalfZ) >= toleranceForError )
  {
    std::ostringstream message;
    message << "Voxels do not fully fill the container: "
            << fContainerSolid->GetName() << G4endl
            << "        DiffX= " << contX-fNoVoxelsX*fVoxelHalfX << G4endl
            << "        DiffY= " << contY-fNoVoxelsY*fVoxelHalfY << G4endl
            << "        DiffZ= " << contZ-fNoVoxelsZ*fVoxelHalfZ << G4endl
            << "        Maximum difference is: " << toleranceForError;
    G4Exception("G4PhantomParameterisation::CheckVoxelsFillContainer()",
                "GeomNav0002", FatalException, message);

  }
  else if( std::fabs(contX-fNoVoxelsX*fVoxelHalfX) >= toleranceForWarning
        || std::fabs(contY-fNoVoxelsY*fVoxelHalfY) >= toleranceForWarning  
        || std::fabs(contZ-fNoVoxelsZ*fVoxelHalfZ) >= toleranceForWarning )
  {
    std::ostringstream message;
    message << "Voxels do not fully fill the container: "
            << fContainerSolid->GetName() << G4endl
            << "          DiffX= " << contX-fNoVoxelsX*fVoxelHalfX << G4endl
            << "          DiffY= " << contY-fNoVoxelsY*fVoxelHalfY << G4endl
            << "          DiffZ= " << contZ-fNoVoxelsZ*fVoxelHalfZ << G4endl
            << "          Maximum difference is: " << toleranceForWarning;
    G4Exception("G4PhantomParameterisation::CheckVoxelsFillContainer()",
                "GeomNav1002", JustWarning, message);
  }
}
  
 
//------------------------------------------------------------------
G4int G4PhantomParameterisation::
GetReplicaNo( const G4ThreeVector& localPoint, const G4ThreeVector& localDir )
{

  // Check first that point is really inside voxels
  //
  if( fContainerSolid->Inside( localPoint ) == kOutside )
  {
    if( std::fabs(localPoint.x()) - fContainerWallX > kCarTolerance
	&& std::fabs(localPoint.y()) - fContainerWallY > kCarTolerance
	&& std::fabs(localPoint.z()) - fContainerWallZ > kCarTolerance )
    {
      std::ostringstream message;
      message << "Point outside voxels!" << G4endl
	      << "        localPoint - " << localPoint
	      << " - is outside container solid: "
	      << fContainerSolid->GetName() << G4endl
	      << "DIFFERENCE WITH PHANTOM WALLS X: "
	      << std::fabs(localPoint.x()) - fContainerWallX
	      << " Y: " << std::fabs(localPoint.y()) - fContainerWallY
	      << " Z: " << std::fabs(localPoint.z()) - fContainerWallZ;
      G4Exception("G4PhantomParameterisation::GetReplicaNo()", "GeomNav0003",
		  FatalErrorInArgument, message);
    }
  }
  
  // Check the voxel numbers corresponding to localPoint
  // When a particle is on a surface, it may be between -kCarTolerance and
  // +kCartolerance. By a simple distance as:
  //   G4int nx = G4int( (localPoint.x()+)/fVoxelHalfX/2.);
  // those between -kCartolerance and 0 will be placed on voxel N-1 and those
  // between 0 and kCarTolerance on voxel N.
  // To avoid precision problems place the tracks that are on the surface on
  // voxel N-1 if they have negative direction and on voxel N if they have
  // positive direction.
  // Add +kCarTolerance so that they are first placed on voxel N, and then
  // if the direction is negative substract 1

  G4double fx = (localPoint.x()+fContainerWallX+kCarTolerance)/(fVoxelHalfX*2.);
  G4int nx = G4int(fx);

  G4double fy = (localPoint.y()+fContainerWallY+kCarTolerance)/(fVoxelHalfY*2.); 
  G4int ny = G4int(fy);

  G4double fz = (localPoint.z()+fContainerWallZ+kCarTolerance)/(fVoxelHalfZ*2.);
  G4int nz = G4int(fz);

  // If it is on the surface side, check the direction: if direction is
  // negative place it in the previous voxel (if direction is positive it is
  // already in the next voxel).
  // Correct also cases where n = -1 or n = fNoVoxels. It is always traced to be
  // due to multiple scattering: track is entering a voxel but multiple
  // scattering changes the angle towards outside
  //
  if( fx - nx < kCarTolerance*fVoxelHalfX )
  {
    if( localDir.x() < 0 )
    {
      if( nx != 0 )
      {
        nx -= 1;
      }
    }
    else
    {
      if( nx == G4int(fNoVoxelsX) )  
      {
        nx -= 1;       
      }
    }
  }
  if( fy - ny < kCarTolerance*fVoxelHalfY )
  {
    if( localDir.y() < 0 )
    {
      if( ny != 0 )
      {
        ny -= 1;
      }
    }
    else
    {
      if( ny == G4int(fNoVoxelsY) )  
      {
        ny -= 1;       
      }
    }
  }
  if( fz - nz < kCarTolerance*fVoxelHalfZ )
  {
    if( localDir.z() < 0 )
    {
      if( nz != 0 )
      {
        nz -= 1;
      }
    }
    else
    {
      if( nz == G4int(fNoVoxelsZ) )  
      {
        nz -= 1;       
      }
    }
  }
  
  G4int copyNo = G4int(nx + fNoVoxelsX*ny + fNoVoxelsXY*nz);

  // Check if there are still errors 
  //
  G4bool isOK = true;
  if( nx < 0 )
  {
    nx = 0;
    isOK = false;
  }
  else if( nx >= G4int(fNoVoxelsX) )
  {
    nx = G4int(fNoVoxelsX)-1;
    isOK = false;
  }
  if( ny < 0 )
  {
    ny = 0;
    isOK = false;
  }
  else if( ny >= G4int(fNoVoxelsY) )
  {
    ny = G4int(fNoVoxelsY)-1;
    isOK = false;
  }
  if( nz < 0 )
  {
    nz = 0;
    isOK = false;
  }
  else if( nz >= G4int(fNoVoxelsZ) )
  {
    nz = G4int(fNoVoxelsZ)-1;
    isOK = false;
  }
  if( !isOK )
  {
    if( std::fabs(localPoint.x()-fContainerWallX) > kCarTolerance &&
	std::fabs(localPoint.y()-fContainerWallY) > kCarTolerance &&
	std::fabs(localPoint.z()-fContainerWallZ) > kCarTolerance ){    // only if difference is big 
      std::ostringstream message;
      message << "Corrected the copy number! It was negative or too big" << G4endl
	      << "          LocalPoint: " << localPoint << G4endl
	      << "          LocalDir: " << localDir << G4endl
	      << "          Voxel container size: " << fContainerWallX
	      << " " << fContainerWallY << " " << fContainerWallZ << G4endl
	      << "          LocalPoint - wall: "
	      << localPoint.x()-fContainerWallX << " "
	      << localPoint.y()-fContainerWallY << " "
	      << localPoint.z()-fContainerWallZ;
      G4Exception("G4PhantomParameterisation::GetReplicaNo()",
		  "GeomNav1002", JustWarning, message);
    }
    
    copyNo = G4int(nx + fNoVoxelsX*ny + fNoVoxelsXY*nz);
  }

  return copyNo;
}


//------------------------------------------------------------------
void G4PhantomParameterisation::CheckCopyNo( const G4long copyNo ) const
{ 
  if( copyNo < 0 || copyNo >= G4int(fNoVoxels) )
  {
    std::ostringstream message;
    message << "Copy number is negative or too big!" << G4endl
            << "        Copy number: " << copyNo << G4endl
            << "        Total number of voxels: " << fNoVoxels;
    G4Exception("G4PhantomParameterisation::CheckCopyNo()",
                "GeomNav0002", FatalErrorInArgument, message);
  }
}
