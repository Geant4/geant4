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
// G4AdjointPosOnPhysVolGenerator class implementation
//
// Author: L. Desorgher, SpaceIT GmbH - 01.06.2006
// Contract: ESA contract 21435/08/NL/AT
// Customer: ESA/ESTEC
// --------------------------------------------------------------------

#include "G4AdjointPosOnPhysVolGenerator.hh"
#include "G4VSolid.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "Randomize.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"

G4ThreadLocal G4AdjointPosOnPhysVolGenerator*
G4AdjointPosOnPhysVolGenerator::theInstance = nullptr;

// --------------------------------------------------------------------
//
G4AdjointPosOnPhysVolGenerator* G4AdjointPosOnPhysVolGenerator::GetInstance()
{
  if(theInstance == nullptr)
  {
    theInstance = new G4AdjointPosOnPhysVolGenerator;
  }
  return theInstance;
}

// --------------------------------------------------------------------
//
G4VPhysicalVolume*
G4AdjointPosOnPhysVolGenerator::DefinePhysicalVolume(const G4String& aName)
{
  thePhysicalVolume = nullptr;
  theSolid = nullptr;
  G4PhysicalVolumeStore* thePhysVolStore = G4PhysicalVolumeStore::GetInstance();
  for ( unsigned int i=0; i< thePhysVolStore->size(); ++i )
  {
    G4String vol_name =(*thePhysVolStore)[i]->GetName();
    if (vol_name.empty())
    {
      vol_name = (*thePhysVolStore)[i]->GetLogicalVolume()->GetName();
    }
    if (vol_name == aName)
    {
      thePhysicalVolume = (*thePhysVolStore)[i];
    }
  }
  if (thePhysicalVolume != nullptr)
  {
    theSolid = thePhysicalVolume->GetLogicalVolume()->GetSolid();
    ComputeTransformationFromPhysVolToWorld();
  }
  else
  {
    G4cout << "The physical volume with name " << aName
           << " does not exist!!" << G4endl;
    G4cout << "Before generating a source on an external surface " << G4endl
           << "of a volume you should select another physical volume."
           << G4endl; 
  }
  return thePhysicalVolume;
}

// --------------------------------------------------------------------
//
void
G4AdjointPosOnPhysVolGenerator::DefinePhysicalVolume1(const G4String& aName)
{
  thePhysicalVolume = DefinePhysicalVolume(aName);
}

// --------------------------------------------------------------------
//
G4double G4AdjointPosOnPhysVolGenerator::ComputeAreaOfExtSurface()
{
  return ComputeAreaOfExtSurface(theSolid); 
}

// --------------------------------------------------------------------
//
G4double G4AdjointPosOnPhysVolGenerator::ComputeAreaOfExtSurface(G4int NStats)
{
  return ComputeAreaOfExtSurface(theSolid,NStats); 
}

// --------------------------------------------------------------------
//
G4double G4AdjointPosOnPhysVolGenerator::ComputeAreaOfExtSurface(G4double eps)
{
  return ComputeAreaOfExtSurface(theSolid,eps); 
}

// --------------------------------------------------------------------
//
G4double
G4AdjointPosOnPhysVolGenerator::ComputeAreaOfExtSurface(G4VSolid* aSolid)
{
  return ComputeAreaOfExtSurface(aSolid,1.e-3); 
}

// --------------------------------------------------------------------
//
G4double
G4AdjointPosOnPhysVolGenerator::ComputeAreaOfExtSurface(G4VSolid* aSolid,
                                                        G4int NStats)
{
  if (ModelOfSurfaceSource == "OnSolid")
  {
    if (UseSphere)
    {
      return ComputeAreaOfExtSurfaceStartingFromSphere(aSolid,NStats);        
    }
    
    return ComputeAreaOfExtSurfaceStartingFromBox(aSolid,NStats);
  }
  
  G4ThreeVector p, dir;
  if (ModelOfSurfaceSource == "ExternalSphere")
  {
    return GenerateAPositionOnASphereBoundary(aSolid, p,dir);
  }
  
  return GenerateAPositionOnABoxBoundary(aSolid, p,dir);
}

// --------------------------------------------------------------------
//
G4double
G4AdjointPosOnPhysVolGenerator::ComputeAreaOfExtSurface(G4VSolid* aSolid,
                                                        G4double eps)
{
  G4int Nstats = G4int(1./(eps*eps));
  return ComputeAreaOfExtSurface(aSolid,Nstats);
}

// --------------------------------------------------------------------
//
void G4AdjointPosOnPhysVolGenerator::
GenerateAPositionOnTheExtSurfaceOfASolid(G4VSolid* aSolid, G4ThreeVector& p,
                                         G4ThreeVector& direction)
{
  if (ModelOfSurfaceSource == "OnSolid")
  {
    GenerateAPositionOnASolidBoundary(aSolid, p,direction);
    return;
  }
  if (ModelOfSurfaceSource == "ExternalSphere")
  {
    GenerateAPositionOnASphereBoundary(aSolid, p, direction);
    return;
  }        
  GenerateAPositionOnABoxBoundary(aSolid, p, direction);
  return;
}

// --------------------------------------------------------------------
//
void G4AdjointPosOnPhysVolGenerator::
GenerateAPositionOnTheExtSurfaceOfTheSolid(G4ThreeVector& p,
                                           G4ThreeVector& direction)
{
  GenerateAPositionOnTheExtSurfaceOfASolid(theSolid,p,direction);
}

// --------------------------------------------------------------------
//
G4double G4AdjointPosOnPhysVolGenerator::
ComputeAreaOfExtSurfaceStartingFromBox(G4VSolid* aSolid, G4int Nstat)
{
  if ( Nstat <= 0 ) { return 0.; }
  G4double area=1.;
  G4int i=0, j=0;
  while (i<Nstat)
  {
    G4ThreeVector p, direction;
    area = GenerateAPositionOnABoxBoundary( aSolid,p, direction);
    G4double dist_to_in = aSolid->DistanceToIn(p,direction);
    if (dist_to_in<kInfinity/2.) { ++i; }
    ++j;
  }
  area=area*G4double(i)/G4double(j);
  return area;
}

// --------------------------------------------------------------------
//
G4double G4AdjointPosOnPhysVolGenerator::
ComputeAreaOfExtSurfaceStartingFromSphere(G4VSolid* aSolid, G4int Nstat)
{
  if ( Nstat <= 0 ) { return 0.; }
  G4double area=1.;
  G4int i=0, j=0;
  while (i<Nstat)
  {
    G4ThreeVector p, direction;
    area = GenerateAPositionOnASphereBoundary( aSolid,p, direction);
    G4double dist_to_in = aSolid->DistanceToIn(p,direction);
    if (dist_to_in<kInfinity/2.)  { ++i; }
    ++j;
  }
  area=area*G4double(i)/G4double(j);
  return area;
}

// --------------------------------------------------------------------
//
void G4AdjointPosOnPhysVolGenerator::
GenerateAPositionOnASolidBoundary(G4VSolid* aSolid, G4ThreeVector& p,
                                  G4ThreeVector& direction)
{ 
  G4bool find_pos = false;
  while (!find_pos)
  {
    if (UseSphere)
    {
      GenerateAPositionOnASphereBoundary( aSolid,p, direction );
    }
    else
    {
      GenerateAPositionOnABoxBoundary( aSolid,p, direction);
    }
    G4double dist_to_in = aSolid->DistanceToIn(p,direction);
    if (dist_to_in<kInfinity/2.)
    {
      find_pos = true;
      p += 0.999999*direction*dist_to_in;
    }
  }
}

// --------------------------------------------------------------------
//
G4double G4AdjointPosOnPhysVolGenerator::
GenerateAPositionOnASphereBoundary(G4VSolid* aSolid, G4ThreeVector& p,
                                   G4ThreeVector& direction)
{
  G4double minX,maxX,minY,maxY,minZ,maxZ;

  // values needed for CalculateExtent signature

  G4VoxelLimits limit;                // Unlimited
  G4AffineTransform origin;

  // min max extents of pSolid along X,Y,Z

  aSolid->CalculateExtent(kXAxis,limit,origin,minX,maxX);
  aSolid->CalculateExtent(kYAxis,limit,origin,minY,maxY);
  aSolid->CalculateExtent(kZAxis,limit,origin,minZ,maxZ);

  G4ThreeVector center = G4ThreeVector((minX+maxX)/2.,
                                       (minY+maxY)/2.,
                                       (minZ+maxZ)/2.);
  G4double dX=(maxX-minX)/2.;
  G4double dY=(maxY-minY)/2.;
  G4double dZ=(maxZ-minZ)/2.;
  G4double scale=1.01;
  G4double r=scale*std::sqrt(dX*dX+dY*dY+dZ*dZ);

  G4double cos_th2 = G4UniformRand();
  G4double theta = std::acos(std::sqrt(cos_th2));
  G4double phi=G4UniformRand()*CLHEP::twopi;
  direction.setRThetaPhi(1.,theta,phi);
  direction=-direction;
  G4double cos_th = (1.-2.*G4UniformRand());
  theta = std::acos(cos_th);
  if (G4UniformRand() < 0.5)  { theta=CLHEP::pi-theta; }
  phi=G4UniformRand()*CLHEP::twopi;
  p.setRThetaPhi(r,theta,phi);
  p+=center;
  direction.rotateY(theta);
  direction.rotateZ(phi);
  return 4.*CLHEP::pi*r*r;;
}

// --------------------------------------------------------------------
//
G4double G4AdjointPosOnPhysVolGenerator::
GenerateAPositionOnABoxBoundary(G4VSolid* aSolid, G4ThreeVector& p,
                                G4ThreeVector& direction)
{

  G4double ran_var,px,py,pz,minX,maxX,minY,maxY,minZ,maxZ;
  
  // values needed for CalculateExtent signature

  G4VoxelLimits limit;                // Unlimited
  G4AffineTransform origin;

  // min max extents of pSolid along X,Y,Z

  aSolid->CalculateExtent(kXAxis,limit,origin,minX,maxX);
  aSolid->CalculateExtent(kYAxis,limit,origin,minY,maxY);
  aSolid->CalculateExtent(kZAxis,limit,origin,minZ,maxZ);
  
  G4double scale=.1;
  minX-=scale*std::abs(minX);
  minY-=scale*std::abs(minY);
  minZ-=scale*std::abs(minZ);
  maxX+=scale*std::abs(maxX);
  maxY+=scale*std::abs(maxY);
  maxZ+=scale*std::abs(maxZ);
  
  G4double dX=(maxX-minX);
  G4double dY=(maxY-minY);
  G4double dZ=(maxZ-minZ);

  G4double XY_prob=2.*dX*dY;
  G4double YZ_prob=2.*dY*dZ;
  G4double ZX_prob=2.*dZ*dX;
  G4double area=XY_prob+YZ_prob+ZX_prob;
  XY_prob/=area;
  YZ_prob/=area;
  ZX_prob/=area;
  
  ran_var=G4UniformRand();
  G4double cos_th2 = G4UniformRand();
  G4double sth = std::sqrt(1.-cos_th2);
  G4double cth = std::sqrt(cos_th2);
  G4double phi = G4UniformRand()*CLHEP::twopi;
  G4double dirX = sth*std::cos(phi);
  G4double dirY = sth*std::sin(phi);
  G4double dirZ = cth;
  if (ran_var <=XY_prob)  // on the XY faces
  {
    G4double ran_var1=ran_var/XY_prob;
    G4double ranX=ran_var1;
    if (ran_var1<=0.5)
    {
      pz=minZ;
      direction=G4ThreeVector(dirX,dirY,dirZ);
      ranX=ran_var1*2.;
    } 
    else
    {
      pz=maxZ;
      direction=-G4ThreeVector(dirX,dirY,dirZ);
      ranX=(ran_var1-0.5)*2.;
    }
    G4double ranY=G4UniformRand();
    px=minX+(maxX-minX)*ranX;
    py=minY+(maxY-minY)*ranY;
  }
  else if (ran_var <=(XY_prob+YZ_prob))  // on the YZ faces
  {
    G4double ran_var1=(ran_var-XY_prob)/YZ_prob;
    G4double ranY=ran_var1;
    if (ran_var1<=0.5)
    {
      px=minX;
      direction=G4ThreeVector(dirZ,dirX,dirY);
      ranY=ran_var1*2.;
    } 
    else
    {
      px=maxX;
      direction=-G4ThreeVector(dirZ,dirX,dirY);
      ranY=(ran_var1-0.5)*2.;
    }
    G4double ranZ=G4UniformRand();
    py=minY+(maxY-minY)*ranY;
    pz=minZ+(maxZ-minZ)*ranZ;
  }
  else  // on the ZX faces
  {
    G4double ran_var1=(ran_var-XY_prob-YZ_prob)/ZX_prob;
    G4double ranZ=ran_var1;
    if (ran_var1<=0.5)
    {
      py=minY;
      direction=G4ThreeVector(dirY,dirZ,dirX);
      ranZ=ran_var1*2.;
    } 
    else
    {
      py=maxY;
      direction=-G4ThreeVector(dirY,dirZ,dirX);
      ranZ=(ran_var1-0.5)*2.;
    }
    G4double ranX=G4UniformRand();
    px=minX+(maxX-minX)*ranX;
    pz=minZ+(maxZ-minZ)*ranZ;
  }
  
  p=G4ThreeVector(px,py,pz);
  return area;
}

// --------------------------------------------------------------------
//
void G4AdjointPosOnPhysVolGenerator::
GenerateAPositionOnTheExtSurfaceOfThePhysicalVolume(G4ThreeVector& p,
                                                    G4ThreeVector& direction)
{
  if (thePhysicalVolume == nullptr)
  {
    G4cout << "Before generating a source on an external surface" << G4endl
           << "of volume you should select a physical volume" << G4endl; 
    return;
  }
  GenerateAPositionOnTheExtSurfaceOfTheSolid(p,direction);
  p = theTransformationFromPhysVolToWorld.TransformPoint(p);
  direction = theTransformationFromPhysVolToWorld.TransformAxis(direction);
}

// --------------------------------------------------------------------
//
void G4AdjointPosOnPhysVolGenerator::
GenerateAPositionOnTheExtSurfaceOfThePhysicalVolume(G4ThreeVector& p,
                                                    G4ThreeVector& direction,
                                                    G4double& costh_to_normal)
{
  GenerateAPositionOnTheExtSurfaceOfThePhysicalVolume(p, direction);
  costh_to_normal = CosThDirComparedToNormal;
}

// --------------------------------------------------------------------
//
void G4AdjointPosOnPhysVolGenerator::ComputeTransformationFromPhysVolToWorld()
{
  G4VPhysicalVolume* daughter = thePhysicalVolume;
  G4LogicalVolume* mother = thePhysicalVolume->GetMotherLogical();
  theTransformationFromPhysVolToWorld = G4AffineTransform();
  G4PhysicalVolumeStore* thePhysVolStore = G4PhysicalVolumeStore::GetInstance();
  while (mother != nullptr)
  {
    theTransformationFromPhysVolToWorld *=
      G4AffineTransform(daughter->GetFrameRotation(),
                        daughter->GetObjectTranslation());
    for ( unsigned int i=0; i<thePhysVolStore->size(); ++i )
    {
      if ((*thePhysVolStore)[i]->GetLogicalVolume() == mother)
      {
        daughter = (*thePhysVolStore)[i];
        mother = daughter->GetMotherLogical();
        break;
      }
    }
  }
}
