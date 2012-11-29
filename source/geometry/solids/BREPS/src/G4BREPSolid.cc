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
// $Id$
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4BREPSolid.cc
//
// ----------------------------------------------------------------------

#include "G4BREPSolid.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBSbox.hh"
#include "G4BoundingBox3D.hh"
#include "G4FPlane.hh"
#include "G4BSplineSurface.hh"
#include "G4ToroidalSurface.hh"
#include "G4SphericalSurface.hh"

G4Ray G4BREPSolid::Track;
G4double G4BREPSolid::ShortestDistance= kInfinity;
G4int G4BREPSolid::NumberOfSolids=0;

G4BREPSolid::G4BREPSolid(const G4String& name)
 : G4VSolid(name),
   Box(0), Convex(0), AxisBox(0), PlaneSolid(0), place(0), bbox(0),
   intersectionDistance(kInfinity), active(1), startInside(0),
   nb_of_surfaces(0), SurfaceVec(0), RealDist(0.), solidname(name), Id(0),
   fStatistics(1000000), fCubVolEpsilon(0.001), fAreaAccuracy(-1.),
   fCubicVolume(0.), fSurfaceArea(0.), fpPolyhedron(0)
{
  static G4bool warn=true;
  if (warn)
  {
    G4cout
         << "--------------------------------------------------------" << G4endl
         << "WARNING: BREPS classes are being dismissed.  They will |" << G4endl
         << "         be removed starting  from next  Geant4  major |" << G4endl
         << "         release.  Please, consider switching to adopt |" << G4endl
         << "         correspondent CSG or specific primitives.     |" << G4endl
         << "--------------------------------------------------------"
         << G4endl;
    warn = false;
  }
}

G4BREPSolid::G4BREPSolid( const G4String&   name        , 
                                G4Surface** srfVec      , 
                                G4int       numberOfSrfs  )
 : G4VSolid(name),
   Box(0), Convex(0), AxisBox(0), PlaneSolid(0), place(0), bbox(0),
   intersectionDistance(kInfinity), active(1), startInside(0),
   nb_of_surfaces(numberOfSrfs), SurfaceVec(srfVec), RealDist(0.),
   solidname(name), Id(0),
   fStatistics(1000000), fCubVolEpsilon(0.001), fAreaAccuracy(-1.),
   fCubicVolume(0.), fSurfaceArea(0.), fpPolyhedron(0)
{
  static G4bool warn=true;
  if (warn)
  {
    G4cout
         << "--------------------------------------------------------" << G4endl
         << "WARNING: BREPS classes are being dismissed.  They will |" << G4endl
         << "         be removed starting  from next  Geant4  major |" << G4endl
         << "         release.  Please, consider switching to adopt |" << G4endl
         << "         correspondent CSG or specific primitives.     |" << G4endl
         << "--------------------------------------------------------"
         << G4endl;
    warn = false;
  }
  Initialize();
}

G4BREPSolid::G4BREPSolid( __void__& a )
  : G4VSolid(a),
   Box(0), Convex(0), AxisBox(0), PlaneSolid(0), place(0), bbox(0),
   intersectionDistance(kInfinity), active(1), startInside(0),
   nb_of_surfaces(0), SurfaceVec(0), RealDist(0.), solidname(""), Id(0),
   fStatistics(1000000), fCubVolEpsilon(0.001), fAreaAccuracy(-1.),
   fCubicVolume(0.), fSurfaceArea(0.), fpPolyhedron(0)
{
}

G4BREPSolid::~G4BREPSolid()
{
  delete place;
  delete bbox;
  
  for(G4int a=0;a<nb_of_surfaces;++a)
    { delete SurfaceVec[a]; }
  
  delete [] SurfaceVec;
  delete fpPolyhedron;
}

G4BREPSolid::G4BREPSolid(const G4BREPSolid& rhs)
  : G4VSolid(rhs), Box(rhs.Box), Convex(rhs.Convex), AxisBox(rhs.AxisBox),
    PlaneSolid(rhs.PlaneSolid), place(0), bbox(0),
    intersectionDistance(rhs.intersectionDistance), active(rhs.active),
    startInside(rhs.startInside), nb_of_surfaces(rhs.nb_of_surfaces),
    intersection_point(rhs.intersection_point), SurfaceVec(0),
    RealDist(rhs.RealDist), solidname(rhs.solidname), Id(rhs.Id),
    fStatistics(rhs.fStatistics), fCubVolEpsilon(rhs.fCubVolEpsilon),
    fAreaAccuracy(rhs.fAreaAccuracy), fCubicVolume(rhs.fCubicVolume),
    fSurfaceArea(rhs.fSurfaceArea), fpPolyhedron(0)
{
  Initialize();
}

G4BREPSolid& G4BREPSolid::operator = (const G4BREPSolid& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VSolid::operator=(rhs);

   // Copy data
   //
   Box= rhs.Box; Convex= rhs.Convex; AxisBox= rhs.AxisBox;
   PlaneSolid= rhs.PlaneSolid;
   intersectionDistance= rhs.intersectionDistance; active= rhs.active;
   startInside= rhs.startInside; nb_of_surfaces= rhs.nb_of_surfaces;
   intersection_point= rhs.intersection_point;
   RealDist= rhs.RealDist; solidname= rhs.solidname; Id= rhs.Id;
   fStatistics= rhs.fStatistics; fCubVolEpsilon= rhs.fCubVolEpsilon;
   fAreaAccuracy= rhs.fAreaAccuracy; fCubicVolume= rhs.fCubicVolume;
   fSurfaceArea= rhs.fSurfaceArea;

   delete place; delete bbox; place= 0; bbox= 0;
   for(G4int a=0;a<nb_of_surfaces;++a)
     { delete SurfaceVec[a]; }
   delete [] SurfaceVec; SurfaceVec= 0;
   delete fpPolyhedron; fpPolyhedron= 0;

   Initialize();

   return *this;
}  

void G4BREPSolid::Initialize()
{
  if(active)
  {
    // Compute bounding box for solids and surfaces
    // Convert concave planes to convex
    //
    ShortestDistance= kInfinity;
    if (!SurfaceVec) { return; }

    IsBox();
    CheckSurfaceNormals();
    
    if(!Box || !AxisBox)
      IsConvex();
    
    CalcBBoxes();
  }
}

G4String G4BREPSolid::GetEntityType() const
{
  return "Closed_Shell";
}

G4VSolid* G4BREPSolid::Clone() const
{
  return new G4BREPSolid(*this);
}

void G4BREPSolid::Reset() const
{
  ((G4BREPSolid*)this)->active=1;
  ((G4BREPSolid*)this)->intersectionDistance=kInfinity;
  ((G4BREPSolid*)this)->startInside=0; 
     
  for(register G4int a=0;a<nb_of_surfaces;a++)
    SurfaceVec[a]->Reset();

  ShortestDistance = kInfinity;
}

void G4BREPSolid::CheckSurfaceNormals()
{
  if(!PlaneSolid)
    return; // All faces must be planar
 
  Convex=1;
 
  // Checks that the normals of the surfaces point outwards.
  // If not, turns the Normal to point out.
  
  // Loop through each face and check the G4Vector3D of the Normal
  //
  G4Surface* srf;
  G4Point3D V;
  
  G4int SrfNum = 0;
  G4double YValue=0;
  G4Point3D Pt;
  
  G4int a, b;
  for(a=0; a<nb_of_surfaces; a++)
  {
    // Find vertex point containing extreme y value
    //
    srf = SurfaceVec[a];
    G4int Points = srf->GetNumberOfPoints();
    
    for(b =0; b<Points; b++)
    {
      Pt = (G4Point3D)srf->GetPoint(b);    
      if(YValue < Pt.y())
      {
        YValue = Pt.y();
        SrfNum = a;   // Save srf number
      }
    }
  }

  // Move the selected face to the first in the List
  //
  srf = SurfaceVec[SrfNum];            
  
  // Start handling the surfaces in order and compare
  // the neighbouring ones and turn their normals if they
  // point inwards
  //
  G4Point3D Pt1;
  G4Point3D Pt2;
  G4Point3D Pt3;
  G4Point3D Pt4;
  
  G4Vector3D N1;
  G4Vector3D N2;    
  G4Vector3D N3;    
  G4Vector3D N4;    

  G4int* ConnectedList = new G4int[nb_of_surfaces];    

  for(a=0; a<nb_of_surfaces; a++)
    ConnectedList[a]=0;
  
  G4Surface* ConnectedSrf;

  for(a=0; a<nb_of_surfaces-1; a++)
  {
    if(ConnectedList[a] == 0) 
      break;
    else
      ConnectedList[a]=1;
    
    srf = SurfaceVec[a];
    G4int SrfPoints = srf->GetNumberOfPoints();
    N1 = (srf->Norm())->GetDir();

    for(b=a+1; b<nb_of_surfaces; b++)
    {
      if(ConnectedList[b] == 1) 
        break;
      else
        ConnectedList[b]=1;
     
      // Get next in List
      // 
      ConnectedSrf = SurfaceVec[b];    

      // Check if it is connected to srf by looping through the points.
      //
      G4int ConnSrfPoints = ConnectedSrf->GetNumberOfPoints();

      for(G4int c=0;c<SrfPoints;c++)
      {
        Pt1 = srf->GetPoint(c);

        for(G4int d=0;d<ConnSrfPoints;d++)
        {
          // Find common points
          //
          Pt2 = (ConnectedSrf)->GetPoint(d);

          if( Pt1 == Pt2 ) 
          {
            // Common point found. Compare normals.
            //
            N2 = ((ConnectedSrf)->Norm())->GetDir();

            // Check cross product.
            //
            G4Vector3D CP1 = G4Vector3D( N1.cross(N2) );
            G4double CrossProd1 = CP1.x()+CP1.y()+CP1.z();

            // Create the other normals
            //
            if(c==0) 
              Pt3 = srf->GetPoint(c+1);
            else
              Pt3 = srf->GetPoint(0);

            N3 = (Pt1-Pt3);

            if(d==0) 
              Pt4 = (ConnectedSrf)->GetPoint(d+1);
            else
              Pt4 = (ConnectedSrf)->GetPoint(0);

            N4 = (Pt1-Pt4);

            G4Vector3D CP2 = G4Vector3D( N3.cross(N4) );
            G4double CrossProd2 = CP2.x()+CP2.y()+CP2.z();

            G4cout << "\nCroosProd2: " << CrossProd2;

            if( (CrossProd1 < 0 && CrossProd2 < 0) ||
                (CrossProd1 > 0 && CrossProd2 > 0)    )
            {
              // Turn Normal
              //
              (ConnectedSrf)->Norm()
                ->SetDir(-1 * (ConnectedSrf)->Norm()->GetDir());

              // Take the CrossProd1 again as the other Normal was turned.
              //
              CP1 = N1.cross(N2);
              CrossProd1 = CP1.x()+CP1.y()+CP1.z();
            }
            if(CrossProd1 > 0) 
              Convex=0;
          }
        }
      }
    }
  }
  
  delete []ConnectedList;
}

G4int G4BREPSolid::IsBox()
{
  // This is done by checking that the solid consists of 6 planes.
  // Then the type is checked to be planar face by face.
  // For each G4Plane the Normal is computed. The dot product
  // of one face Normal and each other face Normal is computed.
  // One result should be 1 and the rest 0 in order to the solid
  // to be a box.

  Box=0;
  G4Surface* srf1, *srf2;
  register G4int a;
  
  // Compute the Normal for the planes
  //
  for(a=0; a < nb_of_surfaces;a++)    
  {
    srf1 = SurfaceVec[a];
    
    if(srf1->MyType()==1)
      (srf1)->Project();     // Compute the projection
    else
    {
      PlaneSolid=0;
      return 0;
    }
  }

  // Check that all faces are planar
  //
  for(a=0; a < nb_of_surfaces;a++)    
  {
    srf1 = SurfaceVec[a];
    
    if (srf1->MyType()!=1)
      return 0;
  }

  PlaneSolid = 1;
  
  // Check that the amount of faces is correct
  //
  if(nb_of_surfaces!=6) return 0;
  
  G4Point3D Pt;
  G4int Points;
  G4int Sides=0;
  G4int Opposite=0;

  srf1 = SurfaceVec[0];
  Points = (srf1)->GetNumberOfPoints();
  
  if(Points!=4)
    return 0;
    
  G4Vector3D Normal1 = (srf1->Norm())->GetDir();
  G4double Result;
  
  for(G4int b=1; b < nb_of_surfaces;b++)
  {
    srf2 = SurfaceVec[b];
    G4Vector3D Normal2 = ((srf2)->Norm())->GetDir();
    Result = std::fabs(Normal1 * Normal2);
    
    if((Result != 0) && (Result != 1))
      return 0;
    else
    {
      if(!(G4int)Result) 
        Sides++;
      else
        if(((G4int)Result) == 1)
          Opposite++;
    }
  }

  if((Opposite != 1) && (Sides != nb_of_surfaces-2))
    return 0;
  
  G4Vector3D x_axis(1,0,0);
  G4Vector3D y_axis(0,1,0);
  
  if(((std::fabs(x_axis * Normal1) == 1) && (std::fabs(y_axis * Normal1) == 0)) ||
     ((std::fabs(x_axis * Normal1) == 0) && (std::fabs(y_axis * Normal1) == 1)) ||
     ((std::fabs(x_axis * Normal1) == 0) && (std::fabs(y_axis * Normal1) == 0)))
    AxisBox=1;
  else 
    Box=1;
  
  return 1;
}

G4bool G4BREPSolid::IsConvex()
{
  if (!PlaneSolid) { return 0; }    // All faces must be planar

  // This is not robust. There can be concave solids
  // where the concavity comes for example from three triangles.

  // Additional checking 20.8. For each face the connecting faces are
  // found and the cross product computed between the face and each
  // connecting face. If the result changes value at any point the
  // solid is concave.
  
  G4Surface* Srf=0;
  G4Surface* ConnectedSrf=0;
  G4int Result;
  Convex = 1;
  
  G4int a, b, c, d;
  for(a=0;a<nb_of_surfaces;a++)
  {
    Srf = SurfaceVec[a];

    // Primary test. Test wether any one of the faces
    // is concave -> solid is concave. This is not enough to
    // distinguish all the cases of concavity.
    //
    Result = Srf->IsConvex();
    
    if (Result != -1)
    {
      Convex = 0;
      return 0;
    }
  }

  Srf = SurfaceVec[0];        

  G4int ConnectingPoints=0;  
  G4Point3D  Pt1, Pt2;
  G4Vector3D N1, N2;    

  // L. Broglia
  // The number of connecting points can be 
  // (nb_of_surfaces-1) * nb_of_surfaces  (loop a & loop b)
  
  // G4int* ConnectedList = new G4int[nb_of_surfaces];
  const G4int maxCNum =  (nb_of_surfaces-1)*nb_of_surfaces; 
  G4int* ConnectedList = new G4int[maxCNum];  

  for(a=0; a<maxCNum; a++)
  {
    ConnectedList[a]=0;
  }
  
  G4int Connections=0;

  for(a=0; a<nb_of_surfaces-1; a++)
  {
    Srf = SurfaceVec[a];
    G4int SrfPoints = Srf->GetNumberOfPoints();
    Result=0;
    
    for(b=0; b<nb_of_surfaces; b++)
    {
      if (b==a) { continue; }

      // Get next in List
      //
      ConnectedSrf = SurfaceVec[b];
      
      // Check if it is connected to Srf by looping through the points.
      //
      G4int ConnSrfPoints = ConnectedSrf->GetNumberOfPoints();
      
      for(c=0; c<SrfPoints; c++)
      {
        const G4Point3D& Pts1 =Srf->GetPoint(c);

        for(d=0; d<ConnSrfPoints; d++)
        {
          // Find common points
          //
          const G4Point3D& Pts2 = ConnectedSrf->GetPoint(d);
          if (Pts1 == Pts2)  { ConnectingPoints++; }
        }
        if (ConnectingPoints > 0)  { break; }
      }
      
      if (ConnectingPoints > 0)
      {
        Connections++;
        ConnectedList[Connections]=b;
      }
      ConnectingPoints=0;
    }
  }
  
  // If connected, check for concavity.
  // Get surfaces from ConnectedList and compare their normals
  //
  for(c=0; c<Connections; c++)
  {
    G4int Left=0; 
    G4int Right =0;
    G4int tmp = ConnectedList[c];
    
    Srf = SurfaceVec[tmp];
    ConnectedSrf = SurfaceVec[tmp+1];
    
    // Get normals
    //
    N1 = Srf->Norm()->GetDir();
    N2 = ConnectedSrf->Norm()->GetDir();

    // Check cross product
    //
    G4Vector3D CP = G4Vector3D( N1.cross(N2) ); 
    G4double CrossProd = CP.x()+CP.y()+CP.z();
    if (CrossProd > 0)  { Left++;  }
    if (CrossProd < 0)  { Right++; }
    if (Left&&Right)
    {
      Convex = 0;
      delete [] ConnectedList;
      return 0;
    }
    Connections=0;
  } 
  
  Convex=1;

  delete [] ConnectedList;

  return 1;
}

G4bool G4BREPSolid::CalculateExtent(const EAxis pAxis,
                                    const G4VoxelLimits& pVoxelLimit,
                                    const G4AffineTransform& pTransform,
                                          G4double& pMin, G4double& pMax) const
{
  G4Point3D Min = bbox->GetBoxMin();
  G4Point3D Max = bbox->GetBoxMax();

  if (!pTransform.IsRotated())
    {
      // Special case handling for unrotated boxes
      // Compute x/y/z mins and maxs respecting limits, with early returns
      // if outside limits. Then switch() on pAxis
      //
      G4double xoffset,xMin,xMax;
      G4double yoffset,yMin,yMax;
      G4double zoffset,zMin,zMax;

      xoffset=pTransform.NetTranslation().x();
      xMin=xoffset+Min.x();
      xMax=xoffset+Max.x();
      if (pVoxelLimit.IsXLimited())
      {
        if (xMin>pVoxelLimit.GetMaxXExtent()
            ||xMax<pVoxelLimit.GetMinXExtent())
        {
          return false;
        }
        else
        {
          if (xMin<pVoxelLimit.GetMinXExtent())
          {
            xMin=pVoxelLimit.GetMinXExtent();
          }
          if (xMax>pVoxelLimit.GetMaxXExtent())
          {
            xMax=pVoxelLimit.GetMaxXExtent();
          }
        }
      }

      yoffset=pTransform.NetTranslation().y();
      yMin=yoffset+Min.y();
      yMax=yoffset+Max.y();
      if (pVoxelLimit.IsYLimited())
      {
        if (yMin>pVoxelLimit.GetMaxYExtent()
            ||yMax<pVoxelLimit.GetMinYExtent())
        {
          return false;
        }
        else
        {
          if (yMin<pVoxelLimit.GetMinYExtent())
          {
            yMin=pVoxelLimit.GetMinYExtent();
          }
          if (yMax>pVoxelLimit.GetMaxYExtent())
          {
            yMax=pVoxelLimit.GetMaxYExtent();
          }
        }
      }

      zoffset=pTransform.NetTranslation().z();
      zMin=zoffset+Min.z();
      zMax=zoffset+Max.z();
      if (pVoxelLimit.IsZLimited())
      {
        if (zMin>pVoxelLimit.GetMaxZExtent()
            ||zMax<pVoxelLimit.GetMinZExtent())
        {
          return false;
        }
        else
        {
          if (zMin<pVoxelLimit.GetMinZExtent())
          {
            zMin=pVoxelLimit.GetMinZExtent();
          }
          if (zMax>pVoxelLimit.GetMaxZExtent())
          {
            zMax=pVoxelLimit.GetMaxZExtent();
          }
        }
      }

      switch (pAxis)
      {
        case kXAxis:
          pMin=xMin;
          pMax=xMax;
          break;
        case kYAxis:
          pMin=yMin;
          pMax=yMax;
          break;
        case kZAxis:
          pMin=zMin;
          pMax=zMax;
          break;
        default:
          break;
      }
      pMin-=kCarTolerance;
      pMax+=kCarTolerance;

      return true;
    }
  else
    {
      // General rotated case - create and clip mesh to boundaries

      G4bool existsAfterClip=false;
      G4ThreeVectorList *vertices;

      pMin=+kInfinity;
      pMax=-kInfinity;

      // Calculate rotated vertex coordinates
      //
      vertices=CreateRotatedVertices(pTransform);
      ClipCrossSection(vertices,0,pVoxelLimit,pAxis,pMin,pMax);
      ClipCrossSection(vertices,4,pVoxelLimit,pAxis,pMin,pMax);
      ClipBetweenSections(vertices,0,pVoxelLimit,pAxis,pMin,pMax);

      if ( (pMin!=kInfinity) || (pMax!=-kInfinity) )
      {
        existsAfterClip=true;
    
        // Add 2*tolerance to avoid precision troubles
        //
        pMin-=kCarTolerance;
        pMax+=kCarTolerance;
      }
      else
      {
        // Check for case where completely enveloping clipping volume.
        // If point inside then we are confident that the solid completely
        // envelopes the clipping volume. Hence set min/max extents according
        // to clipping volume extents along the specified axis.
        //
        G4ThreeVector clipCentre(
                (pVoxelLimit.GetMinXExtent()+pVoxelLimit.GetMaxXExtent())*0.5,
                (pVoxelLimit.GetMinYExtent()+pVoxelLimit.GetMaxYExtent())*0.5,
                (pVoxelLimit.GetMinZExtent()+pVoxelLimit.GetMaxZExtent())*0.5);

        if (Inside(pTransform.Inverse().TransformPoint(clipCentre))!=kOutside)
        {
          existsAfterClip=true;
          pMin=pVoxelLimit.GetMinExtent(pAxis);
          pMax=pVoxelLimit.GetMaxExtent(pAxis);
        }
      }
      delete vertices;
      return existsAfterClip;
    }
}

G4ThreeVectorList*
G4BREPSolid::CreateRotatedVertices(const G4AffineTransform& pTransform) const
{
  G4Point3D Min = bbox->GetBoxMin();
  G4Point3D Max = bbox->GetBoxMax();

  G4ThreeVectorList *vertices;
  vertices=new G4ThreeVectorList();
    
  if (vertices)
  {
    vertices->reserve(8);
    G4ThreeVector vertex0(Min.x(),Min.y(),Min.z());
    G4ThreeVector vertex1(Max.x(),Min.y(),Min.z());
    G4ThreeVector vertex2(Max.x(),Max.y(),Min.z());
    G4ThreeVector vertex3(Min.x(),Max.y(),Min.z());
    G4ThreeVector vertex4(Min.x(),Min.y(),Max.z());
    G4ThreeVector vertex5(Max.x(),Min.y(),Max.z());
    G4ThreeVector vertex6(Max.x(),Max.y(),Max.z());
    G4ThreeVector vertex7(Min.x(),Max.y(),Max.z());

    vertices->push_back(pTransform.TransformPoint(vertex0));
    vertices->push_back(pTransform.TransformPoint(vertex1));
    vertices->push_back(pTransform.TransformPoint(vertex2));
    vertices->push_back(pTransform.TransformPoint(vertex3));
    vertices->push_back(pTransform.TransformPoint(vertex4));
    vertices->push_back(pTransform.TransformPoint(vertex5));
    vertices->push_back(pTransform.TransformPoint(vertex6));
    vertices->push_back(pTransform.TransformPoint(vertex7));
  }
  else
  {
    G4Exception("G4BREPSolid::CreateRotatedVertices()", "GeomSolids0003",
                FatalException, "Out of memory - Cannot allocate vertices!");
  }
  return vertices;
}

EInside G4BREPSolid::Inside(register const G4ThreeVector& Pt) const
{
  // This function finds if the point Pt is inside, 
  // outside or on the surface of the solid

  const G4double sqrHalfTolerance = kCarTolerance*kCarTolerance*0.25;

  G4Vector3D v(1, 0, 0.01);
  G4Vector3D Pttmp(Pt);
  G4Vector3D Vtmp(v);
  G4Ray r(Pttmp, Vtmp);
  
  // Check if point is inside the PCone bounding box
  //
  if( !GetBBox()->Inside(Pttmp) )
    return kOutside;

  // Set the surfaces to active again
  //
  Reset();
  
  // Test if the bounding box of each surface is intersected
  // by the ray. If not, the surface become deactive.
  //
  TestSurfaceBBoxes(r);
  
  G4int hits=0, samehit=0;

  for(G4int a=0; a < nb_of_surfaces; a++)
  {
    if(SurfaceVec[a]->IsActive())
    {
      // Count the number of intersections. If this number is odd,
      // the start of the ray is inside the volume bounded by the surfaces,
      // so increment the number of intersection by 1 if the point is not
      // on the surface and if this intersection was not found before.
      //
      if( (SurfaceVec[a]->Intersect(r)) & 1 )
      {
        // Test if the point is on the surface
        //
        if(SurfaceVec[a]->GetDistance() < sqrHalfTolerance)
          return kSurface;

        // Test if this intersection was found before
        //
        for(G4int i=0; i<a; i++)
          if(SurfaceVec[a]->GetDistance() == SurfaceVec[i]->GetDistance())
          {
            samehit++;
            break;
          }

        // Count the number of surfaces intersected by the ray
        //
        if(!samehit)
          hits++;
      }
    }
  }
   
  // If the number of surfaces intersected is odd,
  // the point is inside the solid
  //
  if(hits&1)
    return kInside;
  else
    return kOutside;
}

G4ThreeVector G4BREPSolid::SurfaceNormal(const G4ThreeVector& Pt) const
{  
  // This function calculates the normal of the surface at a point on the
  // surface. If the point is not on the surface the result is undefined.
  // Note : the sense of the normal depends on the sense of the surface.

  const G4double sqrHalfTolerance = kCarTolerance*kCarTolerance*0.25;
  G4int          iplane;
    
  // Find on which surface the point is
  //
  for(iplane = 0; iplane < nb_of_surfaces; iplane++)
  {
    if(SurfaceVec[iplane]->HowNear(Pt) < sqrHalfTolerance)
      // the point is on this surface
      break;
  }
  
  // Calculate the normal at this point
  //
  G4ThreeVector norm = SurfaceVec[iplane]->SurfaceNormal(Pt);

  return norm.unit();
}

G4double G4BREPSolid::DistanceToIn(const G4ThreeVector& Pt) const
{
  // Calculates the shortest distance ("safety") from a point
  // outside the solid to any boundary of this solid.
  // Return 0 if the point is already inside.

  G4double *dists = new G4double[nb_of_surfaces];
  G4int a;

  // Set the surfaces to active again
  //
  Reset();
  
  // Compute the shortest distance of the point to each surface.
  // Be careful : it's a signed value
  //
  for(a=0; a< nb_of_surfaces; a++)  
    dists[a] = SurfaceVec[a]->HowNear(Pt);
     
  G4double Dist = kInfinity;
  
  // If dists[] is positive, the point is outside, so take the shortest of
  // the shortest positive distances dists[] can be equal to 0 : point on
  // a surface.
  // ( Problem with the G4FPlane : there is no inside and no outside...
  //   So, to test if the point is inside to return 0, utilize the Inside()
  //   function. But I don't know if it is really needed because dToIn is 
  //   called only if the point is outside )
  //
  for(a = 0; a < nb_of_surfaces; a++)
    if( std::fabs(Dist) > std::fabs(dists[a]) ) 
      //if( dists[a] >= 0)
        Dist = dists[a];
  
  delete[] dists;

  if(Dist == kInfinity)
    return 0;  // the point is inside the solid or on a surface
  else
    return std::fabs(Dist);
}

G4double G4BREPSolid::DistanceToIn(register const G4ThreeVector& Pt, 
                                   register const G4ThreeVector& V  ) const
{
  // Calculates the distance from a point outside the solid
  // to the solid's boundary along a specified direction vector.
  //
  // Note : Intersections with boundaries less than the tolerance must be
  //        ignored if the direction is away from the boundary.
  
  G4int a;
  
  // Set the surfaces to active again
  //
  Reset();
  
  const G4double sqrHalfTolerance = kCarTolerance*kCarTolerance*0.25;
  G4Vector3D Pttmp(Pt);
  G4Vector3D Vtmp(V);   
  G4Ray r(Pttmp, Vtmp);

  // Test if the bounding box of each surface is intersected
  // by the ray. If not, the surface become deactive.
  //
  TestSurfaceBBoxes(r);
  
  ShortestDistance = kInfinity;
  
  for(a=0; a< nb_of_surfaces; a++)
  {
    if( SurfaceVec[a]->IsActive() )
    {
      // Test if the ray intersects the surface
      //
      if( SurfaceVec[a]->Intersect(r) )
      {
        G4double surfDistance = SurfaceVec[a]->GetDistance();

        // If more than 1 surface is intersected, take the nearest one
        //
        if( surfDistance < ShortestDistance )
        {
          if( surfDistance > sqrHalfTolerance )
          {
            ShortestDistance = surfDistance;
          }
          else
          {
            // The point is within the boundary. It is ignored it if
            // the direction is away from the boundary
            //
            G4Vector3D Norm = SurfaceVec[a]->SurfaceNormal(Pttmp);

            if( (Norm * Vtmp) < 0 )
            {
              ShortestDistance = surfDistance;
            }
          }
        }
      }
    }
  }
  
  // Be careful !
  // SurfaceVec->Distance is in fact the squared distance
  //
  if(ShortestDistance != kInfinity)
    return std::sqrt(ShortestDistance);
  else
    return kInfinity;  // No intersection
}

G4double G4BREPSolid::DistanceToOut(register const G4ThreeVector& P, 
                                    register const G4ThreeVector& D, 
                                             const G4bool, 
                                                   G4bool *validNorm, 
                                                   G4ThreeVector*    ) const
{
  // Calculates the distance from a point inside the solid to the solid's
  // boundary along a specified direction vector.
  // Returns 0 if the point is already outside.
  //
  // Note : If the shortest distance to a boundary is less than the tolerance,
  //        it is ignored. This allows for a point within a tolerant boundary
  //        to leave immediately.

  // Set the surfaces to active again
  //
  Reset();

  const G4double sqrHalfTolerance = kCarTolerance*kCarTolerance*0.25;
  G4Vector3D Ptv = P;
  G4int a;

  if(validNorm)
    *validNorm=false;

  G4Vector3D Pttmp(Ptv);
  G4Vector3D Vtmp(D);   
  
  G4Ray r(Pttmp, Vtmp);

  // Test if the bounding box of each surface is intersected
  // by the ray. If not, the surface become deactive.
  //
  TestSurfaceBBoxes(r);
  
  ShortestDistance = kInfinity;
 
  for(a=0; a< nb_of_surfaces; a++)
  {
    if(SurfaceVec[a]->IsActive())
    {
      // Test if the ray intersect the surface
      //
      if( (SurfaceVec[a]->Intersect(r)) )
      {
        // If more than 1 surface is intersected, take the nearest one
        //
        G4double surfDistance = SurfaceVec[a]->GetDistance();
        if( surfDistance < ShortestDistance )
        {
          if( surfDistance > sqrHalfTolerance )
          {
            ShortestDistance = surfDistance;
          }
          else
          {
            // The point is within the boundary: ignore it
          }
        }
      }
    }
  }

  // Be careful !
  // SurfaceVec->Distance is in fact the squared distance
  //
  if(ShortestDistance != kInfinity)
    return std::sqrt(ShortestDistance);
  else
    return 0.0;  // No intersection is found, the point is outside
}

G4double G4BREPSolid::DistanceToOut(const G4ThreeVector& Pt)const
{
  // Calculates the shortest distance ("safety") from a point
  // inside the solid to any boundary of this solid.
  // Returns 0 if the point is already outside.

  G4double *dists = new G4double[nb_of_surfaces];
  G4int a;

  // Set the surfaces to active again
  //
  Reset();
  
  // Compute the shortest distance of the point to each surfaces
  // Be careful : it's a signed value
  //
  for(a=0; a< nb_of_surfaces; a++)
    dists[a] = SurfaceVec[a]->HowNear(Pt);  

  G4double Dist = kInfinity;
  
  // If dists[] is negative, the point is inside so take the shortest of the
  // shortest negative distances dists[] can be equal to 0 : point on a
  // surface
  // ( Problem with the G4FPlane : there is no inside and no outside...
  //   So, to test if the point is outside to return 0, utilize the Inside()
  //   function. But I don`t know if it is really needed because dToOut is 
  //   called only if the point is inside )
  //
  for(a = 0; a < nb_of_surfaces; a++)
    if( std::fabs(Dist) > std::fabs(dists[a]) ) 
      //if( dists[a] <= 0)
        Dist = dists[a];
  
  delete[] dists;

  if(Dist == kInfinity)
    return 0;  // The point is ouside the solid or on a surface
  else
    return std::fabs(Dist);
}

void G4BREPSolid::DescribeYourselfTo (G4VGraphicsScene& scene) const 
{
  scene.AddSolid (*this);
}

G4Polyhedron* G4BREPSolid::CreatePolyhedron () const
{
  // Approximate implementation, just a box ...

  G4Point3D Min = bbox->GetBoxMin();
  G4Point3D Max = bbox->GetBoxMax();  

  return new G4PolyhedronBox (Max.x(), Max.y(), Max.z());
}

G4NURBS* G4BREPSolid::CreateNURBS () const 
{
  // Approximate implementation, just a box ...

  G4Point3D Min = bbox->GetBoxMin();
  G4Point3D Max = bbox->GetBoxMax();  

  return new G4NURBSbox (Max.x(), Max.y(), Max.z());
}

void G4BREPSolid::CalcBBoxes()
{
  // First initialization. Calculates the bounding boxes
  // for the surfaces and for the solid.

  G4Surface* srf;
  G4Point3D min, max;
  
  if(active)
  {
    min =  PINFINITY;
    max = -PINFINITY;
    
    for(G4int a = 0;a < nb_of_surfaces;a++)
    {
      // Get first in List
      //
      srf = SurfaceVec[a];
/*
      G4int convex=1; 
      G4int concavepoint=-1;

      if (srf->MyType() == 1) 
      {
        concavepoint = srf->IsConvex();
        convex = srf->GetConvex();
      }

      // Make bbox for face
      //
      if(convex && Concavepoint==-1)
*/
      {
        srf->CalcBBox();
        G4Point3D box_min = srf->GetBBox()->GetBoxMin();
        G4Point3D box_max = srf->GetBBox()->GetBoxMax();

        // Find max and min of face bboxes to make solids bbox.

        // replace by Extend
        // max < box_max
        //
        if(max.x() < box_max.x()) max.setX(box_max.x()); 
        if(max.y() < box_max.y()) max.setY(box_max.y()); 
        if(max.z() < box_max.z()) max.setZ(box_max.z()); 

        // min > box_min
        //
        if(min.x() > box_min.x()) min.setX(box_min.x()); 
        if(min.y() > box_min.y()) min.setY(box_min.y()); 
        if(min.z() > box_min.z()) min.setZ(box_min.z());
      }
    }
    bbox =  new G4BoundingBox3D(min, max);
    return;
  }
  G4Exception("G4BREPSolid::CalcBBoxes()", "GeomSolids1002",
              JustWarning, "No bbox calculated for solid.");
}

void G4BREPSolid::RemoveHiddenFaces(register const G4Ray& rayref,
                                                   G4int In      ) const
{
  // Deactivates the planar faces that are on the "back" side of a solid.
  // B-splines are not handled by this function. Also cases where the ray
  // starting point is Inside the bbox of the solid are ignored as we don't
  // know if the starting point is Inside the actual solid except for
  // axis-oriented box-like solids.
  
  register G4Surface* srf;
  register const G4Vector3D& RayDir = rayref.GetDir();
  register G4double Result;
  G4int a;

  // In all other cases the ray starting point is outside the solid
  //
  if(!In)
    for(a=0; a<nb_of_surfaces; a++)
    {
      // Deactivates the solids faces that are hidden
      //
      srf = SurfaceVec[a];
      if(srf->MyType()==1)
      {
        const G4Vector3D& Normal = (srf->Norm())->GetDir();
        Result = (RayDir * Normal);

        if( Result >= 0 )
          srf->Deactivate();
      }
    }
  else
    for(a=0; a<nb_of_surfaces; a++)
    {
      // Deactivates the AxisBox type solids faces whose normals 
      // point in the G4Vector3D opposite to the rays G4Vector3D
      // i.e. are behind the ray starting point as in this case the
      // ray starts from Inside the solid.
      //
      srf = SurfaceVec[a];
      if(srf->MyType()==1)
      {
        const G4Vector3D& Normal = (srf->Norm())->GetDir();
        Result = (RayDir * Normal);

        if( Result < 0 )
          srf->Deactivate();
      }
    }
}

void G4BREPSolid::TestSurfaceBBoxes(register const G4Ray& rayref) const 
{
  register G4Surface* srf;
  G4int active_srfs = nb_of_surfaces;
  
  // Do the bbox tests to all surfaces in List
  // for planar faces the intersection is instead evaluated.
  //
  G4int intersection=0;
  
  for(G4int a=0;a<nb_of_surfaces;a++)
  {
    // Get first in List
    //
    srf = SurfaceVec[a];
    
    if(srf->IsActive())
    {
      // Get type
      //
      if(srf->MyType() != 1) // 1 == planar face
      {
        if(srf->GetBBox()->Test(rayref))
          srf->SetDistance(bbox->GetDistance());
        else
        {
          // Test failed. Flag as inactive.
          //
          srf->Deactivate();
          active_srfs--;
        }
      }
      else
      {
        // Type was convex planar face
        intersection = srf->Intersect(rayref);

        if(!intersection) 
          active_srfs--;
      }
    }
    else
      active_srfs--;
  }
  
  if(!active_srfs) Active(0);
}


G4int G4BREPSolid::Intersect(register const G4Ray& rayref) const
{
  // Gets the roughly calculated closest intersection point for
  // a b_spline & accurate point for others.

  const G4double sqrHalfTolerance = kCarTolerance*kCarTolerance*0.25;

  register G4Surface* srf;
  G4double HitDistance = -1;
  const G4Point3D& RayStart = rayref.GetStart();
  const G4Point3D& RayDir   = rayref.GetDir();

  G4int result=1;
   
  // Sort List of active surfaces according to
  // bbox distances to ray starting point
  //
  QuickSort(SurfaceVec, 0, nb_of_surfaces-1);
  G4int Number=0;
   
  // Start handling active surfaces in order
  //
  for(register G4int a=0;a<nb_of_surfaces;a++)
  {
    srf = SurfaceVec[a];
    
    if(srf->IsActive())
    {
      result = srf->Intersect(rayref);
      if(result)
      {
        // Get the evaluated point on the surface
        //
        const G4Point3D& closest_point = srf->GetClosestHit();

        // Test for DistanceToIn(pt, vec)
        // if d = 0 and vec.norm > 0, do not see the surface
        //
        if( !( (srf->GetDistance() < sqrHalfTolerance)  || 
             (RayDir.dot(srf->SurfaceNormal(closest_point)) > 0) ) )
        {
  
          if(srf->MyType()==1)
            HitDistance = srf->GetDistance();
          else
          {
            // Check if the evaluated point is in front of the 
            // bbox of the next surface.
            // 
            HitDistance = RayStart.distance2(closest_point);
          }
        }
      }   
      else  // No hit
      {
        srf->Deactivate();
      }
    }
    Number++;
  }
  
  if(HitDistance < 0)
    return 0;
  
  QuickSort(SurfaceVec, 0, nb_of_surfaces-1);
  
  if(!(SurfaceVec[0]->IsActive()))
    return 0;     
  
  ((G4BREPSolid*)this)->intersection_point = SurfaceVec[0]->GetClosestHit();
  bbox->SetDistance(HitDistance);

  return 1;
}

G4int G4BREPSolid::FinalEvaluation(register const G4Ray& rayref, 
                                                  G4int ToIn    ) const
{
  const G4double sqrHalfTolerance = kCarTolerance*kCarTolerance*0.25;
  register G4Surface* srf;
  G4double Dist=0;

  ((G4BREPSolid*)this)->intersectionDistance = kInfinity;
  
  for(register G4int a=0;a<nb_of_surfaces;a++)
  {
    srf = SurfaceVec[a];
    
    if(srf->IsActive())
    {
      const G4Point3D& srf_intersection = srf->Evaluation(rayref);
      
      // Compute hit point distance from ray starting point
      //
      if(srf->MyType() != 1)
      {
        G4Point3D start = rayref.GetStart();
        Dist = srf_intersection.distance2(start);  
      }
      else 
        Dist = srf->GetDistance(); 

      // Skip point wich are on the surface i.e. within tolerance of the
      // surface. Special handling for DistanceToIn & reflections
      //
      if(Dist < sqrHalfTolerance)
      {
        if(ToIn) 
        {
          const G4Vector3D& Dir = rayref.GetDir();
          const G4Point3D& Hit = srf->GetClosestHit();
          const G4Vector3D& Norm = srf->SurfaceNormal(Hit);

          if(( Dir * Norm ) >= 0)
          {
            Dist = kInfinity;
            srf->Deactivate();
          }

          // else continue with the distance, even though < tolerance
        }
        else
        {
          Dist = kInfinity;
          srf->Deactivate();
        }
      }
      
      // If more than one surfaces are evaluated till the
      // final stage, only the closest point is taken
      //
      if(Dist < intersectionDistance)
      {
        // Check that Hit is in the direction of the ray 
        // from the starting point
        //
        const G4Point3D& Pt = rayref.GetStart();
        const G4Vector3D& Dir = rayref.GetDir();

        G4Point3D TestPoint = (0.00001*Dir) + Pt;
        G4double TestDistance = srf_intersection.distance2(TestPoint);

        if(TestDistance > Dist)
        {
          // Hit behind ray starting point, no intersection
          //
          Dist = kInfinity;
          srf->Deactivate();
        }
        else
        {
          ((G4BREPSolid*)this)->intersectionDistance = Dist;
          ((G4BREPSolid*)this)->intersection_point = srf_intersection;
        }

        // Check that the intersection is closer than the
        // next surfaces approximated point
        //
        if(srf->IsActive())
        {
          if(a+1<nb_of_surfaces)
          {
            const G4Point3D& Hit = srf->GetClosestHit();
            const G4Vector3D& Norm = srf->SurfaceNormal(Hit);

            // L. Broglia
            //if(( Dir * Norm ) >= 0)
            if(( Dir * Norm ) < 0)
            {
              Dist = kInfinity;
              srf->Deactivate();
            }

            // else continue with the distance, even though < tolerance

            ShortestDistance = Dist;
          }
          else
          {
            ShortestDistance = Dist;
            return 1;
          }
        }
      }
    }
    else  // if srf NOT active
    {
     /* if(intersectionDistance < kInfinity)
        return 1;
      return 0;*/
    }
  }
  if(intersectionDistance < kInfinity)    
    return 1;
  
  return 0;
}
 
G4Point3D G4BREPSolid::Scope() const
{
  G4Point3D scope;
  G4Point3D Max = bbox->GetBoxMax();
  G4Point3D Min = bbox->GetBoxMin();  
  
  scope.setX(std::fabs(Max.x()) - std::fabs(Min.x()));
  scope.setY(std::fabs(Max.y()) - std::fabs(Min.y()));
  scope.setZ(std::fabs(Max.z()) - std::fabs(Min.z()));
  
  return scope;
}

std::ostream& G4BREPSolid::StreamInfo(std::ostream& os) const
{
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: " << GetEntityType() << "\n"
     << " Parameters: \n"
     << "   Number of solids: " << NumberOfSolids << "\n"
     << "-----------------------------------------------------------\n";

  return os;
}

G4Polyhedron* G4BREPSolid::GetPolyhedron () const
{
  if (!fpPolyhedron ||
      fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fpPolyhedron->GetNumberOfRotationSteps())
    {
      delete fpPolyhedron;
      fpPolyhedron = CreatePolyhedron();
    }
  return fpPolyhedron;
}
