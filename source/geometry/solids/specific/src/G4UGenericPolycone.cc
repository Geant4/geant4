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
// $Id:$
//
// 
// Implementation of G4UGenericPolycone wrapper class
// --------------------------------------------------------------------

#include "G4GenericPolycone.hh"
#include "G4UGenericPolycone.hh"
#include "G4Polyhedron.hh"

using CLHEP::twopi;

////////////////////////////////////////////////////////////////////////
//
// Constructor (generic parameters)
//
G4UGenericPolycone::G4UGenericPolycone(const G4String& name, 
                                             G4double phiStart,
                                             G4double phiTotal,
                                             G4int    numRZ,
                                       const G4double r[],
                                       const G4double z[]   )
  : G4USolid(name, new UGenericPolycone(name, phiStart, phiTotal, numRZ, r, z))
{ 
}


////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UGenericPolycone::G4UGenericPolycone(__void__& a)
  : G4USolid(a)
{
}


//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UGenericPolycone::~G4UGenericPolycone()
{
}


//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UGenericPolycone::G4UGenericPolycone(const G4UGenericPolycone &source)
  : G4USolid(source)
{
}


//////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UGenericPolycone&
G4UGenericPolycone::operator=(const G4UGenericPolycone &source)
{
  if (this == &source) return *this;
  
  G4USolid::operator=( source );
  
  return *this;
}

G4Polyhedron* G4UGenericPolycone::CreatePolyhedron() const
{


 // The following code prepares for:
    // HepPolyhedron::createPolyhedron(int Nnodes, int Nfaces,
    //                                  const double xyz[][3],
    //                                  const int faces_vec[][4])
    // Here is an extract from the header file HepPolyhedron.h:
    /**
     * Creates user defined polyhedron.
     * This function allows to the user to define arbitrary polyhedron.
     * The faces of the polyhedron should be either triangles or planar
     * quadrilateral. Nodes of a face are defined by indexes pointing to
     * the elements in the xyz array. Numeration of the elements in the
     * array starts from 1 (like in fortran). The indexes can be positive
     * or negative. Negative sign means that the corresponding edge is
     * invisible. The normal of the face should be directed to exterior
     * of the polyhedron. 
     * 
     * @param  Nnodes number of nodes
     * @param  Nfaces number of faces
     * @param  xyz    nodes
     * @param  faces_vec  faces (quadrilaterals or triangles)
     * @return status of the operation - is non-zero in case of problem
     */
    const G4int numSide =
          G4int(G4Polyhedron::GetNumberOfRotationSteps()
                * (GetEndPhi() - GetStartPhi()) / twopi) + 1;
    G4int nNodes;
    G4int nFaces;
    typedef G4double double3[3];
    double3* xyz;
    typedef G4int int4[4];
    int4* faces_vec;
    if (IsOpen())
    {
      // Triangulate open ends. Simple ear-chopping algorithm...
      // I'm not sure how robust this algorithm is (J.Allison).
      //
      std::vector<G4bool> chopped(GetNumRZCorner(), false);
      std::vector<G4int*> triQuads;
      G4int remaining = GetNumRZCorner();
      G4int iStarter = 0;
      while (remaining >= 3)
      {
        // Find unchopped corners...
        //
        G4int A = -1, B = -1, C = -1;
        G4int iStepper = iStarter;
        do
        {
          if (A < 0)      { A = iStepper; }
          else if (B < 0) { B = iStepper; }
          else if (C < 0) { C = iStepper; }
          do
          {
            if (++iStepper >= GetNumRZCorner()) { iStepper = 0; }
          }
          while (chopped[iStepper]);
        }
        while (C < 0 && iStepper != iStarter);

        // Check triangle at B is pointing outward (an "ear").
        // Sign of z cross product determines...
        //
        G4double BAr = GetCorner(A).r - GetCorner(B).r;
        G4double BAz = GetCorner(A).z - GetCorner(B).z;
        G4double BCr = GetCorner(C).r - GetCorner(B).r;
        G4double BCz = GetCorner(C).z - GetCorner(B).z;
        if (BAr * BCz - BAz * BCr < kCarTolerance)
        {
          G4int* tq = new G4int[3];
          tq[0] = A + 1;
          tq[1] = B + 1;
          tq[2] = C + 1;
          triQuads.push_back(tq);
          chopped[B] = true;
          --remaining;
        }
        else
        {
          do
          {
            if (++iStarter >= GetNumRZCorner()) { iStarter = 0; }
          }
          while (chopped[iStarter]);
        }
      }
      // Transfer to faces...
      //
      nNodes = (numSide + 1) * GetNumRZCorner();
      nFaces = numSide * GetNumRZCorner() + 2 * triQuads.size();
      faces_vec = new int4[nFaces];
      G4int iface = 0;
      G4int addition = GetNumRZCorner() * numSide;
      G4int d = GetNumRZCorner() - 1;
      for (G4int iEnd = 0; iEnd < 2; ++iEnd)
      {
        for (size_t i = 0; i < triQuads.size(); ++i)
        {
          // Negative for soft/auxiliary/normally invisible edges...
          //
          G4int a, b, c;
          if (iEnd == 0)
          {
            a = triQuads[i][0];
            b = triQuads[i][1];
            c = triQuads[i][2];
          }
          else
          {
            a = triQuads[i][0] + addition;
            b = triQuads[i][2] + addition;
            c = triQuads[i][1] + addition;
          }
          G4int ab = std::abs(b - a);
          G4int bc = std::abs(c - b);
          G4int ca = std::abs(a - c);
          faces_vec[iface][0] = (ab == 1 || ab == d)? a: -a;
          faces_vec[iface][1] = (bc == 1 || bc == d)? b: -b;
          faces_vec[iface][2] = (ca == 1 || ca == d)? c: -c;
          faces_vec[iface][3] = 0;
          ++iface;
        }
      }

      // Continue with sides...

      xyz = new double3[nNodes];
      const G4double dPhi = (GetEndPhi() - GetStartPhi()) / numSide;
      G4double phi = GetStartPhi();
      G4int ixyz = 0;
      for (G4int iSide = 0; iSide < numSide; ++iSide)
      {
        for (G4int iCorner = 0; iCorner < GetNumRZCorner(); ++iCorner)
        {
          xyz[ixyz][0] = GetCorner(iCorner).r * std::cos(phi);
          xyz[ixyz][1] = GetCorner(iCorner).r * std::sin(phi);
          xyz[ixyz][2] = GetCorner(iCorner).z;
          if (iSide == 0)   // startPhi
          {
            if (iCorner < GetNumRZCorner() - 1)
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = -(ixyz + GetNumRZCorner() + 1);
              faces_vec[iface][2] = ixyz + GetNumRZCorner() + 2;
              faces_vec[iface][3] = ixyz + 2;
            }
            else
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = -(ixyz + GetNumRZCorner() + 1);
              faces_vec[iface][2] = ixyz + 2;
              faces_vec[iface][3] = ixyz - GetNumRZCorner() + 2;
            }
          }
          else if (iSide == numSide - 1)   // endPhi
          {
            if (iCorner < GetNumRZCorner() - 1)
              {
                faces_vec[iface][0] = ixyz + 1;
                faces_vec[iface][1] = ixyz + GetNumRZCorner() + 1;
                faces_vec[iface][2] = ixyz + GetNumRZCorner() + 2;
                faces_vec[iface][3] = -(ixyz + 2);
              }
            else
              {
                faces_vec[iface][0] = ixyz + 1;
                faces_vec[iface][1] = ixyz + GetNumRZCorner() + 1;
                faces_vec[iface][2] = ixyz + 2;
                faces_vec[iface][3] = -(ixyz - GetNumRZCorner() + 2);
              }
          }
          else
          {
            if (iCorner < GetNumRZCorner() - 1)
              {
                faces_vec[iface][0] = ixyz + 1;
                faces_vec[iface][1] = -(ixyz + GetNumRZCorner() + 1);
                faces_vec[iface][2] = ixyz + GetNumRZCorner() + 2;
                faces_vec[iface][3] = -(ixyz + 2);
              }
              else
              {
                faces_vec[iface][0] = ixyz + 1;
                faces_vec[iface][1] = -(ixyz + GetNumRZCorner() + 1);
                faces_vec[iface][2] = ixyz + 2;
                faces_vec[iface][3] = -(ixyz - GetNumRZCorner() + 2);
              }
            }
            ++iface;
            ++ixyz;
        }
        phi += dPhi;
      }

      // Last corners...

      for (G4int iCorner = 0; iCorner < GetNumRZCorner(); ++iCorner)
      {
        xyz[ixyz][0] = GetCorner(iCorner).r * std::cos(phi);
        xyz[ixyz][1] = GetCorner(iCorner).r * std::sin(phi);
        xyz[ixyz][2] = GetCorner(iCorner).z;
        ++ixyz;
      }
    }
    else  // !phiIsOpen - i.e., a complete 360 degrees.
    {
      nNodes = numSide * GetNumRZCorner();
      nFaces = numSide * GetNumRZCorner();;
      xyz = new double3[nNodes];
      faces_vec = new int4[nFaces];
      const G4double dPhi = (GetEndPhi() - GetStartPhi()) / numSide;
      G4double phi = GetStartPhi();
      G4int ixyz = 0, iface = 0;
      for (G4int iSide = 0; iSide < numSide; ++iSide)
      {
        for (G4int iCorner = 0; iCorner < GetNumRZCorner(); ++iCorner)
        {
          xyz[ixyz][0] = GetCorner(iCorner).r * std::cos(phi);
          xyz[ixyz][1] = GetCorner(iCorner).r * std::sin(phi);
          xyz[ixyz][2] = GetCorner(iCorner).z;

          if (iSide < numSide - 1)
          {
            if (iCorner < GetNumRZCorner() - 1)
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = -(ixyz + GetNumRZCorner() + 1);
              faces_vec[iface][2] = ixyz + GetNumRZCorner() + 2;
              faces_vec[iface][3] = -(ixyz + 2);
            }
            else
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = -(ixyz + GetNumRZCorner() + 1);
              faces_vec[iface][2] = ixyz + 2;
              faces_vec[iface][3] = -(ixyz - GetNumRZCorner() + 2);
            }
          }
          else   // Last side joins ends...
          {
            if (iCorner < GetNumRZCorner() - 1)
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = -(ixyz + GetNumRZCorner() - nFaces + 1);
              faces_vec[iface][2] = ixyz + GetNumRZCorner() - nFaces + 2;
              faces_vec[iface][3] = -(ixyz + 2);
            }
            else
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = -(ixyz - nFaces + GetNumRZCorner() + 1);
              faces_vec[iface][2] = ixyz - nFaces + 2;
              faces_vec[iface][3] = -(ixyz - GetNumRZCorner() + 2);
            }
          }
          ++ixyz;
          ++iface;
        }
        phi += dPhi;
      }
    }
    G4Polyhedron* polyhedron = new G4Polyhedron;
    G4int problem = polyhedron->createPolyhedron(nNodes, nFaces, xyz, faces_vec);
    delete [] faces_vec;
    delete [] xyz;
    if (problem)
    {
      std::ostringstream message;
      message << "Problem creating G4Polyhedron for: " << GetName();
      G4Exception("G4GenericPolycone::CreatePolyhedron()", "GeomSolids1002",
                  JustWarning, message);
      delete polyhedron;
      return 0;
    }
    else
    {
      return polyhedron;
    }
}

