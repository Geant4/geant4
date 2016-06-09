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
// $Id: G4GDMLReadSolids.hh,v 1.12 2008/11/21 09:32:46 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//
// class G4GDMLReadSolids
//
// Class description:
//
// GDML class for loading solids according to specifications in Geant4.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLREADSOLIDS_INCLUDED_
#define _G4GDMLREADSOLIDS_INCLUDED_

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Ellipsoid.hh"
#include "G4EllipticalCone.hh"
#include "G4EllipticalTube.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Hype.hh"
#include "G4IntersectionSolid.hh"
#include "G4Orb.hh"
#include "G4Para.hh"
#include "G4Paraboloid.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4QuadrangularFacet.hh"
#include "G4ReflectedSolid.hh"
#include "G4Sphere.hh"
#include "G4SolidStore.hh"
#include "G4SubtractionSolid.hh"
#include "G4TessellatedSolid.hh"
#include "G4Tet.hh"
#include "G4Torus.hh"
#include "G4Transform3D.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4TriangularFacet.hh"
#include "G4Tubs.hh"
#include "G4TwistedBox.hh"
#include "G4TwistedTrap.hh"
#include "G4TwistedTrd.hh"
#include "G4TwistedTubs.hh"
#include "G4UnionSolid.hh"
#include "G4OpticalSurface.hh"
#include "G4SurfaceProperty.hh"

#include "G4GDMLReadMaterials.hh"

class G4GDMLReadSolids : public G4GDMLReadMaterials
{
   enum BooleanOp {UNION,SUBTRACTION,INTERSECTION};
   typedef struct { G4double rmin,rmax,z; } zplaneType;

 public:

   G4VSolid* GetSolid(const G4String&) const;
   G4SurfaceProperty* GetSurfaceProperty(const G4String&) const;

   virtual void SolidsRead(const xercesc::DOMElement* const);

 protected:

   void BooleanRead(const xercesc::DOMElement* const,const BooleanOp);
   void BoxRead(const xercesc::DOMElement* const);
   void ConeRead(const xercesc::DOMElement* const);
   void ElconeRead(const xercesc::DOMElement* const);
   void EllipsoidRead(const xercesc::DOMElement* const);
   void EltubeRead(const xercesc::DOMElement* const);
   void XtruRead(const xercesc::DOMElement* const);
   void HypeRead(const xercesc::DOMElement* const);
   void OrbRead(const xercesc::DOMElement* const);
   void ParaRead(const xercesc::DOMElement* const);
   void ParaboloidRead(const xercesc::DOMElement* const);
   void PolyconeRead(const xercesc::DOMElement* const);
   void PolyhedraRead(const xercesc::DOMElement* const);
   G4QuadrangularFacet* QuadrangularRead(const xercesc::DOMElement* const);
   void ReflectedSolidRead(const xercesc::DOMElement* const);
   G4ExtrudedSolid::ZSection SectionRead(const xercesc::DOMElement* const,G4double);
   void SphereRead(const xercesc::DOMElement* const);
   void TessellatedRead(const xercesc::DOMElement* const);
   void TetRead(const xercesc::DOMElement* const);
   void TorusRead(const xercesc::DOMElement* const);
   void TrapRead(const xercesc::DOMElement* const);
   void TrdRead(const xercesc::DOMElement* const);
   void TubeRead(const xercesc::DOMElement* const);
   void TwistedboxRead(const xercesc::DOMElement* const);
   void TwistedtrapRead(const xercesc::DOMElement* const);
   void TwistedtrdRead(const xercesc::DOMElement* const);
   void TwistedtubsRead(const xercesc::DOMElement* const);
   G4TriangularFacet* TriangularRead(const xercesc::DOMElement* const);
   G4TwoVector TwoDimVertexRead(const xercesc::DOMElement* const,G4double);
   zplaneType ZplaneRead(const xercesc::DOMElement* const);
   void OpticalsurfaceRead(const xercesc::DOMElement* const);
};

#endif
