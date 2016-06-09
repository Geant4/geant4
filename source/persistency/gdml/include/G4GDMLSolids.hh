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
// $Id: G4GDMLSolids.hh,v 1.16 2007/11/30 14:51:20 ztorzsok Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
//
// class G4GDMLSolids
//
// Class description:
//
// GDML class for loading solids according to specifications in Geant4.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLSOLIDS_INCLUDED_
#define _G4GDMLSOLIDS_INCLUDED_

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Ellipsoid.hh"
#include "G4EllipticalTube.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Hype.hh"
#include "G4IntersectionSolid.hh"
#include "G4Orb.hh"
#include "G4Para.hh"
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
#include "G4UnionSolid.hh"

#include "G4GDMLMaterials.hh"

class G4GDMLSolids : public G4GDMLMaterials {
private:
   enum BooleanOp {UNION,SUBTRACTION,INTERSECTION};
   typedef struct zplaneType { G4double rmin,rmax,z; };

   void booleanRead(const xercesc::DOMElement* const,const BooleanOp);
   void boxRead(const xercesc::DOMElement* const);
   void coneRead(const xercesc::DOMElement* const);
   void ellipsoidRead(const xercesc::DOMElement* const);
   void eltubeRead(const xercesc::DOMElement* const);
   void hypeRead(const xercesc::DOMElement* const);
   void loopRead(const xercesc::DOMElement* const);
   void orbRead(const xercesc::DOMElement* const);
   void paraRead(const xercesc::DOMElement* const);
   void polyconeRead(const xercesc::DOMElement* const);
   void polyhedraRead(const xercesc::DOMElement* const);
   G4QuadrangularFacet* quadrangularRead(const xercesc::DOMElement* const);
   G4String refRead(const xercesc::DOMElement* const);
   void reflectedSolidRead(const xercesc::DOMElement* const);
   G4ExtrudedSolid::ZSection sectionRead(const xercesc::DOMElement* const,G4double);
   void sphereRead(const xercesc::DOMElement* const);
   void tessellatedRead(const xercesc::DOMElement* const);
   void tetRead(const xercesc::DOMElement* const);
   void torusRead(const xercesc::DOMElement* const);
   void trapRead(const xercesc::DOMElement* const);
   void trdRead(const xercesc::DOMElement* const);
   G4TriangularFacet* triangularRead(const xercesc::DOMElement* const);
   void tubeRead(const xercesc::DOMElement* const);
   G4TwoVector twoDimVertexRead(const xercesc::DOMElement* const,G4double);
   void xtruRead(const xercesc::DOMElement* const);
   zplaneType zplaneRead(const xercesc::DOMElement* const,G4double);
   void solidsRead(const xercesc::DOMElement* const);
protected:
   G4ThreeVector positionRead(const xercesc::DOMElement* const);
   G4ThreeVector rotationRead(const xercesc::DOMElement* const);
   G4ThreeVector scaleRead(const xercesc::DOMElement* const);
   G4VSolid* getSolid(const G4String&) const;
};

#endif
