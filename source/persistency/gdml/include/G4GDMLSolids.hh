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
// $Id: G4GDMLSolids.hh,v 1.10 2007-11-23 10:31:55 ztorzsok Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include <xercesc/dom/DOM.hpp>

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

#include "G4GDMLDefine.hh"
#include "G4GDMLEvaluator.hh"

class G4GDMLSolids {

   G4String file;

   G4GDMLEvaluator* evaluator;

   enum BooleanOp {UNION,SUBTRACTION,INTERSECTION};
   typedef struct zplaneType { double rmin,rmax,z; };

   G4bool booleanRead       (const xercesc::DOMElement* const,const BooleanOp);
   G4bool boxRead           (const xercesc::DOMElement* const);
   G4bool coneRead          (const xercesc::DOMElement* const);
   G4bool ellipsoidRead     (const xercesc::DOMElement* const);
   G4bool eltubeRead        (const xercesc::DOMElement* const);
   G4bool hypeRead          (const xercesc::DOMElement* const);
   G4bool loopRead          (const xercesc::DOMElement* const);
   G4bool orbRead           (const xercesc::DOMElement* const);
   G4bool paraRead          (const xercesc::DOMElement* const);
   G4bool polyconeRead      (const xercesc::DOMElement* const);
   G4bool polyhedraRead     (const xercesc::DOMElement* const);
   G4bool positionRead      (const xercesc::DOMElement* const,G4ThreeVector&);
   G4bool quadrangularRead  (const xercesc::DOMElement* const,G4TessellatedSolid*);
   G4bool refRead           (const xercesc::DOMElement* const,G4String&);
   G4bool reflectedSolidRead(const xercesc::DOMElement* const);
   G4bool rotationRead      (const xercesc::DOMElement* const,G4ThreeVector&);
   G4bool sectionRead       (const xercesc::DOMElement* const,G4ExtrudedSolid::ZSection&,const G4String& lunit);
   G4bool sphereRead        (const xercesc::DOMElement* const);
   G4bool tessellatedRead   (const xercesc::DOMElement* const);
   G4bool tetRead           (const xercesc::DOMElement* const);
   G4bool torusRead         (const xercesc::DOMElement* const);
   G4bool trapRead          (const xercesc::DOMElement* const);
   G4bool trdRead           (const xercesc::DOMElement* const);
   G4bool triangularRead    (const xercesc::DOMElement* const,G4TessellatedSolid*);
   G4bool tubeRead          (const xercesc::DOMElement* const);
   G4bool twoDimVertexRead  (const xercesc::DOMElement* const,G4TwoVector&,const G4String& lunit);
   G4bool xtruRead          (const xercesc::DOMElement* const);
   G4bool zplaneRead        (const xercesc::DOMElement* const,zplaneType&,const G4String& lunit);
public:
   G4GDMLDefine define;   

   G4bool Read(const xercesc::DOMElement* const,G4GDMLEvaluator*,const G4String& file0);
   G4VSolid* getSolid(const G4String&) const;

   std::string nameProcess(const std::string&);
};

#endif
