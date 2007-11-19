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

   G4String prename,postname;

   G4GDMLEvaluator* evaluator;

   enum BooleanOp {UNION,SUBTRACTION,INTERSECTION};
   typedef struct zplaneType { double rmin,rmax,z; };

   bool booleanRead       (const xercesc::DOMElement* const,const BooleanOp);
   bool boxRead           (const xercesc::DOMElement* const);
   bool coneRead          (const xercesc::DOMElement* const);
   bool ellipsoidRead     (const xercesc::DOMElement* const);
   bool eltubeRead        (const xercesc::DOMElement* const);
   bool hypeRead          (const xercesc::DOMElement* const);
   bool loopRead          (const xercesc::DOMElement* const);
   bool orbRead           (const xercesc::DOMElement* const);
   bool paraRead          (const xercesc::DOMElement* const);
   bool polyconeRead      (const xercesc::DOMElement* const);
   bool polyhedraRead     (const xercesc::DOMElement* const);
   bool positionRead      (const xercesc::DOMElement* const,G4ThreeVector&);
   bool quadrangularRead  (const xercesc::DOMElement* const,G4TessellatedSolid*);
   bool refRead           (const xercesc::DOMElement* const,G4String&);
   bool reflectedSolidRead(const xercesc::DOMElement* const);
   bool rotationRead      (const xercesc::DOMElement* const,G4ThreeVector&);
   bool sectionRead       (const xercesc::DOMElement* const,G4ExtrudedSolid::ZSection&,const G4String& lunit);
   bool sphereRead        (const xercesc::DOMElement* const);
   bool tessellatedRead   (const xercesc::DOMElement* const);
   bool tetRead           (const xercesc::DOMElement* const);
   bool torusRead         (const xercesc::DOMElement* const);
   bool trapRead          (const xercesc::DOMElement* const);
   bool trdRead           (const xercesc::DOMElement* const);
   bool triangularRead    (const xercesc::DOMElement* const,G4TessellatedSolid*);
   bool tubeRead          (const xercesc::DOMElement* const);
   bool twoDimVertexRead  (const xercesc::DOMElement* const,G4TwoVector&,const G4String& lunit);
   bool xtruRead          (const xercesc::DOMElement* const);
   bool zplaneRead        (const xercesc::DOMElement* const,zplaneType&,const G4String& lunit);
public:
   G4GDMLDefine define;   

   bool Read(const xercesc::DOMElement* const element,G4GDMLEvaluator*,const G4String& newModule);
   G4VSolid* getSolid(const G4String&) const;
};

#endif
