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
// $Id: G4GDMLWriteSolids.cc 103463 2017-04-11 07:22:55Z gcosmo $
//
// class G4GDMLWriteSolids Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLWriteSolids.hh"

#include "G4SystemOfUnits.hh"
#include "G4BooleanSolid.hh"
#include "G4ScaledSolid.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Ellipsoid.hh"
#include "G4EllipticalCone.hh"
#include "G4EllipticalTube.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Hype.hh"
#include "G4Orb.hh"
#include "G4Para.hh"
#include "G4Paraboloid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Polycone.hh"
#include "G4GenericPolycone.hh"
#include "G4Polyhedra.hh"
#include "G4ReflectedSolid.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4GenericTrap.hh"
#include "G4TessellatedSolid.hh"
#include "G4Tet.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4TwistedBox.hh"
#include "G4TwistedTrap.hh"
#include "G4TwistedTrd.hh"
#include "G4TwistedTubs.hh"
#include "G4UnionSolid.hh"
#include "G4OpticalSurface.hh"
#include "G4SurfaceProperty.hh"
#include "G4MaterialPropertiesTable.hh"

G4GDMLWriteSolids::G4GDMLWriteSolids()
  : G4GDMLWriteMaterials(), solidsElement(0)
{
}

G4GDMLWriteSolids::~G4GDMLWriteSolids()
{
}

void G4GDMLWriteSolids::
MultiUnionWrite(xercesc::DOMElement* solElement,
                const G4MultiUnion* const munionSolid)
{
   G4int numSolids=munionSolid->GetNumberOfSolids();
   G4String tag("multiUnion");
   
   G4VSolid* solid;
   G4Transform3D transform;

   const G4String& name = GenerateName(munionSolid->GetName(),munionSolid);
   xercesc::DOMElement* multiUnionElement = NewElement(tag);
   multiUnionElement->setAttributeNode(NewAttribute("name",name));

   for (G4int i=0; i<numSolids; ++i)
   {
      solid = munionSolid->GetSolid(i);
      transform = munionSolid->GetTransformation(i);

      HepGeom::Rotate3D rot3d;
      HepGeom::Translate3D transl ;
      HepGeom::Scale3D scale;
      transform.getDecomposition(scale,rot3d,transl); 

      G4ThreeVector pos = transl.getTranslation();
      G4RotationMatrix 
        rotm(CLHEP::HepRep3x3(rot3d.xx(), rot3d.xy(), rot3d.xz(),
                              rot3d.yx(), rot3d.yy(), rot3d.yz(),
                              rot3d.zx(), rot3d.zy(), rot3d.zz()));
      G4ThreeVector rot = GetAngles(rotm);

      AddSolid(solid);
      const G4String& solidref = GenerateName(solid->GetName(),solid);
      std::ostringstream os; os << i+1;
      const G4String& nodeName = "Node-" + G4String(os.str());
      xercesc::DOMElement* solidElement = NewElement("solid");
      solidElement->setAttributeNode(NewAttribute("ref",solidref));
      xercesc::DOMElement* multiUnionNodeElement = NewElement("multiUnionNode");
      multiUnionNodeElement->setAttributeNode(NewAttribute("name", nodeName));
      multiUnionNodeElement->appendChild(solidElement); // Append solid to node
      if ( (std::fabs(pos.x()) > kLinearPrecision)
        || (std::fabs(pos.y()) > kLinearPrecision)
        || (std::fabs(pos.z()) > kLinearPrecision) )
      {
        PositionWrite(multiUnionNodeElement,name+"_pos",pos);
      }
      if ( (std::fabs(rot.x()) > kAngularPrecision)
        || (std::fabs(rot.y()) > kAngularPrecision)
        || (std::fabs(rot.z()) > kAngularPrecision) )
      {
        RotationWrite(multiUnionNodeElement,name+"_rot",rot);
      }
      multiUnionElement->appendChild(multiUnionNodeElement); // Append node
   }

   solElement->appendChild(multiUnionElement);
     // Add the multiUnion solid AFTER the constituent nodes!
}

void G4GDMLWriteSolids::
BooleanWrite(xercesc::DOMElement* solElement,
             const G4BooleanSolid* const boolean)
{
   G4int displaced=0;

   G4String tag("undefined");
   if (dynamic_cast<const G4IntersectionSolid*>(boolean))
     { tag = "intersection"; } else
   if (dynamic_cast<const G4SubtractionSolid*>(boolean))
     { tag = "subtraction"; } else
   if (dynamic_cast<const G4UnionSolid*>(boolean))
     { tag = "union"; }
   
   G4VSolid* firstPtr = const_cast<G4VSolid*>(boolean->GetConstituentSolid(0));
   G4VSolid* secondPtr = const_cast<G4VSolid*>(boolean->GetConstituentSolid(1));
   
   G4ThreeVector firstpos,firstrot,pos,rot;

   // Solve possible displacement of referenced solids!
   //
   while (true)
   {
      if ( displaced>8 )
      {
        G4String ErrorMessage = "The referenced solid '"
                              + firstPtr->GetName() +
                              + "in the Boolean shape '" +
                              + boolean->GetName() +
                              + "' was displaced too many times!";
        G4Exception("G4GDMLWriteSolids::BooleanWrite()",
                    "InvalidSetup", FatalException, ErrorMessage);
      }

      if (G4DisplacedSolid* disp = dynamic_cast<G4DisplacedSolid*>(firstPtr))
      {
         firstpos += disp->GetObjectTranslation();
         firstrot += GetAngles(disp->GetObjectRotation());
         firstPtr = disp->GetConstituentMovedSolid();
         displaced++;
         continue;
      }
      break;
   }
   displaced = 0;
   while (true)
   {
      if ( displaced>maxTransforms )
      {
        G4String ErrorMessage = "The referenced solid '"
                              + secondPtr->GetName() +
                              + "in the Boolean shape '" +
                              + boolean->GetName() +
                              + "' was displaced too many times!";
        G4Exception("G4GDMLWriteSolids::BooleanWrite()",
                    "InvalidSetup", FatalException, ErrorMessage);
      }

      if (G4DisplacedSolid* disp = dynamic_cast<G4DisplacedSolid*>(secondPtr))
      {
         pos += disp->GetObjectTranslation();
         rot += GetAngles(disp->GetObjectRotation());
         secondPtr = disp->GetConstituentMovedSolid();
         displaced++;
         continue;
      }
      break;
   }

   AddSolid(firstPtr);   // At first add the constituent solids!
   AddSolid(secondPtr);

   const G4String& name = GenerateName(boolean->GetName(),boolean);
   const G4String& firstref = GenerateName(firstPtr->GetName(),firstPtr);
   const G4String& secondref = GenerateName(secondPtr->GetName(),secondPtr);

   xercesc::DOMElement* booleanElement = NewElement(tag);
   booleanElement->setAttributeNode(NewAttribute("name",name));
   xercesc::DOMElement* firstElement = NewElement("first");
   firstElement->setAttributeNode(NewAttribute("ref",firstref));
   booleanElement->appendChild(firstElement);
   xercesc::DOMElement* secondElement = NewElement("second");
   secondElement->setAttributeNode(NewAttribute("ref",secondref));
   booleanElement->appendChild(secondElement);
   solElement->appendChild(booleanElement);
     // Add the boolean solid AFTER the constituent solids!

   if ( (std::fabs(pos.x()) > kLinearPrecision)
     || (std::fabs(pos.y()) > kLinearPrecision)
     || (std::fabs(pos.z()) > kLinearPrecision) )
   {
     PositionWrite(booleanElement,name+"_pos",pos);
   }

   if ( (std::fabs(rot.x()) > kAngularPrecision)
     || (std::fabs(rot.y()) > kAngularPrecision)
     || (std::fabs(rot.z()) > kAngularPrecision) )
   {
     RotationWrite(booleanElement,name+"_rot",rot);
   }

   if ( (std::fabs(firstpos.x()) > kLinearPrecision)
     || (std::fabs(firstpos.y()) > kLinearPrecision)
     || (std::fabs(firstpos.z()) > kLinearPrecision) )
   {
     FirstpositionWrite(booleanElement,name+"_fpos",firstpos);
   }

   if ( (std::fabs(firstrot.x()) > kAngularPrecision)
     || (std::fabs(firstrot.y()) > kAngularPrecision)
     || (std::fabs(firstrot.z()) > kAngularPrecision) )
   {
     FirstrotationWrite(booleanElement,name+"_frot",firstrot);
   }
}

void G4GDMLWriteSolids::
ScaledWrite(xercesc::DOMElement* solElement,
            const G4ScaledSolid* const scaled)
{
   G4String tag("scaledSolid");
   
   G4VSolid* solid = const_cast<G4VSolid*>(scaled->GetUnscaledSolid());
   G4Scale3D scale = scaled->GetScaleTransform();
   G4ThreeVector sclVector = G4ThreeVector(scale.xx(), scale.yy(), scale.zz());

   AddSolid(solid);   // Add the constituent solid!

   const G4String& name = GenerateName(scaled->GetName(),scaled);
   const G4String& solidref = GenerateName(solid->GetName(),solid);

   xercesc::DOMElement* scaledElement = NewElement(tag);
   scaledElement->setAttributeNode(NewAttribute("name",name));

   xercesc::DOMElement* solidElement = NewElement("solidref");
   solidElement->setAttributeNode(NewAttribute("ref",solidref));
   scaledElement->appendChild(solidElement);

   if ( (std::fabs(scale.xx()) > kLinearPrecision)
     && (std::fabs(scale.yy()) > kLinearPrecision)
     && (std::fabs(scale.zz()) > kLinearPrecision) )
   {
     ScaleWrite(scaledElement, name+"_scl", sclVector);
   }

   solElement->appendChild(scaledElement);
     // Add the scaled solid AFTER its constituent solid!
}

void G4GDMLWriteSolids::
BoxWrite(xercesc::DOMElement* solElement, const G4Box* const box)
{
   const G4String& name = GenerateName(box->GetName(),box);

   xercesc::DOMElement* boxElement = NewElement("box");
   boxElement->setAttributeNode(NewAttribute("name",name));
   boxElement->setAttributeNode(NewAttribute("x",2.0*box->GetXHalfLength()/mm));
   boxElement->setAttributeNode(NewAttribute("y",2.0*box->GetYHalfLength()/mm));
   boxElement->setAttributeNode(NewAttribute("z",2.0*box->GetZHalfLength()/mm));
   boxElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(boxElement);
}

void G4GDMLWriteSolids::
ConeWrite(xercesc::DOMElement* solElement, const G4Cons* const cone)
{
   const G4String& name = GenerateName(cone->GetName(),cone);

   xercesc::DOMElement* coneElement = NewElement("cone");
   coneElement->setAttributeNode(NewAttribute("name",name));
   coneElement->
     setAttributeNode(NewAttribute("rmin1",cone->GetInnerRadiusMinusZ()/mm));
   coneElement->
     setAttributeNode(NewAttribute("rmax1",cone->GetOuterRadiusMinusZ()/mm));
   coneElement->
     setAttributeNode(NewAttribute("rmin2",cone->GetInnerRadiusPlusZ()/mm));
   coneElement->
     setAttributeNode(NewAttribute("rmax2",cone->GetOuterRadiusPlusZ()/mm));
   coneElement->
     setAttributeNode(NewAttribute("z",2.0*cone->GetZHalfLength()/mm));
   coneElement->
     setAttributeNode(NewAttribute("startphi",cone->GetStartPhiAngle()/degree));
   coneElement->
     setAttributeNode(NewAttribute("deltaphi",cone->GetDeltaPhiAngle()/degree));
   coneElement->setAttributeNode(NewAttribute("aunit","deg"));
   coneElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(coneElement);
}

void G4GDMLWriteSolids::
ElconeWrite(xercesc::DOMElement* solElement,
            const G4EllipticalCone* const elcone)
{
   const G4String& name = GenerateName(elcone->GetName(),elcone);

   xercesc::DOMElement* elconeElement = NewElement("elcone");
   elconeElement->setAttributeNode(NewAttribute("name",name));
   elconeElement->setAttributeNode(NewAttribute("dx",elcone->GetSemiAxisX()/mm));
   elconeElement->setAttributeNode(NewAttribute("dy",elcone->GetSemiAxisY()/mm));
   elconeElement->setAttributeNode(NewAttribute("zmax",elcone->GetZMax()/mm));
   elconeElement->setAttributeNode(NewAttribute("zcut",elcone->GetZTopCut()/mm));
   elconeElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(elconeElement);
}

void G4GDMLWriteSolids::
EllipsoidWrite(xercesc::DOMElement* solElement,
               const G4Ellipsoid* const ellipsoid)
{
   const G4String& name = GenerateName(ellipsoid->GetName(),ellipsoid);

   xercesc::DOMElement* ellipsoidElement = NewElement("ellipsoid");
   ellipsoidElement->setAttributeNode(NewAttribute("name",name));
   ellipsoidElement->
     setAttributeNode(NewAttribute("ax",ellipsoid->GetSemiAxisMax(0)/mm));
   ellipsoidElement->
     setAttributeNode(NewAttribute("by",ellipsoid->GetSemiAxisMax(1)/mm));
   ellipsoidElement->
     setAttributeNode(NewAttribute("cz",ellipsoid->GetSemiAxisMax(2)/mm));
   ellipsoidElement->
     setAttributeNode(NewAttribute("zcut1",ellipsoid->GetZBottomCut()/mm));
   ellipsoidElement->
     setAttributeNode(NewAttribute("zcut2",ellipsoid->GetZTopCut()/mm));
   ellipsoidElement->
     setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(ellipsoidElement);
}

void G4GDMLWriteSolids::
EltubeWrite(xercesc::DOMElement* solElement,
            const G4EllipticalTube* const eltube)
{
   const G4String& name = GenerateName(eltube->GetName(),eltube);

   xercesc::DOMElement* eltubeElement = NewElement("eltube");
   eltubeElement->setAttributeNode(NewAttribute("name",name));
   eltubeElement->setAttributeNode(NewAttribute("dx",eltube->GetDx()/mm));
   eltubeElement->setAttributeNode(NewAttribute("dy",eltube->GetDy()/mm));
   eltubeElement->setAttributeNode(NewAttribute("dz",eltube->GetDz()/mm));
   eltubeElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(eltubeElement);
}

void G4GDMLWriteSolids::
XtruWrite(xercesc::DOMElement* solElement,
          const G4ExtrudedSolid* const xtru)
{
   const G4String& name = GenerateName(xtru->GetName(),xtru);

   xercesc::DOMElement* xtruElement = NewElement("xtru");
   xtruElement->setAttributeNode(NewAttribute("name",name));
   xtruElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(xtruElement);

   const G4int NumVertex = xtru->GetNofVertices();
   
   for (G4int i=0;i<NumVertex;i++)
   {
      xercesc::DOMElement* twoDimVertexElement = NewElement("twoDimVertex");
      xtruElement->appendChild(twoDimVertexElement);

      const G4TwoVector& vertex = xtru->GetVertex(i);

      twoDimVertexElement->setAttributeNode(NewAttribute("x",vertex.x()/mm));
      twoDimVertexElement->setAttributeNode(NewAttribute("y",vertex.y()/mm));
   }

   const G4int NumSection = xtru->GetNofZSections();

   for (G4int i=0;i<NumSection;i++)
   {
      xercesc::DOMElement* sectionElement = NewElement("section");
      xtruElement->appendChild(sectionElement);

      const G4ExtrudedSolid::ZSection section = xtru->GetZSection(i);

      sectionElement->setAttributeNode(NewAttribute("zOrder",i));
      sectionElement->setAttributeNode(NewAttribute("zPosition",section.fZ/mm));
      sectionElement->
        setAttributeNode(NewAttribute("xOffset",section.fOffset.x()/mm));
      sectionElement->
        setAttributeNode(NewAttribute("yOffset",section.fOffset.y()/mm));
      sectionElement->
        setAttributeNode(NewAttribute("scalingFactor",section.fScale));
   }
}

void G4GDMLWriteSolids::
HypeWrite(xercesc::DOMElement* solElement, const G4Hype* const hype)
{
   const G4String& name = GenerateName(hype->GetName(),hype);

   xercesc::DOMElement* hypeElement = NewElement("hype");
   hypeElement->setAttributeNode(NewAttribute("name",name));
   hypeElement->setAttributeNode(NewAttribute("rmin",
                hype->GetInnerRadius()/mm));
   hypeElement->setAttributeNode(NewAttribute("rmax",
                hype->GetOuterRadius()/mm));
   hypeElement->setAttributeNode(NewAttribute("inst",
                hype->GetInnerStereo()/degree));
   hypeElement->setAttributeNode(NewAttribute("outst",
                hype->GetOuterStereo()/degree));
   hypeElement->setAttributeNode(NewAttribute("z",
                2.0*hype->GetZHalfLength()/mm));
   hypeElement->setAttributeNode(NewAttribute("aunit","deg"));
   hypeElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(hypeElement);
}

void G4GDMLWriteSolids::
OrbWrite(xercesc::DOMElement* solElement, const G4Orb* const orb)
{
   const G4String& name = GenerateName(orb->GetName(),orb);

   xercesc::DOMElement* orbElement = NewElement("orb");
   orbElement->setAttributeNode(NewAttribute("name",name));
   orbElement->setAttributeNode(NewAttribute("r",orb->GetRadius()/mm));
   orbElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(orbElement);
}

void G4GDMLWriteSolids::
ParaWrite(xercesc::DOMElement* solElement, const G4Para* const para)
{
   const G4String& name = GenerateName(para->GetName(),para);

   const G4ThreeVector simaxis = para->GetSymAxis();
   const G4double alpha = std::atan(para->GetTanAlpha());
   const G4double phi = simaxis.phi();
   const G4double theta = simaxis.theta();

   xercesc::DOMElement* paraElement = NewElement("para");
   paraElement->setAttributeNode(NewAttribute("name",name));
   paraElement->setAttributeNode(NewAttribute("x",
                2.0*para->GetXHalfLength()/mm));
   paraElement->setAttributeNode(NewAttribute("y",
                2.0*para->GetYHalfLength()/mm));
   paraElement->setAttributeNode(NewAttribute("z",
                2.0*para->GetZHalfLength()/mm));
   paraElement->setAttributeNode(NewAttribute("alpha",alpha/degree));
   paraElement->setAttributeNode(NewAttribute("theta",theta/degree));
   paraElement->setAttributeNode(NewAttribute("phi",phi/degree));
   paraElement->setAttributeNode(NewAttribute("aunit","deg"));
   paraElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(paraElement);
}

void G4GDMLWriteSolids::
ParaboloidWrite(xercesc::DOMElement* solElement,
                const G4Paraboloid* const paraboloid)
{
   const G4String& name = GenerateName(paraboloid->GetName(),paraboloid);

   xercesc::DOMElement* paraboloidElement = NewElement("paraboloid");
   paraboloidElement->setAttributeNode(NewAttribute("name",name));
   paraboloidElement->setAttributeNode(NewAttribute("rlo",
                      paraboloid->GetRadiusMinusZ()/mm));
   paraboloidElement->setAttributeNode(NewAttribute("rhi",
                      paraboloid->GetRadiusPlusZ()/mm));
   paraboloidElement->setAttributeNode(NewAttribute("dz",
                      paraboloid->GetZHalfLength()/mm));
   paraboloidElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(paraboloidElement);
}
void G4GDMLWriteSolids::
PolyconeWrite(xercesc::DOMElement* solElement,
              const G4Polycone* const polycone)
{
   const G4String& name = GenerateName(polycone->GetName(),polycone);
  
    xercesc::DOMElement* polyconeElement = NewElement("polycone");
    polyconeElement->setAttributeNode(NewAttribute("name",name));
    polyconeElement->setAttributeNode(NewAttribute("startphi",
                    polycone->GetOriginalParameters()->Start_angle/degree));
    polyconeElement->setAttributeNode(NewAttribute("deltaphi",
                    polycone->GetOriginalParameters()->Opening_angle/degree));
    polyconeElement->setAttributeNode(NewAttribute("aunit","deg"));
    polyconeElement->setAttributeNode(NewAttribute("lunit","mm"));
    solElement->appendChild(polyconeElement);

    const size_t num_zplanes = polycone->GetOriginalParameters()->Num_z_planes;
    const G4double* z_array = polycone->GetOriginalParameters()->Z_values;
    const G4double* rmin_array = polycone->GetOriginalParameters()->Rmin;
    const G4double* rmax_array = polycone->GetOriginalParameters()->Rmax;

    for (size_t i=0; i<num_zplanes; i++)
    {
      ZplaneWrite(polyconeElement,z_array[i],rmin_array[i],rmax_array[i]);
    }
   
  
}

void G4GDMLWriteSolids::
GenericPolyconeWrite(xercesc::DOMElement* solElement,
              const G4GenericPolycone* const polycone)
{
   const G4String& name = GenerateName(polycone->GetName(),polycone);
   xercesc::DOMElement* polyconeElement = NewElement("genericPolycone");
   const G4double startPhi=polycone->GetStartPhi();
    polyconeElement->setAttributeNode(NewAttribute("name",name));
    polyconeElement->setAttributeNode(NewAttribute("startphi",
		    startPhi/degree));
    polyconeElement->setAttributeNode(NewAttribute("deltaphi",
		    (polycone->GetEndPhi()-startPhi)/degree));
    polyconeElement->setAttributeNode(NewAttribute("aunit","deg"));
    polyconeElement->setAttributeNode(NewAttribute("lunit","mm"));
    solElement->appendChild(polyconeElement);

    const size_t num_rzpoints = polycone->GetNumRZCorner();
    for (size_t i=0; i<num_rzpoints; i++)
    {
      const G4double r_point=polycone->GetCorner(i).r;
      const G4double z_point=polycone->GetCorner(i).z;
      RZPointWrite(polyconeElement,r_point,z_point);
    }
   
}

void G4GDMLWriteSolids::
PolyhedraWrite(xercesc::DOMElement* solElement,
               const G4Polyhedra* const polyhedra)
{
   const G4String& name = GenerateName(polyhedra->GetName(),polyhedra);
   if(polyhedra->IsGeneric() == false){
    xercesc::DOMElement* polyhedraElement = NewElement("polyhedra");
    polyhedraElement->setAttributeNode(NewAttribute("name",name));
    polyhedraElement->setAttributeNode(NewAttribute("startphi",
                     polyhedra->GetOriginalParameters()->Start_angle/degree));
    polyhedraElement->setAttributeNode(NewAttribute("deltaphi",
                     polyhedra->GetOriginalParameters()->Opening_angle/degree));
    polyhedraElement->setAttributeNode(NewAttribute("numsides",
                     polyhedra->GetOriginalParameters()->numSide));
    polyhedraElement->setAttributeNode(NewAttribute("aunit","deg"));
    polyhedraElement->setAttributeNode(NewAttribute("lunit","mm"));
    solElement->appendChild(polyhedraElement);

    const size_t num_zplanes = polyhedra->GetOriginalParameters()->Num_z_planes;
    const G4double* z_array = polyhedra->GetOriginalParameters()->Z_values;
    const G4double* rmin_array = polyhedra->GetOriginalParameters()->Rmin;
    const G4double* rmax_array = polyhedra->GetOriginalParameters()->Rmax;

    const G4double convertRad =
         std::cos(0.5*polyhedra->GetOriginalParameters()->Opening_angle
       / polyhedra->GetOriginalParameters()->numSide);

    for (size_t i=0;i<num_zplanes;i++)
    {
      ZplaneWrite(polyhedraElement,z_array[i],
                  rmin_array[i]*convertRad, rmax_array[i]*convertRad);
    }
   }else{//generic polyhedra
   
    xercesc::DOMElement* polyhedraElement = NewElement("genericPolyhedra");
    polyhedraElement->setAttributeNode(NewAttribute("name",name));
    polyhedraElement->setAttributeNode(NewAttribute("startphi",
                     polyhedra->GetOriginalParameters()->Start_angle/degree));
    polyhedraElement->setAttributeNode(NewAttribute("deltaphi",
                     polyhedra->GetOriginalParameters()->Opening_angle/degree));
    polyhedraElement->setAttributeNode(NewAttribute("numsides",
                     polyhedra->GetOriginalParameters()->numSide));
    polyhedraElement->setAttributeNode(NewAttribute("aunit","deg"));
    polyhedraElement->setAttributeNode(NewAttribute("lunit","mm"));
    solElement->appendChild(polyhedraElement);

    const size_t num_rzpoints = polyhedra->GetNumRZCorner();
    
    for (size_t i=0;i<num_rzpoints;i++)
    {
      const G4double r_point = polyhedra->GetCorner(i).r;
      const G4double z_point = polyhedra->GetCorner(i).z;
      RZPointWrite(polyhedraElement,r_point,z_point);
    }
  }
}

void G4GDMLWriteSolids::
SphereWrite(xercesc::DOMElement* solElement, const G4Sphere* const sphere)
{
   const G4String& name = GenerateName(sphere->GetName(),sphere);

   xercesc::DOMElement* sphereElement = NewElement("sphere");
   sphereElement->setAttributeNode(NewAttribute("name",name));
   sphereElement->setAttributeNode(NewAttribute("rmin",
                  sphere->GetInnerRadius()/mm));
   sphereElement->setAttributeNode(NewAttribute("rmax",
                  sphere->GetOuterRadius()/mm));
   sphereElement->setAttributeNode(NewAttribute("startphi",
                  sphere->GetStartPhiAngle()/degree));
   sphereElement->setAttributeNode(NewAttribute("deltaphi",
                  sphere->GetDeltaPhiAngle()/degree));
   sphereElement->setAttributeNode(NewAttribute("starttheta",
                  sphere->GetStartThetaAngle()/degree));
   sphereElement->setAttributeNode(NewAttribute("deltatheta",
                  sphere->GetDeltaThetaAngle()/degree));
   sphereElement->setAttributeNode(NewAttribute("aunit","deg"));
   sphereElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(sphereElement);
}

void G4GDMLWriteSolids::
TessellatedWrite(xercesc::DOMElement* solElement,
                 const G4TessellatedSolid* const tessellated)
{
   const G4String& solid_name = tessellated->GetName();
   const G4String& name = GenerateName(solid_name, tessellated);

   xercesc::DOMElement* tessellatedElement = NewElement("tessellated");
   tessellatedElement->setAttributeNode(NewAttribute("name",name));
   tessellatedElement->setAttributeNode(NewAttribute("aunit","deg"));
   tessellatedElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(tessellatedElement);

   std::map<G4ThreeVector, G4String, G4ThreeVectorCompare> vertexMap;

   const size_t NumFacets = tessellated->GetNumberOfFacets();
   size_t NumVertex = 0;
   
   for (size_t i=0;i<NumFacets;i++)
   {
      const G4VFacet* facet = tessellated->GetFacet(i);
      const size_t NumVertexPerFacet = facet->GetNumberOfVertices();

      G4String FacetTag;
      
      if (NumVertexPerFacet==3) { FacetTag="triangular"; } else
      if (NumVertexPerFacet==4) { FacetTag="quadrangular"; }
      else
      {
        G4Exception("G4GDMLWriteSolids::TessellatedWrite()", "InvalidSetup",
                    FatalException, "Facet should contain 3 or 4 vertices!");
      }

      xercesc::DOMElement* facetElement = NewElement(FacetTag);
      tessellatedElement->appendChild(facetElement);

      for (size_t j=0; j<NumVertexPerFacet; j++)
      {
         std::stringstream name_stream;
         std::stringstream ref_stream;

         name_stream << "vertex" << (j+1);
         ref_stream << solid_name << "_v" << NumVertex;

         const G4String& fname = name_stream.str();  // facet's tag variable
         G4String ref = ref_stream.str();     // vertex tag to be associated

         // Now search for the existance of the current vertex in the
         // map of cached vertices. If existing, do NOT store it as
         // position in the GDML file, so avoiding duplication; otherwise
         // cache it in the local map and add it as position in the
         // "define" section of the GDML file.

         const G4ThreeVector& vertex = facet->GetVertex(j);

         if(vertexMap.find(vertex) != vertexMap.end())  // Vertex is cached
         {
           ref = vertexMap[vertex];     // Set the proper tag for it
         }
         else                                           // Vertex not found
         {
           vertexMap.insert(std::make_pair(vertex,ref)); // Cache vertex and ...
           AddPosition(ref, vertex);    // ... add it to define section!
           NumVertex++;
         }

         // Now create association of the vertex with its facet
         //
         facetElement->setAttributeNode(NewAttribute(fname,ref));
      }
   }
}

void G4GDMLWriteSolids::
TetWrite(xercesc::DOMElement* solElement, const G4Tet* const tet)
{
   const G4String& solid_name = tet->GetName();
   const G4String& name = GenerateName(solid_name, tet);

   std::vector<G4ThreeVector> vertexList = tet->GetVertices();

   xercesc::DOMElement* tetElement = NewElement("tet");
   tetElement->setAttributeNode(NewAttribute("name",name));
   tetElement->setAttributeNode(NewAttribute("vertex1",solid_name+"_v1"));
   tetElement->setAttributeNode(NewAttribute("vertex2",solid_name+"_v2"));
   tetElement->setAttributeNode(NewAttribute("vertex3",solid_name+"_v3"));
   tetElement->setAttributeNode(NewAttribute("vertex4",solid_name+"_v4"));
   tetElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(tetElement);

   AddPosition(solid_name+"_v1",vertexList[0]);
   AddPosition(solid_name+"_v2",vertexList[1]);
   AddPosition(solid_name+"_v3",vertexList[2]);
   AddPosition(solid_name+"_v4",vertexList[3]);
}

void G4GDMLWriteSolids::
TorusWrite(xercesc::DOMElement* solElement, const G4Torus* const torus)
{
   const G4String& name = GenerateName(torus->GetName(),torus);

   xercesc::DOMElement* torusElement = NewElement("torus");
   torusElement->setAttributeNode(NewAttribute("name",name));
   torusElement->setAttributeNode(NewAttribute("rmin",torus->GetRmin()/mm));
   torusElement->setAttributeNode(NewAttribute("rmax",torus->GetRmax()/mm));
   torusElement->setAttributeNode(NewAttribute("rtor",torus->GetRtor()/mm));
   torusElement->
     setAttributeNode(NewAttribute("startphi",torus->GetSPhi()/degree));
   torusElement->
     setAttributeNode(NewAttribute("deltaphi",torus->GetDPhi()/degree));
   torusElement->setAttributeNode(NewAttribute("aunit","deg"));
   torusElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(torusElement);
}

void G4GDMLWriteSolids::
GenTrapWrite(xercesc::DOMElement* solElement,
             const G4GenericTrap* const gtrap)
{
   const G4String& name = GenerateName(gtrap->GetName(),gtrap);

   std::vector<G4TwoVector> vertices = gtrap->GetVertices();

   xercesc::DOMElement* gtrapElement = NewElement("arb8");
   gtrapElement->setAttributeNode(NewAttribute("name",name));
   gtrapElement->setAttributeNode(NewAttribute("dz",
                                           gtrap->GetZHalfLength()/mm));
   gtrapElement->setAttributeNode(NewAttribute("v1x", vertices[0].x()));
   gtrapElement->setAttributeNode(NewAttribute("v1y", vertices[0].y()));
   gtrapElement->setAttributeNode(NewAttribute("v2x", vertices[1].x()));
   gtrapElement->setAttributeNode(NewAttribute("v2y", vertices[1].y()));
   gtrapElement->setAttributeNode(NewAttribute("v3x", vertices[2].x()));
   gtrapElement->setAttributeNode(NewAttribute("v3y", vertices[2].y()));
   gtrapElement->setAttributeNode(NewAttribute("v4x", vertices[3].x()));
   gtrapElement->setAttributeNode(NewAttribute("v4y", vertices[3].y()));
   gtrapElement->setAttributeNode(NewAttribute("v5x", vertices[4].x()));
   gtrapElement->setAttributeNode(NewAttribute("v5y", vertices[4].y()));
   gtrapElement->setAttributeNode(NewAttribute("v6x", vertices[5].x()));
   gtrapElement->setAttributeNode(NewAttribute("v6y", vertices[5].y()));
   gtrapElement->setAttributeNode(NewAttribute("v7x", vertices[6].x()));
   gtrapElement->setAttributeNode(NewAttribute("v7y", vertices[6].y()));
   gtrapElement->setAttributeNode(NewAttribute("v8x", vertices[7].x()));
   gtrapElement->setAttributeNode(NewAttribute("v8y", vertices[7].y()));
   gtrapElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(gtrapElement);
}

void G4GDMLWriteSolids::
TrapWrite(xercesc::DOMElement* solElement, const G4Trap* const trap)
{
   const G4String& name = GenerateName(trap->GetName(),trap);

   const G4ThreeVector& simaxis = trap->GetSymAxis();
   const G4double phi = simaxis.phi();
   const G4double theta = simaxis.theta();
   const G4double alpha1 = std::atan(trap->GetTanAlpha1());
   const G4double alpha2 = std::atan(trap->GetTanAlpha2());

   xercesc::DOMElement* trapElement = NewElement("trap");
   trapElement->setAttributeNode(NewAttribute("name",name));
   trapElement->setAttributeNode(NewAttribute("z",
                2.0*trap->GetZHalfLength()/mm));
   trapElement->setAttributeNode(NewAttribute("theta",theta/degree));
   trapElement->setAttributeNode(NewAttribute("phi",phi/degree));
   trapElement->setAttributeNode(NewAttribute("y1",
                2.0*trap->GetYHalfLength1()/mm));
   trapElement->setAttributeNode(NewAttribute("x1",
                2.0*trap->GetXHalfLength1()/mm));
   trapElement->setAttributeNode(NewAttribute("x2",
                2.0*trap->GetXHalfLength2()/mm));
   trapElement->setAttributeNode(NewAttribute("alpha1",alpha1/degree));
   trapElement->setAttributeNode(NewAttribute("y2",
                2.0*trap->GetYHalfLength2()/mm));
   trapElement->setAttributeNode(NewAttribute("x3",
                2.0*trap->GetXHalfLength3()/mm));
   trapElement->setAttributeNode(NewAttribute("x4",
                2.0*trap->GetXHalfLength4()/mm));
   trapElement->setAttributeNode(NewAttribute("alpha2",alpha2/degree));
   trapElement->setAttributeNode(NewAttribute("aunit","deg"));
   trapElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(trapElement);
}

void G4GDMLWriteSolids::
TrdWrite(xercesc::DOMElement* solElement, const G4Trd* const trd)
{
   const G4String& name = GenerateName(trd->GetName(),trd);

   xercesc::DOMElement* trdElement = NewElement("trd");
   trdElement->setAttributeNode(NewAttribute("name",name));
   trdElement->setAttributeNode(NewAttribute("x1",
               2.0*trd->GetXHalfLength1()/mm));
   trdElement->setAttributeNode(NewAttribute("x2",
               2.0*trd->GetXHalfLength2()/mm));
   trdElement->setAttributeNode(NewAttribute("y1",
               2.0*trd->GetYHalfLength1()/mm));
   trdElement->setAttributeNode(NewAttribute("y2",
               2.0*trd->GetYHalfLength2()/mm));
   trdElement->setAttributeNode(NewAttribute("z",
               2.0*trd->GetZHalfLength()/mm));
   trdElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(trdElement);
}

void G4GDMLWriteSolids::
TubeWrite(xercesc::DOMElement* solElement, const G4Tubs* const tube)
{
   const G4String& name = GenerateName(tube->GetName(),tube);

   xercesc::DOMElement* tubeElement = NewElement("tube");
   tubeElement->setAttributeNode(NewAttribute("name",name));
   tubeElement->setAttributeNode(NewAttribute("rmin",
                tube->GetInnerRadius()/mm));
   tubeElement->setAttributeNode(NewAttribute("rmax",
                tube->GetOuterRadius()/mm));
   tubeElement->setAttributeNode(NewAttribute("z",
                2.0*tube->GetZHalfLength()/mm));
   tubeElement->setAttributeNode(NewAttribute("startphi",
                tube->GetStartPhiAngle()/degree));
   tubeElement->setAttributeNode(NewAttribute("deltaphi",
                tube->GetDeltaPhiAngle()/degree));
   tubeElement->setAttributeNode(NewAttribute("aunit","deg"));
   tubeElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(tubeElement);
}

void G4GDMLWriteSolids::
CutTubeWrite(xercesc::DOMElement* solElement, const G4CutTubs* const cuttube)
{
   const G4String& name = GenerateName(cuttube->GetName(),cuttube);

   xercesc::DOMElement* cuttubeElement = NewElement("cutTube");
   cuttubeElement->setAttributeNode(NewAttribute("name",name));
   cuttubeElement->setAttributeNode(NewAttribute("rmin",
                cuttube->GetInnerRadius()/mm));
   cuttubeElement->setAttributeNode(NewAttribute("rmax",
                cuttube->GetOuterRadius()/mm));
   cuttubeElement->setAttributeNode(NewAttribute("z",
                2.0*cuttube->GetZHalfLength()/mm));
   cuttubeElement->setAttributeNode(NewAttribute("startphi",
                cuttube->GetStartPhiAngle()/degree));
   cuttubeElement->setAttributeNode(NewAttribute("deltaphi",
                cuttube->GetDeltaPhiAngle()/degree));
   cuttubeElement->setAttributeNode(NewAttribute("lowX",
                cuttube->GetLowNorm().getX()/mm));
   cuttubeElement->setAttributeNode(NewAttribute("lowY",
                cuttube->GetLowNorm().getY()/mm));
   cuttubeElement->setAttributeNode(NewAttribute("lowZ",
                cuttube->GetLowNorm().getZ()/mm));
   cuttubeElement->setAttributeNode(NewAttribute("highX",
                cuttube->GetHighNorm().getX()/mm));
   cuttubeElement->setAttributeNode(NewAttribute("highY",
                cuttube->GetHighNorm().getY()/mm));
   cuttubeElement->setAttributeNode(NewAttribute("highZ",
                cuttube->GetHighNorm().getZ()/mm));
   cuttubeElement->setAttributeNode(NewAttribute("aunit","deg"));
   cuttubeElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(cuttubeElement);
}

void G4GDMLWriteSolids::
TwistedboxWrite(xercesc::DOMElement* solElement,
                const G4TwistedBox* const twistedbox)
{
   const G4String& name = GenerateName(twistedbox->GetName(),twistedbox);

   xercesc::DOMElement* twistedboxElement = NewElement("twistedbox");
   twistedboxElement->setAttributeNode(NewAttribute("name",name));
   twistedboxElement->setAttributeNode(NewAttribute("x",
                      2.0*twistedbox->GetXHalfLength()/mm));
   twistedboxElement->setAttributeNode(NewAttribute("y",
                      2.0*twistedbox->GetYHalfLength()/mm));
   twistedboxElement->setAttributeNode(NewAttribute("z",
                      2.0*twistedbox->GetZHalfLength()/mm));
   twistedboxElement->setAttributeNode(NewAttribute("PhiTwist",
                      twistedbox->GetPhiTwist()/degree));
   twistedboxElement->setAttributeNode(NewAttribute("aunit","deg"));
   twistedboxElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(twistedboxElement);
}

void G4GDMLWriteSolids::
TwistedtrapWrite(xercesc::DOMElement* solElement,
                 const G4TwistedTrap* const twistedtrap)
{
   const G4String& name = GenerateName(twistedtrap->GetName(),twistedtrap);

   xercesc::DOMElement* twistedtrapElement = NewElement("twistedtrap");
   twistedtrapElement->setAttributeNode(NewAttribute("name",name));
   twistedtrapElement->setAttributeNode(NewAttribute("y1",
                       2.0*twistedtrap->GetY1HalfLength()/mm));
   twistedtrapElement->setAttributeNode(NewAttribute("x1",
                       2.0*twistedtrap->GetX1HalfLength()/mm));
   twistedtrapElement->setAttributeNode(NewAttribute("x2",
                       2.0*twistedtrap->GetX2HalfLength()/mm));
   twistedtrapElement->setAttributeNode(NewAttribute("y2",
                       2.0*twistedtrap->GetY2HalfLength()/mm));
   twistedtrapElement->setAttributeNode(NewAttribute("x3",
                       2.0*twistedtrap->GetX3HalfLength()/mm));
   twistedtrapElement->setAttributeNode(NewAttribute("x4",
                       2.0*twistedtrap->GetX4HalfLength()/mm));
   twistedtrapElement->setAttributeNode(NewAttribute("z",
                       2.0*twistedtrap->GetZHalfLength()/mm));
   twistedtrapElement->setAttributeNode(NewAttribute("Alph",
                       twistedtrap->GetTiltAngleAlpha()/degree));
   twistedtrapElement->setAttributeNode(NewAttribute("Theta",
                       twistedtrap->GetPolarAngleTheta()/degree));
   twistedtrapElement->setAttributeNode(NewAttribute("Phi",
                       twistedtrap->GetAzimuthalAnglePhi()/degree));
   twistedtrapElement->setAttributeNode(NewAttribute("PhiTwist",
                       twistedtrap->GetPhiTwist()/degree));
   twistedtrapElement->setAttributeNode(NewAttribute("aunit","deg"));
   twistedtrapElement->setAttributeNode(NewAttribute("lunit","mm"));
   
   solElement->appendChild(twistedtrapElement);
}

void G4GDMLWriteSolids::
TwistedtrdWrite(xercesc::DOMElement* solElement,
                const G4TwistedTrd* const twistedtrd)
{
   const G4String& name = GenerateName(twistedtrd->GetName(),twistedtrd);

   xercesc::DOMElement* twistedtrdElement = NewElement("twistedtrd");
   twistedtrdElement->setAttributeNode(NewAttribute("name",name));
   twistedtrdElement->setAttributeNode(NewAttribute("x1",
                      2.0*twistedtrd->GetX1HalfLength()/mm));
   twistedtrdElement->setAttributeNode(NewAttribute("x2",
                      2.0*twistedtrd->GetX2HalfLength()/mm));
   twistedtrdElement->setAttributeNode(NewAttribute("y1",
                      2.0*twistedtrd->GetY1HalfLength()/mm));
   twistedtrdElement->setAttributeNode(NewAttribute("y2",
                      2.0*twistedtrd->GetY2HalfLength()/mm));
   twistedtrdElement->setAttributeNode(NewAttribute("z",
                      2.0*twistedtrd->GetZHalfLength()/mm));
   twistedtrdElement->setAttributeNode(NewAttribute("PhiTwist",
                      twistedtrd->GetPhiTwist()/degree));
   twistedtrdElement->setAttributeNode(NewAttribute("aunit","deg"));
   twistedtrdElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(twistedtrdElement);
}

void G4GDMLWriteSolids::
TwistedtubsWrite(xercesc::DOMElement* solElement,
                 const G4TwistedTubs* const twistedtubs)
{
   const G4String& name = GenerateName(twistedtubs->GetName(),twistedtubs);

   xercesc::DOMElement* twistedtubsElement = NewElement("twistedtubs");
   twistedtubsElement->setAttributeNode(NewAttribute("name",name));
   twistedtubsElement->setAttributeNode(NewAttribute("twistedangle",
                       twistedtubs->GetPhiTwist()/degree));
   twistedtubsElement->setAttributeNode(NewAttribute("endinnerrad",
                       twistedtubs->GetInnerRadius()/mm));
   twistedtubsElement->setAttributeNode(NewAttribute("endouterrad",
                       twistedtubs->GetOuterRadius()/mm));
   twistedtubsElement->setAttributeNode(NewAttribute("zlen",
                       2.0*twistedtubs->GetZHalfLength()/mm));
   twistedtubsElement->setAttributeNode(NewAttribute("phi",
                       twistedtubs->GetDPhi()/degree));
   twistedtubsElement->setAttributeNode(NewAttribute("aunit","deg"));
   twistedtubsElement->setAttributeNode(NewAttribute("lunit","mm"));
   solElement->appendChild(twistedtubsElement);
}

void G4GDMLWriteSolids::
ZplaneWrite(xercesc::DOMElement* element, const G4double& z,
            const G4double& rmin, const G4double& rmax)
{
   xercesc::DOMElement* zplaneElement = NewElement("zplane");
   zplaneElement->setAttributeNode(NewAttribute("z",z/mm));
   zplaneElement->setAttributeNode(NewAttribute("rmin",rmin/mm));
   zplaneElement->setAttributeNode(NewAttribute("rmax",rmax/mm));
   element->appendChild(zplaneElement);
}

void G4GDMLWriteSolids::
RZPointWrite(xercesc::DOMElement* element, const G4double& r,
            const G4double& z)
{
   xercesc::DOMElement* rzpointElement = NewElement("rzpoint");
   rzpointElement->setAttributeNode(NewAttribute("r",r/mm));
   rzpointElement->setAttributeNode(NewAttribute("z",z/mm));
   element->appendChild(rzpointElement);
}

void G4GDMLWriteSolids::
OpticalSurfaceWrite(xercesc::DOMElement* solElement,
                    const G4OpticalSurface* const surf)
{
   xercesc::DOMElement* optElement = NewElement("opticalsurface");
   G4OpticalSurfaceModel smodel = surf->GetModel();
   G4double sval = (smodel==glisur) ? surf->GetPolish() : surf->GetSigmaAlpha();

   optElement->setAttributeNode(NewAttribute("name", surf->GetName()));
   optElement->setAttributeNode(NewAttribute("model", smodel));
   optElement->setAttributeNode(NewAttribute("finish", surf->GetFinish()));
   optElement->setAttributeNode(NewAttribute("type", surf->GetType()));
   optElement->setAttributeNode(NewAttribute("value", sval));

   // Write any property attached to the optical surface...
   //
   if (surf->GetMaterialPropertiesTable())
   {
     PropertyWrite(optElement, surf);
   }

   solElement->appendChild(optElement);
}

void G4GDMLWriteSolids::PropertyWrite(xercesc::DOMElement* optElement,
                                         const G4OpticalSurface* const surf)
{
   xercesc::DOMElement* propElement;
   G4MaterialPropertiesTable* ptable = surf->GetMaterialPropertiesTable();
   const std::map< G4String, G4PhysicsOrderedFreeVector*,
                 std::less<G4String> >* pmap = ptable->GetPropertiesMap();
   const std::map< G4String, G4double,
                 std::less<G4String> >* cmap = ptable->GetPropertiesCMap();
   std::map< G4String, G4PhysicsOrderedFreeVector*,
                 std::less<G4String> >::const_iterator mpos;
   std::map< G4String, G4double,
                 std::less<G4String> >::const_iterator cpos;
   for (mpos=pmap->begin(); mpos!=pmap->end(); mpos++)
   {
      propElement = NewElement("property");
      propElement->setAttributeNode(NewAttribute("name", mpos->first));
      propElement->setAttributeNode(NewAttribute("ref",
                                    GenerateName(mpos->first, mpos->second)));
      if (mpos->second)
      {
         PropertyVectorWrite(mpos->first, mpos->second);
         optElement->appendChild(propElement);
      }
      else
      {
         G4String warn_message = "Null pointer for material property -"
                  + mpos->first + "- of optical surface -" + surf->GetName() + "- !";
         G4Exception("G4GDMLWriteSolids::PropertyWrite()", "NullPointer",
                     JustWarning, warn_message);
         continue;
      }
   }
   for (cpos=cmap->begin(); cpos!=cmap->end(); cpos++)
   {
      propElement = NewElement("property");
      propElement->setAttributeNode(NewAttribute("name", cpos->first));
      propElement->setAttributeNode(NewAttribute("ref", cpos->first));
      xercesc::DOMElement* constElement = NewElement("constant");
      constElement->setAttributeNode(NewAttribute("name", cpos->first));
      constElement->setAttributeNode(NewAttribute("value", cpos->second));
      defineElement->appendChild(constElement);
      optElement->appendChild(propElement);
   }
}

void G4GDMLWriteSolids::SolidsWrite(xercesc::DOMElement* gdmlElement)
{
   G4cout << "G4GDML: Writing solids..." << G4endl;

   solidsElement = NewElement("solids");
   gdmlElement->appendChild(solidsElement);

   solidList.clear();
}

void G4GDMLWriteSolids::AddSolid(const G4VSolid* const solidPtr)
{
   for (size_t i=0; i<solidList.size(); i++)   // Check if solid is
   {                                           // already in the list!
      if (solidList[i] == solidPtr)  { return; }
   }

   solidList.push_back(solidPtr);

   if (const G4BooleanSolid* const booleanPtr
     = dynamic_cast<const G4BooleanSolid*>(solidPtr))
     { BooleanWrite(solidsElement,booleanPtr); } else
   if (const G4ScaledSolid* const scaledPtr
     = dynamic_cast<const G4ScaledSolid*>(solidPtr))
     { ScaledWrite(solidsElement,scaledPtr); } else
   if (solidPtr->GetEntityType()=="G4MultiUnion")
     { const G4MultiUnion* const munionPtr
     = static_cast<const G4MultiUnion*>(solidPtr);
       MultiUnionWrite(solidsElement,munionPtr); } else
   if (solidPtr->GetEntityType()=="G4Box")
     { const G4Box* const boxPtr
     = static_cast<const G4Box*>(solidPtr);
       BoxWrite(solidsElement,boxPtr); } else
   if (solidPtr->GetEntityType()=="G4Cons")
     { const G4Cons* const conePtr
     = static_cast<const G4Cons*>(solidPtr);
       ConeWrite(solidsElement,conePtr); } else
   if (solidPtr->GetEntityType()=="G4EllipticalCone")
     { const G4EllipticalCone* const elconePtr
     = static_cast<const G4EllipticalCone*>(solidPtr);
       ElconeWrite(solidsElement,elconePtr); } else
   if (solidPtr->GetEntityType()=="G4Ellipsoid")
     { const G4Ellipsoid* const ellipsoidPtr
     = static_cast<const G4Ellipsoid*>(solidPtr);
       EllipsoidWrite(solidsElement,ellipsoidPtr); } else
   if (solidPtr->GetEntityType()=="G4EllipticalTube")
     { const G4EllipticalTube* const eltubePtr
     = static_cast<const G4EllipticalTube*>(solidPtr);
       EltubeWrite(solidsElement,eltubePtr); } else
   if (solidPtr->GetEntityType()=="G4ExtrudedSolid")
     { const G4ExtrudedSolid* const xtruPtr
     = static_cast<const G4ExtrudedSolid*>(solidPtr);
       XtruWrite(solidsElement,xtruPtr); } else
   if (solidPtr->GetEntityType()=="G4Hype")
     { const G4Hype* const hypePtr
     = static_cast<const G4Hype*>(solidPtr);
       HypeWrite(solidsElement,hypePtr); } else
   if (solidPtr->GetEntityType()=="G4Orb")
     { const G4Orb* const orbPtr
     = static_cast<const G4Orb*>(solidPtr);
       OrbWrite(solidsElement,orbPtr); } else
   if (solidPtr->GetEntityType()=="G4Para")
     { const G4Para* const paraPtr
     = static_cast<const G4Para*>(solidPtr);
       ParaWrite(solidsElement,paraPtr); } else
   if (solidPtr->GetEntityType()=="G4Paraboloid")
     { const G4Paraboloid* const paraboloidPtr
     = static_cast<const G4Paraboloid*>(solidPtr);
       ParaboloidWrite(solidsElement,paraboloidPtr); } else
   if (solidPtr->GetEntityType()=="G4Polycone")
     { const G4Polycone* const polyconePtr
     = static_cast<const G4Polycone*>(solidPtr);
       PolyconeWrite(solidsElement,polyconePtr); } else
   if (solidPtr->GetEntityType()=="G4GenericPolycone")
     { const G4GenericPolycone* const genpolyconePtr
     = static_cast<const G4GenericPolycone*>(solidPtr);
       GenericPolyconeWrite(solidsElement,genpolyconePtr); } else
   if (solidPtr->GetEntityType()=="G4Polyhedra")
     { const G4Polyhedra* const polyhedraPtr
     = static_cast<const G4Polyhedra*>(solidPtr);
       PolyhedraWrite(solidsElement,polyhedraPtr); } else
   if (solidPtr->GetEntityType()=="G4Sphere")
     { const G4Sphere* const spherePtr
     = static_cast<const G4Sphere*>(solidPtr);
       SphereWrite(solidsElement,spherePtr); } else
   if (solidPtr->GetEntityType()=="G4TessellatedSolid")
     { const G4TessellatedSolid* const tessellatedPtr
     = static_cast<const G4TessellatedSolid*>(solidPtr);
       TessellatedWrite(solidsElement,tessellatedPtr); } else
   if (solidPtr->GetEntityType()=="G4Tet")
     { const G4Tet* const tetPtr
     = static_cast<const G4Tet*>(solidPtr);
       TetWrite(solidsElement,tetPtr); } else
   if (solidPtr->GetEntityType()=="G4Torus")
     { const G4Torus* const torusPtr
     = static_cast<const G4Torus*>(solidPtr);
       TorusWrite(solidsElement,torusPtr); } else
   if (solidPtr->GetEntityType()=="G4GenericTrap")
     { const G4GenericTrap* const gtrapPtr
     = static_cast<const G4GenericTrap*>(solidPtr);
       GenTrapWrite(solidsElement,gtrapPtr); } else
   if (solidPtr->GetEntityType()=="G4Trap")
     { const G4Trap* const trapPtr
     = static_cast<const G4Trap*>(solidPtr);
       TrapWrite(solidsElement,trapPtr); } else
   if (solidPtr->GetEntityType()=="G4Trd")
     { const G4Trd* const trdPtr
     = static_cast<const G4Trd*>(solidPtr);
       TrdWrite(solidsElement,trdPtr); } else
   if (solidPtr->GetEntityType()=="G4Tubs")
     { const G4Tubs* const tubePtr
     = static_cast<const G4Tubs*>(solidPtr);
       TubeWrite(solidsElement,tubePtr); } else
   if (solidPtr->GetEntityType()=="G4CutTubs")
     { const G4CutTubs* const cuttubePtr
     = static_cast<const G4CutTubs*>(solidPtr);
       CutTubeWrite(solidsElement,cuttubePtr); } else
   if (solidPtr->GetEntityType()=="G4TwistedBox")
     { const G4TwistedBox* const twistedboxPtr
     = static_cast<const G4TwistedBox*>(solidPtr);
       TwistedboxWrite(solidsElement,twistedboxPtr); } else
   if (solidPtr->GetEntityType()=="G4TwistedTrap")
     { const G4TwistedTrap* const twistedtrapPtr
     = static_cast<const G4TwistedTrap*>(solidPtr);
       TwistedtrapWrite(solidsElement,twistedtrapPtr); } else
   if (solidPtr->GetEntityType()=="G4TwistedTrd")
     { const G4TwistedTrd* const twistedtrdPtr
     = static_cast<const G4TwistedTrd*>(solidPtr);
       TwistedtrdWrite(solidsElement,twistedtrdPtr); } else
   if (solidPtr->GetEntityType()=="G4TwistedTubs")
     { const G4TwistedTubs* const twistedtubsPtr
     = static_cast<const G4TwistedTubs*>(solidPtr);
       TwistedtubsWrite(solidsElement,twistedtubsPtr); }
   else
   {
     G4String error_msg = "Unknown solid: " + solidPtr->GetName()
                        + "; Type: " + solidPtr->GetEntityType();
     G4Exception("G4GDMLWriteSolids::AddSolid()", "WriteError",
                 FatalException, error_msg);
   }
}
