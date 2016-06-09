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
// $Id: G4GDMLWriteSolids.hh,v 1.32.2.1 2009/08/11 08:27:49 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-03 $
//
//
// class G4GDMLWriteSolids
//
// Class description:
//
// GDML class for writing solids.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLWRITESOLIDS_INCLUDED_
#define _G4GDMLWRITESOLIDS_INCLUDED_

#include "G4BooleanSolid.hh"
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
#include "G4Polyhedra.hh"
#include "G4ReflectedSolid.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4TessellatedSolid.hh"
#include "G4Tet.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4TwistedBox.hh"
#include "G4TwistedTrap.hh"
#include "G4TwistedTrd.hh"
#include "G4TwistedTubs.hh"
#include "G4UnionSolid.hh"

#include "G4GDMLWriteMaterials.hh"

class G4GDMLWriteSolids : public G4GDMLWriteMaterials
{
  public:

   virtual void AddSolid(const G4VSolid* const);
   virtual void SolidsWrite(xercesc::DOMElement*);

  protected:

   G4GDMLWriteSolids();
   virtual ~G4GDMLWriteSolids();

   void BooleanWrite(xercesc::DOMElement*, const G4BooleanSolid* const);
   void BoxWrite(xercesc::DOMElement*, const G4Box* const);
   void ConeWrite(xercesc::DOMElement*, const G4Cons* const);
   void ElconeWrite(xercesc::DOMElement*, const G4EllipticalCone* const);
   void EllipsoidWrite(xercesc::DOMElement*, const G4Ellipsoid* const);
   void EltubeWrite(xercesc::DOMElement*, const G4EllipticalTube* const);
   void XtruWrite(xercesc::DOMElement*, const G4ExtrudedSolid* const);
   void HypeWrite(xercesc::DOMElement*, const G4Hype* const);
   void OrbWrite(xercesc::DOMElement*, const G4Orb* const);
   void ParaWrite(xercesc::DOMElement*, const G4Para* const);
   void ParaboloidWrite(xercesc::DOMElement*, const G4Paraboloid* const);
   void PolyconeWrite(xercesc::DOMElement*, const G4Polycone* const);
   void PolyhedraWrite(xercesc::DOMElement*, const G4Polyhedra* const);
   void SphereWrite(xercesc::DOMElement*, const G4Sphere* const);
   void TessellatedWrite(xercesc::DOMElement*, const G4TessellatedSolid* const);
   void TetWrite(xercesc::DOMElement*, const G4Tet* const);
   void TorusWrite(xercesc::DOMElement*, const G4Torus* const);
   void TrapWrite(xercesc::DOMElement*, const G4Trap* const);
   void TrdWrite(xercesc::DOMElement*, const G4Trd* const);
   void TubeWrite(xercesc::DOMElement*, const G4Tubs* const);
   void TwistedboxWrite(xercesc::DOMElement*, const G4TwistedBox* const);
   void TwistedtrapWrite(xercesc::DOMElement*, const G4TwistedTrap* const);
   void TwistedtrdWrite(xercesc::DOMElement*, const G4TwistedTrd* const);
   void TwistedtubsWrite(xercesc::DOMElement*, const G4TwistedTubs* const);
   void ZplaneWrite(xercesc::DOMElement*, const G4double&,
                    const G4double&, const G4double&);

  protected:

   std::vector<const G4VSolid*> solidList;
   xercesc::DOMElement* solidsElement;
   static const G4int maxTransforms = 8; // Constant for limiting the number
                                         // of displacements/reflections
                                         // applied to a single solid
};

#endif
