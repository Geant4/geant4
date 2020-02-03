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

#include "G4Types.hh"
#include "G4GDMLReadMaterials.hh"
#include "G4ExtrudedSolid.hh"
#include "G4MultiUnion.hh"
#include "G4MaterialPropertiesTable.hh"

class G4VSolid;
class G4QuadrangularFacet;
class G4TriangularFacet;
class G4SurfaceProperty;
class G4OpticalSurface;

class G4GDMLReadSolids : public G4GDMLReadMaterials
{
   enum BooleanOp {UNION,SUBTRACTION,INTERSECTION};
  // typedef struct { G4double rmin,rmax,z; } zplaneType;
   typedef struct { G4double r,z; } rzPointType;
  
 public:

   G4VSolid* GetSolid(const G4String&) const;
   G4SurfaceProperty* GetSurfaceProperty(const G4String&) const;

   virtual void SolidsRead(const xercesc::DOMElement* const);

 protected:
   typedef struct { G4double rmin,rmax,z; } zplaneType;
   G4GDMLReadSolids();
   virtual ~G4GDMLReadSolids();

   void BooleanRead(const xercesc::DOMElement* const,const BooleanOp);
   void BoxRead(const xercesc::DOMElement* const);
   void ConeRead(const xercesc::DOMElement* const);
   void ElconeRead(const xercesc::DOMElement* const);
   void EllipsoidRead(const xercesc::DOMElement* const);
   void EltubeRead(const xercesc::DOMElement* const);
   void XtruRead(const xercesc::DOMElement* const);
   void HypeRead(const xercesc::DOMElement* const);
   void MultiUnionNodeRead(const xercesc::DOMElement* const, G4MultiUnion* const);
   void MultiUnionRead(const xercesc::DOMElement* const);
   void OrbRead(const xercesc::DOMElement* const);
   void ParaRead(const xercesc::DOMElement* const);
   void ParaboloidRead(const xercesc::DOMElement* const);
   void PolyconeRead(const xercesc::DOMElement* const);
   void GenericPolyconeRead(const xercesc::DOMElement* const);
   void PolyhedraRead(const xercesc::DOMElement* const);
   void GenericPolyhedraRead(const xercesc::DOMElement* const);
   G4QuadrangularFacet* QuadrangularRead(const xercesc::DOMElement* const);
   void ReflectedSolidRead(const xercesc::DOMElement* const);
   void ScaledSolidRead(const xercesc::DOMElement* const);
   G4ExtrudedSolid::ZSection SectionRead(const xercesc::DOMElement* const,G4double);
   void SphereRead(const xercesc::DOMElement* const);
   void TessellatedRead(const xercesc::DOMElement* const);
   void TetRead(const xercesc::DOMElement* const);
   void TorusRead(const xercesc::DOMElement* const);
   void GenTrapRead(const xercesc::DOMElement* const);
   void TrapRead(const xercesc::DOMElement* const);
   void TrdRead(const xercesc::DOMElement* const);
   void TubeRead(const xercesc::DOMElement* const);
   void CutTubeRead(const xercesc::DOMElement* const);
   void TwistedboxRead(const xercesc::DOMElement* const);
   void TwistedtrapRead(const xercesc::DOMElement* const);
   void TwistedtrdRead(const xercesc::DOMElement* const);
   void TwistedtubsRead(const xercesc::DOMElement* const);
   G4TriangularFacet* TriangularRead(const xercesc::DOMElement* const);
   G4TwoVector TwoDimVertexRead(const xercesc::DOMElement* const,G4double);
   zplaneType ZplaneRead(const xercesc::DOMElement* const);
   rzPointType RZPointRead(const xercesc::DOMElement* const);
   void OpticalSurfaceRead(const xercesc::DOMElement* const);
   void PropertyRead(const xercesc::DOMElement* const,G4OpticalSurface*);

private:
  std::map<G4String, G4MaterialPropertyVector*> mapOfMatPropVects;
  
};

#endif
