// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AdvancedBrepShapeRepresentationCreator.cc,v 1.4 2000-11-20 18:17:26 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4AdvancedBrepShapeRepresentationCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include <stdio.h>
#include "G4AdvancedBrepShapeRepresentationCreator.hh"
#include "G4GeometryTable.hh"

G4AdvancedBrepShapeRepresentationCreator
   G4AdvancedBrepShapeRepresentationCreator::csc;


G4AdvancedBrepShapeRepresentationCreator::
   G4AdvancedBrepShapeRepresentationCreator()
{
  G4GeometryTable::RegisterObject(this);
}


G4AdvancedBrepShapeRepresentationCreator::
   ~G4AdvancedBrepShapeRepresentationCreator(){}


void G4AdvancedBrepShapeRepresentationCreator::CreateG4Geometry(STEPentity& sEnt)
{
  G4String attrName("items");
  STEPattribute *Attr = GetNamedAttribute(attrName, sEnt);
  EntityAggregate *ea = (EntityAggregate*)Attr->ptr.a;
  
  // Get parts of this shape  
  G4int partCount = ea->EntryCount();
  EntityNode* en = (EntityNode*)ea->GetHead();

  // An Advaced_BREP_Shape_Representation is made of
  // 0 or 1 G4Axis2Placement3D and 
  // 1 or more Manifold_Solid_BREP 

  G4int placement = 0, nbofsolids = 0;
  G4PlacedSolidVector* placedSldV = new G4PlacedSolidVector(); 
  G4PlacedSolid* placedSld; 
  G4BREPSolid* sld; 
  STEPentity* ent;
  G4Axis2Placement3D* place = 0;

  // default place
  G4Axis2Placement3D* place0 = new (G4Axis2Placement3D)( G4Vector3D(1, 0, 0),
							 G4Vector3D(0, 1, 0),
							 G4Point3D(0, 0, 0)  );
  

  for(G4int a = 0; a < partCount; a++) 
  {
   ent = en->node; 
   void *tmp = G4GeometryTable::CreateObject(*ent); 

   if (tmp)
   {
     if( !strcmp(ent->EntityName(), "Axis2_Placement_3d") )
     {
       // Get placement 
       place = (G4Axis2Placement3D*)tmp;  
       placement = 1;
     }
     else
       if( !strcmp(ent->EntityName(), "Manifold_Solid_Brep") )
       {
         // Get solid    
         sld = (G4BREPSolid*)tmp; 
         nbofsolids++;

         G4int id = ent->GetFileId();
         sld->SetId(id);

         // Create the placed solid 
         if(placement)
	   placedSld = new G4PlacedSolid(sld, place);
         else
	   placedSld = new G4PlacedSolid(sld, place0);
       
         placedSldV->append(placedSld);
       }
   }
   else
     G4cerr << "WARNING - G4AdvancedBrepShapeRepresentationCreator::CreateG4Geometry" << G4endl
            << "\tEntity part for " << ent->EntityName() << " NOT found." << G4endl;

   en = (EntityNode*)en->NextNode();
  }

  // Get geometric_representation_context, but it is not used
  attrName = "context_of_items";
  Attr = GetNamedAttribute(attrName, sEnt);
  
  createdObject = placedSldV;
}


void G4AdvancedBrepShapeRepresentationCreator::CreateSTEPGeometry(void* G4obj)
{
}
