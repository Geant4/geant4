// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RepresentationRelationshipCreator.cc,v 1.5 2000-11-20 18:17:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4RepresentationRelationshipCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include <STEPcomplex.h>

#include "G4RepresentationRelationshipCreator.hh"
#include "G4GeometryTable.hh"
#include "G4PlacementVector.hh"

G4RepresentationRelationshipCreator G4RepresentationRelationshipCreator::csc;

G4RepresentationRelationshipCreator::G4RepresentationRelationshipCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4RepresentationRelationshipCreator::~G4RepresentationRelationshipCreator() {}

void G4RepresentationRelationshipCreator::CreateG4Geometry(STEPentity& Ent)
{

  STEPcomplex* complexEnt = (STEPcomplex*)&Ent;
  STEPentity* subEnt=0;
  STEPentity* subEnt1=0;
  STEPentity* subEnt2=0;  
  G4String attrName;
  STEPattribute *Attr=0;

  G4PlacementVector*   itemDefPlaces=0;
  G4PlacementVector*   places=0;
  G4PlacedSolid*       sld=0;
  G4PlacedSolid*       tmpsld=0; 
  G4PlacedSolidVector* sldV = new G4PlacedSolidVector();
  G4PlacedSolidVector* tmpsldV=0;

  if(complexEnt->EntityExists("Representation_Relationship"))
  {
    subEnt = complexEnt->EntityPart("Representation_Relationship");

    // get subparts    
    attrName = "rep_1";
    Attr = GetNamedAttribute(attrName, *subEnt);
    subEnt1 = *Attr->ptr.c;
    places = (G4PlacementVector*)G4GeometryTable::CreateObject(*subEnt1);
#ifdef G4_STEPINTERFACE_DEBUG
    G4cout << "Rep_1 entries: " << places->entries()
           << " - Entity: " << subEnt1->EntityName() << G4endl;
#endif

    attrName = "rep_2";
    Attr = GetNamedAttribute(attrName, *subEnt);
    subEnt2 = *Attr->ptr.c;
    tmpsldV = (G4PlacedSolidVector*)G4GeometryTable::CreateObject(*subEnt2);
#ifdef G4_STEPINTERFACE_DEBUG
    G4cout << "Rep_2 entries: " << tmpsldV->entries()
           << " - Entity: " << subEnt2->EntityName() << G4endl;
#endif
  }
  G4String e1Name = subEnt1->EntityName();
  G4String e2Name = subEnt2->EntityName();
  if ((e1Name != "Shape_Representation") ||
      (e2Name != "Advanced_Brep_Shape_Representation") )
  {
    G4cerr << "WARNING - G4RepresentationRelationshipCreator::CreateG4Geometry" << G4endl
	   << "\tInconsistent association for representations:" << G4endl
	   << "\t" << e1Name << " - Rep_1 entries: " << places->entries() << G4endl
	   << "\t" << e2Name << " - Rep_2 entries: " << tmpsldV->entries() << G4endl
	   << "\tRepresentation Relationship NOT created." << G4endl;
    createdObject=0;
    return;
  }

  if( complexEnt->EntityExists
      ("Representation_Relationship_With_Transformation") )
  {
    subEnt = complexEnt->EntityPart
      ("Representation_Relationship_With_Transformation");
      
    // get subparts
    attrName = "transformation_operator";
    Attr = GetNamedAttribute(attrName, *subEnt);
    SdaiSelect* sel = Attr->ptr.sh;
    
    G4String underlyingTypeName(G4String(sel->UnderlyingTypeName()));
    G4String itemDefined("Item_Defined_Transformation");
    
    if(underlyingTypeName == itemDefined)
    {
      G4String koe;
      SdaiTransformation* transf = (SdaiTransformation*)sel;
      
      SdaiItem_defined_transformation* iTransf = 
	transf->operator SdaiItem_defined_transformationH();
	  
      SdaiRepresentation_itemH repItem1 = iTransf->transform_item_1_();
      SdaiRepresentation_itemH repItem2 = iTransf->transform_item_2_();
      
      subEnt = (STEPentity*)repItem1;
      itemDefPlaces = new G4PlacementVector();
      itemDefPlaces->append( ( (G4Axis2Placement3D*)G4GeometryTable::
			       CreateObject(*subEnt)                  ) );
      
      subEnt  = (STEPentity*)repItem2;
      itemDefPlaces->append( ( (G4Axis2Placement3D*)G4GeometryTable::
			       CreateObject(*subEnt)                  ) );
    }
    else
    {  
    }
      
    //      subEnt2 = *Attr->ptr.c;
    //      G4GeometryTable::CreateObject(*subEnt2);
  }

  // Place the advanced BREP, so place each manifold solids.
  // At this point, one cannot do anything else that guessing (!) the
  // logic behind the association of placements with each manifold solid.
  // See also creator for Advanced BREP shape Representations.
  // All this needs to be reviewed ! - GC
  //
  G4Axis2Placement3D* place = (G4Axis2Placement3D*)(places->at(0));
  G4int entr = tmpsldV->entries();
  for(G4int a = 0; a < entr; a++)
  {
    if (0 < a < places->entries())
      place = (G4Axis2Placement3D*)(places->at(a));
    tmpsld = tmpsldV->at(a);
    if (tmpsld)
    {
      sld = new G4PlacedSolid((G4BREPSolid*)(tmpsld->GetSolid()), place);
      sldV->append(sld);
    }
    else
    {
      G4cerr << "WARNING - G4RepresentationRelationshipCreator::CreateG4Geometry" << G4endl
             << "\tUnexpected NULL manifold solids vector !" << G4endl
	     << "\tPlacement failed." << G4endl;
    }
  }
  createdObject = sldV;
}

void G4RepresentationRelationshipCreator::CreateSTEPGeometry(void* G4obj)
{
}
