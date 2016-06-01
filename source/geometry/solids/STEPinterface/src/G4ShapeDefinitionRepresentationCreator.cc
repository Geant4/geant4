#include "G4ShapeDefinitionRepresentationCreator.hh"

G4ShapeDefinitionRepresentationCreator G4ShapeDefinitionRepresentationCreator::csc;

G4ShapeDefinitionRepresentationCreator::G4ShapeDefinitionRepresentationCreator(){G4GeometryTable::RegisterObject(this);}

G4ShapeDefinitionRepresentationCreator::~G4ShapeDefinitionRepresentationCreator(){}

void G4ShapeDefinitionRepresentationCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4String attrName("used_representation");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  G4PlacedSolid* ps = (G4PlacedSolid*)tmp;
  createdObject = ps;

/*
  // Output messages
  G4cout<<"Shape definition representation created with "
	<<"the Advanced Brep shape :"<<endl;

  G4cout<<" - BBox of the solid :"<<endl;
  G4cout<< "    box min: "
	<< ((G4BREPSolid*)ps->GetSolid())->GetBBox()->GetBoxMin().x()<< " "
	<< ((G4BREPSolid*)ps->GetSolid())->GetBBox()->GetBoxMin().y()<< " "
	<< ((G4BREPSolid*)ps->GetSolid())->GetBBox()->GetBoxMin().z()<< endl;

  G4cout<< "    box max: "
	<< ((G4BREPSolid*)ps->GetSolid())->GetBBox()->GetBoxMax().x()<< " "
	<< ((G4BREPSolid*)ps->GetSolid())->GetBBox()->GetBoxMax().y()<< " "
	<< ((G4BREPSolid*)ps->GetSolid())->GetBBox()->GetBoxMax().z()<< endl;

*/
}

void G4ShapeDefinitionRepresentationCreator::CreateSTEPGeometry(void* G4obj)
{

}
