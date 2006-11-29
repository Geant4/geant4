#include "G4tgbMaterialMixtureByNoAtoms.hh"
#include "G4tgbMaterialMixture.hh"
#include "G4tgbMaterial.hh"
#include "G4tgbMaterialMgr.hh"
#include "G4tgrMessenger.hh"


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4tgbMaterialMixtureByNoAtoms::G4tgbMaterialMixtureByNoAtoms( G4tgrMaterial* hg)
{
  theTgrMate = hg;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4Material* G4tgbMaterialMixtureByNoAtoms::BuildG4Material()
{ 
  //----- construct new G4Material with components materials (a Mixture)
  G4Material* mate = new G4Material( theTgrMate->GetName(), theTgrMate->GetDensity(), theTgrMate->GetNumberOfComponents() );
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbMaterialMixtureByNoAtoms::BuildG4Material: constructing  new G4Material " <<  " " << theTgrMate->GetName() << " " << " " <<   theTgrMate->GetDensity() << endl;
#endif
  //--- add components
  G4Element* compElem;
  G4Material* compMate;
  compAreElements = 0;
  compAreMaterials = 0;
  G4tgbMaterialMgr* mf = G4tgbMaterialMgr::GetInstance();
  for( int ii = 0; ii < theTgrMate->GetNumberOfComponents(); ii++) {
    //look if this component is an element
    G4tgbElement* hselem = mf->FindG4tgbElement( GetComponent(ii) );
    if( hselem != 0 ) {
      compElem = mf->FindOrBuildG4Element( GetComponent(ii) );
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
	cout << " adding g4 element " << endl;
#endif
      if( compAreMaterials ) { 
	G4Exception("G4tgbMaterialMixtureByNoAtoms::BuildG4Material: material with some components elements and some materials "  + theTgrMate->GetName());
      }
      compAreElements = 1;
      //add it by number of atoms
      mate->AddElement( compElem, int(GetFraction(ii)) );
    //if it is not an element look if it is a material
    } else { 
      compMate = mf->FindOrBuildG4Material( GetComponent(ii) );
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
	cout << " compMAte " << GetFraction(ii)<<endl;
#endif
      if( compMate != 0 ) { 
	if( compAreElements ) { 
	  G4Exception("G4tgbMaterialMixtureByNoAtoms::BuildG4Material: material with some components materials and some elements " + theTgrMate->GetName());
	}
	compAreMaterials = 1;

	G4Exception("G4tgbMaterialMixtureByNoAtoms::buildG4Material: adding material by atoms is not supported: " + theTgrMate->GetName());

	// if it is a material add it by No Atoms
        double fr = GetFraction(ii);
        if( fr > 1.0 ) fr = 1.;
	mate->AddMaterial( compMate, GetFraction( ii ) );
	//-	cout << "material added" << endl;
      } else {
	G4Exception("G4tgbMaterialMixtureByWeight::buildG4Material. componenent " + GetComponent(ii) + " of material " +  theTgrMate->GetName() + " is not an element nor a material ");
      }
    } 
  }

  return mate;

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void G4tgbMaterialMixtureByNoAtoms::TransformToFractionsByWeight() 
{

}
