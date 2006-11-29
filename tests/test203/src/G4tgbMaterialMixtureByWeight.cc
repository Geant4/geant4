#include "G4tgbMaterialMixtureByWeight.hh"
#include "G4tgbMaterialMixture.hh"
#include "G4tgbMaterial.hh"
#include "G4tgbMaterialMgr.hh"
#include "G4tgrMessenger.hh"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4tgbMaterialMixtureByWeight::G4tgbMaterialMixtureByWeight( G4tgrMaterial* hg)
{
  theTgrMate = hg;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4Material* G4tgbMaterialMixtureByWeight::BuildG4Material()
{
  
  //----- construct new G4Material with component materials (a Mixture)
  G4Material* mate = new G4Material( theTgrMate->GetName(), theTgrMate->GetDensity(), theTgrMate->GetNumberOfComponents() );
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbMaterialMixtureByWeight::BuildG4Material: constructing  new G4Material " <<  " " << theTgrMate->GetName() << " " << " " << theTgrMate->GetDensity() << endl;
#endif
  
  //--- add components
  G4Element* compElem;
  G4Material* compMate;
  G4tgbMaterialMgr* mf = G4tgbMaterialMgr::GetInstance();
  G4bool bCompElem = FALSE;
  for( int ii = 0; ii < theTgrMate->GetNumberOfComponents(); ii++) {
    //look if this component is an element
    G4tgbElement* hselem = mf->FindG4tgbElement( GetComponent(ii) );
    if( hselem != 0 ) {
      compElem = mf->FindOrBuildG4Element( GetComponent(ii) );
      //-      cout << " adding g4 element " << endl;
      mate->AddElement( compElem, GetFraction( ii ) );
    //if it is not an element, it must be a material
    } else { 
      compMate = mf->FindOrBuildG4Material( GetComponent(ii) );
      //-      cout << " compMAte " << fraction(ii)<<endl;
      if( compMate != 0 ) { 
	// if it is a material add it by weight fraction 
	mate->AddMaterial( compMate, GetFraction( ii ) );
	//-	cout << "material added" << endl;
      } else {
	G4Exception("!!! EXITING G4tgbMaterialMixtureByWeight::BuildG4Material. componenent " + GetComponent(ii) + " of material " + theTgrMate->GetName() + " is not an element nor a material ");
      }
    } 
  }

  return mate;

}

