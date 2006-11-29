#include "G4tgbMaterialMixtureByVolume.hh"
#include "G4tgbMaterialMixture.hh"
#include "G4tgbMaterial.hh"
#include "G4tgbMaterialMgr.hh"
#include "G4tgrMessenger.hh"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4tgbMaterialMixtureByVolume::G4tgbMaterialMixtureByVolume( G4tgrMaterial* hg)
{
  theTgrMate = hg;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4Material* G4tgbMaterialMixtureByVolume::BuildG4Material()
{
 
  //----- construct new G4Material with components materials (a Mixture)
  G4Material* mate = new G4Material( theTgrMate->GetName(), theTgrMate->GetDensity(), theTgrMate->GetNumberOfComponents() );
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbMaterialMixtureByVolume::buildG4Material: constructing new G4Material " <<  " " << theTgrMate->GetName() << endl;
#endif
 
  //----- Transform fractions by volume to fractions by weight
  TransformToFractionsByWeight();

  //----- add components
  G4Element* compElem;
  G4Material* compMate;
  G4tgbMaterialMgr* mf = G4tgbMaterialMgr::GetInstance();
  for( uint ii = 0; ii < theTgrMate->GetNumberOfComponents(); ii++) {
    //Look if this component is a material
    compMate = mf->FindOrBuildG4Material( GetComponent(ii) );
    //-      cout << " compMAte " <<  theFractionsByWeight[ii]<<endl;
    if( compMate != 0 ) { 
      // if it is a material add it by weight fraction 
      mate->AddMaterial( compMate,  theFractionsByWeight[ii] );
      //-	cout << "material added" << endl;
    } else {
      G4Exception("G4tgbMaterialMixtureByVolume::buildG4Material. componenent " + GetComponent(ii) + " of material " + theTgrMate->GetName() + " is not an element nor a material ");
    }
  } 

  return mate;

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void G4tgbMaterialMixtureByVolume::TransformToFractionsByWeight() 
{ 
  theFractionsByWeight = new double( theTgrMate->GetNumberOfComponents() );

  G4tgbMaterialMgr* mf = G4tgbMaterialMgr::GetInstance();
  G4Material* compMate;
  double totalfd = 0.;
  for( uint ii = 0; ii < theTgrMate->GetNumberOfComponents(); ii++) {
    compMate = mf->FindOrBuildG4Material( GetComponent(ii) );
    //-      cout << " compMAte " <<  theFractionsByWeight[ii]<<endl;
    if( compMate != 0 ) { 
      // if it is a material add it by weight fraction 
      theFractionsByWeight[ii] = GetFraction(ii) * compMate->GetDensity();
      totalfd += theFractionsByWeight[ii];
      //      cout << compMate->GetName() << " adding mate " << fraction(ii) << " " << compMate->GetDensity() << " " << totalfd << endl;
      //-	cout << "material added" << endl;
    } else {
      G4Exception("G4tgbMaterialMixtureByVolume::BuildG4Material. componenent " + GetComponent(ii) + " of material " + theTgrMate->GetName() + " is not a material ");
    }   
  }
  for( uint ii = 0; ii < theTgrMate->GetNumberOfComponents(); ii++) {
    theFractionsByWeight[ii] /= totalfd;
    //    cout << ii << "  theFractiosByWeight " <<  theFractionsByWeight[ii] << endl;
  }
}

