#include "G4tgbElement.hh"
#include "G4tgbMaterialMgr.hh"
#include "G4tgrElementSimple.hh"
#include "G4tgrElementFromIsotopes.hh"
#include "globals.hh"
#include "G4tgrMessenger.hh"

//----------------------------------------------------------------------
G4tgbElement::G4tgbElement( G4tgrElement* hg )
{
  theTgrElem = hg;
  theG4Elem = 0;
}


//----------------------------------------------------------------------
G4Element* G4tgbElement::BuildG4ElementSimple()
{
  G4Element* elem = 0;

  //-------- if G4Element not found, construct it 
  if( theG4Elem == 0 ) { 
  //----- construct new G4Element 
    G4tgrElementSimple* tgrElem = static_cast<G4tgrElementSimple*>(theTgrElem);

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " constructing  new G4Element " <<  " " << tgrElem->GetName() << " " << tgrElem->GetSymbol() << " " << tgrElem->GetZ() << " " << tgrElem->GetA() << endl;
#endif
    elem = new G4Element(tgrElem->GetName(), tgrElem->GetSymbol(), tgrElem->GetZ(), tgrElem->GetA() );
    theG4Elem = elem;
  } else {
    elem = theG4Elem; 
  }

  return elem;
}


//----------------------------------------------------------------------
G4Element* G4tgbElement::BuildG4ElementFromIsotopes()
{
  G4Element* elem = 0;

  //-------- if G4Element not found, construct it 
  if( theG4Elem == 0 ) { 
    //----- construct new G4Element 
    G4tgrElementFromIsotopes* tgrElem = static_cast<G4tgrElementFromIsotopes*>(theTgrElem);
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " constructing  new G4Element " <<  " " << tgrElem->GetName() << " " << tgrElem->GetSymbol() << " " << tgrElem->GetNumberOfIsotopes() << endl;
#endif
    elem = new G4Element(tgrElem->GetName(), tgrElem->GetSymbol(), tgrElem->GetNumberOfIsotopes() );

    //----- add isotopes
    G4Isotope* compIsot;
    G4tgbMaterialMgr* mf = G4tgbMaterialMgr::GetInstance();
    for( int ii = 0; ii < tgrElem->GetNumberOfIsotopes(); ii++) {
      //Look if this component is a material
      compIsot = mf->FindOrBuildG4Isotope( tgrElem->GetComponent(ii) );
    //-      cout << " compMAte " <<  theFractionsByWeight[ii]<<endl;
      if( compIsot != 0 ) { 
	elem->AddIsotope( compIsot, tgrElem->GetAbundance(ii) );
	//-	cout << "isotrial added" << endl;
      } else {
	G4Exception("!!! EXITING G4tgbElement::BuildG4ElementFromIsotopes(). componenent " + tgrElem->GetComponent(ii) + " of element " + tgrElem->GetName() + " is not an isotope ");
      }
    } 
    
    theG4Elem = elem;
  } else {
    elem = theG4Elem; 
  }

  return elem;
}
