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
// $Id: G4tgbElement.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgbElement

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgbElement.hh"
#include "G4tgbMaterialMgr.hh"
#include "G4tgrElementSimple.hh"
#include "G4tgrElementFromIsotopes.hh"
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
  if( theG4Elem == 0 )
  { 
    //----- construct new G4Element 
    G4tgrElementSimple* tgrElem = static_cast<G4tgrElementSimple*>(theTgrElem);

    elem = new G4Element(tgrElem->GetName(), tgrElem->GetSymbol(),
                         tgrElem->GetZ(), tgrElem->GetA() );
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 1 )
    {
      G4cout << " Constructing new G4Element: " 
             << *elem << G4endl;
    }
#endif
    theG4Elem = elem;
  }
  else
  {
    elem = theG4Elem; 
  }

  return elem;
}


//----------------------------------------------------------------------
G4Element* G4tgbElement::BuildG4ElementFromIsotopes()
{
  G4Element* elem = 0;

  //-------- if G4Element not found, construct it 
  if( theG4Elem == 0 )
  { 
    //----- construct new G4Element 
    G4tgrElementFromIsotopes* tgrElem
      = static_cast<G4tgrElementFromIsotopes*>(theTgrElem);

    elem = new G4Element(tgrElem->GetName(), tgrElem->GetSymbol(),
                         tgrElem->GetNumberOfIsotopes() );

    //----- add isotopes
    G4Isotope* compIsot;
    G4tgbMaterialMgr* mf = G4tgbMaterialMgr::GetInstance();
    for( G4int ii = 0; ii < tgrElem->GetNumberOfIsotopes(); ii++)
    {
      // Look if this component is a material

      compIsot = mf->FindOrBuildG4Isotope( tgrElem->GetComponent(ii) );
      if( compIsot != 0 )
      { 
        elem->AddIsotope( compIsot, tgrElem->GetAbundance(ii) );
      }
      else
      {
        G4String ErrMessage = "Component " + tgrElem->GetComponent(ii)
                            + " of element " + tgrElem->GetName()
                            + " is not an isotope !";
        G4Exception("G4tgbElement::BuildG4ElementFromIsotopes()",
                    "InvalidSetup", FatalException, ErrMessage );
      }
    } 
    theG4Elem = elem;
  }
  else
  {
    elem = theG4Elem; 
  }


#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 1 )
    {
      G4cout << " Constructing  new G4Element from isotopes: "
             << *elem << G4endl;
    }
#endif

  return elem;
}
