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
// $Id: Tst33MaterialFactory.hh,v 1.5 2006-06-29 21:59:39 gunter Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class Tst33MaterialFactory
//
// Class description:
//
// Creates the material used in this test.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33MaterialFactory_hh
#define Tst33MaterialFactory_hh Tst33MaterialFactory_hh

#include "globals.hh"
#include <map>

class G4Material;
class G4Element;

typedef std::map< G4String , G4Element* > Tst33MapSymbolElement;
typedef std::map< G4Element* , G4double > Tst33MapElementFraction;

class Tst33MaterialFactory{
public:
  Tst33MaterialFactory();
  ~Tst33MaterialFactory();
  
  G4Material *CreateConcrete();
  G4Material *CreateLightConcrete();
  G4Material *CreateGalactic();
  
private:
  Tst33MaterialFactory(const Tst33MaterialFactory &);

  void FillElementMap(const G4String &name, 
		      const G4String &symbol,
		      G4int Z,
		      G4double A);

  Tst33MaterialFactory &operator=(const Tst33MaterialFactory &);

  Tst33MapSymbolElement fMapSymbolElement;
  Tst33MapElementFraction fConcreteFractions;
};


#endif
