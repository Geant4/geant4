//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: Tst33MaterialFactory.hh,v 1.2 2002-10-29 16:37:09 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class Tst33MaterialFactory
//
// Class description:
//
// ...

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33MaterialFactory_hh
#define Tst33MaterialFactory_hh Tst33MaterialFactory_hh

#include "globals.hh"
#include "g4std/map"

class G4Material;
class G4Element;

typedef G4std::map< G4String , G4Element* > Tst33MapSymbolElement;
typedef G4std::map< G4Element* , G4double > Tst33MapElementFraction;

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
