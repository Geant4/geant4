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
// $Id: B08ScorePrinter.hh,v 1.1 2002/06/04 11:14:52 dressel Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

#ifndef B08ScorePrinter_hh
#define B08ScorePrinter_hh B08ScorePrinter_hh

#include "g4std/iostream"
#include "g4std/iomanip"

#include "globals.hh"
#include "G4PMapPtkTallys.hh"

class G4VIStore;

class B08ScorePrinter
{
public:
  B08ScorePrinter(const G4VIStore *aIStore = 0);
  ~B08ScorePrinter(){}
  void PrintHeader(G4std::ostream *out);
  void PrintTable(const G4PMapPtkTallys &aMapPtkTallys, 
		  G4std::ostream *out);
  G4std::string FillString(const G4std::string &name, 
			   char c, G4int n, G4bool back = true);
private:
  const G4VIStore *fIStore;
  G4int FieldName;
  G4int FieldValue;

};

#endif
