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
// $Id: G4ScorePrinter.hh,v 1.1 2002-06-13 07:12:33 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4ScorePrinter_hh
#define G4ScorePrinter_hh G4ScorePrinter_hh

#include "g4std/iostream"
#include "g4std/iomanip"
#include "g4std/vector"

#include "globals.hh"
#include "G4PMapPtkTallys.hh"

class G4VIStore;
class G4VCombiTally;
class G4VImpCombiTally;

typedef G4std::vector<G4VCombiTally *> G4VecCombiTally;
typedef G4std::vector<G4VImpCombiTally *> G4VecImpCombiTally;

struct G4TallyImp {
  G4PMapNameTally tally;
  G4double imp;
};

class G4ScorePrinter
{
public:
  G4ScorePrinter(const G4VIStore *aIStore = 0);
  ~G4ScorePrinter(){}
  void PrintHeader(G4std::ostream *out);
  void PrintTable(const G4PMapPtkTallys &aMapPtkTallys, 
		  G4std::ostream *out);

  G4std::string FillString(const G4std::string &name, 
			   char c, G4int n, G4bool back = true);
private:
  G4String CreateName(G4PTouchableKey ptk);
  void PrintLine(G4String &name,
		 G4std::ostream *out);
  void CalculateTalliesForAVolume(G4PMapNameTally &mtallies,
				  G4double importance);
  
  const G4VIStore *fIStore;
  G4int FieldName;
  G4int FieldValue;
  G4VecCombiTally fVecCombiTally;
  G4VecImpCombiTally fVecImpCombiTally;
};

#endif
