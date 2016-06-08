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
// $Id: B08ImportanceMessenger.hh,v 1.1 2002/06/04 11:14:51 dressel Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

#ifndef B08ImportanceMessenger_hh
#define B08ImportanceMessenger_hh B08ImportanceMessenger_hh

#include "G4UImessenger.hh"
#include "globals.hh"
#include "g4std/map"

class G4VIStore;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcommand;
class G4VPhysicalVolume;

class B08BaseExp {
public:
  B08BaseExp(G4double b = 1, G4double e = 1) : fBase(b), fExp(e){}
  void SetBase(G4double b){fBase =  b;}
  void SetExponent(G4double e){fExp = e;}
  G4double GetImportance() const {return pow(fBase, fExp);}
private:
  G4double fBase;
  G4double fExp;
};

class B08ImportanceMessenger: public G4UImessenger{
public: 
  B08ImportanceMessenger(G4VIStore &istore);
  ~B08ImportanceMessenger(){}
  void SetNewValue(G4UIcommand* pCmd,G4String szValue);
  void AddCell(const G4String &cellname, G4VPhysicalVolume *p);
  void SetImportanceBase(const G4String &cellname, G4double base);
  void SetImportanceExponent(const G4String &cellname, G4double exp);
  
  typedef G4std::map<G4UIcmdWithADouble *, G4VPhysicalVolume *> B08BaseMap;
  typedef G4std::map<G4UIcmdWithADouble *, G4VPhysicalVolume *> B08ExpMap;
  typedef G4std::map<G4VPhysicalVolume *, B08BaseExp> B08BaseExpMap;
  typedef G4std::map<G4String, G4VPhysicalVolume *> B08CellVolMap;
  
private:
  G4VIStore &fIStore;
  
  B08BaseMap fBaseMap;
  B08ExpMap fExpMap;
  B08BaseExpMap fBaseExpMap;
  B08CellVolMap fCellVolMap;
};


#endif
