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
// $Id: VisActionWithAttDefs.hh,v 1.1 2005-03-23 17:43:25 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef VISACTIONWITHATTDEFS_HH
#define VISACTIONWITHATTDEFS_HH

#include "G4VUserVisAction.hh"

#include <map>
#include <vector>
#include "G4Transform3D.hh"

class G4string;
class G4AttDef;
class G4AttValue;
class G4Box;
class G4VisAttributes;

class VisActionWithAttDefs: public G4VUserVisAction {
public:
  VisActionWithAttDefs();
  ~VisActionWithAttDefs();
  void Draw();
private:
  G4Box* fpBox;
  G4VisAttributes* fpVisAtts;
  std::map<G4String,G4AttDef>* fpAttDefs;
  std::vector<G4AttValue>* fpAttValues;
  G4Transform3D* fpTransform;
};

#endif

