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
// $Id: G4AttCheck.cc,v 1.1 2005-03-22 16:51:34 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4AttCheck.hh"

#include "G4AttDef.hh"
#include "G4AttValue.hh"

G4AttCheck::G4AttCheck
(const std::vector<G4AttValue>* values,
 const std::map<G4String,G4AttDef>* definitions):
  fpValues(values),
  fpDefinitions(definitions)
{}


G4AttCheck::~G4AttCheck() {}

std::ostream& operator<< (std::ostream& os, const G4AttCheck& ac) {
  using namespace std;
  vector<G4AttValue>::const_iterator i;
  for (i = ac.fpValues->begin(); i != ac.fpValues->end(); ++i) {
    map<G4String,G4AttDef>::const_iterator j =
      ac.fpDefinitions->find(i->GetName());
    if (j != ac.fpDefinitions->end()) {
      os << j->second.GetDesc() << ": " << i->GetValue() << endl;
    } else {
      os << "WARNING: No G4AttDef for G4AttValue \""
	 <<  i->GetName() << "\": " << i->GetValue() << endl;
    }
  }
  return os;
}
