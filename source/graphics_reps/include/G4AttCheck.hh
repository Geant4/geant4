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
// $Id: G4AttCheck.hh,v 1.1 2005-03-22 16:51:34 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4ATTCHECK_HH
#define G4ATTCHECK_HH

// Class Description:
//
// @class G4AttCheck
//
// @brief Checks G4AttValue's and their corresponding G4AttDef map.
//
// Usage (a): G4AttCheck(values,definitions);
//    or (b): G4cout << G4AttCheck(values,definitions) << G4endl;
//
// For further details, see the HepRep home page at http://heprep.freehep.org
//  
// @author J.Allison
// @author J.Perl
// Class Description - End:

#include "globals.hh"

#include <vector>
#include<map>
#include<iostream>

class G4AttValue;
class G4AttDef;

class G4AttCheck {
public:
  G4AttCheck
  (const std::vector<G4AttValue>* values,
   const std::map<G4String,G4AttDef>* definitions);
  ~G4AttCheck();
  friend std::ostream& operator<< (std::ostream&, const G4AttCheck&);
private:
  const std::vector<G4AttValue>* fpValues;
  const std::map<G4String,G4AttDef>* fpDefinitions;
};

#endif //G4ATTCHECK_HH
