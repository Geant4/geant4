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
// $Id$
//
// Base type key class
//
// Jane Tinslay, September 2006
//
#ifndef G4TYPEKEY_HH
#define G4TYPEKEY_HH

#include "globals.hh"
#include <ostream>

class G4TypeKey {

public:

  typedef unsigned long Key;

  // Constructor
  G4TypeKey():fMyKey(0) {}

  // Destructor
  virtual ~G4TypeKey() {}

  G4bool IsValid() {
    return (0 == fMyKey ? false : true);
  }
  
  // Operators
  Key operator()() const {return fMyKey;}  
  bool operator==(const G4TypeKey& rhs) const {return fMyKey == rhs.fMyKey;}
  bool operator!=(const G4TypeKey& rhs) const {return !operator==(rhs);}
  bool operator<(const G4TypeKey& rhs) const {return fMyKey < rhs.fMyKey;}
  bool operator>(const G4TypeKey& rhs) const {return fMyKey > rhs.fMyKey;}

  friend std::ostream& operator<<(std::ostream& out, const G4TypeKey& key){ 
    return out<< key.fMyKey;
  }

protected:

  Key NextKey() const {
    static Key nKey = 0;
    return ++nKey;
  }

  Key fMyKey;

};

#endif
