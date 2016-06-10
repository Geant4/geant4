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
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VParticlePropertyReporter.hh 67971 2013-03-13 10:13:24Z gcosmo $
//
// 
// ---------------------------------------------------------------
#ifndef G4VParticlePropertyReporter_h
#define G4VParticlePropertyReporter_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "G4ParticlePropertyData.hh"
#include "G4ParticlePropertyTable.hh"

class G4VParticlePropertyReporter
{
 public:
  //constructors
  G4VParticlePropertyReporter();
  
  //destructor
  virtual ~G4VParticlePropertyReporter();
  
 public:
  // equality operators
  G4int operator==(const G4VParticlePropertyReporter &right) const 
  {   return (this == &right);    }
  
  G4int operator!=(const G4VParticlePropertyReporter &right) const 
  {   return (this != &right);    }
  
 public:
  typedef std::vector<G4ParticlePropertyData*> G4PPDContainer;

 public:
  virtual G4bool FillList(G4String name); 
  // fill out properties for specified particle type/name  
  // (return false if specified name/type does not exist) 

  virtual void Clear();
  // clear the list 

  virtual void Print(const G4String& option) = 0;
  // print out particle properties in the list

  const G4PPDContainer&     GetList() const {return pList;}
  
 protected:
  G4PPDContainer            pList;
  G4ParticlePropertyTable*  pPropertyTable;

};

#endif

