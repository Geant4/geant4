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
// ------------------------------------------------------------
//	GEANT 4 class header file 
// Class Description:
//      This class is a ordering parameter only used by G4PhysicsListHelper
// ------------------------------------------- 
//	History
//        first version                   29 Apr. 2011 by H.Kurashige 
// ------------------------------------------------------------

#ifndef G4PhysicsListOrderingParameter_h
#define G4PhysicsListOrderingParameter_h 1
#include "globals.hh"
#include "G4ios.hh"

class G4PhysicsListHelper; 
class G4PhysicsListOrderingParameter
{
  friend class G4PhysicsListHelper;  
  
public:
  // Hide constructor and destructor 
  G4PhysicsListOrderingParameter();
  virtual ~G4PhysicsListOrderingParameter();

  G4String  GetTypeName() const {return processTypeName;}
  G4int     GetType() const { return processType;}
  G4int     GetSubType() const { return processSubType;}
  G4int     GetOrdering(int idx) const;
  G4bool    GetDuplicable() const {return isDuplicable;}
  
private:
  G4String  processTypeName;
  G4int     processType;
  G4int     processSubType;
  G4int     ordering[3];
  G4bool    isDuplicable;
};

inline
 G4int     G4PhysicsListOrderingParameter::GetOrdering(int idx) const
{
  if ((idx<-1)||(idx>2)) return -1;
  else return ordering[idx]; 
}

#endif

