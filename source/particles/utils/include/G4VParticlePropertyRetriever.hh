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
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VParticlePropertyRetriever.hh,v 1.1 2004/03/11 09:47:44 kurasige Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// ---------------------------------------------------------------
#ifndef G4VParticlePropertyRetriever_h
#define G4VParticlePropertyRetriever_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "G4ParticlePropertyData.hh"
#include "G4ParticlePropertyTable.hh"

class G4VParticlePropertyRetriever
{
 public:
  //constructors
  G4VParticlePropertyRetriever();
  
  //destructor
  virtual ~G4VParticlePropertyRetriever();
  
 public:
  // equality operators
  G4int operator==(const G4VParticlePropertyRetriever &right) const 
  {   return (this == &right);    }
  
  G4int operator!=(const G4VParticlePropertyRetriever &right) const 
  {   return (this != &right);    }
  
 public:
  virtual void Retrieve(const G4String& option) = 0;
  // print out particle properties in the list

 protected:
  G4ParticlePropertyTable*  pPropertyTable;

};

inline
 G4VParticlePropertyRetriever::G4VParticlePropertyRetriever()
{
  pPropertyTable =   G4ParticlePropertyTable::GetParticlePropertyTable();
}

/////////////////////////////
inline
 G4VParticlePropertyRetriever::~G4VParticlePropertyRetriever()
{
  pPropertyTable->Clear();
}    


#endif

