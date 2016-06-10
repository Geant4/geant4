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
// $Id: G4MCCIndexConversionTable.hh 70369 2013-05-29 14:59:24Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
// Class description:
//
// G4MCCIndexConversionTable is used by G4ProductionTable 
// when the cut table is retrieved from the file. 
// An index pointing a Material-Cut-Couple can be different 
// from the index pointing the same MCC in the file. This class
// has a map between them.
//
// ------------------------------------------------------------
//
// History:
// -------
// - First implementation   20th August 2004  by H.Kurashige
//-------------------------------------

#ifndef G4MCCIndexConversionTable_h
#define G4MCCIndexConversionTable_h 1

#include <vector>
#include "globals.hh"
#include "G4ios.hh"

class G4MCCIndexConversionTable
{
 public: // with description

  G4MCCIndexConversionTable();
    // Default constructor.

  virtual ~G4MCCIndexConversionTable();
    // Destructor.

  void Reset(size_t size);
    // reset conversion table 
 
  G4bool IsUsed(size_t index) const;
    // returns 'true' if the indicated MCC in the file 
    // is used in the current production cut table
  
  void SetNewIndex(size_t index, size_t new_value);   
    // set the index in the current production cut table
    // for the indicated MCC in the file

  G4int GetIndex(size_t index) const;
    // get the index in the current production cut table
    // for the indicated MCC in the file
 
  size_t size() const;

  protected:
   typedef std::vector<G4int> G4IntVector;
   G4IntVector vecNewIndex;
};

inline
 G4bool G4MCCIndexConversionTable::IsUsed(size_t index) const
{
  if (index >= vecNewIndex.size()) return false;

  // returns 'true' if the indicated MCC in the file 
  // is used in the current production cut table
  return (vecNewIndex[index] >= 0); 
}

inline
 void G4MCCIndexConversionTable::SetNewIndex(size_t index, size_t new_value)
{
  if (index >= vecNewIndex.size()) return;
  // set the index in the current production cut table
  // for the indicated MCC in the file
  vecNewIndex[index]=new_value;  
}  

inline
  G4int G4MCCIndexConversionTable::GetIndex(size_t index) const
{
  if (index >= vecNewIndex.size()) return -1;
  // get the index in the current production cut table
  // for the indicated MCC in the file
  return (vecNewIndex[index]);  
}

inline
  size_t G4MCCIndexConversionTable::size() const
{
  return vecNewIndex.size();
}

#endif
