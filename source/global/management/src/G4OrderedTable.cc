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
// $Id: G4OrderedTable.cc 67970 2013-03-13 10:10:06Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation
//
//      G4OrderedTable
//
// ------------------------------------------------------------

#include "G4DataVector.hh"
#include "G4OrderedTable.hh"
#include <iostream>
#include <fstream>
#include <iomanip>

G4OrderedTable::G4OrderedTable()
  : std::vector<G4DataVector*>()
{
}

G4OrderedTable::G4OrderedTable(size_t cap)
  : std::vector<G4DataVector*>(cap, (G4DataVector*)(0) )
{
}

G4OrderedTable::~G4OrderedTable()
{
}

G4bool G4OrderedTable::Store(const G4String& fileName,
                             G4bool          ascii)
{
  std::ofstream fOut;  
  
  // open output file //
  if (!ascii)
    { fOut.open(fileName, std::ios::out|std::ios::binary); }
  else
    { fOut.open(fileName, std::ios::out); }

  // check if the file has been opened successfully 
  if (!fOut)
  {
#ifdef G4VERBOSE  
    G4cerr << "G4OrderedTable::::Store():";
    G4cerr << " Cannot open file: " << fileName << G4endl;
#endif
    fOut.close();
    return false;
  }

 // Number of elements
  size_t tableSize = size(); 
  if (!ascii)
  {
    fOut.write( (char*)(&tableSize), sizeof tableSize); 
  }
  else
  {
    fOut << tableSize << G4endl;
  }

  // Data Vector
  G4int vType = G4DataVector::T_G4DataVector;
  for (G4OrderedTableIterator itr=begin(); itr!=end(); ++itr)
  {
    if (!ascii)
    {
      fOut.write( (char*)(&vType), sizeof vType); 
    }
    else
    {
      fOut << vType << G4endl;
    }
    (*itr)->Store(fOut,ascii);
  }
  fOut.close();
  return true;
}



G4bool G4OrderedTable::Retrieve(const G4String& fileName,
                                G4bool          ascii)
{
  std::ifstream fIn;  
  // open input file //
  if (ascii)
    { fIn.open(fileName,std::ios::in|std::ios::binary); }
  else
    { fIn.open(fileName,std::ios::in); }

  // check if the file has been opened successfully 
  if (!fIn)
  {
#ifdef G4VERBOSE  
    G4cerr << "G4OrderedTable::Retrieve():";
    G4cerr << " Cannot open file: " << fileName << G4endl;
#endif
    fIn.close();
    return false;
  }

  // clear 
  clearAndDestroy();
  
  // Number of elements
  G4int tableSize=0; 
  if (!ascii)
  {
    fIn.read((char*)(&tableSize), sizeof tableSize); 
  }
  else
  {
    fIn >> tableSize;
  }
  if (tableSize<=0)
  {
#ifdef G4VERBOSE  
    G4cerr << "G4OrderedTable::Retrieve():";
    G4cerr << " Invalid table size: " << tableSize << G4endl;
#endif
    return false;
  }
  reserve(tableSize); 

  // Physics Vector
  for (G4int idx=0; idx<tableSize; ++idx)
  {
    G4int vType=0;
    if (!ascii)
    {
      fIn.read( (char*)(&vType), sizeof vType); 
    }
    else
    {
      fIn >>  vType;
    }
    if (vType != G4DataVector::T_G4DataVector)
    {
#ifdef G4VERBOSE  
      G4cerr << "G4OrderedTable::Retrieve():";
      G4cerr << " Illegal Data Vector type: " << vType << " in  ";
      G4cerr << fileName << G4endl;
#endif          
      fIn.close();
      return false;
    }

    G4DataVector* pVec = new G4DataVector;

    if (! (pVec->Retrieve(fIn,ascii)) )
    {
#ifdef G4VERBOSE  
      G4cerr << "G4OrderedTable::Retrieve(): ";
      G4cerr << " Rrror in retreiving " << idx
             << "-th Physics Vector from file: ";
      G4cerr << fileName << G4endl;
#endif          
      fIn.close();
      delete pVec;
      return false;
    }

    // add a PhysicsVector to this OrderedTable
    push_back(pVec);
  } 
  fIn.close();
  return true;
}

std::ostream& operator<<(std::ostream& out, 
                         G4OrderedTable& right)
{
  // Printout Data Vector
  size_t i=0;
  for (G4OrderedTableIterator itr=right.begin(); itr!=right.end(); ++itr)
  {
    out << std::setw(8) << i << "-th Vector   ";
    out << ": Type    " << G4DataVector::T_G4DataVector << G4endl;
    out << *(*itr);
    i +=1;
  }
  out << G4endl;
  return out; 
}
