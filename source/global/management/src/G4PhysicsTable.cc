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
// $Id: G4PhysicsTable.cc 98864 2016-08-15 11:53:26Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation
//
//      G4PhysicsTable
//
// ------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <iomanip>

#include "G4PhysicsVector.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVectorType.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4PhysicsLnVector.hh"
 
G4PhysicsTable::G4PhysicsTable()
  : G4PhysCollection()
{
}

G4PhysicsTable::G4PhysicsTable(size_t cap)
  : G4PhysCollection()
{
  reserve(cap);
  vecFlag.reserve(cap);
}

/*
G4PhysicsTable::G4PhysicsTable(const G4PhysicsTable& right)
  : G4PhysCollection()
{
  *this = right;
}

G4PhysicsTable& G4PhysicsTable::operator=(const G4PhysicsTable& right)
{
  if (this != &right)
  {
    size_t idx = 0;
    for (G4PhysCollection::const_iterator itr=right.begin();
         itr!=right.end(); ++itr )
    {
      G4PhysCollection::push_back(*itr);
      vecFlag.push_back(right.GetFlag(idx));
      idx +=1;
    }
  }
  return *this;
}
*/
G4PhysicsTable::~G4PhysicsTable()
{
  G4PhysCollection::clear();
  vecFlag.clear();
}
 
void   G4PhysicsTable::resize(size_t siz, G4PhysicsVector* vec)
{
  G4PhysCollection::resize(siz, vec);
  vecFlag.resize(siz, true);
}

G4bool G4PhysicsTable::StorePhysicsTable(const G4String& fileName,
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
    G4cerr << "G4PhysicsTable::StorePhysicsTable():";
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

  // Physics Vector
  for (G4PhysicsTableIterator itr=begin(); itr!=end(); ++itr)
  {
    G4int vType = (*itr)->GetType();
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


G4bool G4PhysicsTable::ExistPhysicsTable(const G4String& fileName) const
{
  std::ifstream fIn;  
  G4bool value=true;
  // open input file
  fIn.open(fileName,std::ios::in);

  // check if the file has been opened successfully 
  if (!fIn)
  {
    value = false;
  }
  fIn.close();
  return value;
}
    
G4bool G4PhysicsTable::RetrievePhysicsTable(const G4String& fileName,
                                            G4bool          ascii)
{
  std::ifstream fIn;  
  // open input file
  if (ascii)
    { fIn.open(fileName,std::ios::in|std::ios::binary); }
  else
    { fIn.open(fileName,std::ios::in);} 

  // check if the file has been opened successfully 
  if (!fIn)
  {
#ifdef G4VERBOSE  
    G4cerr << "G4PhysicsTable::RetrievePhysicsTable():";
    G4cerr << " Cannot open file: " << fileName << G4endl;
#endif
    fIn.close();
    return false;
  }

  // clear 
  clearAndDestroy();
  
  // Number of elements
  size_t tableSize=0; 
  if (!ascii)
  {
    fIn.read((char*)(&tableSize), sizeof tableSize); 
  }
  else
  {
    fIn >> tableSize;
  }
  reserve(tableSize); 
  vecFlag.clear();

  // Physics Vector
  for (size_t idx=0; idx<tableSize; ++idx)
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
    G4PhysicsVector* pVec = CreatePhysicsVector(vType);
    if (pVec==nullptr)
    {
#ifdef G4VERBOSE  
      G4cerr << "G4PhysicsTable::RetrievePhysicsTable():";
      G4cerr << " Illegal Physics Vector type: " << vType << " in: ";
      G4cerr << fileName << G4endl;
#endif          
      fIn.close();
      return false;
    }

    if (! (pVec->Retrieve(fIn,ascii)) )
    {
#ifdef G4VERBOSE  
      G4cerr << "G4PhysicsTable::RetrievePhysicsTable():";
      G4cerr << " Rrror in retreiving " << idx
             << "-th Physics Vector from file: ";
      G4cerr << fileName << G4endl;
#endif          
      fIn.close();
      return false;
    }

    // add a PhysicsVector to this PhysicsTable
    G4PhysCollection::push_back(pVec);
    vecFlag.push_back(true);
    
  } 
  fIn.close();
  return true;
}

std::ostream& operator<<(std::ostream& out, 
                         G4PhysicsTable& right)
{
  // Printout Physics Vector
  size_t i=0;
  for (G4PhysicsTableIterator itr=right.begin(); itr!=right.end(); ++itr)
  {
    out << std::setw(8) << i << "-th Vector   ";
    out << ": Type    " << G4int((*itr)->GetType()) ;
    out << ": Flag    ";
    if (right.GetFlag(i))
    {
      out << " T";
    } 
    else
    {
      out << " F";
    } 
    out << G4endl;
    out << *(*itr);
    i +=1;
  }
  out << G4endl;
  return out; 
}

void G4PhysicsTable::ResetFlagArray()
{
  size_t tableSize = G4PhysCollection::size(); 
  vecFlag.clear();
  for (size_t idx=0; idx<tableSize; idx++)
  {
    vecFlag.push_back(true);
  }
}

G4PhysicsVector* G4PhysicsTable::CreatePhysicsVector(G4int type)
{
  G4PhysicsVector* pVector = nullptr;
  switch (type)
  {
  case T_G4PhysicsLinearVector: 
    pVector = new G4PhysicsLinearVector();
    break;

  case T_G4PhysicsLogVector: 
    pVector = new G4PhysicsLogVector();
    break;

  case T_G4PhysicsLnVector: 
    pVector = new G4PhysicsLogVector();
    break;

  case T_G4PhysicsFreeVector: 
    pVector = new G4PhysicsFreeVector();
    break;

  case T_G4PhysicsOrderedFreeVector: 
    pVector = new G4PhysicsOrderedFreeVector();
    break;

  case T_G4LPhysicsFreeVector: 
    pVector = new G4PhysicsFreeVector();
    break;
  
  default:
    break;
  }
  return pVector;
}
