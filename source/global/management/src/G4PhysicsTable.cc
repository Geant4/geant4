// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsTable.cc,v 1.3 2001-03-09 03:39:31 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation
//
//	G4PhysicsTable
//
// ------------------------------------------------------------

#include "G4PhysicsVector.hh"
#include "G4PhysicsTable.hh"
#include "g4std/iostream"
#include "g4std/fstream"
#include "g4std/iomanip"

G4PhysicsTable::G4PhysicsTable()
  : G4std::vector<G4PhysicsVector*>()
{
}

G4PhysicsTable::G4PhysicsTable(size_t capacity)
  : G4std::vector<G4PhysicsVector*>()
{
  reserve(capacity);
}

G4PhysicsTable::G4PhysicsTable(const G4PhysicsTable& right)
{
  *this = right;
}

G4PhysicsTable& G4PhysicsTable::operator=(const G4PhysicsTable& right)
{
  if (this != &right)
  {
    G4std::vector<G4PhysicsVector*>::const_iterator itr;
    for (itr=right.begin(); itr!=right.end(); ++itr)
    {
      G4std::vector<G4PhysicsVector*>::push_back(*itr);
    }
  }
  return *this;
}

G4PhysicsTable::~G4PhysicsTable()
{
  clear();
}
 

G4bool G4PhysicsTable::StorePhysicsTable(const G4String& fileName,
					 G4bool          ascii)
{
  G4std::ofstream fOut;  
  
  // open output file //
  if (!ascii){
    fOut.open(fileName, G4std::ios::out|G4std::ios::binary);
  }else{
    fOut.open(fileName, G4std::ios::out);
  }
  // check if the file has been opened successfully 
  if (!fOut) {
#ifdef G4VERBOSE  
    G4cerr << "G4PhysicsTable::::StorePhysicsTable  ";
    G4cerr << " Can not open file " << fileName << G4endl;
#endif
    fOut.close();
    return false;
  }

  // Number of elements
  size_t tableSize = size(); 
  if (!ascii){
    fOut.write( (char*)(&tableSize), sizeof tableSize); 
  } else {
    fOut << tableSize << G4endl;
  }

  // Physics Vector
  G4PhysicsTableIterator itr;
  for (itr=begin(); itr!=end(); ++itr) {
    G4int vType = (*itr)->GetType();
    if (!ascii){
      fOut.write( (char*)(&vType), sizeof vType); 
    } else {
      fOut << vType << G4endl;
    }
    (*itr)->Store(fOut,ascii);
  }
  fOut.close();
  return true;
}



G4bool G4PhysicsTable::RetrievePhysicsTable(const G4String& fileName,
					    G4bool          ascii)
{
  G4std::ifstream fIn;  
  // open input file //
  if (ascii) {
    fIn.open(fileName,G4std::ios::in|G4std::ios::binary);
  } else {
    fIn.open(fileName,G4std::ios::in);
  }
  // check if the file has been opened successfully 
  if (!fIn) {
#ifdef G4VERBOSE  
    G4cerr << "G4PhysicsTable::RetrievePhysicsTable  ";
    G4cerr << " Can not open file " << fileName << G4endl;
#endif
    fIn.close();
    return false;
  }

  // clear 
  clearAndDestroy();
  
  // Number of elements
  size_t tableSize; 
  if (!ascii){
    fIn.read((char*)(&tableSize), sizeof tableSize); 
  } else {
    fIn >> tableSize;
  }
  reserve(tableSize); 

  // Physics Vector
  for (size_t idx=0; idx<tableSize; ++idx) {
    G4int vType;
    if (!ascii){
      fIn.read( (char*)(&vType), sizeof vType); 
    } else {
      fIn >>  vType;
    }
    G4PhysicsVector* pVec = CreatePhysicsVector(vType);
    if (pVec==0) {
#ifdef G4VERBOSE  
      G4cerr << "G4PhysicsTable::RetrievePhysicsTable  ";
      G4cerr << " illegal Physics Vector type " << vType << " in  ";
      G4cerr << fileName << G4endl;
#endif          
      fIn.close();
      return false;
    }

    if (! (pVec->Retrieve(fIn,ascii)) ){
#ifdef G4VERBOSE  
      G4cerr << "G4PhysicsTable::RetrievePhysicsTable  ";
      G4cerr << " error in retreiving " << idx << "-th Physics Vector from file ";
      G4cerr << fileName << G4endl;
#endif          
      fIn.close();
      return false;
    }

    // add a PhysicsVector to this PhysicsTable
    push_back(pVec);
  } 
  fIn.close();
  return true;
}

G4std::ostream& operator<<(G4std::ostream& out, 
			   G4PhysicsTable& right)
{
  // Printout Physics Vector
  G4PhysicsTableIterator itr;
  size_t i=0;
  for (itr=right.begin(); itr!=right.end(); ++itr) {
    out << G4std::setw(8) << i << "-th Vector   ";
    out << ": Type    " << (*itr)->GetType() << G4endl;
    out << *(*itr);
    i +=1;
  }
  out << G4endl;
  return out; 
}

#include "G4PhysicsVectorType.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4PhysicsLnVector.hh"
 
G4PhysicsVector* G4PhysicsTable::CreatePhysicsVector(G4int type)
{
  G4PhysicsVector* pVector=0;
  switch (type) {

  case T_G4PhysicsLinearVector: 
    pVector = new G4PhysicsLinearVector();
    break;

  case T_G4PhysicsLogVector: 
    pVector = new G4PhysicsLogVector();
    break;

  case T_G4PhysicsLnVector: 
    pVector = new G4PhysicsLnVector();
    break;

  case T_G4PhysicsFreeVector: 
    pVector = new G4PhysicsFreeVector();
    break;

  case T_G4PhysicsOrderedFreeVector: 
    pVector = new G4PhysicsOrderedFreeVector();
    break;

  case T_G4LPhysicsFreeVector: 
    pVector = new G4LPhysicsFreeVector();
    break;
  
  default:
    break;
 
  }
  return pVector;
}
 
