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
// $Id: G4OrderedTable.cc,v 1.1 2001-09-17 08:17:58 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation
//
//	G4OrderedTable
//
// ------------------------------------------------------------

#include "G4DataVector.hh"
#include "G4OrderedTable.hh"
#include "g4std/iostream"
#include "g4std/fstream"
#include "g4std/iomanip"


G4bool G4OrderedTable::Store(const G4String& fileName,
			     G4bool          ascii)
{
  G4std::ofstream fOut;  
  
  // open output file //
#ifdef G4USE_STD_NAMESPACE
  if (!ascii)
    fOut.open(fileName, G4std::ios::out|G4std::ios::binary);
  else
#endif
    fOut.open(fileName, G4std::ios::out);

  // check if the file has been opened successfully 
  if (!fOut) {
#ifdef G4VERBOSE  
    G4cerr << "G4OrderedTable::::Store  ";
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

  // Data Vector
  G4OrderedTableIterator itr;
  G4int vType = G4DataVector::T_G4DataVector;
  for (itr=begin(); itr!=end(); ++itr) {
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



G4bool G4OrderedTable::Retrieve(const G4String& fileName,
				G4bool          ascii)
{
  G4std::ifstream fIn;  
  // open input file //
#ifdef G4USE_STD_NAMESPACE
  if (ascii)
    fIn.open(fileName,G4std::ios::in|G4std::ios::binary);
  else
#endif
    fIn.open(fileName,G4std::ios::in);

  // check if the file has been opened successfully 
  if (!fIn) {
#ifdef G4VERBOSE  
    G4cerr << "G4OrderedTable::Retrieve  ";
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
    if (vType != G4DataVector::T_G4DataVector) {
#ifdef G4VERBOSE  
      G4cerr << "G4OrderedTable::Retrieve  ";
      G4cerr << " illegal Data Vector type " << vType << " in  ";
      G4cerr << fileName << G4endl;
#endif          
      fIn.close();
      return false;
    }

    G4DataVector* pVec = new G4DataVector;

    if (! (pVec->Retrieve(fIn,ascii)) ){
#ifdef G4VERBOSE  
      G4cerr << "G4OrderedTable::Retrieve  ";
      G4cerr << " error in retreiving " << idx << "-th Physics Vector from file ";
      G4cerr << fileName << G4endl;
#endif          
      fIn.close();
      return false;
    }

    // add a PhysicsVector to this OrderedTable
    push_back(pVec);
  } 
  fIn.close();
  return true;
}

G4std::ostream& operator<<(G4std::ostream& out, 
			   G4OrderedTable& right)
{
  // Printout Data Vector
  G4OrderedTableIterator itr;
  size_t i=0;
  for (itr=right.begin(); itr!=right.end(); ++itr) {
    out << G4std::setw(8) << i << "-th Vector   ";
    out << ": Type    " << G4DataVector::T_G4DataVector << G4endl;
    out << *(*itr);
    i +=1;
  }
  out << G4endl;
  return out; 
}


