// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProcessVector.cc,v 1.2 2000-11-03 06:15:34 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//  This class is a container for pointers to physics process objects.
// ------------------------------------------------------------

#include "G4VProcess.hh"
#include "G4ProcessVector.hh"


/////// Constructors
G4ProcessVector::G4ProcessVector(size_t):pProcVector(0)
{
   pProcVector = new G4ProcVector();
}

G4ProcessVector::G4ProcessVector(const G4ProcessVector& right)
{
  *this == right;
}

/////// destructor
G4ProcessVector::~G4ProcessVector()
{
  // delete pProcVector
  if (pProcVector != 0 ) {
    pProcVector->clear();
    delete pProcVector;
  }
}

////// assignment oeprator
G4ProcessVector & G4ProcessVector::operator=(G4ProcessVector &right)
{
  if (this != &right) {
    // delete pProcVector
    if (pProcVector != 0 ) {
      pProcVector->clear();
      delete pProcVector;
    }
    
    pProcVector = new G4ProcVector();
    
    // copy all contents in  pProcVector
    if (pProcVector != 0) {
      pProcVector->clear();
      delete pProcVector;
    }
    pProcVector = new G4ProcVector();
    G4ProcVector::iterator i;
    for (i = pProcVector->begin(); i!= pProcVector->end(); ++i) {
      pProcVector->push_back(*i);
    }
  }
  return *this;
}


