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
// $Id: G4ProcessVector.cc,v 1.3 2001-07-11 10:08:21 gunter Exp $
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


