// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VBayon.cc,v 1.1 1999-10-03 09:13:38 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
// --------------------------------------------------------------

#include "G4VBaryon.hh"

const G4VBaryon & G4VBaryon::operator=(const G4VBaryon &right)
{
  if (this != &right) {
  } return right;
}
