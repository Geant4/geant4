// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPPolynomExpansion.hh,v 1.1 1999-01-07 16:13:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPPolynomExpansion_h
#define G4NeutronHPPolynomExpansion_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <fstream.h>

class G4NeutronHPPolynomExpansion
{
  public:
  G4NeutronHPPolynomExpansion()
  {
    theCoeff = NULL;
    nPoly=0;
  }
  ~G4NeutronHPPolynomExpansion()
  {
    if(theCoeff!=NULL) delete [] theCoeff;
  }
  
  inline void Init(ifstream & theData)
  {
    theData >> nPoly;
    theCoeff = new G4double[nPoly];
    G4int i;
    for(i=0;i<nPoly;i++)
    {
      theData >> theCoeff[i];
    }
  }
  
  inline G4double GetValue(G4double anEnergy) 
  {
    G4int i;
    G4double result=0;
    G4double base = anEnergy/eV;
    G4double running = 1;
    for(i=0; i<nPoly; i++)
    {
      result+=theCoeff[i]*running;
      running *= base;
    }
    return result;
  }
  private:
  
  G4int nPoly;
  G4double * theCoeff;
};

#endif
