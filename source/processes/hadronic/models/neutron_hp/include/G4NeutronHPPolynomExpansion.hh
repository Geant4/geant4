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
// $Id: G4NeutronHPPolynomExpansion.hh,v 1.9 2005/06/04 13:44:43 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
#ifndef G4NeutronHPPolynomExpansion_h
#define G4NeutronHPPolynomExpansion_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>

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
  
  inline void Init(std::ifstream & theData)
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
