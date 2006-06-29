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
// $Id: G4NeutronHPPolynomExpansion.hh,v 1.10 2006-06-29 20:49:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
