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
///////////////////////////////////////////////////////////////////////////////
// File: CCalAMaterial.hh
// Description: CCalAMaterial holds the basic information needed to make a
//              G4Material from atomic proportions.
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalAMaterial_h
#define CCalAMaterial_h 1

#include "CCalMaterial.hh"

#include "G4Element.hh"

class CCalAMaterial : public CCalMaterial
{
  friend std::ostream& operator<<(std::ostream&, const CCalAMaterial&);

public:
  //Construct from list of constituents
  CCalAMaterial(G4String mat, G4double dens, G4int nelem, 
                CCalAMaterial** constituents, G4double* weights);
  //Construct from one element
  CCalAMaterial(G4String elemat, G4double Aeff, G4double dens);
  //Copy constructor
  CCalAMaterial(const CCalAMaterial&);
  virtual ~CCalAMaterial();

  G4double Aeff() const {return aEff;}

  CCalAMaterial& operator= (const CCalAMaterial&);       //Assignment

protected:
  void computeAeff(G4int nconst, CCalAMaterial** constituents,
                   G4double* weights);

protected:
  double aEff;  //Effective mass

};

#endif
