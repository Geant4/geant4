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
//////////////////////////////////////////////////////////////////////////////
// File: CCalMaterial.hh
// Description: CCalMaterial holds the basic information needed to make a
//              G4Material. 
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalMaterial_h
#define CCalMaterial_h 1
#include <iostream>
#include "globals.hh"

class CCalMaterial
{

friend std::ostream& operator<<(std::ostream&, const CCalMaterial&);
  
public:
  enum FractionType {FTWeight, FTVolume};

  //Constructors and destructors
  CCalMaterial(G4String mat, G4double dens, G4int nelem, 
               CCalMaterial** constituents, G4double* weights,
               FractionType=FTWeight);
  CCalMaterial(const CCalMaterial&);
  virtual ~CCalMaterial();

  //Get methods
  G4String Name() const         {return name;}           //Material name
  G4double Density() const      {return density;}        //Density in g/cm3
  G4int    NElements() const    {return nElem;}          //Number of Elements
  G4String Element(G4int i) const {return theElements[i];} //Should be protected
  G4double Weight(G4int i) const {return theWeights[i];} //Should be protected

  //Operators
  G4bool        operator==(const CCalMaterial&) const; //Compares ONLY names
  G4bool        operator!=(const CCalMaterial&) const; //Compares ONLY names
  CCalMaterial& operator= (const CCalMaterial&);       //Assignment

protected:
  CCalMaterial(){} //Default constructor
  void computeDensity(int nconst,
                      CCalMaterial** constituents, double* weights,
                      FractionType ft);
  void closeMaterial(); //Closes material construction.

protected:
  G4String  name;            //Material name
  G4double    density;       //Density in g/cm3
  G4int       nElem;         //Number of constituents.
  G4String* theElements;     //Basic constituents
  G4double*   theWeights;    //Elements' weight fractions
};
#endif
