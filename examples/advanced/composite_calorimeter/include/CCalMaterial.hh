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
//////////////////////////////////////////////////////////////////////////////
// File: CCalMaterial.hh
// Description: CCalMaterial holds the basic information needed to make a
//              G4Material. 
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalMaterial_h
#define CCalMaterial_h 1
#include "g4std/iostream"
#include "globals.hh"

class CCalMaterial {

friend G4std::ostream& operator<<(G4std::ostream&, const CCalMaterial&);
  
public:
  enum FractionType {FTWeight, FTVolume};

  //Constructors and destructors
  CCalMaterial(G4String mat, double dens, int nelem, 
	       CCalMaterial** constituents, double* weights,
	       FractionType=FTWeight);
  CCalMaterial(const CCalMaterial&);
  virtual ~CCalMaterial();

  //Get methods
  G4String Name() const         {return name;}           //Material name.
  double   Density() const      {return density;}        //Density in g/cm3.
  int      NElements() const    {return nElem;}          //Number of Elements.
  G4String Element(int i) const {return theElements[i];} //Should be protected.
  double   Weight(int i) const  {return theWeights[i];}  //Should be protected.

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
  double    density;         //Density in g/cm3
  int       nElem;           //Number of constituents.
  G4String* theElements;     //Basic constituents
  double*   theWeights;      //Elements' weight fractions
};
#endif
