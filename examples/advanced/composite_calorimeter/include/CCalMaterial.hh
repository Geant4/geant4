//////////////////////////////////////////////////////////////////////////////
// File: CCalMaterial.hh
// Description: CCalMaterial holds the basic information needed to make a
//              G4Material. 
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalMaterial_h
#define CCalMaterial_h 1
#include <iostream>
#include "globals.hh"

class CCalMaterial {

friend ostream& operator<<(ostream&, const CCalMaterial&);
  
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
