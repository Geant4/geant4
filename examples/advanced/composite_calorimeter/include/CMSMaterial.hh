//////////////////////////////////////////////////////////////////////////////
// File: CMSMaterial.hh
// Description: CMSMaterial holds the basic information needed to make a
//              G4Material. Temporary solution?... We'll see.
// Date: 12/03/98 
// Modifications: 31/08/98 I.G. -> G4type moved to type for int, double.
//                                 NElements modified to return int, not double
///////////////////////////////////////////////////////////////////////////////
#ifndef CMSMaterial_h
#define CMSMaterial_h 1
#include <iostream>
#include "globals.hh"

class CMSMaterial {

friend ostream& operator<<(ostream&, const CMSMaterial&);
  
public:
  enum FractionType {FTWeight, FTVolume};

  //Constructors and destructors
  CMSMaterial(G4String mat, double dens, int nelem, 
	      CMSMaterial** constituents, double* weights,
	      FractionType=FTWeight);
  CMSMaterial(const CMSMaterial&);
  virtual ~CMSMaterial();

  //Get methods
  G4String Name() const         {return name;}           //Material name.
  double   Density() const      {return density;}        //Density in g/cm3.
  int      NElements() const    {return nElem;}          //Number of Elements.
  G4String Element(int i) const {return theElements[i];} //Should be protected.
  double   Weight(int i) const  {return theWeights[i];}  //Should be protected.

  //Operators
  G4bool       operator==(const CMSMaterial&) const; //Compares ONLY names
  G4bool       operator!=(const CMSMaterial&) const; //Compares ONLY names
  CMSMaterial& operator= (const CMSMaterial&);       //Assignment

protected:
  CMSMaterial(){} //Default constructor
  void computeDensity(int nconst,
		      CMSMaterial** constituents, double* weights,
		      FractionType ft);
  void closeMaterial(); //Closes material construction.

protected:
  G4String  name;            //Material name
  double    density;         //Density in g/cms3
  int       nElem;           //Number of constituents.
  G4String* theElements;     //Basic constituents
  double*   theWeights;      //Elements' weight fractions
};
#endif
