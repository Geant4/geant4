///////////////////////////////////////////////////////////////////////////////
// File: CCalAMaterial.hh
// Description: CCalAMaterial holds the basic information needed to make a
//              G4Material from atomic proportions.
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalAMaterial_h
#define CCalAMaterial_h 1

#include "CCalMaterial.hh"

#include "G4Element.hh"

class CCalAMaterial : public CCalMaterial {
  friend G4std::ostream& operator<<(G4std::ostream&, const CCalAMaterial&);

public:
  //Construct from list of constituents
  CCalAMaterial(G4String mat, double dens, int nelem, 
		CCalAMaterial** constituents, double* weights);
  //Construct from one element
  CCalAMaterial(G4String elemat, double Aeff, double dens);
  //Copy constructor
  CCalAMaterial(const CCalAMaterial&);
  virtual ~CCalAMaterial();

  G4double Aeff() const {return aEff;}

  CCalAMaterial& operator= (const CCalAMaterial&);       //Assignment

protected:
  void computeAeff(G4int nconst, CCalAMaterial** constituents, double* weights);

protected:
  double aEff;  //Effective mass

};

#endif
