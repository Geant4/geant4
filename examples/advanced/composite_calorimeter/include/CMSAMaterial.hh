///////////////////////////////////////////////////////////////////////////////
// File: CMSAMaterial.hh
// Date: 12/03/98 I. Gonzalez
// Description: CMSAMaterial holds the basic information needed to make a
//              G4Material from atomic proportions.
// Modifications: 31/08/98 I.G. -> G4type moved to type
///////////////////////////////////////////////////////////////////////////////
#ifndef CMSAMaterial_h
#define CMSAMaterial_h 1

#include "CMSMaterial.hh"

#include "G4Element.hh"

class CMSAMaterial : public CMSMaterial {
  friend ostream& operator<<(ostream&, const CMSAMaterial&);

public:
  //Construct from list of constituents
  CMSAMaterial(G4String mat, double dens, int nelem, 
	       CMSAMaterial** constituents, double* weights);
  //Construct from one element
  CMSAMaterial(G4String elemat, double Aeff, double dens);
  //Copy constructor
  CMSAMaterial(const CMSAMaterial&);
  virtual ~CMSAMaterial();

  G4double Aeff() const {return aEff;}

  CMSAMaterial& operator= (const CMSAMaterial&);       //Assignment

protected:
  void computeAeff(G4int nconst, CMSAMaterial** constituents, double* weights);

protected:
  double aEff;  //Effective mass

};

#endif
