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
