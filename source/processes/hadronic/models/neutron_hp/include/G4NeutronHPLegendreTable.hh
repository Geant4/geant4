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
// $Id: G4NeutronHPLegendreTable.hh,v 1.7 2002-12-12 19:18:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPLegendreTable_h
#define G4NeutronHPLegendreTable_h 1

#include "globals.hh"
#include "g4std/fstream"
#include "G4ios.hh"
#include "G4InterpolationManager.hh"

class G4NeutronHPLegendreTable
{
  public:
  G4NeutronHPLegendreTable()
  {
    nCoeff=0; 
    theCoeff = NULL;
  }
  ~G4NeutronHPLegendreTable(){if(theCoeff!=NULL) delete [] theCoeff;}
  
  void operator= (const G4NeutronHPLegendreTable & aSet)
  {
    if(&aSet!=this)
    {
      theRep = aSet.theRep;
      theEnergy = aSet.theEnergy;
      theTemp = aSet.theTemp;
      theManager = aSet.theManager;
      nCoeff = aSet.nCoeff;
      if(theCoeff!=NULL) delete [] theCoeff;
      theCoeff = new G4double[nCoeff];
      for(G4int i=0; i<nCoeff; i++)
      {
        theCoeff[i] = aSet.theCoeff[i];
      }
    }
  }
  
  inline void Init(G4std::ifstream & aDataFile) 
  {
    G4double eNeu, coeff;
    G4int nPoly;
    aDataFile >> eNeu >> nPoly;
    eNeu *= eV;
    Init(eNeu, nPoly);
    for(G4int l=0; l<nPoly; l++)
    {
      aDataFile >> coeff;
      SetCoeff(l+1, coeff);
    }
  }
  
  inline void Init(G4double e, G4int n)
  {
    theCoeff = new G4double[n+1];
    theCoeff[0]=1.;
    theEnergy = e;
    nCoeff = n+1;
//    G4cout << "G4NeutronHPLegendreTable::Init called "<<e<<" "<<n<<G4endl;
  }
  inline void SetEnergy(G4double energy){ theEnergy = energy; }
  inline void SetTemperature(G4double temp){ theTemp = temp; }
  inline void SetCoeff(G4int l, G4double coeff) {theCoeff[l]=coeff;}
  inline void SetRepresentation(G4int aRep) {theRep = aRep;}
  
  inline G4double GetCoeff(G4int l) {return theCoeff[l];}
  inline G4double GetEnergy(){return theEnergy;}
  inline G4double GetTemperature(){return theTemp;}
  inline G4int GetNumberOfPoly() {return nCoeff;}
  inline G4int GetRepresentation() {return theRep;}
  inline const G4InterpolationManager & GetManager() {return theManager;}
  private:
  
  G4int theRep;
  G4double theEnergy;
  G4double theTemp;
  G4int nCoeff;
  G4InterpolationManager theManager; // knows the interpolation between stores
  G4double * theCoeff;
};

#endif
