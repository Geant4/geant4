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
// R&D: Vladimir.Grichine@cern.ch


#ifndef G4HadDataReading_HH
#define G4HadDataReading_HH 1

#include "globals.hh"
#include "g4std/vector"
#include "g4std/map"

class G4PhysicsTable;
class G4PhysicsVector;
class G4DataVector;
class G4Isotope;
class G4Element;
class G4Material;



class G4HadDataReading
{
public:

  G4HadDataReading();

  virtual ~G4HadDataReading();

  G4double GetTkinBin() const {return fTkinBin;}
  G4int    GetNumber() const {return fNo;}

  G4DataVector*  GetTkinVector() const {return fTkinVector;}
  G4DataVector*  GetXscVector() const {return fXscVector;}
  G4DataVector*  GetMultiplicityVector() const {return fMultiplicityVector;}

  G4DataVector*  GetMomentumVector() const {return fMomentumCVector;}
  G4DataVector*  GetAngleVector() const {return fAngleVector;}

  G4PhysicsTable*  GetAngleTable() const {return fAngleTable;}
  G4PhysicsTable*  GetMomentumCTable() const {return fMomentumCTable;}

  G4std::vector<G4PhysicsTable*>* GetDoubleDiffXscBank() const 
  {return fDoubleDiffXscBank;};


protected:

  void LoadIntegralXSC(G4String&);
  void LoadMultiplicity(G4String&);

  void LoadDifferentialXSC(G4String&,G4bool);
  void LoadDoubleDiffXSC(G4String&);

  void SetEnergyUnit(G4String&);
  void SetAngleUnit(G4String&);
  void SetXscUnit(G4String&);
  void SetDdXscUnit(G4String&);

private:

  G4double fEnergyUnit, fAngleUnit, fXscUnit;
  G4double fXscPeAngleUnit, fXscPerMomCUnit, fDdXscUnit;

  G4double fTkinBin;
  G4int    fNo;

  G4DataVector*  fTkinVector;
  G4DataVector*  fTkinBinVector;

  G4DataVector*  fXscVector;
  G4DataVector*  fMultiplicityVector;

  G4DataVector*  fMomentumCVector;
  G4DataVector*  fMomentumCBinVector;

  G4DataVector*  fAngleVector;
  G4DataVector*  fAngleBinVector;

  G4PhysicsTable*  fAngleTable;
  G4PhysicsTable*  fMomentumCTable;

  G4std::vector<G4PhysicsTable*>* fDoubleDiffXscBank;


  /*
  typedef G4std::map<G4int,G4std::
                 vector<G4double>,G4std::less<G4int> > xsc_Table;

  xsc_Table xscTable;

  G4std::vector<G4int> nXSC;
  G4std::vector<G4int> numberOfSecondaries;
  */  
};

#endif





