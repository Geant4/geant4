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
#include <vector>
#include <map> 
#include "G4DataVector.hh"
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

  G4String GetAnyNumber() const { return fAnyNumber;}
  G4String GetAnyEmptySpace() const { return fAnyEmptySpace;}
  G4String GetAnyHidden() const { return fAnyHidden;}

  std::vector<G4int> GetEnergyNoVector() const { return fEnergyNoVector;}
  std::vector<G4int> GetAngleNoVector() const { return fAngleNoVector;}
  std::vector<G4int> GetOmegaNoVector() const { return fOmegaNoVector;}

  G4DataVector GetEnergyUnitVector() const { return fEnergyUnitVector;}
  G4DataVector GetAngleUnitVector() const { return fAngleUnitVector;}
  G4DataVector GetXscUnitVector() const { return fXscUnitVector;}
  G4DataVector GetXscPerAngleUnitVector() const { return fXscPerAngleUnitVector;}
  G4DataVector GetXscPerMomCUnitVector() const { return fXscPerMomCUnitVector;}
  G4DataVector GetDdXscUnitVector() const { return fDdXscUnitVector;}

  G4DataVector*  GetTkinVector() const {return fTkinVector;}
  G4DataVector*  GetTkinBinVector() const {return fTkinBinVector;}

  G4DataVector*  GetXscVector() const {return fXscVector;}
  G4DataVector*  GetDeltaXscVector() const {return fDeltaXscVector;}

  G4DataVector*  GetMultiplicityVector() const {return fMultiplicityVector;}

  G4DataVector*  GetMomentumVector() const {return fMomentumCVector;}
  G4DataVector*  GetMomentumCBinVector() const {return fMomentumCBinVector;}
  G4DataVector*  GetDeltaMomCVector() const {return fDeltaMomCVector;}

  G4DataVector*  GetAngleVector() const {return fAngleVector;}
  G4DataVector*  GetAngleBinVector() const {return fAngleBinVector;}

  std::vector<G4DataVector*>* GetAngleDdTable() const {return fAngleDdTable;} 

  G4PhysicsTable*  GetAngleTable() const {return fAngleTable;}
  G4PhysicsTable*  GetMomentumCTable() const {return fMomentumCTable;}

  std::vector<G4PhysicsTable*>* GetDoubleDiffXscBank() const 
  {return fDoubleDiffXscBank;}

  std::vector<G4PhysicsTable*>* GetDoubleDiffXscErrorBank() const 
  {return fDoubleDiffXscErrorBank;};


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
  G4double fXscPerAngleUnit, fXscPerMomCUnit, fDdXscUnit;
 
  G4DataVector fEnergyUnitVector, fAngleUnitVector, fXscUnitVector;
  G4DataVector fXscPerAngleUnitVector, fXscPerMomCUnitVector, fDdXscUnitVector;
 

  G4double fTkinBin;
  G4int    fNo;

  static const G4String fAnyNumber;
  static const G4String fAnyEmptySpace;
  static const G4String fAnyHidden;

  std::vector<G4int> fEnergyNoVector;
  std::vector<G4int> fAngleNoVector;
  std::vector<G4int> fOmegaNoVector;

  G4DataVector*  fTkinVector;
  G4DataVector*  fTkinBinVector;

  G4DataVector*  fXscVector;
  G4DataVector*  fDeltaXscVector;
  G4DataVector*  fMultiplicityVector;

  G4DataVector*  fMomentumCVector;
  G4DataVector*  fDeltaMomCVector;
  G4DataVector*  fMomentumCBinVector;

  G4DataVector*  fAngleVector;
  G4DataVector*  fAngleBinVector;

  std::vector<G4DataVector*>* fAngleDdTable; 
  G4PhysicsTable*  fAngleTable;
  G4PhysicsTable*  fMomentumCTable;

  std::vector<G4PhysicsTable*>* fDoubleDiffXscBank;
  std::vector<G4PhysicsTable*>* fDoubleDiffXscErrorBank;
  std::vector<G4String>* fCommentVector;

};

#endif





