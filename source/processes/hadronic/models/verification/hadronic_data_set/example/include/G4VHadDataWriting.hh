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
//
// R&D: Vladimir.Grichine@cern.ch

#ifndef G4VHadDataWriting_HH
#define G4VHadDataWriting_HH 1


#include "globals.hh"
#include <vector>
#include <map>

#include "G4HadFileSpec.hh"



class G4PhysicsTable;
class G4PhysicsVector;
class G4DataVector;

class G4Isotope;
class G4Element;
class G4Material;

class G4VHadDataWriting 
{ 
public:

  G4VHadDataWriting();

  virtual ~G4VHadDataWriting();

  G4double GetTkinBin() const {return fTkinBin;}
  G4int    GetNumber() const {return fNo;}

  G4DataVector*  GetTkinVector() const {return fTkinVector;}
  G4DataVector*  GetXscVector() const {return fXscVector;}
  G4DataVector*  GetMultiplicityVector() const {return fMultiplicityVector;}

  G4DataVector*  GetMomentumVector() const {return fMomentumCVector;}
  G4DataVector*  GetAngleVector() const {return fAngleVector;}

  G4PhysicsTable*  GetAngleTable() const {return fAngleTable;}
  G4PhysicsTable*  GetMomentumCTable() const {return fMomentumCTable;}

  std::vector<G4PhysicsTable*>* GetDoubleDiffXscBank() const 
  {return fDoubleDiffXscBank;}

  std::vector<G4String>* GetCommentVector() const {return fCommentVector;}
 
protected:

  virtual void WriteDataFile(G4HadFileSpec&) = 0;

  void FillIntegralXSC(G4String&);
  void FillMultiplicity(G4String&);

  void FillDifferentialXSC(G4String&,G4bool);
  virtual  void FillDoubleDiffXSC(G4String&);

  void SetEnergyUnit( G4String&);
  void SetAngleUnit( G4String&);
  void SetXscUnit( G4String&);
  void SetDdXscUnit( G4String&);

  G4String fInputFileName;
  void SetInputFileName(G4String fn) {fInputFileName=fn;}

protected:

  G4String fEnergyUnit, fAngleUnit, fXscUnit;
  G4String fXscPeAngleUnit, fXscPerMomCUnit, fDdXscUnit;

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

  G4PhysicsTable*  fAngleTable;
  G4PhysicsTable*  fMomentumCTable;

  std::vector<G4PhysicsTable*>* fDoubleDiffXscBank;
  std::vector<G4String>* fCommentVector;
};
 
#endif
 
