//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4SandiaTable.hh 96794 2016-05-09 10:09:30Z gcosmo $

// class description
//
// This class is an interface to G4StaticSandiaData.
// it provides - Sandia coeff for an element, given its Z
//             - sandia coeff for a material, given a pointer to it
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....
//
// History:
//
// 10.06.97 created. V. Grichine
// 18.11.98 simplified public interface; new methods for materials.  mma
// 30.01.01 major bug in the computation of AoverAvo and in the units (/g!)
//          in GetSandiaCofPerAtom(). mma
// 03.04.01 fnulcof[4] added; returned if energy < emin
// 05.03.04 V.Grichine, new methods for old sorting algorithm for PAI model
// 21.21.13 V.Ivanchenko, changed signature of methods, reduced number of 
//                        static variables, methods
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

#ifndef G4SANDIATABLE_HH
#define G4SANDIATABLE_HH

#include "G4OrderedTable.hh"      
#include "G4ios.hh"
#include "globals.hh"
#include <assert.h>
#include <vector>

#include <CLHEP/Units/PhysicalConstants.h>

class G4Material;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

class G4SandiaTable
{
public:  // with description

  G4SandiaTable(G4Material*);	         
	         
  ~G4SandiaTable();

  //main computation per atom:
  void GetSandiaCofPerAtom(G4int Z, G4double energy, 
			   std::vector<G4double>& coeff) const;

  void GetSandiaCofWater(G4double energy, 
			 std::vector<G4double>& coeff) const;

  G4double GetWaterEnergyLimit() const;
  G4double GetWaterCofForMaterial(G4int,G4int) const;

  static G4double GetZtoA(G4int Z);

  //per volume of a material:
  G4int GetMatNbOfIntervals() const;
  G4double  GetSandiaCofForMaterial(G4int,G4int) const;
  G4double  GetSandiaMatTable(G4int,G4int) const;
  const G4double* GetSandiaCofForMaterial(G4double energy) const;

  G4double  GetSandiaMatTablePAI(G4int,G4int) const;
  const G4double* GetSandiaCofForMaterialPAI(G4double energy) const;

  inline void SetVerbose(G4int ver) { fVerbose = ver; };

public:  // without description

  G4SandiaTable(__void__&);
  // Fake default constructor for usage restricted to direct object
  // persistency for clients requiring preallocation of memory for
  // persistifiable objects.

private:

  void ComputeMatSandiaMatrix();
  void ComputeMatSandiaMatrixPAI();

  // methods per atom
  G4double  GetSandiaPerAtom(G4int Z, G4int, G4int) const;

#ifdef G4VERBOSE
  static G4int PrintErrorZ(G4int Z, const G4String&);
  static void PrintErrorV(const G4String&);
#endif

  // computed once
  static G4int    fCumulInterval[101];
  static const G4double funitc[5];

  // used at initialisation
  std::vector<G4double>   fSandiaCofPerAtom;
  
  // members of the class
  G4Material*     fMaterial;
  G4int           fMatNbOfIntervals;
  G4OrderedTable* fMatSandiaMatrix;
  G4OrderedTable* fMatSandiaMatrixPAI;
  		   
/////////////////////////////////////////////////////////////////////
//
// Methods for implementation of PAI model
//
/////////////////////////////////////////////////////////////////////

public:  // without description

  G4SandiaTable(G4int matIndex);
        
  G4SandiaTable();

  void Initialize(G4Material*);	         

  G4int SandiaIntervals(G4int Z[], G4int el);

  G4int SandiaMixing(G4int Z[], const G4double* fractionW,
		     G4int el, G4int mi);

  G4double GetPhotoAbsorpCof(G4int i , G4int j) const;

  G4int GetMaxInterval() const; 

  inline G4bool GetLowerI1() {return fLowerI1;};
  inline void   SetLowerI1(G4bool flag) {fLowerI1=flag;};

private:

  void ComputeMatTable();

  void SandiaSwap(G4double** da, G4int i, G4int j);

  void SandiaSort(G4double** da, G4int sz);

  G4double** GetPointerToCof(); 

  // copy constructor and hide assignment operator
  G4SandiaTable(G4SandiaTable &) = delete;
  G4SandiaTable & operator=(const G4SandiaTable &right) = delete;

  static const G4double fSandiaTable[981][5];
  static const G4int fNumberOfElements;
  static const G4int fIntervalLimit;
  static const G4int fNumberOfIntervals;
  static const G4int fH2OlowerInt;

  // data members for PAI model
  G4double** fPhotoAbsorptionCof;  // SandiaTable  for mixture

  G4int fMaxInterval ;
  G4int fVerbose;
  G4bool fLowerI1;  
};
    
#endif 
