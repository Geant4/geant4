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
// $Id: G4NistMaterialBuilder.hh,v 1.3 2005/05/12 17:29:08 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $

#ifndef G4NistMaterialBuilder_h
#define G4NistMaterialBuilder_h 1

//---------------------------------------------------------------------------
//
// ClassName:   G4NistMaterialBuilder
//
// Description: Utility class to hold and manipulate G4Materials
//
// Author:      V.Ivanchenko 21-11-2004
//
// Modifications:
//
//
//----------------------------------------------------------------------------
//
// Class Description:
//
// Element data from the NIST DB on Atomic Weights and Isotope Compositions
// http://physics.nist.gov/PhysRefData/Compositions/index.html
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "G4Material.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4NistManager;
class G4NistElementBuilder;

class G4NistMaterialBuilder
{
public:

  G4NistMaterialBuilder(G4NistManager* mm=0, G4NistElementBuilder* eb=0,
                        G4int verb=0);
			
 ~G4NistMaterialBuilder();
 
  // Find or build a G4Material by name, from dataBase
  //
  G4Material* FindOrBuildMaterial (const G4String& name, G4bool isotopes=true);
					    
  // construct a G4Material from scratch by atome count
  // 
  G4Material* ConstructNewMaterial (const G4String& name,
                                      const std::vector<G4String>& elm,
                                      const std::vector<G4int>& nbAtoms,
				      G4double dens, G4bool isotopes=true);
				      
  // construct a G4Material from scratch by fraction mass
  //
  G4Material* ConstructNewMaterial (const G4String& name,
                                      const std::vector<G4String>& elm,
                                      const std::vector<G4double>& weight,
				      G4double dens, G4bool isotopes=true);

  void SetVerbose(G4int val);
  void ListMaterials(const G4String&);
  void ListNistSimpleMaterials();
  void ListNistCompoundMaterials();
  void ListHepMaterials();

private:

  void Initialise();
  void NistSimpleMaterials();
  void NistCompoundMaterials();
  void HepAndNuclearMaterials();

  void AddMaterial(const G4String& nameMat, G4double dens, G4int Z=0,
                         G4double pot=0.0, G4int ncomp=1,
                         G4State=kStateSolid);

  void AddElementByWeightFraction(G4int Z, G4double);
  void AddElementByAtomCount     (G4int Z, G4int);

  void AddElementByWeightFraction(const G4String& name, G4double);
  void AddElementByAtomCount     (const G4String& name, G4int);

  // build a G4Material from dataBase
  G4Material* BuildMaterial(const G4String& name, G4bool isotopes);
  
  void DumpElm(G4int);
  void DumpMix(G4int);
  
private:

  G4NistManager*         matManager;
  G4NistElementBuilder*  elmBuilder;

  G4int                  verbose;
  G4int                  nMaterials;
  G4int                  nComponents;
  G4int                  nCurrent;
  G4int                  nElementary;
  G4int                  nNIST;

  std::vector<G4String>  names;
  std::vector<G4String>  chFormulas;

  std::vector<G4double>  densities;
  std::vector<G4double>  ionPotentials;
  std::vector<G4State>   states;
  std::vector<G4double>  fractions;
  std::vector<G4int>     components;
  std::vector<G4int>     indexes;
  std::vector<G4int>     elements;

  G4bool                 first;

};

#endif
