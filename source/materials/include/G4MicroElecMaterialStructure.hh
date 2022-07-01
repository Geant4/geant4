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
// G4MicroElecMaterialStructure.hh, 2011/08/29 A.Valentin, M. Raine are with CEA [a]
//                   	    	2020/05/20 P. Caron, C. Inguimbert are with ONERA [b] 
//				       	   Q. Gibaru is with CEA [a], ONERA [b] and CNES [c]
//				           M. Raine and D. Lambert are with CEA [a]
//
// A part of this work has been funded by the French space agency(CNES[c])
// [a] CEA, DAM, DIF - 91297 ARPAJON, France
// [b] ONERA - DPHY, 2 avenue E.Belin, 31055 Toulouse, France
// [c] CNES, 18 av.E.Belin, 31401 Toulouse CEDEX, France
//
// Based on the following publications
//	- A.Valentin, M. Raine, 
//		Inelastic cross-sections of low energy electrons in silicon
//	      for the simulation of heavy ion tracks with the Geant4-DNA toolkit,
//	      NSS Conf. Record 2010, pp. 80-85
//             https://doi.org/10.1109/NSSMIC.2010.5873720
//
//      - A.Valentin, M. Raine, M.Gaillardin, P.Paillet
//	      Geant4 physics processes for microdosimetry simulation:
//	      very low energy electromagnetic models for electrons in Silicon,
//             https://doi.org/10.1016/j.nimb.2012.06.007
//	      NIM B, vol. 288, pp. 66-73, 2012, part A
//	      heavy ions in Si, NIM B, vol. 287, pp. 124-129, 2012, part B
//             https://doi.org/10.1016/j.nimb.2012.07.028
//
//	- M. Raine, M. Gaillardin, P. Paillet
//	      Geant4 physics processes for silicon microdosimetry simulation: 
//	      Improvements and extension of the energy-range validity up to 10 GeV/nucleon
//	      NIM B, vol. 325, pp. 97-100, 2014
//             https://doi.org/10.1016/j.nimb.2014.01.014
//
//      - J. Pierron, C. Inguimbert, M. Belhaj, T. Gineste, J. Puech, M. Raine
//	      Electron emission yield for low energy electrons: 
//	      Monte Carlo simulation and experimental comparison for Al, Ag, and Si
//	      Journal of Applied Physics 121 (2017) 215107. 
//               https://doi.org/10.1063/1.4984761
//
//      - P. Caron,
//	      Study of Electron-Induced Single-Event Upset in Integrated Memory Devices
//	      PHD, 16th October 2019
//
//	- Q.Gibaru, C.Inguimbert, P.Caron, M.Raine, D.Lambert, J.Puech, 
//	      Geant4 physics processes for microdosimetry and secondary electron emission simulation : 
//	      Extension of MicroElec to very low energies and new materials
//	      NIM B, 2020, in review.
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef G4MICROELECMATERIALSTRUCTURE_HH 
#define G4MICROELECMATERIALSTRUCTURE_HH 1 

#include "globals.hh" 
#include "G4Material.hh"
#include <vector>

class G4MicroElecMaterialStructure
{
public:
  G4MicroElecMaterialStructure(const G4String& matName = "");
  virtual ~G4MicroElecMaterialStructure() = default;
	
  void ReadMaterialFile();
  G4double Energy(G4int level);
  G4int NumberOfLevels() { return nLevels; }
  G4double GetZ(G4int Shell);
  G4double ConvertUnit(const G4String& unitName);
  G4double GetEnergyGap() { return energyGap; }
  G4double GetInitialEnergy() { return initialEnergy; }
  G4int GetEADL_Enumerator(G4int shell) { return EADL_Enumerator[shell]; };
  G4double GetWorkFunction() { return workFunction; };
  G4String GetMaterialName() { return materialName; };
  G4double GetLimitEnergy(G4int level);
  G4double GetElasticModelLowLimit() {return limitElastic[0];}
  G4double GetElasticModelHighLimit() { return limitElastic[1]; }
  G4double GetInelasticModelLowLimit(G4int pdg);
  G4double GetInelasticModelHighLimit(G4int pdg);
  G4bool IsShellWeaklyBound(G4int level);
  
private:
  // private elements     
  G4int nLevels = 3;   	// Number of levels of material
  G4bool isCompound = false;
  G4String materialName = "";
  std::vector<G4bool> isShellWeaklyBoundVector;
  std::vector<G4double> energyConstant;
  std::vector<G4double> LimitEnergy;
  std::vector<G4int> EADL_Enumerator;
  G4double workFunction = 0.0;
  G4double initialEnergy = 0.0;
  std::vector<G4double> compoundShellZ;
  G4double Z = 0.0;
  G4double energyGap = 0.0;
  G4double limitElastic[2] = { 0,0 };
  G4double limitInelastic[4] = { 0,0,0,0 };
};

#endif 
