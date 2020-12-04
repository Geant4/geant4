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
// G4MicroElecInelastic_new.hh, 2011/08/29 A.Valentin, M. Raine are with CEA [a]
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

#ifndef G4MICROELEINELASTIC_HH
#define G4MICROELEINELASTIC_HH 1

#include "G4VEmProcess.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MicroElecInelastic_new : public G4VEmProcess
{
public: 
  G4MicroElecInelastic_new(const G4String& processName ="MicroElecIonisation",
			   G4ProcessType type = fElectromagnetic);
  ~G4MicroElecInelastic_new() override;
  
  G4bool IsApplicable(const G4ParticleDefinition&) override;
  void PrintInfo() override;
  
protected:
  void InitialiseProcess(const G4ParticleDefinition*) override;
  
private:
  G4bool isInitialised = false;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
#endif
