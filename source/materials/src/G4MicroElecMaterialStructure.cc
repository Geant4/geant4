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
// G4MicroElecMaterialStructure.cc, 2011/08/29 A.Valentin, M. Raine are with CEA [a]
//                   	    	2020/05/20 P. Caron, C. Inguimbert are with ONERA [b] 
//				       	   Q. Gibaru is with CEA [a], ONERA [b] and CNES [c]
//				            M. Raine and D. Lambert are with CEA [a]
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

#include "G4MicroElecMaterialStructure.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MicroElecMaterialStructure::G4MicroElecMaterialStructure(const G4String& matName)
{
  materialName = matName;
  if (matName == "Vacuum" || matName == "uum") {
    workFunction = 0;
    initialEnergy = 0;
  }
  else {
    ReadMaterialFile();
  }
  nLevels = (G4int)energyConstant.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MicroElecMaterialStructure::ReadMaterialFile() 
{
  const char* path = G4FindDataDir("G4LEDATA");

  if (materialName[0] == 'G' && materialName[1] == '4') {
    //in the case the NIST database is used
    materialName.erase(0, 1);
    materialName.erase(0, 1);
    materialName.erase(0, 1);
  }
  
  std::ostringstream fileName;
  fileName << path << "/microelec/Structure/Data_" + materialName + ".dat";
  std::ifstream fichier(fileName.str().c_str());
  
  int varLength = 0;
  G4String nameParameter;
  
  G4String unitName;	
  G4double unitValue;
  G4double data;
  G4String filler;	
  G4String type;
  
  if (fichier)
    {
      fichier >> filler >> type;
      materialName = filler;
      if (type == "Compound") {isCompound = true; Z = 0; }
      else { isCompound = false; Z = std::stoi(type); }
      while(!fichier.eof()) {
	
	getline(fichier, filler);
	std::stringstream line(filler);
	
	if (filler[0] == '#' || filler.empty()) {continue;}
	
	line >> varLength;
	line >> nameParameter;
	line >> unitName;
	unitValue = ConvertUnit(unitName);
	
	for (int i = 0; i < varLength; i++)
	  {
	    line >> data;	data = data*unitValue;

      if(nameParameter == "WorkFunction")
      {
        workFunction = data;
      }
      if(nameParameter == "EnergyGap")
      {
        energyGap = data;
      }

      if(nameParameter == "EnergyPeak")
      {
        energyConstant.push_back(data);
      }
      if(nameParameter == "EnergyLimit")
      {
        LimitEnergy.push_back(data);
      }
      if(nameParameter == "EADL")
      {
        EADL_Enumerator.push_back(data);
      }

      if (nameParameter == "WeaklyBoundShell")
	      {if (data == 0) { isShellWeaklyBoundVector.push_back(false); }
		else {isShellWeaklyBoundVector.push_back(true);}}

        if(nameParameter == "WeaklyBoundInitialEnergy")
        {
          initialEnergy = data;
        }

        if(nameParameter == "ShellAtomicNumber")
        {
          compoundShellZ.push_back(data);
        }

        if(nameParameter == "DielectricModelLowEnergyLimit_e")
        {
          limitInelastic[0] = data;
        }
        if(nameParameter == "DielectricModelHighEnergyLimit_e")
        {
          limitInelastic[1] = data;
        }
        if(nameParameter == "DielectricModelLowEnergyLimit_p")
        {
          limitInelastic[2] = data;
        }
        if(nameParameter == "DielectricModelHighEnergyLimit_p")
        {
          limitInelastic[3] = data;
        }

        if(nameParameter == "ElasticModelLowEnergyLimit")
        {
          limitElastic[0] = data;
        }
        if(nameParameter == "ElasticModelHighEnergyLimit")
        {
          limitElastic[1] = data;
        }
    }
      }
      fichier.close();  // on ferme le fichier
    }
  else {
    G4String str = "file ";
    str += fileName.str() + " not found!";
    G4Exception("G4MicroElecMaterialStructure::ReadMaterialFile", "em0002", FatalException, str);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecMaterialStructure::Energy(G4int level)
{
  return (level >= 0 && level < nLevels) ? energyConstant[level] : 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecMaterialStructure::GetZ(G4int Shell)
{
  if (Shell >= 0 && Shell < nLevels) {
    if(!isCompound)
    {
      return Z;
    }
    else
    {
      return compoundShellZ[Shell];
    }
  }
  else
  {
    return 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecMaterialStructure::ConvertUnit(const G4String& unitName)
{
  G4double unitValue = 0;
  if(unitName == "meV")
  {
    unitValue = 1e-3 * CLHEP::eV;
  }
  else if(unitName == "eV")
  {
    unitValue = CLHEP::eV;
  }
  else if(unitName == "keV")
  {
    unitValue = CLHEP::keV;
  }
  else if(unitName == "MeV")
  {
    unitValue = CLHEP::MeV;
  }
  else if(unitName == "noUnit")
  {
    unitValue = 1;
  }

  return unitValue;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecMaterialStructure::GetLimitEnergy(G4int level)
{
  G4double E = LimitEnergy[level];
  if (IsShellWeaklyBound(level)) { E = energyGap+ initialEnergy; }
  return E;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecMaterialStructure::GetInelasticModelLowLimit(G4int pdg)
{
  G4double res = 0.0;
  if(pdg == 11)
  {
    res = limitInelastic[0];
  }
  else if(pdg == 2212)
  {
    res = limitInelastic[2];
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecMaterialStructure::GetInelasticModelHighLimit(G4int pdg)
{
  G4double res = 0.0;
  if(pdg == 11)
  {
    res = limitInelastic[1];
  }
  else if(pdg == 2212)
  {
    res = limitInelastic[3];
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4MicroElecMaterialStructure::IsShellWeaklyBound(G4int level)
{
  return isShellWeaklyBoundVector[level];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

