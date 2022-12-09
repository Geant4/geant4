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
// G4MicroElecInelasticModel_new.cc, 2011/08/29 A.Valentin, M. Raine are with CEA [a]
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


#include "globals.hh"
#include "G4MicroElecInelasticModel_new.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4ionEffectiveCharge.hh"
#include "G4MicroElecMaterialStructure.hh"
#include "G4DeltaAngle.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MicroElecInelasticModel_new::G4MicroElecInelasticModel_new(
       const G4ParticleDefinition*, const G4String& nam)
  :G4VEmModel(nam),isInitialised(false)
{
  
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods
  
  if( verboseLevel>0 )
  {
    G4cout << "MicroElec inelastic model is constructed " << G4endl;
  }
  
  //Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);
  fAtomDeexcitation = nullptr;
  fParticleChangeForGamma = nullptr;
  
  // default generator
  SetAngularDistribution(new G4DeltaAngle());

  // Selection of computation method
  fasterCode = true;
  SEFromFermiLevel = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MicroElecInelasticModel_new::~G4MicroElecInelasticModel_new()
{
  // Cross section  
  // (0)
  TCSMap::iterator pos2;
  for (pos2 = tableTCS.begin(); pos2 != tableTCS.end(); ++pos2) {
    MapData* tableData = pos2->second;
    for (auto pos = tableData->begin(); pos != tableData->end(); ++pos)
      {
	G4MicroElecCrossSectionDataSet_new* table = pos->second;
	delete table;
      }
    delete tableData;
  }
  tableTCS.clear();
    
  dataDiffCSMap::iterator iterator_proba;  
  // (1)
  for (iterator_proba = eNrjTransStorage.begin(); iterator_proba != eNrjTransStorage.end(); ++iterator_proba) {
    vector<TriDimensionMap>* eNrjTransfData = iterator_proba->second;
    eNrjTransfData->clear();
    delete eNrjTransfData;
  }
  eNrjTransStorage.clear();
  
  for (iterator_proba = pNrjTransStorage.begin(); iterator_proba != pNrjTransStorage.end(); ++iterator_proba) {
    vector<TriDimensionMap>* pNrjTransfData = iterator_proba->second;
    pNrjTransfData->clear();
    delete pNrjTransfData;
  }
  pNrjTransStorage.clear();
  
  // (2)
  for (iterator_proba = eDiffDatatable.begin(); iterator_proba != eDiffDatatable.end(); ++iterator_proba) {
    vector<TriDimensionMap>* eDiffCrossSectionData = iterator_proba->second;
    eDiffCrossSectionData->clear();
    delete eDiffCrossSectionData;
  }
  eDiffDatatable.clear();

  for (iterator_proba = pDiffDatatable.begin(); iterator_proba != pDiffDatatable.end(); ++iterator_proba) {
    vector<TriDimensionMap>* pDiffCrossSectionData = iterator_proba->second;
    pDiffCrossSectionData->clear();
    delete pDiffCrossSectionData;
  }
  pDiffDatatable.clear();
  
  // (3)
  dataProbaShellMap::iterator iterator_probaShell;
  
  for (iterator_probaShell = eProbaShellStorage.begin(); iterator_probaShell != eProbaShellStorage.end(); ++iterator_probaShell) {
    vector<VecMap>* eProbaShellMap = iterator_probaShell->second;
    eProbaShellMap->clear();
    delete eProbaShellMap;
  }
  eProbaShellStorage.clear();
  
  for (iterator_probaShell = pProbaShellStorage.begin(); iterator_probaShell != pProbaShellStorage.end(); ++iterator_probaShell) {
    vector<VecMap>* pProbaShellMap = iterator_probaShell->second;
    pProbaShellMap->clear();
    delete pProbaShellMap;
  }
  pProbaShellStorage.clear();
  
  // (4)
  TranfEnergyMap::iterator iterator_nrjtransf;
  for (iterator_nrjtransf = eVecmStorage.begin(); iterator_nrjtransf != eVecmStorage.end(); ++iterator_nrjtransf) {
    VecMap* eVecm = iterator_nrjtransf->second;
    eVecm->clear();
    delete eVecm;
  }
  eVecmStorage.clear();
  for (iterator_nrjtransf = pVecmStorage.begin(); iterator_nrjtransf != pVecmStorage.end(); ++iterator_nrjtransf) {
    VecMap* pVecm = iterator_nrjtransf->second;
    pVecm->clear();
    delete pVecm;
  }
  pVecmStorage.clear();
  
  // (5)
  incidentEnergyMap::iterator iterator_energy;
  for (iterator_energy = eIncidentEnergyStorage.begin(); iterator_energy != eIncidentEnergyStorage.end(); ++iterator_energy) {
    std::vector<G4double>* eTdummyVec = iterator_energy->second;
    eTdummyVec->clear();
    delete eTdummyVec;
  }
  eIncidentEnergyStorage.clear();
  
  for (iterator_energy = pIncidentEnergyStorage.begin(); iterator_energy != pIncidentEnergyStorage.end(); ++iterator_energy) {
    std::vector<G4double>* pTdummyVec = iterator_energy->second;
    pTdummyVec->clear();
    delete pTdummyVec;
  }
  pIncidentEnergyStorage.clear();
  
  // (6)
  MapStructure::iterator iterator_matStructure;
  for (iterator_matStructure = tableMaterialsStructures.begin(); 
       iterator_matStructure != tableMaterialsStructures.end(); ++iterator_matStructure) {
    currentMaterialStructure = iterator_matStructure->second;
    delete currentMaterialStructure;
  }
  tableMaterialsStructures.clear();
  currentMaterialStructure = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MicroElecInelasticModel_new::Initialise(const G4ParticleDefinition* particle,
				    const G4DataVector& /*cuts*/)
{
  if (isInitialised) { return; }

  if (verboseLevel > 3)
    G4cout << "Calling G4MicroElecInelasticModel_new::Initialise()" << G4endl;
  
  const char* path = G4FindDataDir("G4LEDATA");
  if (!path)
    {
      G4Exception("G4MicroElecElasticModel_new::Initialise","em0006",FatalException,"G4LEDATA environment variable not set.");
      return;
    }
  
  G4String modelName = "mermin";
  G4cout << "****************************" << G4endl;
  G4cout << modelName << " model loaded !" << G4endl;
  
  // Energy limits 
  G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();
  G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();
  G4String electron = electronDef->GetParticleName();
  G4String proton = protonDef->GetParticleName();
  
  G4double scaleFactor = 1.0;
  
    // *** ELECTRON 
  lowEnergyLimit[electron] = 2 * eV;
  highEnergyLimit[electron] = 10.0 * MeV;
  
  // Cross section
  G4ProductionCutsTable* theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();
  
  for (G4int i = 0; i < numOfCouples; ++i) {
    const G4Material* material = theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();
    G4cout << "Material " << i + 1 << " / " << numOfCouples << " : " << material->GetName() << G4endl;
    if (material->GetName() == "Vacuum") continue;
    G4String mat = material->GetName().substr(3, material->GetName().size());	
    MapData* tableData = new MapData;
    currentMaterialStructure = new G4MicroElecMaterialStructure(mat);
    
    tableMaterialsStructures[mat] = currentMaterialStructure;
    if (particle == electronDef) {
      //TCS
      G4String fileElectron("Inelastic/" + modelName + "_sigma_inelastic_e-_" + mat);
      G4cout << fileElectron << G4endl;
      G4MicroElecCrossSectionDataSet_new* tableE = new G4MicroElecCrossSectionDataSet_new(new G4LogLogInterpolation, MeV, scaleFactor);
      tableE->LoadData(fileElectron);
      tableData->insert(make_pair(electron, tableE));
      
      // DCS
      std::ostringstream eFullFileName;
      if (fasterCode) {
	eFullFileName << path << "/microelec/Inelastic/cumulated_" + modelName + "_sigmadiff_inelastic_e-_" + mat + ".dat";
	G4cout << "Faster code = true" << G4endl;
	G4cout << "Inelastic/cumulated_" + modelName + "_sigmadiff_inelastic_e-_" + mat + ".dat" << G4endl;
      }
      else {
	eFullFileName << path << "/microelec/Inelastic/" + modelName + "_sigmadiff_inelastic_e-_" + mat + ".dat";
	G4cout << "Faster code = false" << G4endl;
	G4cout << "Inelastic/" + modelName + "_sigmadiff_inelastic_e-_" + mat + ".dat" << G4endl;
      }
      
      std::ifstream eDiffCrossSection(eFullFileName.str().c_str());
      if (!eDiffCrossSection)
	{
	  std::stringstream ss;
	  ss << "Missing data " << eFullFileName.str().c_str();
	  std::string sortieString = ss.str();
	 	  
	  if (fasterCode) G4Exception("G4MicroElecInelasticModel_new::Initialise", "em0003",
				      FatalException, sortieString.c_str());
				      				      
	  else {
	    G4Exception("G4MicroElecInelasticModel_new::Initialise", "em0003",
			FatalException, "Missing data file:/microelec/sigmadiff_inelastic_e_Si.dat");
	  }	  
	}
      
      // Clear the arrays for re-initialization case (MT mode)
      // Octobre 22nd, 2014 - Melanie Raine      
      //Creating vectors of maps for DCS and Cumulated DCS for the current material.
      //Each vector is storing one map for each shell.      
      vector<TriDimensionMap>* eDiffCrossSectionData = 
	new vector<TriDimensionMap>; //Storage of [IncidentEnergy, TransfEnergy, DCS values], used in slower code
      vector<TriDimensionMap>* eNrjTransfData = 
	new vector<TriDimensionMap>; //Storage of possible transfer energies by shell
      vector<VecMap>* eProbaShellMap = new vector<VecMap>; //Storage of the vectors containing all cumulated DCS values for an initial energy, by shell
      vector<G4double>* eTdummyVec = new vector<G4double>; //Storage of incident energies for interpolation
      VecMap* eVecm = new VecMap; //Transfered energy map for slower code
      
      for (int j = 0; j < currentMaterialStructure->NumberOfLevels(); j++) //Filling the map vectors with an empty map for each shell
	{
	  eDiffCrossSectionData->push_back(TriDimensionMap());
	  eNrjTransfData->push_back(TriDimensionMap());
	  eProbaShellMap->push_back(VecMap());
	}
      
      eTdummyVec->push_back(0.);
      while (!eDiffCrossSection.eof())
	{
	  G4double tDummy; //incident energy
	  G4double eDummy; //transfered energy
	  eDiffCrossSection >> tDummy >> eDummy;
	  if (tDummy != eTdummyVec->back()) eTdummyVec->push_back(tDummy);
	  
	  G4double tmp; //probability
	  for (int j = 0; j < currentMaterialStructure->NumberOfLevels(); j++)
	    {
	      eDiffCrossSection >> tmp;	      
	      (*eDiffCrossSectionData)[j][tDummy][eDummy] = tmp;
	      
	      if (fasterCode)
		{
		  (*eNrjTransfData)[j][tDummy][(*eDiffCrossSectionData)[j][tDummy][eDummy]] = eDummy;
		  (*eProbaShellMap)[j][tDummy].push_back((*eDiffCrossSectionData)[j][tDummy][eDummy]);
		}
	      else {  // SI - only if eof is not reached !
		if (!eDiffCrossSection.eof()) (*eDiffCrossSectionData)[j][tDummy][eDummy] *= scaleFactor;
		(*eVecm)[tDummy].push_back(eDummy);
	      }
	    }
	}
      //
      G4cout << "add to material vector" << G4endl;
      
      //Filing maps for the current material into the master maps
      if (fasterCode) {
	eNrjTransStorage[mat] = eNrjTransfData;
	eProbaShellStorage[mat] = eProbaShellMap;
      }
      else {
	eDiffDatatable[mat] = eDiffCrossSectionData;
	eVecmStorage[mat] = eVecm;
      }
      eIncidentEnergyStorage[mat] = eTdummyVec;

      //Cleanup support vectors
      // delete eProbaShellMap;
      // delete eDiffCrossSectionData;
      // delete eNrjTransfData;
    }
    
    // *** PROTON
    if (particle == protonDef)
      {
	// Cross section
	G4String fileProton("Inelastic/" + modelName + "_sigma_inelastic_p_" + mat); G4cout << fileProton << G4endl;
	G4MicroElecCrossSectionDataSet_new* tableP = new G4MicroElecCrossSectionDataSet_new(new G4LogLogInterpolation, MeV, scaleFactor);
	tableP->LoadData(fileProton);
	tableData->insert(make_pair(proton, tableP));

	// DCS
	std::ostringstream pFullFileName;	
	if (fasterCode) {
	  pFullFileName << path << "/microelec/Inelastic/cumulated_" + modelName + "_sigmadiff_inelastic_p_" + mat + ".dat";
	  G4cout << "Faster code = true" << G4endl;
	  G4cout << "Inelastic/cumulated_" + modelName + "_sigmadiff_inelastic_p_" + mat + ".dat" << G4endl;
	}
	else {
	  pFullFileName << path << "/microelec/Inelastic/" + modelName + "_sigmadiff_inelastic_p_" + mat + ".dat";
	  G4cout << "Faster code = false" << G4endl;
	  G4cout << "Inelastic/" + modelName + "_sigmadiff_inelastic_e-_" + mat + ".dat" << G4endl;
	}
	
	std::ifstream pDiffCrossSection(pFullFileName.str().c_str());
	if (!pDiffCrossSection)
	  {
	    if (fasterCode) G4Exception("G4MicroElecInelasticModel_new::Initialise", "em0003",
					FatalException, "Missing data file:/microelec/sigmadiff_cumulated_inelastic_p_Si.dat");
	    else {	      
	      G4Exception("G4MicroElecInelasticModel_new::Initialise", "em0003",
			  FatalException, "Missing data file:/microelec/sigmadiff_inelastic_p_Si.dat");
	    }
	  }
	
	//
	// Clear the arrays for re-initialization case (MT mode)
	// Octobre 22nd, 2014 - Melanie Raine
	//Creating vectors of maps for DCS and Cumulated DCS for the current material.
	//Each vector is storing one map for each shell.
	
	vector<TriDimensionMap>* pDiffCrossSectionData = 
	  new vector<TriDimensionMap>; //Storage of [IncidentEnergy, TransfEnergy, DCS values], used in slower code
	vector<TriDimensionMap>* pNrjTransfData = 
	  new vector<TriDimensionMap>; //Storage of possible transfer energies by shell
	vector<VecMap>* pProbaShellMap = 
	  new vector<VecMap>; //Storage of the vectors containing all cumulated DCS values for an initial energy, by shell
	vector<G4double>* pTdummyVec = 
	  new vector<G4double>; //Storage of incident energies for interpolation
	VecMap* eVecm = new VecMap; //Transfered energy map for slower code

	for (int j = 0; j < currentMaterialStructure->NumberOfLevels(); ++j) 
	  //Filling the map vectors with an empty map for each shell
	  {
	    pDiffCrossSectionData->push_back(TriDimensionMap());
	    pNrjTransfData->push_back(TriDimensionMap());
	    pProbaShellMap->push_back(VecMap());
	  }
	
	pTdummyVec->push_back(0.);
	while (!pDiffCrossSection.eof())
	  {
	    G4double tDummy; //incident energy
	    G4double eDummy; //transfered energy
	    pDiffCrossSection >> tDummy >> eDummy;
	    if (tDummy != pTdummyVec->back()) pTdummyVec->push_back(tDummy);
	    
	    G4double tmp; //probability
	    for (int j = 0; j < currentMaterialStructure->NumberOfLevels(); j++)
	      {
		pDiffCrossSection >> tmp;	
		(*pDiffCrossSectionData)[j][tDummy][eDummy] = tmp; 
		// ArrayofMaps[j] -> fill with 3DMap(incidentEnergy, 
		// 2Dmap (transferedEnergy,proba=tmp) ) -> fill map for shell j 
		// with proba for transfered energy eDummy
		
		if (fasterCode)
		  {
		    (*pNrjTransfData)[j][tDummy][(*pDiffCrossSectionData)[j][tDummy][eDummy]] = eDummy;
		    (*pProbaShellMap)[j][tDummy].push_back((*pDiffCrossSectionData)[j][tDummy][eDummy]);
		  }
		else {  // SI - only if eof is not reached !
		  if (!pDiffCrossSection.eof()) (*pDiffCrossSectionData)[j][tDummy][eDummy] *= scaleFactor;
		  (*eVecm)[tDummy].push_back(eDummy);
		}
	      }
	  }
	
	//Filing maps for the current material into the master maps
	if (fasterCode) {
	  pNrjTransStorage[mat] = pNrjTransfData;
	  pProbaShellStorage[mat] = pProbaShellMap;
	}
	else {
	  pDiffDatatable[mat] = pDiffCrossSectionData;
	  pVecmStorage[mat] = eVecm;
	}
	pIncidentEnergyStorage[mat] = pTdummyVec;
         
	//Cleanup support vectors
	// delete pNrjTransfData;
	// delete eVecm;
	// delete pDiffCrossSectionData;
	// delete pProbaShellMap;
      }
    tableTCS[mat] = tableData;
} 
  if (particle==electronDef)
    {
      SetLowEnergyLimit(lowEnergyLimit[electron]);
      SetHighEnergyLimit(highEnergyLimit[electron]);
    }
  
  if (particle==protonDef)
    {
      SetLowEnergyLimit(100*eV);
      SetHighEnergyLimit(10*MeV);
    }
  
  if( verboseLevel>1 )
    {
      G4cout << "MicroElec Inelastic model is initialized " << G4endl
	     << "Energy range: "
	     << LowEnergyLimit() / keV << " keV - "
	     << HighEnergyLimit() / MeV << " MeV for "
	     << particle->GetParticleName()
	     << " with mass (amu) " << particle->GetPDGMass()/proton_mass_c2
	     << " and charge " << particle->GetPDGCharge()
	     << G4endl << G4endl ;
    }
  
  fAtomDeexcitation  = G4LossTableManager::Instance()->AtomDeexcitation();
  
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecInelasticModel_new::CrossSectionPerVolume(const G4Material* material,
						   const G4ParticleDefinition* particleDefinition,
						   G4double ekin,
						   G4double,
						   G4double)
{
  if (verboseLevel > 3) G4cout << "Calling CrossSectionPerVolume() of G4MicroElecInelasticModel" << G4endl;
  
  G4double density = material->GetTotNbOfAtomsPerVolume();
  currentMaterial = material->GetName().substr(3, material->GetName().size());
  
  MapStructure::iterator structPos;
  structPos = tableMaterialsStructures.find(currentMaterial);
  
  // Calculate total cross section for model
  TCSMap::iterator tablepos;
  tablepos = tableTCS.find(currentMaterial);
  
  if (tablepos == tableTCS.end() )
    {
      G4String str = "Material ";
      str += currentMaterial + " TCS Table not found!";
      G4Exception("G4MicroElecInelasticModel_new::ComputeCrossSectionPerVolume", "em0002", FatalException, str);
      return 0;
    }
  else if(structPos == tableMaterialsStructures.end())
    {
      G4String str = "Material ";
      str += currentMaterial + " Structure not found!";
      G4Exception("G4MicroElecInelasticModel_new::ComputeCrossSectionPerVolume", "em0002", FatalException, str);
      return 0;
    }
  else {
    MapData* tableData = tablepos->second;
    currentMaterialStructure = structPos->second;
    
    G4double sigma = 0;
   
    const G4String& particleName = particleDefinition->GetParticleName();
    G4String nameLocal = particleName;
    G4int pdg = particleDefinition->GetPDGEncoding();
    G4int Z = particleDefinition->GetAtomicNumber();
    
    G4double Zeff = 1.0, Zeff2 = Zeff*Zeff;
    G4double Mion_c2 = particleDefinition->GetPDGMass();
    
    if (Mion_c2 > proton_mass_c2)
      {
	ekin *= proton_mass_c2 / Mion_c2;
        nameLocal = "proton";
      }

    G4double lowLim = currentMaterialStructure->GetInelasticModelLowLimit(pdg);
    G4double highLim = currentMaterialStructure->GetInelasticModelHighLimit(pdg);
    
    if (ekin >= lowLim && ekin < highLim)
      {
	std::map< G4String, G4MicroElecCrossSectionDataSet_new*, std::less<G4String> >::iterator pos;
	pos = tableData->find(nameLocal); //find particle type
	
	if (pos != tableData->end())
	  {
	    G4MicroElecCrossSectionDataSet_new* table = pos->second;
	    if (table != 0)
	      {
		sigma = table->FindValue(ekin);
				
		if (Mion_c2 > proton_mass_c2) {
		  sigma = 0.;
		  for (G4int i = 0; i < currentMaterialStructure->NumberOfLevels(); i++) {
		    Zeff = BKZ(ekin / (proton_mass_c2 / Mion_c2), Mion_c2 / c_squared, Z, currentMaterialStructure->Energy(i)); // il faut garder le vrai ekin car le calcul à l'interieur de la methode convertie l'énergie en vitesse
		    Zeff2 = Zeff*Zeff;
		    sigma += Zeff2*table->FindShellValue(ekin, i); 
// il faut utiliser le ekin mis à l'echelle pour chercher la bonne 
// valeur dans les tables proton
		
		  }
		}
		else {
		  sigma = table->FindValue(ekin);
		}
	      }
	  }
	else
	  {
	    G4Exception("G4MicroElecInelasticModel_new::CrossSectionPerVolume", 
                        "em0002", FatalException, 
                        "Model not applicable to particle type.");
	  }
      }
    else
      {
	return 1 / DBL_MAX;
      }
    
    if (verboseLevel > 3)
      {
	G4cout << "---> Kinetic energy (eV)=" << ekin / eV << G4endl;
	G4cout << " - Cross section per Si atom (cm^2)=" << sigma / cm2 << G4endl;
	G4cout << " - Cross section per Si atom (cm^-1)=" << sigma*density / (1. / cm) << G4endl;
      }
        
    return (sigma)*density;} 
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MicroElecInelasticModel_new::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                  const G4MaterialCutsCouple* couple,
                                                  const G4DynamicParticle* particle,
                                                  G4double,
                                                  G4double)
{
  
  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4MicroElecInelasticModel" << G4endl;
  
  G4int pdg = particle->GetParticleDefinition()->GetPDGEncoding();
  G4double lowLim = currentMaterialStructure->GetInelasticModelLowLimit(pdg);
  G4double highLim = currentMaterialStructure->GetInelasticModelHighLimit(pdg);
  
  G4double ekin = particle->GetKineticEnergy();
  G4double k = ekin ;
  
  G4ParticleDefinition* PartDef = particle->GetDefinition();
  const G4String& particleName = PartDef->GetParticleName();
  G4String nameLocal2 = particleName ;
  G4double particleMass = particle->GetDefinition()->GetPDGMass();
  G4double originalMass = particleMass; // a passer en argument dans samplesecondaryenergy pour évaluer correctement Qmax
  G4int originalZ = particle->GetDefinition()->GetAtomicNumber();
  
  if (particleMass > proton_mass_c2)
    {
      k *= proton_mass_c2/particleMass ;
      PartDef = G4Proton::ProtonDefinition();
      nameLocal2 = "proton" ;
    }
   
  if (k >= lowLim && k < highLim)
    {
      G4ParticleMomentum primaryDirection = particle->GetMomentumDirection();
      G4double totalEnergy = ekin + particleMass;
      G4double pSquare = ekin * (totalEnergy + particleMass);
      G4double totalMomentum = std::sqrt(pSquare);
      
      G4int Shell = 1;
      
      Shell = RandomSelect(k,nameLocal2,originalMass, originalZ);
            
      G4double bindingEnergy = currentMaterialStructure->Energy(Shell);
      G4double limitEnergy = currentMaterialStructure->GetLimitEnergy(Shell);
      
      if (verboseLevel > 3)
	{
	  G4cout << "---> Kinetic energy (eV)=" << k/eV << G4endl ;
	  G4cout << "Shell: " << Shell << ", energy: " << bindingEnergy/eV << G4endl;
	}
      
      // sample deexcitation
      
      std::size_t secNumberInit = 0;  // need to know at a certain point the energy of secondaries
      std::size_t secNumberFinal = 0; // So I'll make the difference and then sum the energies
      
      //SI: additional protection if tcs interpolation method is modified
      //if (k<bindingEnergy) return;
      if (k<limitEnergy) return;
      //	G4cout << currentMaterial << G4endl;
      G4int Z = currentMaterialStructure->GetZ(Shell);     
      G4int shellEnum = currentMaterialStructure->GetEADL_Enumerator(Shell);
      if (currentMaterialStructure->IsShellWeaklyBound(Shell)) { shellEnum = -1; }
      
      if(fAtomDeexcitation && shellEnum >=0) 
	{
	  //		G4cout << "enter if deex and shell 0" << G4endl;
	  G4AtomicShellEnumerator as = G4AtomicShellEnumerator(shellEnum);
	  const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);
	  secNumberInit = fvect->size();
	  fAtomDeexcitation->GenerateParticles(fvect, shell, Z, 0, 0);
	  secNumberFinal = fvect->size();
	}
            
      G4double secondaryKinetic=-1000*eV;
      SEFromFermiLevel = false;
      if (!fasterCode)
	{
	  secondaryKinetic = RandomizeEjectedElectronEnergy(PartDef, k, Shell, originalMass, originalZ);	  
	}
      else 
	{
	  secondaryKinetic = RandomizeEjectedElectronEnergyFromCumulatedDcs(PartDef, k, Shell) ;
	}

      if (verboseLevel > 3)
	{
	  G4cout << "Ionisation process" << G4endl;
	  G4cout << "Shell: " << Shell << " Kin. energy (eV)=" << k/eV
		 << " Sec. energy (eV)=" << secondaryKinetic/eV << G4endl;
	}
      G4ThreeVector deltaDirection =
	GetAngularDistribution()->SampleDirectionForShell(particle, secondaryKinetic,
							  Z, Shell,
							  couple->GetMaterial());
      
      if (particle->GetDefinition() == G4Electron::ElectronDefinition())
	{
	  G4double deltaTotalMomentum = std::sqrt(secondaryKinetic*(secondaryKinetic + 2.*electron_mass_c2 ));	  
	  
	  G4double finalPx = totalMomentum*primaryDirection.x() - deltaTotalMomentum*deltaDirection.x();
	  G4double finalPy = totalMomentum*primaryDirection.y() - deltaTotalMomentum*deltaDirection.y();
	  G4double finalPz = totalMomentum*primaryDirection.z() - deltaTotalMomentum*deltaDirection.z();
	  G4double finalMomentum = std::sqrt(finalPx*finalPx + finalPy*finalPy + finalPz*finalPz);
	  finalPx /= finalMomentum;
	  finalPy /= finalMomentum;
	  finalPz /= finalMomentum;
	  
	  G4ThreeVector direction;
	  direction.set(finalPx,finalPy,finalPz);
	  
	  fParticleChangeForGamma->ProposeMomentumDirection(direction.unit());
	}
      else fParticleChangeForGamma->ProposeMomentumDirection(primaryDirection);
      
      // note that secondaryKinetic is the energy of the delta ray, not of all secondaries.
      G4double deexSecEnergy = 0;
      for (std::size_t j=secNumberInit; j < secNumberFinal; ++j) {
        deexSecEnergy = deexSecEnergy + (*fvect)[j]->GetKineticEnergy();
      }      
      if (SEFromFermiLevel) limitEnergy = currentMaterialStructure->GetEnergyGap();
      fParticleChangeForGamma->SetProposedKineticEnergy(ekin - secondaryKinetic - limitEnergy); //Ef = Ei-(Q-El)-El = Ei-Q
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(limitEnergy - deexSecEnergy);
           
      if (secondaryKinetic>0)
	{  
	  G4DynamicParticle* dp = new G4DynamicParticle(G4Electron::Electron(), deltaDirection, secondaryKinetic); //Esec = Q-El
	  fvect->push_back(dp);
	}      
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MicroElecInelasticModel_new::RandomizeEjectedElectronEnergy(
         const G4ParticleDefinition* particleDefinition,
	 G4double k, G4int shell, G4double originalMass, G4int)
{
  G4double secondaryElectronKineticEnergy=0.;
  if (particleDefinition == G4Electron::ElectronDefinition())
    {
      G4double maximumEnergyTransfer=k;      
      G4double crossSectionMaximum = 0.;
      G4double minEnergy = currentMaterialStructure->GetLimitEnergy(shell);
      G4double maxEnergy = maximumEnergyTransfer;
      G4int nEnergySteps = 100;
      
      G4double value(minEnergy);
      G4double stpEnergy(std::pow(maxEnergy/value, 1./static_cast<G4double>(nEnergySteps-1)));
      G4int step(nEnergySteps);
      while (step>0)
	{
	  --step;
	  G4double differentialCrossSection = 
	    DifferentialCrossSection(particleDefinition, k, value, shell);
	  crossSectionMaximum = std::max(crossSectionMaximum, differentialCrossSection);
	  value*=stpEnergy;
	}
      
      do
	{
	  secondaryElectronKineticEnergy = G4UniformRand() * 
	    (maximumEnergyTransfer-currentMaterialStructure->GetLimitEnergy(shell));
	} while(G4UniformRand()*crossSectionMaximum > 
		DifferentialCrossSection(particleDefinition, k,
		  (secondaryElectronKineticEnergy+currentMaterialStructure->GetLimitEnergy(shell)),shell));
    }
  else if (particleDefinition == G4Proton::ProtonDefinition())
    {      
      G4double maximumEnergyTransfer =
	ComputeElasticQmax(k/(proton_mass_c2/originalMass), 
			   currentMaterialStructure->Energy(shell), 
			   originalMass/c_squared, electron_mass_c2/c_squared);
      
      G4double crossSectionMaximum = 0.;
      
      G4double minEnergy = currentMaterialStructure->GetLimitEnergy(shell);
      G4double maxEnergy = maximumEnergyTransfer;
      G4int nEnergySteps = 100;
      
      G4double value(minEnergy);
      G4double stpEnergy(std::pow(maxEnergy/value, 1./static_cast<G4double>(nEnergySteps-1)));
      G4int step(nEnergySteps);
      
      while (step>0)
	{
	  --step;
	  G4double differentialCrossSection = 
	    DifferentialCrossSection(particleDefinition, k, value, shell);
	  crossSectionMaximum = std::max(crossSectionMaximum, differentialCrossSection);
	  value*=stpEnergy;
	}
      
      G4double energyTransfer = 0.;      
      do
	{
	  energyTransfer = G4UniformRand() * maximumEnergyTransfer;
	} while(G4UniformRand()*crossSectionMaximum > 
		DifferentialCrossSection(particleDefinition, k,energyTransfer,shell));
      
      secondaryElectronKineticEnergy = 
	energyTransfer-currentMaterialStructure->GetLimitEnergy(shell);
      
    }
  return std::max(secondaryElectronKineticEnergy, 0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MicroElecInelasticModel_new::RandomizeEjectedElectronEnergyFromCumulatedDcs(
         const G4ParticleDefinition* particleDefinition, G4double k, G4int shell)
{
  G4double secondaryElectronKineticEnergy = 0.;
  G4double random = G4UniformRand();

  secondaryElectronKineticEnergy = TransferedEnergy(particleDefinition, k, shell, random)
    - currentMaterialStructure->GetLimitEnergy(shell) ;

  if (isnan(secondaryElectronKineticEnergy)) { secondaryElectronKineticEnergy = k - currentMaterialStructure->GetLimitEnergy(shell); }

  if (secondaryElectronKineticEnergy < 0.) {
    secondaryElectronKineticEnergy = k - currentMaterialStructure->GetEnergyGap();
    SEFromFermiLevel = true;
  }
  return secondaryElectronKineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MicroElecInelasticModel_new::TransferedEnergy(
	 const G4ParticleDefinition* particleDefinition,
	 G4double k,
	 G4int ionizationLevelIndex,
	 G4double random)
{
  G4double nrj = 0.;
  G4double valueK1 = 0;
  G4double valueK2 = 0;
  G4double valuePROB21 = 0;
  G4double valuePROB22 = 0;
  G4double valuePROB12 = 0;
  G4double valuePROB11 = 0;
  G4double nrjTransf11 = 0;
  G4double nrjTransf12 = 0;
  G4double nrjTransf21 = 0;
  G4double nrjTransf22 = 0;
  
  G4double maximumEnergyTransfer1 = 0;
  G4double maximumEnergyTransfer2 = 0;
  G4double maximumEnergyTransferP = 4.* (electron_mass_c2 / proton_mass_c2) * k;
  G4double bindingEnergy = currentMaterialStructure->GetLimitEnergy(ionizationLevelIndex);
  
  if (particleDefinition == G4Electron::ElectronDefinition())
    {      
      dataDiffCSMap::iterator iterator_Nrj;
      iterator_Nrj = eNrjTransStorage.find(currentMaterial);
      
      dataProbaShellMap::iterator iterator_Proba;
      iterator_Proba = eProbaShellStorage.find(currentMaterial);
      
      incidentEnergyMap::iterator iterator_Tdummy;
      iterator_Tdummy = eIncidentEnergyStorage.find(currentMaterial);
      
      if(iterator_Nrj == eNrjTransStorage.end() || iterator_Proba == eProbaShellStorage.end() || 
	 iterator_Tdummy == eIncidentEnergyStorage.end())
	{
	  G4String str = "Material ";
	  str += currentMaterial + " not found!";
	  G4Exception("G4MicroElecInelasticModel_new::TransferedEnergy", "em0002", 
		      FatalException, str);
	}
      else {
	vector<TriDimensionMap>* eNrjTransfData = iterator_Nrj->second; //Storage of possible transfer energies
	vector<VecMap>* eProbaShellMap = iterator_Proba->second; //Storage of probabilities for energy transfer
	vector<G4double>* eTdummyVec = iterator_Tdummy->second; //Incident energies for interpolation
	
	// k should be in eV
	auto k2 = std::upper_bound(eTdummyVec->begin(),
							    eTdummyVec->end(),
							    k);
	auto k1 = k2 - 1;

	// SI : the following condition avoids situations where random >last vector element
	if (random <= (*eProbaShellMap)[ionizationLevelIndex][(*k1)].back()
	    && random <= (*eProbaShellMap)[ionizationLevelIndex][(*k2)].back())
	  {
	    auto prob12 =
	      std::upper_bound((*eProbaShellMap)[ionizationLevelIndex][(*k1)].begin(),
			       (*eProbaShellMap)[ionizationLevelIndex][(*k1)].end(),
			       random);
	    
	    auto prob11 = prob12 - 1;
	    
	    auto prob22 =
	      std::upper_bound((*eProbaShellMap)[ionizationLevelIndex][(*k2)].begin(),
			       (*eProbaShellMap)[ionizationLevelIndex][(*k2)].end(),
			       random);
	    
	    auto prob21 = prob22 - 1;
	    
	    valueK1 = *k1;
	    valueK2 = *k2;
	    valuePROB21 = *prob21;
	    valuePROB22 = *prob22;
	    valuePROB12 = *prob12;
	    valuePROB11 = *prob11;
	    
	    // The following condition avoid getting transfered energy < binding energy and forces cumxs = 1 for maximum energy transfer.
	    if (valuePROB11 == 0) nrjTransf11 = bindingEnergy;
	    else nrjTransf11 = (*eNrjTransfData)[ionizationLevelIndex][valueK1][valuePROB11];
	    if (valuePROB12 == 1)
	      {
		if ((valueK1 + bindingEnergy) / 2. > valueK1) 
		  maximumEnergyTransfer1 = valueK1;
		else 
		  maximumEnergyTransfer1 = (valueK1 + bindingEnergy) / 2.;
		
		nrjTransf12 = maximumEnergyTransfer1;
	      }
	    else 
	      nrjTransf12 = (*eNrjTransfData)[ionizationLevelIndex][valueK1][valuePROB12];
	    
	    if (valuePROB21 == 0) nrjTransf21 = bindingEnergy;
	    else nrjTransf21 = (*eNrjTransfData)[ionizationLevelIndex][valueK2][valuePROB21];
	    if (valuePROB22 == 1)
	      {
		if ((valueK2 + bindingEnergy) / 2. > valueK2) maximumEnergyTransfer2 = valueK2;
		else maximumEnergyTransfer2 = (valueK2 + bindingEnergy) / 2.;
				
		nrjTransf22 = maximumEnergyTransfer2;
	      }
	    else nrjTransf22 = (*eNrjTransfData)[ionizationLevelIndex][valueK2][valuePROB22];
	    
	  }
	// Avoids cases where cum xs is zero for k1 and is not for k2 (with always k1<k2)
	if (random > (*eProbaShellMap)[ionizationLevelIndex][(*k1)].back())
	  {
	    auto prob22 =
	      std::upper_bound((*eProbaShellMap)[ionizationLevelIndex][(*k2)].begin(),
			       (*eProbaShellMap)[ionizationLevelIndex][(*k2)].end(),
			       random);	    
	    auto prob21 = prob22 - 1;
	    
	    valueK1 = *k1;
	    valueK2 = *k2;
	    valuePROB21 = *prob21;
	    valuePROB22 = *prob22;
	    	    
	    nrjTransf21 = (*eNrjTransfData)[ionizationLevelIndex][valueK2][valuePROB21];
	    nrjTransf22 = (*eNrjTransfData)[ionizationLevelIndex][valueK2][valuePROB22];
	    
	    G4double interpolatedvalue2 = Interpolate(valuePROB21,
						      valuePROB22,
						      random,
						      nrjTransf21,
						      nrjTransf22);
	    
	    // zeros are explicitly set
	    G4double value = Interpolate(valueK1, valueK2, k, 0., interpolatedvalue2);
	    
	    return value;
	  }
      }
    }
  else if (particleDefinition == G4Proton::ProtonDefinition())
    {
      // k should be in eV      
      dataDiffCSMap::iterator iterator_Nrj;
      iterator_Nrj = pNrjTransStorage.find(currentMaterial);
      
      dataProbaShellMap::iterator iterator_Proba;
      iterator_Proba = pProbaShellStorage.find(currentMaterial);
      
      incidentEnergyMap::iterator iterator_Tdummy;
      iterator_Tdummy = pIncidentEnergyStorage.find(currentMaterial);
      
      if (iterator_Nrj == pNrjTransStorage.end() || iterator_Proba == pProbaShellStorage.end() || 
	  iterator_Tdummy == pIncidentEnergyStorage.end())
	{
	  G4String str = "Material ";
	  str += currentMaterial + " not found!";
	  G4Exception("G4MicroElecInelasticModel_new::TransferedEnergy", "em0002", 
		      FatalException, str);
	}
      else 
	{
	  vector<TriDimensionMap>* pNrjTransfData = iterator_Nrj->second; //Storage of possible transfer energies
	  vector<VecMap>* pProbaShellMap = iterator_Proba->second; //Storage of probabilities for energy transfer
	  vector<G4double>* pTdummyVec = iterator_Tdummy->second; //Incident energies for interpolation
	  
	  auto k2 = std::upper_bound(pTdummyVec->begin(),
			 				      pTdummyVec->end(),
							      k);
	  
	  auto k1 = k2 - 1;
	  
	  // SI : the following condition avoids situations where random > last vector element,
	  //      for eg. when the last element is zero
	  if (random <= (*pProbaShellMap)[ionizationLevelIndex][(*k1)].back()
	      && random <= (*pProbaShellMap)[ionizationLevelIndex][(*k2)].back())
	    {
	      auto prob12 =
		std::upper_bound((*pProbaShellMap)[ionizationLevelIndex][(*k1)].begin(),
				 (*pProbaShellMap)[ionizationLevelIndex][(*k1)].end(),
				 random);
	      auto prob11 = prob12 - 1;	      
	      auto prob22 =
		std::upper_bound((*pProbaShellMap)[ionizationLevelIndex][(*k2)].begin(),
				 (*pProbaShellMap)[ionizationLevelIndex][(*k2)].end(),
				 random);
	      auto prob21 = prob22 - 1;
	      
	      valueK1 = *k1;
	      valueK2 = *k2;
	      valuePROB21 = *prob21;
	      valuePROB22 = *prob22;
	      valuePROB12 = *prob12;
	      valuePROB11 = *prob11;
	      
	      // The following condition avoid getting transfered energy < binding energy 
	      // and forces cumxs = 1 for maximum energy transfer.
	      if (valuePROB11 == 0) nrjTransf11 = bindingEnergy;
	      else nrjTransf11 = (*pNrjTransfData)[ionizationLevelIndex][valueK1][valuePROB11];
	      
	      if (valuePROB12 == 1) nrjTransf12 = maximumEnergyTransferP;
	      else nrjTransf12 = (*pNrjTransfData)[ionizationLevelIndex][valueK1][valuePROB12];
	      
	      if (valuePROB21 == 0) nrjTransf21 = bindingEnergy;
	      else nrjTransf21 = (*pNrjTransfData)[ionizationLevelIndex][valueK2][valuePROB21];
	      
	      if (valuePROB22 == 1) nrjTransf22 = maximumEnergyTransferP;
	      else nrjTransf22 = (*pNrjTransfData)[ionizationLevelIndex][valueK2][valuePROB22];
	      
	    }
	  
	  // Avoids cases where cum xs is zero for k1 and is not for k2 (with always k1<k2)	  
	  if (random > (*pProbaShellMap)[ionizationLevelIndex][(*k1)].back())
	    {
	      auto prob22 =
		std::upper_bound((*pProbaShellMap)[ionizationLevelIndex][(*k2)].begin(),
				 (*pProbaShellMap)[ionizationLevelIndex][(*k2)].end(),
				 random);
	      
	      auto prob21 = prob22 - 1;
	      
	      valueK1 = *k1;
	      valueK2 = *k2;
	      valuePROB21 = *prob21;
	      valuePROB22 = *prob22;
	      
	      nrjTransf21 = (*pNrjTransfData)[ionizationLevelIndex][valueK2][valuePROB21];
	      nrjTransf22 = (*pNrjTransfData)[ionizationLevelIndex][valueK2][valuePROB22];
	      
	      G4double interpolatedvalue2 = Interpolate(valuePROB21,
							valuePROB22,
							random,
							nrjTransf21,
							nrjTransf22);
	      
	      // zeros are explicitly set	      
	      G4double value = Interpolate(valueK1, valueK2, k, 0., interpolatedvalue2);	      
	      return value;
	    }
	}
    }
  // End electron and proton cases
  
  G4double nrjTransfProduct = nrjTransf11 * nrjTransf12 * nrjTransf21 * nrjTransf22;
  
  if (nrjTransfProduct != 0.)
    {
      nrj = QuadInterpolator(valuePROB11,
			     valuePROB12,
			     valuePROB21,
			     valuePROB22,
			     nrjTransf11,
			     nrjTransf12,
			     nrjTransf21,
			     nrjTransf22,
			     valueK1,
			     valueK2,
			     k,
			     random);
    }
  
  return nrj;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MicroElecInelasticModel_new::DifferentialCrossSection(
         const G4ParticleDefinition * particleDefinition,
	 G4double k,
	 G4double energyTransfer,
	 G4int LevelIndex)
{
  G4double sigma = 0.;
  
  if (energyTransfer >= currentMaterialStructure->GetLimitEnergy(LevelIndex))
    {
      G4double valueT1 = 0;
      G4double valueT2 = 0;
      G4double valueE21 = 0;
      G4double valueE22 = 0;
      G4double valueE12 = 0;
      G4double valueE11 = 0;
      
      G4double xs11 = 0;
      G4double xs12 = 0;
      G4double xs21 = 0;
      G4double xs22 = 0;
      
      if (particleDefinition == G4Electron::ElectronDefinition())
	{
	  
	  dataDiffCSMap::iterator iterator_Proba;
	  iterator_Proba = eDiffDatatable.find(currentMaterial);
	  
	  incidentEnergyMap::iterator iterator_Nrj;
	  iterator_Nrj = eIncidentEnergyStorage.find(currentMaterial);
	  
	  TranfEnergyMap::iterator iterator_TransfNrj;
	  iterator_TransfNrj = eVecmStorage.find(currentMaterial);
	  
	  if (iterator_Proba != eDiffDatatable.end() && iterator_Nrj != eIncidentEnergyStorage.end() 
	      && iterator_TransfNrj!= eVecmStorage.end())
	    {	
	      vector<TriDimensionMap>* eDiffCrossSectionData = (iterator_Proba->second);
	      vector<G4double>* eTdummyVec = iterator_Nrj->second; //Incident energies for interpolation
	      VecMap* eVecm = iterator_TransfNrj->second;
	      
	      // k should be in eV and energy transfer eV also
	      auto t2 = std::upper_bound(eTdummyVec->begin(), eTdummyVec->end(), k);
	      auto t1 = t2 - 1;
	      // SI : the following condition avoids situations where energyTransfer >last vector element
	      if (energyTransfer <= (*eVecm)[(*t1)].back() && energyTransfer <= (*eVecm)[(*t2)].back())
		{
		  auto e12 = std::upper_bound((*eVecm)[(*t1)].begin(), (*eVecm)[(*t1)].end(), energyTransfer);
		  auto e11 = e12 - 1;		  
		  auto e22 = std::upper_bound((*eVecm)[(*t2)].begin(), (*eVecm)[(*t2)].end(), energyTransfer);
		  auto e21 = e22 - 1;
		  
		  valueT1 = *t1;
		  valueT2 = *t2;
		  valueE21 = *e21;
		  valueE22 = *e22;
		  valueE12 = *e12;
		  valueE11 = *e11;

		  xs11 = (*eDiffCrossSectionData)[LevelIndex][valueT1][valueE11];
		  xs12 = (*eDiffCrossSectionData)[LevelIndex][valueT1][valueE12];
		  xs21 = (*eDiffCrossSectionData)[LevelIndex][valueT2][valueE21];
		  xs22 = (*eDiffCrossSectionData)[LevelIndex][valueT2][valueE22];
		}
	    }
	  else {
	    G4String str = "Material ";
	    str += currentMaterial + " not found!";
	    G4Exception("G4MicroElecDielectricModels::DifferentialCrossSection", "em0002", FatalException, str);
	  }
	}
      
      if (particleDefinition == G4Proton::ProtonDefinition())
	{	  
	  dataDiffCSMap::iterator iterator_Proba;
	  iterator_Proba = pDiffDatatable.find(currentMaterial);
	  
	  incidentEnergyMap::iterator iterator_Nrj;
	  iterator_Nrj = pIncidentEnergyStorage.find(currentMaterial);
	  
	  TranfEnergyMap::iterator iterator_TransfNrj;
	  iterator_TransfNrj = pVecmStorage.find(currentMaterial);
	  
	  if (iterator_Proba != pDiffDatatable.end() && iterator_Nrj != pIncidentEnergyStorage.end() 
	      && iterator_TransfNrj != pVecmStorage.end())
	    {
	      vector<TriDimensionMap>* pDiffCrossSectionData = (iterator_Proba->second);
	      vector<G4double>* pTdummyVec = iterator_Nrj->second; //Incident energies for interpolation
	      VecMap* pVecm = iterator_TransfNrj->second;
	      	      
	      // k should be in eV and energy transfer eV also
	      auto t2 = 
		std::upper_bound(pTdummyVec->begin(), pTdummyVec->end(), k);
	      auto t1 = t2 - 1;
	      if (energyTransfer <= (*pVecm)[(*t1)].back() && energyTransfer <= (*pVecm)[(*t2)].back())
		{
		  auto e12 = std::upper_bound((*pVecm)[(*t1)].begin(), (*pVecm)[(*t1)].end(), energyTransfer);
		  auto e11 = e12 - 1;		  
		  auto e22 = std::upper_bound((*pVecm)[(*t2)].begin(), (*pVecm)[(*t2)].end(), energyTransfer);
		  auto e21 = e22 - 1;
		  
		  valueT1 = *t1;
		  valueT2 = *t2;
		  valueE21 = *e21;
		  valueE22 = *e22;
		  valueE12 = *e12;
		  valueE11 = *e11;
		  
		  xs11 = (*pDiffCrossSectionData)[LevelIndex][valueT1][valueE11];
		  xs12 = (*pDiffCrossSectionData)[LevelIndex][valueT1][valueE12];
		  xs21 = (*pDiffCrossSectionData)[LevelIndex][valueT2][valueE21];
		  xs22 = (*pDiffCrossSectionData)[LevelIndex][valueT2][valueE22];
		}
	    }
	  else {
	    G4String str = "Material ";
	    str += currentMaterial + " not found!";
	    G4Exception("G4MicroElecDielectricModels::DifferentialCrossSection", "em0002", FatalException, str);
	  }
	}
      
      G4double xsProduct = xs11 * xs12 * xs21 * xs22;
      if (xsProduct != 0.)
	{
	  sigma = QuadInterpolator(     valueE11, valueE12,
					valueE21, valueE22,
					xs11, xs12,
					xs21, xs22,
					valueT1, valueT2,
					k, energyTransfer);
	}
      
    }
  
  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double G4MicroElecInelasticModel_new::Interpolate(G4double e1,
					 G4double e2,
					 G4double e,
					 G4double xs1,
					 G4double xs2)
{
  G4double value = 0.;

  // Log-log interpolation by default
  if (e1 != 0 && e2 != 0 && (e2-e1) != 0 && !fasterCode)
    {
      G4double a = std::log(xs2/xs1)/ std::log(e2/e1);
      G4double b = std::log(xs2) - a * std::log(e2);
      G4double sigma = a * std::log(e) + b;
      value = (std::exp(sigma));
    }

  // Switch to log-lin interpolation for faster code
  if ((e2 - e1) != 0 && xs1 != 0 && xs2 != 0 && fasterCode)
    {
      G4double d1 = std::log(xs1);
      G4double d2 = std::log(xs2);
      value = std::exp((d1 + (d2 - d1) * (e - e1) / (e2 - e1)));
    }

  // Switch to lin-lin interpolation for faster code
  // in case one of xs1 or xs2 (=cum proba) value is zero
  if ((e2 - e1) != 0 && (xs1 == 0 || xs2 == 0) && fasterCode)
    {
      G4double d1 = xs1;
      G4double d2 = xs2;
      value = (d1 + (d2 - d1) * (e - e1) / (e2 - e1));
    }

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MicroElecInelasticModel_new::QuadInterpolator(G4double e11, G4double e12,
					      G4double e21, G4double e22,
					      G4double xs11, G4double xs12,
					      G4double xs21, G4double xs22,
					      G4double t1, G4double t2,
					      G4double t, G4double e)
{
  G4double interpolatedvalue1 = Interpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = Interpolate(e21, e22, e, xs21, xs22);
  G4double value = Interpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4MicroElecInelasticModel_new::RandomSelect(G4double k, const G4String& particle, G4double originalMass, G4int originalZ )
{
  G4int level = 0;
  
  TCSMap::iterator tablepos;
  tablepos = tableTCS.find(currentMaterial);
  MapData* tableData = tablepos->second;
  
  std::map< G4String,G4MicroElecCrossSectionDataSet_new*,std::less<G4String> >::iterator pos;
  pos = tableData->find(particle);
  
  std::vector<G4double> Zeff(currentMaterialStructure->NumberOfLevels(), 1.0);
  if(originalMass>proton_mass_c2) {
    for(G4int nl=0;nl<currentMaterialStructure->NumberOfLevels();nl++) {
      Zeff[nl] = BKZ(k/(proton_mass_c2/originalMass), originalMass/c_squared, originalZ, currentMaterialStructure->Energy(nl));
    }
  }
  
  if (pos != tableData->end())
    {
      G4MicroElecCrossSectionDataSet_new* table = pos->second;
      
      if (table != 0)
	{
	  G4double* valuesBuffer = new G4double[table->NumberOfComponents()];
	  const G4int n = (G4int)table->NumberOfComponents();
	  G4int i = (G4int)n;
	  G4double value = 0.;
	  
	  while (i>0)
	    {
	      --i;
	      valuesBuffer[i] = table->GetComponent(i)->FindValue(k)*Zeff[i]*Zeff[i];
	      value += valuesBuffer[i];
	    }	  
	  value *= G4UniformRand();
	  
	  i = n;
	  
	  while (i > 0)
	    {
	      --i;
	      
	      if (valuesBuffer[i] > value)
		{
		  delete[] valuesBuffer;
		  return i;
		}
	      value -= valuesBuffer[i];
	    }
	  
	  if (valuesBuffer) delete[] valuesBuffer;
	  
	}
    }
  else
    {
      G4Exception("G4MicroElecInelasticModel_new::RandomSelect","em0002",FatalException,"Model not applicable to particle type.");
    }
  
  return level;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecInelasticModel_new::ComputeRelativistVelocity(G4double E, G4double mass) {
  G4double x = E/mass;
  return c_light*std::sqrt(x*(x + 2.0))/(x + 1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecInelasticModel_new::ComputeElasticQmax(G4double T1i, G4double T2i, G4double M1, G4double M2) {
  G4double v1i = ComputeRelativistVelocity(T1i, M1);
  G4double v2i = ComputeRelativistVelocity(T2i, M2);
  
  G4double v2f = 2*M1/(M1+M2)*v1i + (M2-M1)/(M1+M2)*-1*v2i;
  G4double vtransfer2a = v2f*v2f-v2i*v2i;

  v2f = 2*M1/(M1+M2)*v1i + (M2-M1)/(M1+M2)*v2i;
  G4double vtransfer2b = v2f*v2f-v2i*v2i;

  G4double vtransfer2 = std::max(vtransfer2a, vtransfer2b);
  return 0.5*M2*vtransfer2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecInelasticModel_new::stepFunc(G4double x) {
  return (x < 0.) ? 1.0 : 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecInelasticModel_new::vrkreussler(G4double v, G4double vF) 
{
  G4double r = vF*( std::pow(v/vF+1., 3.) - fabs(std::pow(v/vF-1., 3.)) 
		    + 4.*(v/vF)*(v/vF) ) + stepFunc(v/vF-1.) * (3./2.*v/vF - 
								4.*(v/vF)*(v/vF) + 3.*std::pow(v/vF, 3.) 
								- 0.5*std::pow(v/vF, 5.));
  return r/(10.*v/vF);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecInelasticModel_new::BKZ(G4double Ep, G4double mp, G4int Zp, G4double Eplasmon) 
{
  // need atomic unit conversion
  G4double hbar = hbar_Planck, hbar2 = hbar*hbar, me = electron_mass_c2/c_squared, Ry = me*elm_coupling*elm_coupling/(2*hbar2);
  G4double hartree = 2*Ry,	a0 = Bohr_radius,	velocity = a0*hartree/hbar;
  G4double vp = ComputeRelativistVelocity(Ep,  mp);

  vp /= velocity;

  G4double wp = Eplasmon/hartree;
  G4double a = std::pow(4./9./CLHEP::pi, 1./3.);
  G4double vF = std::pow(wp*wp/(3.*a*a*a), 1./3.);
  G4double c = 0.9;
  G4double vr = vrkreussler(vp /*in u.a*/, vF /*in u.a*/);
  G4double yr = vr/std::pow(Zp, 2./3.);
  G4double q = 0.;
  if(Zp==2) q = 1-exp(-c*vr/(Zp-5./16.));
  else q = 1.-exp(-c*(yr-0.07));
  G4double Neq = Zp*(1.-q);
  G4double l0 = 0.;
  if(Neq<=2) l0 = 3./(Zp-0.3*(Neq-1))/2.;
  else l0 = 0.48*std::pow(Neq, 2./3.)/(Zp-Neq/7.);
  if(Zp==2)  c = 1.0;
  else c = 3./2.;
  return Zp*(q + c*(1.-q)/vF/vF/2.0 * log(1.+std::pow(2.*l0*vF,2.)));	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
