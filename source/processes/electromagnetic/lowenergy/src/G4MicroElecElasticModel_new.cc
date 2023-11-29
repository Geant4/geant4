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
// G4MicroElecElasticModel_new.cc, 2011/08/29 A.Valentin, M. Raine are with CEA [a]
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
#include "G4MicroElecElasticModel_new.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exp.hh"
#include "G4Material.hh"
#include "G4String.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MicroElecElasticModel_new::G4MicroElecElasticModel_new(const G4ParticleDefinition*,
						 const G4String& nam)
  :G4VEmModel(nam), isInitialised(false)
{ 
  killBelowEnergy = 0.1*eV;     // Minimum e- energy for energy loss by excitation
  lowEnergyLimit = 0.1 * eV;
  lowEnergyLimitOfModel = 10 * eV; // The model lower energy is 10 eV
  highEnergyLimit = 500. * keV;
  SetLowEnergyLimit(lowEnergyLimit);
  SetHighEnergyLimit(highEnergyLimit);
  
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods
  
  if( verboseLevel>0 )
    {
      G4cout << "MicroElec Elastic model is constructed " << G4endl
	     << "Energy range: "
	     << lowEnergyLimit / eV << " eV - "
	     << highEnergyLimit / MeV << " MeV"
	     << G4endl;
    }
  fParticleChangeForGamma = 0;
  
  killElectron = false;
  acousticModelEnabled = false;
  currentMaterialName = "";  
  isOkToBeInitialised = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MicroElecElasticModel_new::~G4MicroElecElasticModel_new()
{
  // For total cross section
  TCSMap::iterator pos2;
  for (pos2 = tableTCS.begin(); pos2 != tableTCS.end(); ++pos2) {
    MapData* tableData = pos2->second;
    std::map< G4String, G4MicroElecCrossSectionDataSet_new*, std::less<G4String> >::iterator pos;
    for (pos = tableData->begin(); pos != tableData->end(); ++pos)
      {
	G4MicroElecCrossSectionDataSet_new* table = pos->second;
	delete table;
      }
    delete tableData;
  }
  
  //Clearing DCS maps
  
  ThetaMap::iterator iterator_angle;
  for (iterator_angle = thetaDataStorage.begin(); iterator_angle != thetaDataStorage.end(); ++iterator_angle) {
    TriDimensionMap* eDiffCrossSectionData = iterator_angle->second;
    eDiffCrossSectionData->clear();
    delete eDiffCrossSectionData;
  }
  
  energyMap::iterator iterator_energy;
  for (iterator_energy = eIncidentEnergyStorage.begin(); iterator_energy != eIncidentEnergyStorage.end(); ++iterator_energy) {
    std::vector<G4double>* eTdummyVec = iterator_energy->second;
    eTdummyVec->clear();
    delete eTdummyVec;
  }
  
  ProbaMap::iterator iterator_proba;
  for (iterator_proba = eProbaStorage.begin(); iterator_proba != eProbaStorage.end(); ++iterator_proba) {
    VecMap* eVecm = iterator_proba->second;
    eVecm->clear();
    delete eVecm;
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MicroElecElasticModel_new::Initialise(const G4ParticleDefinition* /*particle*/,
					 const G4DataVector& /*cuts*/)
{
  if (isOkToBeInitialised == true && isInitialised == false) {
      
    if (verboseLevel > -1)
      G4cout << "Calling G4MicroElecElasticModel_new::Initialise()" << G4endl;
    // Energy limits
    // Reading of data files
    
    G4double scaleFactor = 1e-18 * cm * cm;
    
    G4ProductionCutsTable* theCoupleTable =
      G4ProductionCutsTable::GetProductionCutsTable();
    G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();
    
    for (G4int i = 0; i < numOfCouples; ++i) {
      const G4Material* material =
	theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();
      
      //theCoupleTable->GetMaterialCutsCouple(i)->;
      
      G4cout << "MicroElasticModel, Material " << i + 1 << " / " << numOfCouples << " : " << material->GetName() << G4endl;
      if (material->GetName() == "Vacuum") continue;
      
      G4String matName = material->GetName().substr(3, material->GetName().size());
      G4cout<< matName<< G4endl;
      
      currentMaterialStructure = new G4MicroElecMaterialStructure(matName);
      lowEnergyLimitTable[matName]=currentMaterialStructure->GetElasticModelLowLimit();
      highEnergyLimitTable[matName]=currentMaterialStructure->GetElasticModelHighLimit();
      workFunctionTable[matName] = currentMaterialStructure->GetWorkFunction();
      
      delete currentMaterialStructure;
      
      G4cout << "Reading TCS file" << G4endl;       
      G4String fileElectron = "Elastic/elsepa_elastic_cross_e_" + matName;
      G4cout << "Elastic Total Cross file : " << fileElectron << G4endl;
      
      G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();
      G4String electron = electronDef->GetParticleName();
      
      // For total cross section      
      MapData* tableData = new MapData();
      
      G4MicroElecCrossSectionDataSet_new* tableE = new G4MicroElecCrossSectionDataSet_new(new G4LogLogInterpolation, eV, scaleFactor);
      tableE->LoadData(fileElectron);
      tableData->insert(make_pair(electron, tableE));
      tableTCS[matName] = tableData; //Storage of TCS
      
      // For final state
      const char* path = G4FindDataDir("G4LEDATA");
      if (!path)
	{
	  G4Exception("G4MicroElecElasticModel_new::Initialise","em0006",FatalException,"G4LEDATA environment variable not set.");
	  return;
	}
      
      //Reading DCS file
      std::ostringstream eFullFileName;
      eFullFileName << path << "/microelec/Elastic/elsepa_elastic_cumulated_diffcross_e_" + matName + ".dat";
      G4cout << "Elastic Cumulated Diff Cross : " << eFullFileName.str().c_str() << G4endl;
      std::ifstream eDiffCrossSection(eFullFileName.str().c_str());
      
      if (!eDiffCrossSection)
	G4Exception("G4MicroElecElasticModel_new::Initialise", "em0003", FatalException, "Missing data file: /microelec/sigmadiff_cumulated_elastic_e_Si.dat");
     
      // October 21th, 2014 - Melanie Raine
      // Added clear for MT
      // Diff Cross Sections in cumulated mode      
      TriDimensionMap* eDiffCrossSectionData = new TriDimensionMap(); //Angles 
      std::vector<G4double>* eTdummyVec = new std::vector<G4double>; //Incident energy vector
      VecMap* eProbVec = new VecMap; //Probabilities
      
      eTdummyVec->push_back(0.);
      
      while (!eDiffCrossSection.eof())
	{
	  G4double tDummy; //incident energy
	  G4double eProb; //Proba
	  eDiffCrossSection >> tDummy >> eProb;
	  
	  // SI : mandatory eVecm initialization	  
	  if (tDummy != eTdummyVec->back())
	    {
	      eTdummyVec->push_back(tDummy); //adding values for incident energy points
	      (*eProbVec)[tDummy].push_back(0.); //adding probability for the first angle, equal to 0 		
	    }
	  
	  eDiffCrossSection >> (*eDiffCrossSectionData)[tDummy][eProb]; //adding Angle Value to map
	  
	  if (eProb != (*eProbVec)[tDummy].back()) {
	    (*eProbVec)[tDummy].push_back(eProb); //Adding cumulated proba to map
	  }
	  
	}
      
      //Filling maps for the material
      thetaDataStorage[matName] = eDiffCrossSectionData;
      eIncidentEnergyStorage[matName] = eTdummyVec;
      eProbaStorage[matName] = eProbVec;
    }
    // End final state
    
    if (verboseLevel > 2)
      G4cout << "Loaded cross section files for MicroElec Elastic model" << G4endl;
    
    if (verboseLevel > 0)
      {
	G4cout << "MicroElec Elastic model is initialized " << G4endl
	       << "Energy range: "
	       << LowEnergyLimit() / eV << " eV - "
	       << HighEnergyLimit() / MeV << " MeV"
	       << G4endl; // system("pause"); linux doesn't like
      }
    
    if (isInitialised) { return; }
    fParticleChangeForGamma = GetParticleChangeForGamma();
    isInitialised = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecElasticModel_new::CrossSectionPerVolume(const G4Material* material,
							const G4ParticleDefinition* p,
							G4double ekin,
							G4double,
							G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling CrossSectionPerVolume() of G4MicroElecElasticModel" << G4endl;
  
  isOkToBeInitialised = true;
  currentMaterialName = material->GetName().substr(3, material->GetName().size());
  const G4DataVector cuts;
  Initialise(p, cuts);
  // Calculate total cross section for model
  MapEnergy::iterator lowEPos;
  lowEPos = lowEnergyLimitTable.find(currentMaterialName);
  
  MapEnergy::iterator highEPos;
  highEPos = highEnergyLimitTable.find(currentMaterialName);
  
  MapEnergy::iterator killEPos;
  killEPos = workFunctionTable.find(currentMaterialName);
  
  if (lowEPos == lowEnergyLimitTable.end() || highEPos == highEnergyLimitTable.end() || killEPos == workFunctionTable.end())
    {
      G4String str = "Material ";
      str += currentMaterialName + " not found!";
      G4Exception("G4MicroElecElasticModel_new::EnergyLimits", "em0002", FatalException, str);
      return 0;
    }
  else {
    //	 G4cout << "normal elastic " << G4endl;
    lowEnergyLimit = lowEPos->second;
    highEnergyLimit = highEPos->second;
    killBelowEnergy = killEPos->second;
    
  }
  
  if (ekin < killBelowEnergy) {

    return DBL_MAX; }
  
  G4double sigma=0;
  
  //Phonon for SiO2
  if (currentMaterialName == "SILICON_DIOXIDE" && ekin < 100 * eV) {
    acousticModelEnabled = true;
    
    //Values for SiO2
    G4double kbz = 11.54e9,
      rho = 2.2 * 1000, // [g/cm3] * 1000
      cs = 3560, //Sound speed
      Ebz = 5.1 * 1.6e-19,
      Aac = 17 * Ebz, //A screening parameter
      Eac = 3.5 * 1.6e-19, //C deformation potential
      prefactor = 2.2;// Facteur pour modifier les MFP  
    
    return AcousticCrossSectionPerVolume(ekin, kbz, rho, cs, Aac, Eac, prefactor);
  }
  
  //Elastic
  else {
    acousticModelEnabled = false;
    
    G4double density = material->GetTotNbOfAtomsPerVolume();
    const G4String& particleName = p->GetParticleName();    
    
    TCSMap::iterator tablepos;
    tablepos = tableTCS.find(currentMaterialName);
    
    if (tablepos != tableTCS.end())
      {
	MapData* tableData = tablepos->second;
	
	if (ekin >= lowEnergyLimit && ekin < highEnergyLimit)
	  { 
	    std::map< G4String, G4MicroElecCrossSectionDataSet_new*, std::less<G4String> >::iterator pos;
	    pos = tableData->find(particleName);
	    
	    if (pos != tableData->end())
	      {
		G4MicroElecCrossSectionDataSet_new* table = pos->second;
		if (table != 0)
		  {
		    sigma = table->FindValue(ekin);
		  }
	      }
	    else
	      {
		G4Exception("G4MicroElecElasticModel_new::ComputeCrossSectionPerVolume", "em0002", FatalException, "Model not applicable to particle type.");
	      }
	  }
	else return 1 / DBL_MAX;
      }
    else
      {
	G4String str = "Material ";
	str += currentMaterialName + " TCS table not found!";
	G4Exception("G4MicroElecElasticModel_new::ComputeCrossSectionPerVolume", "em0002", FatalException, str);
      }
    
    if (verboseLevel > 3)
      {
	G4cout << "---> Kinetic energy(eV)=" << ekin / eV << G4endl;
	G4cout << " - Cross section per Si atom (cm^2)=" << sigma / cm / cm << G4endl;
	G4cout << " - Cross section per Si atom (cm^-1)=" << sigma*density / (1. / cm) << G4endl;
      }
    return sigma*density;   
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecElasticModel_new::AcousticCrossSectionPerVolume(G4double ekin, 
								    G4double kbz,
								    G4double rho, 
								    G4double cs, 
								    G4double Aac, 
								    G4double Eac, 
								    G4double prefactor)
{
  
  G4double e = 1.6e-19,
    m0 = 9.10938356e-31,
    h = 1.0546e-34,
    kb = 1.38e-23;
  
  G4double E = (ekin / eV) * e;
  G4double D = (2 / (std::sqrt(2) * std::pow(pi, 2) * std::pow(h, 3))) * (1 + 2 * E) * std::pow(m0, 1.5) * std::sqrt(E);
  
  // Parametres SiO2
  G4double T = 300,
    Ebz = (std::pow(h, 2) * std::pow(kbz, 2)) / (2 * m0),
    hwbz = cs * kbz * h,
    nbz = 1.0 / (exp(hwbz / (kb * T)) - 1),
    Pac;
  
  if (E < Ebz / 4.0)
    {
      Pac = ((pi * kb * T) / (h * std::pow(cs, 2) * rho)) * (std::pow(Eac, 2) * D) / (1 + (E / Aac));
    }
  
  else if (E > Ebz) //Screened relationship
    {
      Pac = ((2 * pi * m0 * (2 * nbz + 1)) / (h * rho * hwbz)) * std::pow(Eac, 2) * D * E * 2 * std::pow((Aac / E), 2) * (((-E / Aac) / (1 + (E / Aac))) + log(1 + (E / Aac)));
    }
  else //Linear interpolation
    {
      G4double fEbz = ((2 * pi * m0 * (2 * nbz + 1)) / (h * rho * hwbz)) * std::pow(Eac, 2) * D * Ebz * 2 * std::pow((Aac / Ebz), 2) * (((-Ebz / Aac) / (1 + (Ebz / Aac))) + log(1 + (Ebz / Aac)));
      G4double fEbz4 = ((pi * kb * T) / (h * std::pow(cs, 2) * rho)) * (std::pow(Eac, 2) * D) / (1 + ((Ebz / 4) / Aac));
      G4double alpha = ((fEbz - fEbz4) / (Ebz - (Ebz / 4)));
      Pac = alpha * E + (fEbz - alpha * Ebz);
    }
  
  G4double MFP = (std::sqrt(2 * E / m0) / (prefactor * Pac)) * m;
  
  return  1 / MFP;  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MicroElecElasticModel_new::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
						const G4MaterialCutsCouple* /*couple*/,
						const G4DynamicParticle* aDynamicElectron,
						G4double,
						G4double)
{

  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4MicroElecElasticModel" << G4endl;
  
  G4double electronEnergy0 = aDynamicElectron->GetKineticEnergy();
  
 if (electronEnergy0 < killBelowEnergy)
   {
     fParticleChangeForGamma->SetProposedKineticEnergy(0.);
     fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
     fParticleChangeForGamma->ProposeLocalEnergyDeposit(electronEnergy0);
     return;
   }
 
 if (electronEnergy0 < highEnergyLimit)
   {
     G4double cosTheta = 0;
     if (acousticModelEnabled)
       {
	 cosTheta = 1 - 2 * G4UniformRand(); //Isotrope
       }
     else if (electronEnergy0 >= lowEnergyLimit)
       {
	 cosTheta = RandomizeCosTheta(electronEnergy0);
       }
     
     G4double phi = 2. * pi * G4UniformRand();
     
     G4ThreeVector zVers = aDynamicElectron->GetMomentumDirection();
     G4ThreeVector xVers = zVers.orthogonal();
     G4ThreeVector yVers = zVers.cross(xVers);
     
     G4double xDir = std::sqrt(1. - cosTheta*cosTheta);
     G4double yDir = xDir;
     xDir *= std::cos(phi);
     yDir *= std::sin(phi);
     
     G4ThreeVector zPrimeVers((xDir*xVers + yDir*yVers + cosTheta*zVers));
     
     fParticleChangeForGamma->ProposeMomentumDirection(zPrimeVers.unit());
     fParticleChangeForGamma->SetProposedKineticEnergy(electronEnergy0);
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double  G4MicroElecElasticModel_new::DamageEnergy(G4double T,G4double A, G4double Z)
{
  //.................. T in  eV!!!!!!!!!!!!!
  G4double Z2= Z;
  G4double M2= A;
  G4double k_d;
  G4double epsilon_d;
  G4double g_epsilon_d;
  G4double E_nu;

  k_d=0.1334*std::pow(Z2,(2./3.))*std::pow(M2,(-1./2.));
  epsilon_d=0.01014*std::pow(Z2,(-7./3.))*(T/eV);
  g_epsilon_d= epsilon_d+0.40244*std::pow(epsilon_d,(3./4.))+3.4008*std::pow(epsilon_d,(1./6.));
  
  E_nu=1./(1.+ k_d*g_epsilon_d);
  
  return E_nu;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecElasticModel_new::Theta
  (G4ParticleDefinition * particleDefinition, G4double k, G4double integrDiff)
{
  
  G4double theta = 0.;
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
    ThetaMap::iterator iterator_angle;
    iterator_angle = thetaDataStorage.find(currentMaterialName);
    
    energyMap::iterator iterator_energy;
    iterator_energy = eIncidentEnergyStorage.find(currentMaterialName);
    
    ProbaMap::iterator iterator_proba;
    iterator_proba = eProbaStorage.find(currentMaterialName);
    
    if (iterator_angle != thetaDataStorage.end() && iterator_energy != eIncidentEnergyStorage.end() && iterator_proba != eProbaStorage.end())
      {
	TriDimensionMap* eDiffCrossSectionData = iterator_angle->second; //Theta points
	std::vector<G4double>* eTdummyVec = iterator_energy->second;
	VecMap* eVecm = iterator_proba->second;
	
	auto t2 = std::upper_bound(eTdummyVec->begin(), eTdummyVec->end(), k);
	auto t1 = t2 - 1;	
        auto e12 = std::upper_bound((*eVecm)[(*t1)].begin(), (*eVecm)[(*t1)].end(), integrDiff);
	auto e11 = e12 - 1;	
	auto e22 = std::upper_bound((*eVecm)[(*t2)].begin(), (*eVecm)[(*t2)].end(), integrDiff);
	auto e21 = e22 - 1;
	
	valueT1 = *t1;
	valueT2 = *t2;
	valueE21 = *e21;
	valueE22 = *e22;
	valueE12 = *e12;
	valueE11 = *e11;
	
	xs11 = (*eDiffCrossSectionData)[valueT1][valueE11]; 
	xs12 = (*eDiffCrossSectionData)[valueT1][valueE12];
	xs21 = (*eDiffCrossSectionData)[valueT2][valueE21];
	xs22 = (*eDiffCrossSectionData)[valueT2][valueE22];
      }
    else
      {
	G4String str = "Material ";
	str += currentMaterialName + " not found!";
	G4Exception("G4MicroElecElasticModel_new::ComputeCrossSectionPerVolume", "em0002", FatalException, str);
      }
    
  }
  
  if (xs11==0 || xs12==0 ||xs21==0 ||xs22==0) return (0.);
  
  theta = QuadInterpolator(  valueE11, valueE12,
    			     valueE21, valueE22,
			     xs11, xs12,
			     xs21, xs22,
			     valueT1, valueT2,
			     k, integrDiff );
  
  return theta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecElasticModel_new::LinLogInterpolate(G4double e1,
						    G4double e2,
						    G4double e,
						    G4double xs1,
						    G4double xs2)
{
  G4double d1 = std::log(xs1);
  G4double d2 = std::log(xs2);
  G4double value = G4Exp(d1 + (d2 - d1)*(e - e1)/ (e2 - e1));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecElasticModel_new::LinLinInterpolate(G4double e1,
						    G4double e2,
						    G4double e,
						    G4double xs1,
						    G4double xs2)
{
  G4double d1 = xs1;
  G4double d2 = xs2;
  G4double value = (d1 + (d2 - d1)*(e - e1)/ (e2 - e1));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecElasticModel_new::LogLogInterpolate(G4double e1,
						    G4double e2,
						    G4double e,
						    G4double xs1,
						    G4double xs2)
{
  G4double a = (std::log10(xs2)-std::log10(xs1)) / (std::log10(e2)-std::log10(e1));
  G4double b = std::log10(xs2) - a*std::log10(e2);
  G4double sigma = a*std::log10(e) + b;
  G4double value = (std::pow(10.,sigma));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecElasticModel_new::QuadInterpolator(G4double e11, G4double e12,
						   G4double e21, G4double e22,
						   G4double xs11, G4double xs12,
						   G4double xs21, G4double xs22,
						   G4double t1, G4double t2,
						   G4double t, G4double e)
{
  
  
  // Lin-Lin
  G4double interpolatedvalue1 = LinLinInterpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = LinLinInterpolate(e21, e22, e, xs21, xs22);
  G4double value = LinLinInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
  
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecElasticModel_new::RandomizeCosTheta(G4double k)
{
  G4double integrdiff=0;
  G4double uniformRand=G4UniformRand();
  integrdiff = uniformRand;

  G4double theta=0.;
  G4double cosTheta=0.;
  theta = Theta(G4Electron::ElectronDefinition(),k/eV,integrdiff);
  
  cosTheta= std::cos(theta*pi/180.);
  
  return cosTheta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void G4MicroElecElasticModel_new::SetKillBelowThreshold (G4double threshold) 
{ 
  killBelowEnergy = threshold; 
  
  if (threshold < 5*CLHEP::eV)
    {
      G4Exception ("*** WARNING : the G4MicroElecElasticModel class is not validated below 5 eV !","",JustWarning,"") ;
      threshold = 5*CLHEP::eV;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
