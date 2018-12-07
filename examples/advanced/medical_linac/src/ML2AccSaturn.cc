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
// The code was written by :
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#include "ML2AccSaturn.hh"
#include "ML2AccSaturnMessenger.hh"
#include "ML2Accelerator.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"




using namespace std;

CML2AccSaturn::CML2AccSaturn()
{
	PVWorld = 0;
	accSaturnMessenger = new CML2AccSaturnMessenger(this);
	buildMaterial_SSteel();
	buildMaterial_XC10();
	buildMaterial_WNICU();
	buildMaterial_Kapton();
}

CML2AccSaturn::~CML2AccSaturn(void)
{
}

//void CML2AccSaturn::reset()
//{
//	vec_leavesA.clear();
//	vec_leavesB.clear();
//}

CML2AccSaturn* CML2AccSaturn::instance = 0;

CML2AccSaturn* CML2AccSaturn::GetInstance(void)
{
	if (instance == 0)
	{
//		G4cout << "in CML2AccSaturn::GetInstance" << G4endl;
		instance = new CML2AccSaturn();
	}
	return instance;
}

void CML2AccSaturn::writeInfo()
{
	G4cout << "----------------------------------------------------------------" << G4endl;
	G4cout << "Accelerator VARIAN SATURN 43" << G4endl;
	G4cout <<"\n\n\tnominal beam energy: "<<idEnergy << G4endl;
	G4cout << "\tdistance isocentre [mm]:"<< isoCentre/mm << G4endl;
	G4cout <<"\tJaw X aperture: 1) "<< jaw1XAperture/mm<<"[mm]\t2) " << jaw2XAperture/mm<< " [mm]"<< G4endl;
	G4cout <<"\tJaw Y aperture: 1) "<< jaw1YAperture/mm<<"[mm]\t2) " << jaw2YAperture/mm<< " [mm]"<< G4endl;
	if (vec_leavesA.size()>0)
	{
		G4cout << "\tleaves A aperture [mm]" << G4endl;
		for (int i=0; i< (int)vec_leavesA.size(); i++)
		{
			G4cout<<"\t" <<i <<") "<< vec_leavesA[i]/mm << G4endl;
		}
	}
	else
	{
		G4cout << "\tNo leaves A" << G4endl;
	}
	if (vec_leavesB.size()>0)
	{
		G4cout << "\tleaves B aperture [mm]" << G4endl;
		for (int i=0; i< (int)vec_leavesB.size(); i++)
		{
			G4cout<<"\t" <<i <<") "<< vec_leavesB[i]/mm << G4endl;
		}
	}
	else
	{
		G4cout << "\tNo leaves B" << G4endl;
	}
	G4cout << "______________________________________________________________" << G4endl;
}


void CML2AccSaturn::buildMaterial_SSteel()
{

	G4Material *elFe = G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");
	G4Material *elCr = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cr");
	G4Material *elNi = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ni");
	G4Material *elMn = G4NistManager::Instance()->FindOrBuildMaterial("G4_Mn");
	G4Material *elSi = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
	G4Material *elC  = G4NistManager::Instance()->FindOrBuildMaterial("G4_C");

	mat_ssteel = new G4Material("StainlessSteel", 7.8 *g/cm3, 6);
	mat_ssteel -> AddMaterial(elFe, 0.6898);
	mat_ssteel -> AddMaterial(elCr, 0.18);
	mat_ssteel -> AddMaterial(elNi, 0.10);
	mat_ssteel -> AddMaterial(elMn, 0.02);
	mat_ssteel -> AddMaterial(elSi, 0.01);
	mat_ssteel -> AddMaterial(elC,  0.0002);

}

void CML2AccSaturn::buildMaterial_XC10()
{
	G4Material *elFe = G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");
	G4Material *elMn = G4NistManager::Instance()->FindOrBuildMaterial("G4_Mn");
	G4Material *elC  = G4NistManager::Instance()->FindOrBuildMaterial("G4_C");

	mat_XC10 = new G4Material("CARBON_STEEL", 7.8 *g/cm3, 3);
	mat_XC10 -> AddMaterial(elFe, 0.993);
	mat_XC10 -> AddMaterial(elMn, 0.006);
	mat_XC10 -> AddMaterial(elC,  0.001);

}

void CML2AccSaturn::buildMaterial_WNICU()
{
	G4Material *elW  = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
	G4Material *elNi = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ni");
	G4Material *elCu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");

	mat_WNICU = new G4Material("Denal(WNICU)", 16.8*g/cm3, 3);
	mat_WNICU -> AddMaterial(elW,  0.905);
	mat_WNICU -> AddMaterial(elNi, 0.07);
	mat_WNICU -> AddMaterial(elCu, 0.025);

}

void CML2AccSaturn::buildMaterial_Kapton()
{
	G4Material *elH  = G4NistManager::Instance()->FindOrBuildMaterial("G4_H");
	G4Material *elC  = G4NistManager::Instance()->FindOrBuildMaterial("G4_C");
	G4Material *elN  = G4NistManager::Instance()->FindOrBuildMaterial("G4_N");
	G4Material *elO  = G4NistManager::Instance()->FindOrBuildMaterial("G4_O");

	mat_Kapton = new G4Material("kapton", 1.42*g/cm3, 4);
	mat_Kapton -> AddMaterial(elC, 0.69113);
	mat_Kapton -> AddMaterial(elO, 0.20924);
	mat_Kapton -> AddMaterial(elN, 0.07327);
	mat_Kapton -> AddMaterial(elH, 0.026362);

}

G4Material* CML2AccSaturn::getMaterial(const G4String materialName)
{
	G4Material *material = 0;

	if (materialName == "steel")
	{
		material = mat_ssteel;
	}
	else if (materialName == "xc10")
	{
		material = mat_XC10;
	}
	else if (materialName == "wnicu")
	{
		material = mat_WNICU;
	}
	else if (materialName == "kapton")
	{
		material = mat_Kapton;
	}
	else
	{
		G4cout << "Sbajato" << G4endl;
	}
	return material;
}

void CML2AccSaturn::Construct(G4VPhysicalVolume *PWorld, G4double iso)
{
	PVWorld = PWorld;
    setIsoCentre(iso);
	target();
	vacuumWindow();
	ionizationChamber();
	flatteningFilter();
	primaryCollimator();
	Jaw1X();
	Jaw2X();
	Jaw1Y();
	Jaw2Y();
}

bool CML2AccSaturn::target()
{
	bool bCreated = false;

	G4ThreeVector targetCentre, boxHalfSize, tubeCentre;
	G4double tubeRadius, tubeHeight;
	targetCentre.set(0,0,-3.5*mm);
	boxHalfSize.set(5.*mm,5.*mm, 7.5*mm);
	tubeCentre.set(0,0,-5.5*mm);
	tubeRadius = 3.*mm;
	tubeHeight = 5.5*mm;

	G4Material *diskMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ti");
	G4Material *targetMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");

	// Region for cuts
	G4Region *regVol;
	regVol = new G4Region("TargetR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts -> SetProductionCut(1.*cm);
	regVol -> SetProductionCuts(cuts);

	// Physical volumes
	G4Box *box = new G4Box("targetBox", boxHalfSize.getX(), boxHalfSize.getY(), boxHalfSize.getZ());
	G4Tubs *tube = new G4Tubs("targetTube", 0.*mm, tubeRadius, tubeHeight, 0.*deg, 360.*deg);
	G4SubtractionSolid* BoxMinusTube = new G4SubtractionSolid("TargetSolid", box, tube,0, G4ThreeVector(0.,0.,-4.*mm));
	G4LogicalVolume* BoxMinusTubeLV = new G4LogicalVolume(BoxMinusTube, targetMaterial, "BoxMinusTubeLV",0,0,0);
	new G4PVPlacement(0, targetCentre, "TargetPV", BoxMinusTubeLV, this->PVWorld, false, 0);

	G4Tubs *disk = new G4Tubs("targetDisk", 0.*mm, 4.*mm, 0.025*mm, 0.*deg, 360.*deg);
	G4LogicalVolume *diskLV = new G4LogicalVolume(disk,diskMaterial,"targetDiskLV",0,0,0);

	new G4PVPlacement(0, G4ThreeVector(0.,0.,-15.*mm), "diskPV", diskLV, this->PVWorld, false, 0);

	// Visualization
	G4VisAttributes* simpleAlSVisAtt;
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Grey());
	simpleAlSVisAtt -> SetVisibility(true);
	simpleAlSVisAtt -> SetForceSolid(true);
	BoxMinusTubeLV -> SetVisAttributes(simpleAlSVisAtt);
	diskLV -> SetVisAttributes(simpleAlSVisAtt);

	// Region for cuts
	BoxMinusTubeLV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(BoxMinusTubeLV);
	diskLV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(diskLV);

	bCreated = true;
	return bCreated;

}
bool CML2AccSaturn::primaryCollimator()
{
	bool bCreated = false;

	G4double tube1Radius, tube1Height, tube1Z;
	G4double tube2Radius, tube2Height, tube2Z;
	G4double tube3Radius, tube3Height;
	G4double tube4Radius, tube4Height, tube4Z;
	G4double cone1RadiusMin, cone1RadiusMax, cone1Height;
	G4double cone2RadiusMin, cone2RadiusMax, cone2Height;
	G4double cone3RadiusMin, cone3RadiusMax, cone3Height;

	tube1Radius = 141./2.*mm; // upper principal tube
	tube1Height = 79./2.*mm;
	tube1Z = 5.*mm+tube1Height;

	tube2Radius = 141./2.*mm; // middle principal tube
	tube2Height = 57.5/2.*mm;
	tube2Z = 1.35*mm+tube1Z+tube1Height+tube2Height;

	tube3Radius = 34./2.*mm; // tube of air to be subtracted
	tube3Height = 15.86/2.*mm;

	tube4Radius = 241./2.*mm; // lower principal tube
	tube4Height = 40.50/2.*mm;
	tube4Z = 161.*mm+tube4Height;

	cone1RadiusMin=34./2.*mm; // upper cone
	cone1RadiusMax=54./2.*mm;
	cone1Height=(79.-15.86)/2.*mm;

	cone2RadiusMin=53./2.*mm; // middle cone
	cone2RadiusMax=81./2.*mm;
	cone2Height=57.50/2.*mm;

	cone3RadiusMin=88.06/2.*mm; // lower cone
	cone3RadiusMax=109.76/2.*mm;
	cone3Height=40.50/2.*mm;

//	G4Material *upperTubeMaterial = SetMaterials("XC10");
//	G4Material *middleTubeMaterial = SetMaterials("Denal(WNICU)");
	G4Material *upperTubeMaterial = getMaterial("xc10");
	G4Material *middleTubeMaterial = getMaterial("wnicu");

	G4Material *lowerTubeMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
	G4Region *regVol;
	G4VisAttributes* simpleAlSVisAtt;

	// Region for cuts
	regVol = new G4Region("primaryCollimatorR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(1.*cm);
	regVol->SetProductionCuts(cuts);

	// Physical volumes
	G4Tubs *tube1 = new G4Tubs(
			"PriCollTube1",
			0.*mm,
			tube1Radius,
			tube1Height,
			0.*deg,
			360.*deg);
	G4Cons *cone1 = new G4Cons(
			"PriCollCone1",
			0.,
			cone1RadiusMin,
			0.,
			cone1RadiusMax,
			cone1Height,
			0,
			360.*deg);
	G4Tubs *tube2 = new G4Tubs("PriCollTube2", 0.*mm, tube2Radius, tube2Height, 0.*deg, 360.*deg);
	G4Cons *cone2 = new G4Cons("PriCollCone2", 0., cone2RadiusMin, 0.,cone2RadiusMax,cone2Height, 0, 360.*deg);
	G4Tubs *tube3 = new G4Tubs("PriCollTube3", 0.*mm, tube3Radius, tube3Height, 0.*deg, 360.*deg);
	G4Cons *cone3 = new G4Cons("PriCollCone3", 0., cone3RadiusMin,0.,cone3RadiusMax, cone3Height, 0, 360.*deg);
	G4Tubs *tube4 = new G4Tubs("PriCollTube4", 0.*mm, tube4Radius, tube4Height, 0.*deg, 360.*deg);

	G4UnionSolid *tube3AndCone1 = new G4UnionSolid("PriCollTube3AndCone1",tube3,cone1,0, G4ThreeVector(0.,0.,tube3Height+cone1Height));
	G4SubtractionSolid* tube1NotTube3AndCone1 = new G4SubtractionSolid("PriColltube1NotTube3AndCone1", tube1, tube3AndCone1,0, G4ThreeVector(0.,0.,-tube1Height+tube3Height));
	G4LogicalVolume* tube1NotTube3AndCone1LV = new G4LogicalVolume(tube1NotTube3AndCone1, upperTubeMaterial, "PriColltube1NotTube3AndCone1LV",0,0,0);

	G4SubtractionSolid* tube2NotCone2 = new G4SubtractionSolid("PriCollTube2NotCone2", tube2, cone2);
	G4LogicalVolume* tube2NotCone2LV  = new G4LogicalVolume(tube2NotCone2, middleTubeMaterial, "PriCollTube2NotCone2LV",0,0,0);

	G4SubtractionSolid* tube4NotCone3 = new G4SubtractionSolid("PriCollTube4NotCone3", tube4, cone3);
	G4LogicalVolume* tube4NotCone3LV  = new G4LogicalVolume(tube4NotCone3, lowerTubeMaterial, "PriCollTube4NotCone3LV",0,0,0);

	new G4PVPlacement(0, G4ThreeVector(0,0,tube1Z), "PriCollUpperPV", tube1NotTube3AndCone1LV, PVWorld, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0,0,tube2Z), "PriCollMiddlePV", tube2NotCone2LV, PVWorld, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0,0,tube4Z), "PriCollLowerPV", tube4NotCone3LV, PVWorld, false, 0);

	// Visualization
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAtt -> SetVisibility(true);
	simpleAlSVisAtt -> SetForceSolid(true);
	tube1NotTube3AndCone1LV -> SetVisAttributes(simpleAlSVisAtt);
	simpleAlSVisAtt -> SetColour(G4Colour::Red());
	tube2NotCone2LV -> SetVisAttributes(simpleAlSVisAtt);
	simpleAlSVisAtt -> SetColour(G4Colour::Blue());
	tube4NotCone3LV -> SetVisAttributes(simpleAlSVisAtt);


	// Region for cuts
	tube1NotTube3AndCone1LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(tube1NotTube3AndCone1LV);
	tube2NotCone2LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(tube2NotCone2LV);
	tube4NotCone3LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(tube4NotCone3LV);

	bCreated = true;
	return bCreated;
}

bool CML2AccSaturn::vacuumWindow()
{
	bool bCreated = false;

	G4double tubeRadius, tubeHeight, tubeZ;

	tubeRadius = 116.53/2.*mm; // upper principal tube
	tubeHeight = 2./2.*mm;
	tubeZ = 215.75*mm+tubeHeight;

	G4Material *elAl= G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

	G4Region *regVol;
	G4VisAttributes* simpleAlSVisAtt;
	// Region for cuts
	regVol= new G4Region("VacWindowR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*cm);
	regVol->SetProductionCuts(cuts);

	// Physical volumes
	G4Tubs  *tube = new G4Tubs("VacWindowTube", 0.*mm, tubeRadius, tubeHeight, 0.*deg, 360.*deg);
	G4LogicalVolume* tubeLV = new G4LogicalVolume(tube, elAl, "VacWindowTubeLV",0,0,0);

	new G4PVPlacement(0, G4ThreeVector(0,0,tubeZ), "VacWindowTubePV", tubeLV, PVWorld, false, 0);

	// Visualization
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAtt -> SetVisibility(true);
	simpleAlSVisAtt -> SetForceSolid(true);
	tubeLV -> SetVisAttributes(simpleAlSVisAtt);

	// Region for cuts
	tubeLV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(tubeLV);

	bCreated = true;
	return bCreated;
}
bool CML2AccSaturn::flatteningFilter()
{
	bool bCreated = false;

	G4double tube1Radius, tube1Height;
	G4double cone1RadiusMin, cone1RadiusMax, cone1Height;
	G4double cone2RadiusMin, cone2RadiusMax, cone2Height;
	G4double cone3RadiusMin, cone3RadiusMax, cone3Height, ffZ;

	tube1Radius = 108./2.*mm; // upper principal tube
	tube1Height = 7.5/2.*mm;

	cone1RadiusMin = 54./2.*mm; // upper cone
	cone1RadiusMax = 76./2.*mm;
	cone1Height = 13.7/2.*mm;

	cone2RadiusMin = 8./2.*mm; // middle cone
	cone2RadiusMax = 54./2.*mm;
	cone2Height = (44.3-13.7)/2.*mm;

	cone3RadiusMin = 0.000001*mm; // lower cone
	cone3RadiusMax = 8./2.*mm;
	cone3Height = (46.8-44.3)/2.*mm;
	ffZ=149.50+tube1Height; // the half point is located ad the centre of the tube1 because of the solids unions and translations

//	G4Material *ffMaterial = SetMaterials("StainlessSteel");
	G4Material *ffMaterial = getMaterial("steel");

	G4Region *regVol;
	G4VisAttributes* simpleAlSVisAtt;
	// Region for cuts
	regVol= new G4Region("ffR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*cm);
	regVol->SetProductionCuts(cuts);

	// Physical volumes

	G4Tubs *tube1 = new G4Tubs("ffTube1", 0.*mm, tube1Radius, tube1Height, 0.*deg, 360.*deg);
	G4Cons *cone1 = new G4Cons("ffCone1",0.,cone1RadiusMax, 0., cone1RadiusMin, cone1Height, 0, 360.*deg);
	G4Cons *cone2 = new G4Cons("ffCone2", 0.,cone2RadiusMax, 0., cone2RadiusMin,cone2Height, 0, 360.*deg);
	G4Cons *cone3 = new G4Cons("ffCone3", 0.,cone3RadiusMax, 0., cone3RadiusMin,cone3Height, 0, 360.*deg);

	G4double pos = (cone1Height-tube1Height-1.)*mm;
	G4UnionSolid *tube1AndCone1 = new G4UnionSolid(
			"ffTube1AndCone1",
			tube1,
			cone1,
			0,
			G4ThreeVector(0., 0., pos));

	pos+=cone1Height+cone2Height;
	G4UnionSolid *tubeCone1AndCone2 = new G4UnionSolid(
			"fftubeCone1AndCone2",
			tube1AndCone1,
			cone2,
			0,
			G4ThreeVector(0., 0., pos));
	pos+=cone2Height+cone3Height;
	G4UnionSolid *tubeCone12AndCone3 = new G4UnionSolid(
			"fftubeCone12AndCone3",
			tubeCone1AndCone2,
			cone3,
			0,
			G4ThreeVector(0., 0., pos));

	G4LogicalVolume* ffLV = new G4LogicalVolume(tubeCone12AndCone3, ffMaterial, "ffLV",0,0,0);

	new G4PVPlacement(0, G4ThreeVector(0,0,ffZ), "ffPV", ffLV, PVWorld, false, 0);

	// Visualization
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Magenta());
	simpleAlSVisAtt->SetVisibility(true);
	simpleAlSVisAtt->SetForceSolid(true);
	ffLV->SetVisAttributes(simpleAlSVisAtt);

    // Region for cuts
	ffLV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(ffLV);

	bCreated = true;
	return bCreated;
}
bool CML2AccSaturn::ionizationChamber()
{
	bool bCreated = false;

	G4double tubeRadius, tubeA1Z, tubeA2Z, tubeA3Z, tubeA4Z, tubeA5Z, tubeA6Z;
	G4double kaptonThickness, AlThickness1, AlThickness8;

	kaptonThickness = 0.025/2.*mm;
	AlThickness1 = 1.6e-5/2.*cm;
	AlThickness8 = 8e-6/2.*cm;

    tubeRadius = 110./2.*mm; // upper principal tube
    tubeA1Z=(202.5+AlThickness8)*mm;
	tubeA2Z=(204.5+AlThickness1)*mm;
	tubeA3Z=(206.5+AlThickness8)*mm;
	tubeA4Z=(207.5+AlThickness8)*mm;
	tubeA5Z=(209.5+AlThickness1)*mm;
	tubeA6Z=(211.5+AlThickness8)*mm;

	G4double pos1=(AlThickness1+kaptonThickness)*mm;
	G4double pos8=(AlThickness8+kaptonThickness)*mm;


	G4Material *elAL = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
//	G4Material *kapton = SetMaterials("kapton");
	G4Material *kapton = getMaterial("kapton");

	G4Region *regVol;
	G4VisAttributes* simpleAlSVisAttAl1;
	G4VisAttributes* simpleAlSVisAttAl8;
	G4VisAttributes* simpleAlSVisAttK;
		// Region for cuts
	regVol = new G4Region("IonChamberR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts -> SetProductionCut(0.1*mm);
	regVol -> SetProductionCuts(cuts);

		// Physical volumes

	G4Tubs  *tubeAl8 = new G4Tubs("ICTube1A", 0.*mm, tubeRadius, AlThickness8, 0.*deg, 360.*deg);
	G4Tubs  *tubeK = new G4Tubs("ICTube1B", 0.*mm, tubeRadius, kaptonThickness, 0.*deg, 360.*deg);
	G4Tubs  *tubeAl1 = new G4Tubs("ICTube2A", 0.*mm, tubeRadius, AlThickness1, 0.*deg, 360.*deg);


	G4LogicalVolume* IC_Al8_LV = new G4LogicalVolume(tubeAl8, elAL, "IC_Al8_LV", 0, 0, 0);
	G4LogicalVolume* IC_Al1_LV = new G4LogicalVolume(tubeAl1, elAL, "IC_Al1_LV", 0, 0, 0);
	G4LogicalVolume* IC_K_LV   = new G4LogicalVolume(tubeK, kapton, "IC_Al_LV",  0, 0, 0);

	new G4PVPlacement(0, G4ThreeVector(0,0,tubeA1Z), "IC_AL8_PV1", IC_Al8_LV, PVWorld, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0,0,tubeA1Z+pos8), "IC_K_PV1", IC_K_LV, PVWorld, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0,0,tubeA2Z), "IC_AL1_PV2", IC_Al1_LV, PVWorld, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0,0,tubeA2Z+pos1), "IC_K_PV2", IC_K_LV, PVWorld, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0,0,tubeA3Z), "IC_AL8_PV3", IC_Al8_LV, PVWorld, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0,0,tubeA3Z+pos8), "IC_K_PV3", IC_K_LV, PVWorld, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0,0,tubeA4Z), "IC_AL8_PV4", IC_Al8_LV, PVWorld, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0,0,tubeA4Z+pos8), "IC_K_PV4", IC_K_LV, PVWorld, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0,0,tubeA5Z), "IC_AL1_PV5", IC_Al1_LV, PVWorld, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0,0,tubeA5Z+pos1), "IC_K_PV5", IC_K_LV, PVWorld, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0,0,tubeA6Z), "IC_AL8_PV6", IC_Al8_LV, PVWorld, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0,0,tubeA6Z+pos8), "IC_K_PV6", IC_K_LV, PVWorld, false, 0);


	// Visualization
	simpleAlSVisAttAl8 = new G4VisAttributes(G4Colour::Magenta());
	simpleAlSVisAttAl8 -> SetVisibility(true);
	simpleAlSVisAttAl8 -> SetForceSolid(true);
	IC_Al8_LV -> SetVisAttributes(simpleAlSVisAttAl8);

	simpleAlSVisAttAl1 = new G4VisAttributes(G4Colour::Red());
	simpleAlSVisAttAl1 -> SetVisibility(true);
	simpleAlSVisAttAl1 -> SetForceSolid(true);
	IC_Al1_LV -> SetVisAttributes(simpleAlSVisAttAl1);

	simpleAlSVisAttK = new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAttK -> SetVisibility(true);
	simpleAlSVisAttK -> SetForceSolid(true);
	IC_K_LV -> SetVisAttributes(simpleAlSVisAttK);

    // Region for cuts
	IC_Al1_LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(IC_Al1_LV);
	IC_Al8_LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(IC_Al8_LV);
	IC_K_LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(IC_K_LV);

	bCreated = true;
	return bCreated;
}

void CML2AccSaturn::SetJawAperture(G4int idJaw, G4ThreeVector &centre, G4ThreeVector halfSize,  G4RotationMatrix *cRotation)
{
	using namespace std;
	G4double theta, x, y, z, dx, dy;
	x=centre.getX();
	y=centre.getY();
	z=centre.getZ();
	dx=halfSize.getX();
	dy=halfSize.getY();
	
	
	switch (idJaw)
	{
	case 1: //idJaw1XV2100:
		theta=fabs(atan(jaw1XAperture/isoCentre));
		centre.set(z*sin(theta)+dx*cos(theta), y, z*cos(theta)-dx*sin(theta));
		cRotation->rotateY(-theta);
		break;
	case 2: //idJaw2XV2100:
		theta=fabs(atan(jaw2XAperture/isoCentre));
		centre.set(-(z*sin(theta)+dx*cos(theta)), y, z*cos(theta)-dx*sin(theta));
		cRotation->rotateY(theta);
		break;
	case 3: //idJaw1YV2100:
		theta=fabs(atan(jaw1YAperture/isoCentre));
		centre.set(x, z*sin(theta)+dy*cos(theta), z*cos(theta)-dy*sin(theta));
		cRotation->rotateX(theta);
		break;
	case 4: //idJaw2YV2100:
		theta=fabs(atan(jaw2YAperture/isoCentre));
		centre.set(x, -(z*sin(theta)+dy*cos(theta)), z*cos(theta)-dy*sin(theta));
		cRotation->rotateX(-theta);
		break;
	}
}



bool CML2AccSaturn::Jaw1X()
{
	bool bCreated = false;

	G4double boxSide, box1Height, box1Z, box2Height, box2Z, box3Height, box3Z, box4Height, box4Z, box5Height, box5Z, box6Height, box6Z;

	boxSide=101./2.*mm; // upper principal box

	box1Height=3./2.*mm;
	box1Z=275.50*mm+box1Height;

	box2Height=35./2.*mm;
	box2Z=box1Z+box1Height+box2Height;

	box3Height=35./2.*mm;
	box3Z=box2Z+box2Height+box3Height;

	box4Height=27./2.*mm;
	box4Z=box3Z+box3Height+box4Height;

	box5Height=10./2.*mm;
	box5Z=472.50*mm+box5Height;

	box6Height=5./2.*mm;
	box6Z=box5Z+box5Height+box6Height;

	G4Material *elPb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
//	G4Material *XC10 = SetMaterials("XC10");
//	G4Material *WNICU = SetMaterials("Denal(WNICU)");
	G4Material *XC10 = getMaterial("xc10");
	G4Material *WNICU = getMaterial("wnicu");

	G4Region *regVol;
	G4VisAttributes* simpleAlSVisAttPb;
	G4VisAttributes* simpleAlSVisAttXC10;
	G4VisAttributes* simpleAlSVisAttWNICU;
		// Region for cuts
	regVol = new G4Region("Jaws1XR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts -> SetProductionCut(2.*cm);
	regVol -> SetProductionCuts(cuts);

		// Physical volumes

	G4Box *box1 = new G4Box("Jaws1XBox1", boxSide, boxSide, box1Height);
	G4Box *box2 = new G4Box("Jaws1XBox2", boxSide, boxSide, box2Height);
	G4Box *box3 = new G4Box("Jaws1XBox3", boxSide, boxSide, box3Height);
	G4Box *box4 = new G4Box("Jaws1XBox4", boxSide, boxSide, box4Height);
	G4Box *box5 = new G4Box("Jaws1XBox5", boxSide, boxSide, box5Height);
	G4Box *box6 = new G4Box("Jaws1XBox6", boxSide, boxSide, box6Height);
	G4LogicalVolume* box1LV = new G4LogicalVolume(box1, XC10,  "Jaws1XLV1", 0, 0, 0);
	G4LogicalVolume* box2LV = new G4LogicalVolume(box2, elPb,  "Jaws1XLV2", 0, 0, 0);
	G4LogicalVolume* box3LV = new G4LogicalVolume(box3, WNICU, "Jaws1XLV3", 0, 0, 0);
	G4LogicalVolume* box4LV = new G4LogicalVolume(box4, elPb,  "Jaws1XLV4", 0, 0, 0);
	G4LogicalVolume* box5LV = new G4LogicalVolume(box5, elPb,  "Jaws1XLV5", 0, 0, 0);
	G4LogicalVolume* box6LV = new G4LogicalVolume(box6, WNICU, "Jaws1XLV6", 0, 0, 0);

	G4ThreeVector centre;
	G4RotationMatrix *cRotationId = new G4RotationMatrix();
	G4RotationMatrix *cRotation = new G4RotationMatrix();

	*cRotation = *cRotationId;
	centre.set(boxSide,0.,box1Z);
	SetJawAperture(1, centre, G4ThreeVector(boxSide, boxSide, box1Height), cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws1XPV1", box1LV, PVWorld, false, 0);

	*cRotation = *cRotationId;
	centre.set(boxSide,0.,box2Z);
	SetJawAperture(1, centre, G4ThreeVector(boxSide, boxSide, box2Height), cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws1XPV2", box2LV, PVWorld, false, 0);

	*cRotation = *cRotationId;
	centre.set(boxSide,0.,box3Z);
	SetJawAperture(1, centre, G4ThreeVector(boxSide, boxSide, box3Height), cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws1XPV3", box3LV, PVWorld, false, 0);

	*cRotation = *cRotationId;
	centre.set(boxSide,0.,box4Z);
	SetJawAperture(1, centre, G4ThreeVector(boxSide, boxSide, box4Height), cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws1XPV4", box4LV, PVWorld, false, 0);

	*cRotation = *cRotationId;
	centre.set(boxSide,0.,box5Z);
	SetJawAperture(1, centre, G4ThreeVector(boxSide, boxSide, box5Height), cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws1XPV5", box5LV, PVWorld, false, 0);

	*cRotation = *cRotationId;
	centre.set(boxSide,0.,box6Z);
	SetJawAperture(1, centre, G4ThreeVector(boxSide, boxSide, box6Height), cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws1XPV6", box6LV, PVWorld, false, 0);

	// Visualization
	simpleAlSVisAttPb = new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAttPb -> SetVisibility(true);
	simpleAlSVisAttPb -> SetForceSolid(true);
	simpleAlSVisAttXC10 = new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAttXC10 -> SetVisibility(true);
	simpleAlSVisAttXC10 ->SetForceSolid(true);
	simpleAlSVisAttWNICU = new G4VisAttributes(G4Colour::Red());
	simpleAlSVisAttWNICU ->SetVisibility(true);
	simpleAlSVisAttWNICU ->SetForceSolid(true);

	box1LV -> SetVisAttributes(simpleAlSVisAttXC10);
	box2LV -> SetVisAttributes(simpleAlSVisAttPb);
	box3LV -> SetVisAttributes(simpleAlSVisAttWNICU);
	box4LV -> SetVisAttributes(simpleAlSVisAttPb);
	box5LV -> SetVisAttributes(simpleAlSVisAttPb);
	box6LV -> SetVisAttributes(simpleAlSVisAttWNICU);

	// Region for cuts
	box1LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(box1LV);
	box2LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(box2LV);
	box3LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(box3LV);
	box4LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(box4LV);
	box5LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(box5LV);
	box6LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(box6LV);

	bCreated = true;
	return bCreated;
}
bool CML2AccSaturn::Jaw2X()
{
	bool bCreated = false;

	G4double boxSide;
	G4double box1Height, box1Z;
	G4double box2Height, box2Z;
	G4double box3Height, box3Z;
	G4double box4Height, box4Z;
	G4double box5Height, box5Z;
	G4double box6Height, box6Z;

	boxSide=101./2.*mm; // upper principal box

	box1Height=3./2.*mm;
	box1Z=275.50*mm+box1Height;

	box2Height=35./2.*mm;
	box2Z=box1Z+box1Height+box2Height;

	box3Height=35./2.*mm;
	box3Z=box2Z+box2Height+box3Height;

	box4Height=27./2.*mm;
	box4Z=box3Z+box3Height+box4Height;

	box5Height=10./2.*mm;
	box5Z=472.50*mm+box5Height;

	box6Height=5./2.*mm;
	box6Z=box5Z+box5Height+box6Height;

	G4Material *elPb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
//	G4Material *XC10 = SetMaterials("XC10");
//	G4Material *WNICU = SetMaterials("Denal(WNICU)");
	G4Material *XC10 = getMaterial("xc10");
	G4Material *WNICU = getMaterial("wnicu");

	G4Region *regVol;
	G4VisAttributes* simpleAlSVisAttPb;
	G4VisAttributes* simpleAlSVisAttXC10;
	G4VisAttributes* simpleAlSVisAttWNICU;
	// Region for cuts
	regVol= new G4Region("Jaws2XR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(2.*cm);
	regVol->SetProductionCuts(cuts);

	// Physical volumes

	G4Box *box1 = new G4Box("Jaws2XBox1", boxSide, boxSide, box1Height);
	G4Box *box2 = new G4Box("Jaws2XBox2", boxSide, boxSide, box2Height);
	G4Box *box3 = new G4Box("Jaws2XBox3", boxSide, boxSide, box3Height);
	G4Box *box4 = new G4Box("Jaws2XBox4", boxSide, boxSide, box4Height);
	G4Box *box5 = new G4Box("Jaws2XBox5", boxSide, boxSide, box5Height);
	G4Box *box6 = new G4Box("Jaws2XBox6", boxSide, boxSide, box6Height);
	G4LogicalVolume* box1LV = new G4LogicalVolume(box1, XC10, "Jaws2XLV1",0,0,0);
	G4LogicalVolume* box2LV = new G4LogicalVolume(box2, elPb, "Jaws2XLV2",0,0,0);
	G4LogicalVolume* box3LV = new G4LogicalVolume(box3, WNICU, "Jaws2XLV3",0,0,0);
	G4LogicalVolume* box4LV = new G4LogicalVolume(box4, elPb, "Jaws2XLV4",0,0,0);
	G4LogicalVolume* box5LV = new G4LogicalVolume(box5, elPb, "Jaws2XLV5",0,0,0);
	G4LogicalVolume* box6LV = new G4LogicalVolume(box6, WNICU, "Jaws2XLV6",0,0,0);

	G4ThreeVector centre;
	G4RotationMatrix *cRotationId = new G4RotationMatrix();
	G4RotationMatrix *cRotation = new G4RotationMatrix();

	*cRotation=*cRotationId;
	centre.set(boxSide,0.,box1Z);
	SetJawAperture(2,centre,G4ThreeVector(boxSide ,boxSide ,box1Height),cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws2XPV1", box1LV, PVWorld, false, 0);

	*cRotation=*cRotationId;
	centre.set(boxSide,0.,box2Z);
	SetJawAperture(2,centre,G4ThreeVector(boxSide ,boxSide ,box2Height),cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws2XPV2", box2LV,  PVWorld, false, 0);

	*cRotation=*cRotationId;
	centre.set(boxSide,0.,box3Z);
	SetJawAperture(2,centre,G4ThreeVector(boxSide ,boxSide ,box3Height),cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws2XPV3", box3LV, PVWorld, false, 0);

	*cRotation=*cRotationId;
	centre.set(boxSide,0.,box4Z);
	SetJawAperture(2,centre,G4ThreeVector(boxSide ,boxSide ,box4Height),cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws2XPV4", box4LV, PVWorld, false, 0);

	*cRotation=*cRotationId;
	centre.set(boxSide,0.,box5Z);
	SetJawAperture(2,centre,G4ThreeVector(boxSide ,boxSide ,box5Height),cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws2XPV5", box5LV, PVWorld, false, 0);

	*cRotation=*cRotationId;
	centre.set(boxSide,0.,box6Z);
	SetJawAperture(2,centre,G4ThreeVector(boxSide ,boxSide ,box6Height),cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws2XPV6", box6LV, PVWorld, false, 0);

	// Visualization
	simpleAlSVisAttPb = new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAttPb -> SetVisibility(true);
	simpleAlSVisAttPb -> SetForceSolid(true);
	simpleAlSVisAttXC10 = new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAttXC10 -> SetVisibility(true);
	simpleAlSVisAttXC10 -> SetForceSolid(true);
	simpleAlSVisAttWNICU = new G4VisAttributes(G4Colour::Red());
	simpleAlSVisAttWNICU -> SetVisibility(true);
	simpleAlSVisAttWNICU -> SetForceSolid(true);

	box1LV -> SetVisAttributes(simpleAlSVisAttXC10);
	box2LV -> SetVisAttributes(simpleAlSVisAttPb);
	box3LV -> SetVisAttributes(simpleAlSVisAttWNICU);
	box4LV -> SetVisAttributes(simpleAlSVisAttPb);
	box5LV -> SetVisAttributes(simpleAlSVisAttPb);
	box6LV -> SetVisAttributes(simpleAlSVisAttWNICU);

	// Region for cuts
	box1LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box1LV);
	box2LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box2LV);
	box3LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box3LV);
	box4LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box4LV);
	box5LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box5LV);
	box6LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box6LV);

	bCreated = true;
	return bCreated;
}
bool CML2AccSaturn::Jaw1Y()
{
	bool bCreated = false;

	G4double boxSide;
	G4double box1Height, box1Z;
	G4double box2Height, box2Z;
	G4double box3Height, box3Z;
	G4double box4Height, box4Z;

	boxSide=101./2.*mm; // upper principal box

	box1Height = 15./2.*mm;
	box1Z = 248.50*mm+box1Height;

	box2Height = 31./2.*mm;
	box2Z = 380.5+box2Height;

	box3Height = 35./2.*mm;
	box3Z = box2Z+box2Height+box3Height;

	box4Height = 21./2.*mm;
	box4Z = box3Z+box3Height+box4Height;

	G4Material *elPb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
//	G4Material *XC10 = SetMaterials("XC10");
//	G4Material *WNICU = SetMaterials("Denal(WNICU)");
	G4Material *XC10 = getMaterial("xc10");
	G4Material *WNICU = getMaterial("wnicu");

	G4Region *regVol;
	G4VisAttributes* simpleAlSVisAttPb;
	G4VisAttributes* simpleAlSVisAttXC10;
	G4VisAttributes* simpleAlSVisAttWNICU;
	// Region for cuts
	regVol= new G4Region("Jaws1YR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(2.*cm);
	regVol->SetProductionCuts(cuts);

	// Physical volumes

	G4Box *box1 = new G4Box("Jaws1YBox1", boxSide, boxSide, box1Height);
	G4Box *box2 = new G4Box("Jaws1YBox2", boxSide, boxSide, box2Height);
	G4Box *box3 = new G4Box("Jaws1YBox3", boxSide, boxSide, box3Height);
	G4Box *box4 = new G4Box("Jaws1YBox4", boxSide, boxSide, box4Height);
	G4LogicalVolume* box1LV = new G4LogicalVolume(box1, XC10, "Jaws1YLV1",0,0,0);
	G4LogicalVolume* box2LV = new G4LogicalVolume(box2, elPb, "Jaws1YLV2",0,0,0);
	G4LogicalVolume* box3LV = new G4LogicalVolume(box3, WNICU, "Jaws1YLV3",0,0,0);
	G4LogicalVolume* box4LV = new G4LogicalVolume(box4, elPb, "Jaws1YLV4",0,0,0);

	G4ThreeVector centre;
	G4RotationMatrix *cRotationId = new G4RotationMatrix();
	G4RotationMatrix *cRotation = new G4RotationMatrix();

	*cRotation = *cRotationId;
	centre.set(0.,boxSide,box1Z);
	SetJawAperture(3,centre,G4ThreeVector(boxSide ,boxSide ,box1Height),cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws1YPV1", box1LV, PVWorld, false, 0);

	*cRotation=*cRotationId;
	centre.set(0.,boxSide,box2Z);
	SetJawAperture(3,centre,G4ThreeVector(boxSide ,boxSide ,box2Height),cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws1YPV2", box2LV, PVWorld, false, 0);

	*cRotation=*cRotationId;
	centre.set(0.,boxSide,box3Z);
	SetJawAperture(3,centre,G4ThreeVector(boxSide ,boxSide ,box3Height),cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws1YPV3", box3LV, PVWorld, false, 0);

	*cRotation=*cRotationId;
	centre.set(0.,boxSide,box4Z);
	SetJawAperture(3,centre,G4ThreeVector(boxSide ,boxSide ,box4Height),cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws1YPV4", box4LV, PVWorld, false, 0);


	// Visualization
	simpleAlSVisAttPb= new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAttPb->SetVisibility(true);
	simpleAlSVisAttPb->SetForceSolid(true);
	simpleAlSVisAttXC10= new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAttXC10->SetVisibility(true);
	simpleAlSVisAttXC10->SetForceSolid(true);
	simpleAlSVisAttWNICU= new G4VisAttributes(G4Colour::Red());
	simpleAlSVisAttWNICU->SetVisibility(true);
	simpleAlSVisAttWNICU->SetForceSolid(true);

	box1LV->SetVisAttributes(simpleAlSVisAttWNICU);
	box2LV->SetVisAttributes(simpleAlSVisAttPb);
	box3LV->SetVisAttributes(simpleAlSVisAttWNICU);
	box4LV->SetVisAttributes(simpleAlSVisAttPb);

	// Region for cuts
	box1LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(box1LV);
	box2LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(box2LV);
	box3LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(box3LV);
	box4LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(box4LV);

	bCreated = true;
	return bCreated;
}
bool CML2AccSaturn::Jaw2Y()
{
	bool bCreated = false;

	G4double boxSide, box1Height, box1Z, box2Height, box2Z, box3Height, box3Z, box4Height, box4Z;

	boxSide = 101./2.*mm; // upper principal box

	box1Height = 15./2.*mm;
	box1Z = 248.50*mm+box1Height;

	box2Height = 31./2.*mm;
	box2Z = 380.5+box2Height;

	box3Height = 35./2.*mm;
	box3Z = box2Z+box2Height+box3Height;

	box4Height = 21./2.*mm;
	box4Z = box3Z+box3Height+box4Height;

	G4Material *elPb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
//	G4Material *XC10 = SetMaterials("XC10");
//	G4Material *WNICU = SetMaterials("Denal(WNICU)");
	G4Material *XC10 = getMaterial("xc10");
	G4Material *WNICU = getMaterial("wnicu");

	G4Region *regVol;
	G4VisAttributes* simpleAlSVisAttPb;
	G4VisAttributes* simpleAlSVisAttXC10;
	G4VisAttributes* simpleAlSVisAttWNICU;
	// Region for cuts
	regVol= new G4Region("Jaws2YR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts -> SetProductionCut(2.*cm);
	regVol -> SetProductionCuts(cuts);

	// Physical volumes
	G4Box *box1 = new G4Box("Jaws2YBox1", boxSide, boxSide, box1Height);
	G4Box *box2 = new G4Box("Jaws2YBox2", boxSide, boxSide, box2Height);
	G4Box *box3 = new G4Box("Jaws2YBox3",boxSide, boxSide, box3Height);
	G4Box *box4 = new G4Box("Jaws2YBox4", boxSide, boxSide, box4Height);
	G4LogicalVolume* box1LV = new G4LogicalVolume(box1, XC10, "Jaws2YLV1",0,0,0);
	G4LogicalVolume* box2LV = new G4LogicalVolume(box2, elPb, "Jaws2YLV2",0,0,0);
	G4LogicalVolume* box3LV = new G4LogicalVolume(box3, WNICU, "Jaws2YLV3",0,0,0);
	G4LogicalVolume* box4LV = new G4LogicalVolume(box4, elPb, "Jaws2YLV4",0,0,0);

	G4ThreeVector centre;
	G4RotationMatrix *cRotationId = new G4RotationMatrix();
	G4RotationMatrix *cRotation = new G4RotationMatrix();

	*cRotation = *cRotationId;
	centre.set(0.,boxSide,box1Z);
	SetJawAperture(4,centre,G4ThreeVector(boxSide ,boxSide ,box1Height),cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws2YPV1", box1LV, PVWorld, false, 0);

	*cRotation = *cRotationId;
	centre.set(0.,boxSide,box2Z);
	SetJawAperture(4,centre,G4ThreeVector(boxSide ,boxSide ,box2Height),cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws2YPV2", box2LV, PVWorld, false, 0);

	*cRotation = *cRotationId;
	centre.set(0.,boxSide,box3Z);
	SetJawAperture(4,centre,G4ThreeVector(boxSide ,boxSide ,box3Height),cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws2YPV3", box3LV, PVWorld, false, 0);

	*cRotation = *cRotationId;
	centre.set(0.,boxSide,box4Z);
	SetJawAperture(4,centre,G4ThreeVector(boxSide ,boxSide ,box4Height),cRotation);
	new G4PVPlacement(cRotation, centre, "Jaws2YPV4", box4LV, PVWorld, false, 0);


	// Visualization
	simpleAlSVisAttPb = new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAttPb -> SetVisibility(true);
	simpleAlSVisAttPb -> SetForceSolid(true);
	simpleAlSVisAttXC10 = new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAttXC10 -> SetVisibility(true);
	simpleAlSVisAttXC10 -> SetForceSolid(true);
	simpleAlSVisAttWNICU = new G4VisAttributes(G4Colour::Red());
	simpleAlSVisAttWNICU -> SetVisibility(true);
	simpleAlSVisAttWNICU -> SetForceSolid(true);

	box1LV -> SetVisAttributes(simpleAlSVisAttWNICU);
	box2LV -> SetVisAttributes(simpleAlSVisAttPb);
	box3LV -> SetVisAttributes(simpleAlSVisAttWNICU);
	box4LV -> SetVisAttributes(simpleAlSVisAttPb);

	// Region for cuts
	box1LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box1LV);
	box2LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box2LV);
	box3LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box3LV);
	box4LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box4LV);

	bCreated = true;
	return bCreated;
}

