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


#include "ML2Acc1.hh"
#include "ML2Acc1Messenger.hh"
#include "ML2Accelerator.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

using namespace std;

CML2Acc1::CML2Acc1()
{
	PVWorld = 0;
	acc1Messenger = new CML2Acc1Messenger(this);
	buildMaterial_SSteel1();
}

CML2Acc1::~CML2Acc1(void)
{
}
CML2Acc1* CML2Acc1::instance = 0;

CML2Acc1* CML2Acc1::GetInstance(void)
{
	if (instance == 0)
	{
		instance = new CML2Acc1();
	}
return instance;
}

void CML2Acc1::writeInfo()
{
	G4cout <<"\n\n\tnominal beam energy: "<<idEnergy << G4endl;
	G4cout <<"\tJaw X aperture: 1) "<< jaw1XAperture/mm<<"[mm]\t2) " << jaw2XAperture/mm<< " [mm]"<< G4endl;
	G4cout <<"\tJaw Y aperture: 1) "<< jaw1YAperture/mm<<"[mm]\t2) " << jaw2YAperture/mm<< " [mm]\n"<< G4endl;
}

void CML2Acc1::buildMaterial_SSteel1()
{

	G4Element* elFe = new G4Element("Iron", "Fe", 26., 55.85*g/mole);
	G4Element* elS  = new G4Element("Sulfur", "S", 16., 32.064*g/mole);
	G4Element* elMn = new G4Element("Manganese", "Mn", 25., 54.94*g/mole);
	G4Element* elC  = new G4Element("Carbon", "C", 6., 12.011*g/mole);

	steel1 = new G4Material("steel1", 7.76 *g/cm3, 4);
	steel1 -> AddElement(elFe, 0.935);
        steel1 -> AddElement(elS, 0.01);
	steel1 -> AddElement(elMn, 0.05);
	steel1 -> AddElement(elC, 0.005);
}

void CML2Acc1::Construct(G4VPhysicalVolume *PWorld, G4double iso)
{
	PVWorld = PWorld;
	setIsoCentre(iso);
	target();
	vacuumWindow();
	ionizationChamber();
	flatteningFilter();
	mirror();
	primaryCollimator();
	MLC();
	Jaw1X();
	Jaw2X();
	Jaw1Y();
	Jaw2Y();
}

bool CML2Acc1::target()
{
	bool bCreated = false;
	switch (idEnergy)
	{
		case 6:
			//    materials
			G4Material* Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
			G4Material* W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");

			//    volumes
			//    beam line along z axis
			//------------------------target 6MV------------------------
			G4double targetADim_x = 0.6*cm;
			G4double targetADim_y = 0.6*cm;
			G4double targetADim_z = 0.04445*cm;
			G4Box* targetA_box = new G4Box("targetA_box",targetADim_x,targetADim_y,targetADim_z);
			G4LogicalVolume *targetA_log = new G4LogicalVolume(targetA_box,W,"targetA_log",0,0,0);
			G4double targetAPos_x = 0.0*m;
			G4double targetAPos_y = 0.0*m;
			G4double targetAPos_z = 0.20055*cm;
			new G4PVPlacement(0,
					G4ThreeVector(targetAPos_x,targetAPos_y,targetAPos_z),
					"targetA",targetA_log,PVWorld,false,0);

			G4double targetBDim_x = 0.6*cm;
			G4double targetBDim_y = 0.6*cm;
			G4double targetBDim_z = 0.07874*cm;
			G4Box* targetB_box = new G4Box("targetB_box",targetBDim_x,targetBDim_y,targetBDim_z);
			G4LogicalVolume *targetB_log = new G4LogicalVolume(targetB_box,Cu,"targetB_log",0,0,0);
			G4double targetBPos_x = 0.0*m;
			G4double targetBPos_y = 0.0*m;
			G4double targetBPos_z = 0.07736*cm;
			new G4PVPlacement(0,
					G4ThreeVector(targetBPos_x,targetBPos_y,targetBPos_z),
					"targetB",targetB_log,PVWorld,false,0);

			// ***********  REGIONS for CUTS
			G4Region *regVol;
			regVol= new G4Region("targetR");
			G4ProductionCuts* cuts = new G4ProductionCuts;
			cuts->SetProductionCut(0.1*cm);
			regVol->SetProductionCuts(cuts);

			targetA_log -> SetRegion(regVol);
			regVol -> AddRootLogicalVolume(targetA_log);
			targetB_log -> SetRegion(regVol);
			regVol -> AddRootLogicalVolume(targetB_log);

			// Visualization attributes
			G4VisAttributes* simpleWSVisAtt, *simpleCuSVisAtt;
			simpleWSVisAtt = new G4VisAttributes(G4Colour::Magenta());
			simpleWSVisAtt -> SetVisibility(true);
			simpleCuSVisAtt = new G4VisAttributes(G4Colour::Cyan());
			simpleCuSVisAtt -> SetVisibility(true);
			targetA_log -> SetVisAttributes(simpleWSVisAtt);
			targetB_log -> SetVisAttributes(simpleCuSVisAtt);

			bCreated = true;
			return bCreated;
			break;
	}
	return false;
}
bool CML2Acc1::primaryCollimator()
{
	//    materials
	G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	G4Material* W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
 
	//---------rotation matrix first collimator --------
	G4RotationMatrix*  rotateMatrix=new G4RotationMatrix();
	rotateMatrix->rotateX(180.0*deg);

	//--------------------upper collimator----------------
	G4double innerRadiusOfTheTubeEx = 1.0*cm;
	G4double outerRadiusOfTheTubeEx = 8.*cm;
	G4double hightOfTheTubeEx = 3.0*cm;
	G4double startAngleOfTheTubeEx = 0.*deg;
	G4double spanningAngleOfTheTubeEx = 360.*deg;
	G4Tubs* UpperCollimator = new G4Tubs("UpperCollimator",
			innerRadiusOfTheTubeEx,
			outerRadiusOfTheTubeEx,
			hightOfTheTubeEx,
			startAngleOfTheTubeEx,
			spanningAngleOfTheTubeEx);
	G4LogicalVolume *UpperCollimator_log = new G4LogicalVolume(UpperCollimator, W, "UpperCollimator_log", 0, 0, 0);

	G4double UpperCollimatorPosX = 0.*cm;
	G4double UpperCollimatorPosY = 0.*cm;
	G4double UpperCollimatorPosZ = -1.*cm;
	new G4PVPlacement(0,
			G4ThreeVector(UpperCollimatorPosX, UpperCollimatorPosY, UpperCollimatorPosZ),
			"UpperCollimator",
			UpperCollimator_log,
			PVWorld,
			false,
			0);

 
	//--------------------lower collimator----------------

	G4double  pRmin1 = 0.*cm;
	G4double  pRmax1 = 0.5*cm;
	G4double  pRmin2 = 0.*cm;
	G4double  pRmax2 = 1.7658592*cm;
	G4double  hightOfTheCone =3.2*cm;
	G4double  startAngleOfTheCone = 0.*deg;
	G4double  spanningAngleOfTheCone = 360.*deg;

	G4Cons* collim_cone = new G4Cons("collim_cone",pRmin1,pRmax1,pRmin2,
		  pRmax2,hightOfTheCone,startAngleOfTheCone,
		  spanningAngleOfTheCone);
	G4LogicalVolume *collim_log = new G4LogicalVolume(collim_cone,Vacuum,"collim_log",0,0,0);


	G4double innerRadiusOfTheTube = 0.*cm;
	G4double outerRadiusOfTheTube = 8.*cm;
	G4double hightOfTheTube = 3.1*cm;
	G4double startAngleOfTheTube = 0.*deg;
	G4double spanningAngleOfTheTube = 360.*deg;
	G4Tubs* tracker_tube = new G4Tubs("tracker_tube",innerRadiusOfTheTube,
			outerRadiusOfTheTube,hightOfTheTube,
			startAngleOfTheTube,spanningAngleOfTheTube);

	G4SubtractionSolid* CylMinusCone = new G4SubtractionSolid("Cyl-Cone",
			tracker_tube,collim_cone);
	G4LogicalVolume *CylMinusCone_log = new G4LogicalVolume(CylMinusCone,W,"CylminusCone_log",0,0,0);
	G4double CminusCPos_x = 0.*cm;
	G4double CminusCPos_y = 0.*cm;
	G4double CminusCPos_z = +6.2*cm;
	new G4PVPlacement(rotateMatrix,
			G4ThreeVector(CminusCPos_x, CminusCPos_y, CminusCPos_z),
			"CylMinusCone",
			CylMinusCone_log,
			PVWorld,
			false,
			0);

	//--------- Visualization attributes -------------------------------
	G4VisAttributes* simpleTungstenWVisAtt= new G4VisAttributes(G4Colour::Magenta());
	simpleTungstenWVisAtt->SetVisibility(true);
	collim_log->SetVisAttributes(simpleTungstenWVisAtt);

	CylMinusCone_log->SetVisAttributes(simpleTungstenWVisAtt);
	UpperCollimator_log->SetVisAttributes(simpleTungstenWVisAtt);

	// ***********  REGIONS for CUTS
	G4Region *regVol;
	regVol= new G4Region("PrymCollR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*cm);
	regVol->SetProductionCuts(cuts);

	collim_log->SetRegion(regVol);
	regVol->AddRootLogicalVolume(collim_log);

	CylMinusCone_log->SetRegion(regVol);
	regVol->AddRootLogicalVolume(CylMinusCone_log);

	UpperCollimator_log->SetRegion(regVol);
	regVol->AddRootLogicalVolume(UpperCollimator_log);

	return true;
}
bool CML2Acc1::vacuumWindow()
{
	bool bCreated = false;
	G4Material *Be = G4NistManager::Instance()->FindOrBuildMaterial("G4_Be");
	G4Region *regVol;
	G4VisAttributes* simpleAlSVisAtt;
	// Region for cuts
	regVol= new G4Region("BeWindow");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*cm);
	regVol->SetProductionCuts(cuts);

	G4Tubs* BeWTube = new G4Tubs("BeWindowTube", 0., 36.*mm, 0.2*mm, 0.*deg, 360.*deg);
	G4LogicalVolume *BeWTubeLV = new G4LogicalVolume(BeWTube, Be, "BeWTubeLV", 0, 0, 0);
	new G4PVPlacement(0, G4ThreeVector(0.,0.,100.*mm), "BeWTubePV", BeWTubeLV,
			PVWorld, false, 0);

	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Yellow());
	simpleAlSVisAtt -> SetVisibility(true);
	BeWTubeLV -> SetVisAttributes(simpleAlSVisAtt);
	BeWTubeLV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(BeWTubeLV);

	bCreated = true;
	return bCreated;
}
bool CML2Acc1::flatteningFilter()
{
	bool bCreated = false;
	switch (idEnergy)
	{
	case 6:
		G4double z0, h0;
		G4ThreeVector centre, halSize;
		G4Material *Cu=G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
		// Region for cuts
		G4Region *regVol;
		regVol= new G4Region("flatfilterR");
		G4ProductionCuts* cuts = new G4ProductionCuts;
		cuts->SetProductionCut(0.5*cm);
		regVol->SetProductionCuts(cuts);

		G4VisAttributes* simpleAlSVisAtt;

		// one
		z0=130.0*mm;
		h0=5.0/2.*cm;
		centre.set(0.,0.,z0);
		G4Cons *FFL1A_1Cone = new G4Cons("FFL1A_1", 0.*cm, 0.3*cm, 0.*cm, 5.*cm, h0, 0.*deg, 360.*deg);
		G4LogicalVolume *FFL1A_1LV = new G4LogicalVolume(FFL1A_1Cone, Cu, "FFL1A_1LV", 0, 0, 0);
		new G4PVPlacement(0, centre, "FFL1A_1PV", FFL1A_1LV, PVWorld, false, 0);

		// two
		z0+=h0;
		h0=0.081/2.*cm;
		z0+=h0;
		centre.setZ(z0);
		z0+=h0;
		G4Tubs *FFL2_1Tube = new G4Tubs("FFL6_1", 0.*cm, 2.5*cm, h0, 0.*deg, 360.*deg);
		G4LogicalVolume *FFL2_1LV = new G4LogicalVolume(FFL2_1Tube, Cu, "FFL2_1LV", 0, 0, 0);
		new G4PVPlacement(0, centre, "FFL2_1PV", FFL2_1LV, PVWorld, false, 0);

		simpleAlSVisAtt= new G4VisAttributes(G4Colour::Red());
		simpleAlSVisAtt->SetVisibility(true);
		FFL1A_1LV->SetVisAttributes(simpleAlSVisAtt);
		FFL2_1LV->SetVisAttributes(simpleAlSVisAtt);

		FFL1A_1LV->SetRegion(regVol);
		FFL2_1LV->SetRegion(regVol);

		regVol->AddRootLogicalVolume(FFL1A_1LV);
		regVol->AddRootLogicalVolume(FFL2_1LV);
		bCreated = true;
		break;
		}

	return bCreated;
}
bool CML2Acc1::ionizationChamber()
{
	bool bCreated = false;

	G4Material *material=G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
	G4VisAttributes* simpleAlSVisAtt;
	// Region for cuts
	G4Region *regVol;
	regVol= new G4Region("ionizationChamber");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*cm);
	regVol->SetProductionCuts(cuts);

	G4Tubs* ICTubeW = new G4Tubs("ionizationChamberTube", 0., 2.*2.54*10.*mm, 0.016*25.4*mm, 0.*deg, 360.*deg);
	G4Tubs* ICTubeP = new G4Tubs("ionizationChamberTube", 0., 2.*2.54*10.*mm, 0.010*25.4*mm, 0.*deg, 360.*deg);

	G4ThreeVector centre;
	// W1
	centre.set(0.,0.,157.*mm);
	G4LogicalVolume *PCUTubeW1LV = new G4LogicalVolume(ICTubeW, material, "ionizationChamberTubeW1LV", 0, 0, 0);
	new G4PVPlacement(0, centre, "ionizationChamberTubeW1PV", PCUTubeW1LV, PVWorld, false, 0);
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAtt->SetVisibility(true);
	PCUTubeW1LV->SetVisAttributes(simpleAlSVisAtt);
	PCUTubeW1LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(PCUTubeW1LV);

	// P1
	centre.set(0.,0.,158.*mm);
	G4LogicalVolume *PCUTubeP1LV = new G4LogicalVolume(ICTubeP, material, "ionizationChamberTubeP1LV", 0, 0, 0);
	new G4PVPlacement(0, centre, "ionizationChamberTubeP1PV", PCUTubeP1LV, PVWorld, false, 0);
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Yellow());
	simpleAlSVisAtt->SetVisibility(true);
	PCUTubeP1LV->SetVisAttributes(simpleAlSVisAtt);
	PCUTubeP1LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(PCUTubeP1LV);

	// W2
	centre.set(0.,0.,159.*mm);
	G4LogicalVolume *PCUTubeW2LV = new G4LogicalVolume(ICTubeW, material, "ionizationChamberTubeW2LV", 0, 0, 0);
	new G4PVPlacement(0, centre, "ionizationChamberTubeW2PV", PCUTubeW2LV, PVWorld, false, 0);
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAtt->SetVisibility(true);
	PCUTubeW2LV->SetVisAttributes(simpleAlSVisAtt);
	PCUTubeW2LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(PCUTubeW2LV);

	// P2
	centre.set(0.,0.,160.*mm);
	G4LogicalVolume *PCUTubeP2LV = new G4LogicalVolume(ICTubeP, material, "ionizationChamberTubeP2LV", 0, 0, 0);
	new G4PVPlacement(0, centre, "ionizationChamberTubeP2PV", PCUTubeP2LV, PVWorld, false, 0);
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Yellow());
	simpleAlSVisAtt->SetVisibility(true);
	PCUTubeP2LV->SetVisAttributes(simpleAlSVisAtt);
	PCUTubeP2LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(PCUTubeP2LV);

	// W3
	centre.set(0.,0.,161.*mm);
	G4LogicalVolume *PCUTubeW3LV = new G4LogicalVolume(ICTubeW, material, "ionizationChamberTubeW3LV", 0, 0, 0);
	new G4PVPlacement(0, centre, "ionizationChamberTubeW3PV", PCUTubeW3LV, PVWorld, false, 0);
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAtt->SetVisibility(true);
	PCUTubeW3LV->SetVisAttributes(simpleAlSVisAtt);
	PCUTubeW3LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(PCUTubeW3LV);

	// P3
	centre.set(0.,0.,162.*mm);
	G4LogicalVolume *PCUTubeP3LV = new G4LogicalVolume(ICTubeP, material, "ionizationChamberTubeP3LV", 0, 0, 0);
	new G4PVPlacement(0, centre, "ionizationChamberTubeP3PV", PCUTubeP3LV, PVWorld, false, 0);
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Yellow());
	simpleAlSVisAtt->SetVisibility(true);
	PCUTubeP3LV->SetVisAttributes(simpleAlSVisAtt);
	PCUTubeP3LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(PCUTubeP3LV);

	bCreated = true;
	return bCreated;
}
bool CML2Acc1::mirror()
{
	bool bCreated=false;
	G4Material *MYLAR = G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR");
	G4VisAttributes* simpleAlSVisAtt;
	// Region for cuts
	G4Region *regVol;
	regVol = new G4Region("Mirror");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts -> SetProductionCut(0.1*cm);
	regVol -> SetProductionCuts(cuts);

	G4Tubs* MirrorTube = new G4Tubs("MirrorTube", 0., 63.*mm, .5*mm, 0.*deg, 360.*deg);
	G4LogicalVolume *MirrorTubeLV = new G4LogicalVolume(MirrorTube, MYLAR, "MirrorTubeLV", 0, 0, 0);
	G4RotationMatrix *cRotation = new G4RotationMatrix();
	cRotation -> rotateY(12.0*deg);
	new G4PVPlacement(cRotation, G4ThreeVector(0., 0., 175.*mm), "MirrorTubePV", MirrorTubeLV,PVWorld, false, 0);

	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAtt -> SetVisibility(true);
	MirrorTubeLV -> SetVisAttributes(simpleAlSVisAtt);
	MirrorTubeLV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(MirrorTubeLV);
	bCreated = true;
	return bCreated;
}

void CML2Acc1::SetJawAperture(G4int idJaw, G4ThreeVector &centre, G4ThreeVector halfSize, G4RotationMatrix *cRotation)
{
	G4double theta, x, y, z, dx, dy, dz;
	x=centre.getX();
	y=centre.getY();
	z=centre.getZ();

	dx=halfSize.getX();
	dy=halfSize.getY();
	dz=halfSize.getZ();
	
	switch (idJaw)
	{
	case 1: //idJaw1XV2100:
		theta = fabs(atan(jaw1XAperture/isoCentre));
		centre.set(z*sin(theta)+dx*cos(theta), y, z*cos(theta)-dx*sin(theta));
		cRotation->rotateY(-theta);
		halfSize.set(fabs(dx*cos(theta)+dz*sin(theta)), fabs(dy), fabs(dz*cos(theta)+dx*sin(theta)));
		break;

	case 2: //idJaw2XV2100:
		theta = fabs(atan(jaw2XAperture/isoCentre));
		centre.set(-(z*sin(theta)+dx*cos(theta)), y, z*cos(theta)-dx*sin(theta));
		cRotation->rotateY(theta);
		halfSize.set(fabs(dx*cos(theta)+dz*sin(theta)), fabs(dy), fabs(dz*cos(theta)+dx*sin(theta)));
		break;

	case 3: //idJaw1YV2100:
		theta = fabs(atan(jaw1YAperture/isoCentre));
		centre.set(x, z*sin(theta)+dy*cos(theta), z*cos(theta)-dy*sin(theta));
		cRotation->rotateX(theta);
		halfSize.set(fabs(dx), fabs(dy*cos(theta)+dz*sin(theta)), fabs(dz*cos(theta)+dy*sin(theta)));
		break;

	case 4: //idJaw2YV2100:
		theta = fabs(atan(jaw2YAperture/isoCentre));
		centre.set(x, -(z*sin(theta)+dy*cos(theta)), z*cos(theta)-dy*sin(theta));
		cRotation->rotateX(-theta);
		halfSize.set(fabs(dx), fabs(dy*cos(theta)+dz*sin(theta)), fabs(dz*cos(theta)+dy*sin(theta)));
		break;
	}
}



bool CML2Acc1::Jaw1X()
{
	bool bCreated = false;

	G4String name = "Jaws1X";
	G4Box *box;
	G4LogicalVolume *logVol;
	G4VisAttributes* simpleAlSVisAtt;

	G4ThreeVector centre, halfSize;
	G4RotationMatrix *cRotation = new G4RotationMatrix();

	centre.set(0.,0.,(320.+80./2.)*mm);
	halfSize.set(45.*mm, 93.*mm, 78./2.*mm);
	box = new G4Box(name+"Box", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	logVol = new G4LogicalVolume(box, steel1, name+"LV", 0, 0, 0);
	SetJawAperture(1, centre, halfSize, cRotation);
	new G4PVPlacement(
			cRotation,
			centre,
			name+"PV",
			logVol,
			PVWorld,
			false,
			0);

	// Region for cuts
	G4Region *regVol;
	regVol= new G4Region(name+"R");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(2.*cm);
	regVol->SetProductionCuts(cuts);
	logVol->SetRegion(regVol);
	regVol->AddRootLogicalVolume(logVol);

	// Visibility
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAtt -> SetVisibility(true);
	logVol -> SetVisAttributes(simpleAlSVisAtt);

	bCreated = true;
	return bCreated;
}
bool CML2Acc1::Jaw2X()
{
	bool bCreated = false;

	G4String name = "Jaws2X";
	G4Box *box;
	G4LogicalVolume *logVol;
	G4VisAttributes* simpleAlSVisAtt;

	G4ThreeVector centre, halfSize;
	G4RotationMatrix *cRotation=new G4RotationMatrix();

	centre.set(0.,0.,(320.+80./2.)*mm);
	halfSize.set(45.*mm, 93.*mm, 78./2.*mm);
	box = new G4Box(name+"Box", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	logVol = new G4LogicalVolume(box, steel1, name+"LV", 0, 0, 0);
	SetJawAperture(2, centre, halfSize, cRotation);
	new G4PVPlacement(cRotation, centre, name+"PV", logVol, PVWorld, false, 0);

	// Region for cuts
	G4Region *regVol;
	regVol = new G4Region(name+"R");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts -> SetProductionCut(2.*cm);
	regVol -> SetProductionCuts(cuts);
	logVol -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(logVol);

	// Visibility
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Cyan());
	simpleAlSVisAtt -> SetVisibility(true);
	logVol->SetVisAttributes(simpleAlSVisAtt);

	bCreated = true;
	return bCreated;
}
bool CML2Acc1::Jaw1Y()
{
	bool bCreated = false;

	G4String name = "Jaws1Y";
	G4Box *box;
	G4LogicalVolume *logVol;
	G4VisAttributes* simpleAlSVisAtt;

	G4ThreeVector centre, halfSize;
	G4RotationMatrix *cRotation=new G4RotationMatrix();

	centre.set(0.,0.,(230.+80./2.)*mm);
	halfSize.set(93.*mm, 35.*mm, 78./2.*mm);
	box = new G4Box(name+"Box", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	logVol = new G4LogicalVolume(box, steel1, name+"LV", 0, 0, 0);
	SetJawAperture(3, centre, halfSize, cRotation);
	new G4PVPlacement(cRotation, centre, name+"PV", logVol, PVWorld, false, 0);

	// Region for cuts
	G4Region *regVol;
	regVol = new G4Region(name+"R");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts -> SetProductionCut(2.*cm);
	regVol -> SetProductionCuts(cuts);
	logVol -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(logVol);

	// Visibility
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Red());
	simpleAlSVisAtt -> SetVisibility(true);
	logVol -> SetVisAttributes(simpleAlSVisAtt);

	bCreated = true;
	return bCreated;
}
bool CML2Acc1::Jaw2Y()
{
	bool bCreated = false;

	G4String name = "Jaws2Y";
	G4Box *box;
	G4LogicalVolume *logVol;
	G4VisAttributes* simpleAlSVisAtt;

	G4ThreeVector centre, halfSize;
	G4RotationMatrix *cRotation=new G4RotationMatrix();

	centre.set(0.,0.,(230.+80./2.)*mm);
	halfSize.set(93.*mm, 35.*mm, 78./2.*mm);
	box = new G4Box(name+"Box", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	logVol = new G4LogicalVolume(box, steel1, name+"LV", 0, 0, 0);
	SetJawAperture(4, centre, halfSize, cRotation);
	new G4PVPlacement(cRotation, centre, name+"PV", logVol, PVWorld, false, 0);

	// Region for cuts
	G4Region *regVol;
	regVol = new G4Region(name+"R");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts -> SetProductionCut(2.*cm);
	regVol -> SetProductionCuts(cuts);
	logVol -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(logVol);

	// Visibility
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Magenta());
	simpleAlSVisAtt -> SetVisibility(true);
	logVol -> SetVisAttributes(simpleAlSVisAtt);

	bCreated = true;
	return bCreated;
}
bool CML2Acc1::MLC() //MultiLeaf Collimator
{
	bool bCreated = false;
	//    material
	G4Material *Fe = G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");
	G4VisAttributes* simpleAlSVisAtt;
	// Region for cuts
	G4Region *regVol;
	regVol = new G4Region("MLCR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts -> SetProductionCut(1.0*cm);
	regVol -> SetProductionCuts(cuts);

	G4ThreeVector boxSize;
	G4ThreeVector centreStart;
	centreStart.set(0.,0.,(330.+600.)/2.*mm);
	boxSize.set(6./2.*mm, 180./2.*mm, 50./2.*mm);

	// single leaf
	G4Box* boxLeaf = new G4Box("LeafBox", boxSize.getX(), boxSize.getY(), boxSize.getZ());

	G4LogicalVolume *leafLVA = new G4LogicalVolume(boxLeaf, Fe, "leafSolidALV", 0, 0, 0);
	G4LogicalVolume *leafLVB = new G4LogicalVolume(boxLeaf, Fe, "leafSolidBLV", 0, 0, 0);

	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Cyan());
	simpleAlSVisAtt -> SetVisibility(true);
	leafLVA -> SetVisAttributes(simpleAlSVisAtt);
	leafLVA -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(leafLVA);
	
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAtt -> SetVisibility(true);
	leafLVB -> SetVisAttributes(simpleAlSVisAtt);
	leafLVB -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(leafLVB);

	int i;
	G4String PVname;
	G4ThreeVector centre;
	int nhalfLeaves = (int)(vec_leavesA.size()/2.);
	centre = centreStart + G4ThreeVector(-nhalfLeaves*boxSize.getX(), 0.,0.);
	for (i = 1; i < (int)vec_leavesA.size(); i++)
	{
		G4String str = std::to_string(i);
		PVname = "leafA"+str;
		centre.setX(centre.getX()+boxSize.getX()*2.);
		centre.setY(-boxSize.getY()-vec_leavesA[i]);
		new G4PVPlacement(0, centre, PVname, leafLVA, PVWorld, false, i);
	}
	nhalfLeaves = (int)(vec_leavesB.size()/2.);
	centre = centreStart+G4ThreeVector(-nhalfLeaves*boxSize.getX(), 0.,0.);
	for (i = 1; i < (int)vec_leavesB.size(); i++)
	{
		G4String str = std::to_string(i);
		PVname = "leafB"+str;
		centre.setX(centre.getX()+boxSize.getX()*2.);
		centre.setY(+boxSize.getY()+vec_leavesB[i]);
		new G4PVPlacement(0, centre, PVname, leafLVB, PVWorld, false, i);
	}
	bCreated = true;
	return bCreated;
}

