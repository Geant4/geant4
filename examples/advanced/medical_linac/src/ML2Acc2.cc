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


#include "ML2Acc2.hh"
#include "ML2Acc2Messenger.hh"
#include "ML2Accelerator.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

using namespace std;

CML2Acc2::CML2Acc2()
{
	PVWorld = 0;
	acc2Messenger = new CML2Acc2Messenger(this);
}

CML2Acc2::~CML2Acc2(void)
{
}

CML2Acc2* CML2Acc2::instance = 0;

CML2Acc2* CML2Acc2::GetInstance(void)
{
	if (instance == 0)
	{
		instance = new CML2Acc2();
	}
	return instance;
}
void CML2Acc2::writeInfo()
{
	G4cout << "----------------------------------------------------------------" << G4endl;;
	G4cout << "Accelerator VARIAN LINAC 2100 " << G4endl;
	G4cout <<"\n\n\tnominal beam energy: "<<idEnergy << G4endl;
	G4cout << "\tdistance isocentre [mm]:"<< isoCentre/mm << G4endl;
	G4cout <<"\tJaw X aperture: 1) "<< jaw1XAperture/mm<<"[mm]\t2) " << jaw2XAperture/mm<< " [mm]"<< G4endl;
	G4cout <<"\tJaw Y aperture: 1) "<< jaw1YAperture/mm<<"[mm]\t2) " << jaw2YAperture/mm<< " [mm]\n"<< G4endl;
	if (vec_leavesA.size()>0)
	{
		G4cout << "\tvec_leaves A aperture [mm]" << G4endl;;
		for (int i=0; i< (int)vec_leavesA.size(); i++)
		{
			G4cout<<"\t" <<i <<") "<< vec_leavesA[i]/mm << G4endl;
		}
	}
	else
	{
		G4cout << "\tNo vec_leaves A" << G4endl;
	}
	if (vec_leavesB.size()>0)
	{
		G4cout << "\tvec_leaves B aperture [mm]" << G4endl;
		for (int i=0; i< (int)vec_leavesB.size(); i++)
		{
			G4cout<<"\t" <<i <<") "<< vec_leavesB[i]/mm << G4endl;
		}
	}
	else
	{
		G4cout << "\tNo vec_leaves B" << G4endl;
	}
	G4cout << "______________________________________________________________" << G4endl;
}

void CML2Acc2::Construct(G4VPhysicalVolume *PWorld, G4double iso)
{
	setIsoCentre(iso);
	PVWorld = PWorld;
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
void CML2Acc2::SetJawAperture(G4int idJaw, G4ThreeVector &centre, G4ThreeVector halfSize, G4RotationMatrix *cRotation)
{
	using namespace std;
	G4double theta, x, y, z, dx, dy; //, dz, top;
//	G4double beta, R;
	x=centre.getX();
	y=centre.getY();
	z=centre.getZ();
//	top=z-78./2.;
	dx=halfSize.getX();
	dy=halfSize.getY();
//	dz=halfSize.getZ();
	
//	G4double p1x, p1y, p2x, p2y, dist;
	
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


bool CML2Acc2::target()
{
	bool bCreated = false;
	switch (idEnergy)
	{
	case 6:
	{
		// Physical and logical volumes
		G4Material *W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
		G4Box *targetABox = new G4Box("targetABox", 5.*mm, 5.*mm, (.035*25.4/2.)*mm);
		G4LogicalVolume *targetALV = new G4LogicalVolume(targetABox, W, "targetALV", 0, 0, 0);

		G4Material *Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
		G4Box *targetBBox = new G4Box("targetABox", 5.*mm, 5.*mm, (0.062*25.4/2.)*mm);
		G4LogicalVolume *targetBLV = new G4LogicalVolume(targetBBox, Cu, "targetBLV", 0, 0, 0);

		// specific translations for the various parts of the component
		new G4PVPlacement(
				0,
				G4ThreeVector(0.,0.,targetABox->GetZHalfLength()),
				"targetAPV",
				targetALV,
				PVWorld,
				false,
				0);

		new G4PVPlacement(
				0,
				G4ThreeVector(0.,0.,targetABox->GetZHalfLength()*2.+targetBBox->GetZHalfLength()),
				"targetBPV",
				targetBLV,
				PVWorld,
				false,
				0);

		// Region for cuts
		G4Region *regVol= new G4Region("targetR");
		G4ProductionCuts* cuts = new G4ProductionCuts;
		cuts->SetProductionCut(0.2*cm);
		regVol->SetProductionCuts(cuts);

		targetALV->SetRegion(regVol);
		regVol->AddRootLogicalVolume(targetALV);

		targetBLV->SetRegion(regVol);
		regVol->AddRootLogicalVolume(targetBLV);

		// Visibility
		G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::Red());
		simpleAlSVisAtt->SetVisibility(true);
		targetALV->SetVisAttributes(simpleAlSVisAtt);

		G4VisAttributes* simpleBlSVisAtt= new G4VisAttributes(G4Colour::Yellow());
		simpleBlSVisAtt->SetVisibility(true);
		targetBLV->SetVisAttributes(simpleBlSVisAtt);

		bCreated = true;
		break;
	}
	case 15:
	{
		// Physical and logical volumes
		G4Material *W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
		G4Box *targetABox = new G4Box("targetABox", 5.*mm, 5.*mm, (.025*25.4/2.)*mm);
		G4LogicalVolume *targetALV = new G4LogicalVolume(targetABox, W, "targetALV", 0, 0, 0);

		G4Material *Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
		G4Box *targetBBox = new G4Box("targetABox", 5.*mm, 5.*mm, (0.312*25.4/2.)*mm);
		G4LogicalVolume *targetBLV = new G4LogicalVolume(targetBBox, Cu, "targetBLV", 0, 0, 0);

		// specific translations for the various parts of the component
		new G4PVPlacement(
				0,
				G4ThreeVector(0.,0.,targetABox->GetZHalfLength()),
				"targetAPV",
				targetALV,
				PVWorld,
				false,
				0);

		new G4PVPlacement(
				0,
				G4ThreeVector(0.,0.,targetABox->GetZHalfLength()*2.+targetBBox->GetZHalfLength()),
				"targetBPV",
				targetBLV,
				PVWorld,
				false,
				0);

		// Region for cuts
		G4Region *regVol= new G4Region("targetR");
		G4ProductionCuts* cuts = new G4ProductionCuts;
		cuts->SetProductionCut(0.2*cm);
		regVol->SetProductionCuts(cuts);

		targetALV->SetRegion(regVol);
		regVol->AddRootLogicalVolume(targetALV);

		targetBLV->SetRegion(regVol);
		regVol->AddRootLogicalVolume(targetBLV);

		// Visibility
		G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::Red());
		simpleAlSVisAtt->SetVisibility(true);
		targetALV->SetVisAttributes(simpleAlSVisAtt);

		G4VisAttributes* simpleBlSVisAtt= new G4VisAttributes(G4Colour::Yellow());
		simpleBlSVisAtt->SetVisibility(true);
		targetBLV->SetVisAttributes(simpleBlSVisAtt);

		bCreated=true;
		break;
	}
	}
	return bCreated;
}
bool CML2Acc2::primaryCollimator()
{
	bool bCreated = false;
	// the component as a whole
	G4double totalHeight = 80.0*mm;

	G4Material *Pb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
	G4Material *Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

	// Region for cuts
	G4Region *regVol = new G4Region("primaryCollimator");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(1.*cm);
	regVol->SetProductionCuts(cuts);

	// //-------------------- upper part----------------
	G4ThreeVector centre=G4ThreeVector(0.,0.,16.)-G4ThreeVector(0.,0.,(16. + totalHeight-76.)/2.);
	G4Tubs *PCUTube = new G4Tubs("PrimaryCollimatorUTube", 10.*mm, 40.0*mm, 10.*mm, 0.*deg, 360.*deg);
	G4LogicalVolume *PCUTubeLV = new G4LogicalVolume(PCUTube, Pb, "PrimaryCollimatorUTubeLV", 0, 0, 0);

	new G4PVPlacement(0, centre, "PrimaryCollimatorUTubePV", PCUTubeLV, PVWorld, false, 0);
	G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::White());
	simpleAlSVisAtt->SetVisibility(true);
	PCUTubeLV->SetVisAttributes(simpleAlSVisAtt);
	PCUTubeLV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(PCUTubeLV);

	// //-------------------- lower part----------------
	// Tube
	G4Tubs* PCLTube = new G4Tubs("PrimaryCollimatorLTube", 0., 40.*mm, 60./2.*mm, 0.*deg, 360.*deg);
	// Cone
	G4double coneAperture = 14.*deg;
	G4Cons* collimCone = new G4Cons("PrimaryCollimatorLCone", 0., (16.*std::tan(coneAperture))*mm, 0., (76.*std::tan(coneAperture))*mm, 30.*mm,  0.*deg, 360.*deg);

	G4LogicalVolume* PCLTubeLV = new G4LogicalVolume(PCLTube, Pb, "PCLTubeLV",0,0,0);
	G4LogicalVolume* collimConeLV = new G4LogicalVolume(collimCone, Vacuum, "collimConeLV",0,0,0);
	centre = G4ThreeVector(0.,0.,16.+60./2.);
	G4VPhysicalVolume *PCLTubePV = new G4PVPlacement(0, centre, "PCLTubePV", PCLTubeLV, PVWorld, false, 0);

	centre = G4ThreeVector(0.,0.,0.);
	new G4PVPlacement(0, centre, "TubeMinusConeLPV", collimConeLV, PCLTubePV, false, 0);

	// Visualization
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAtt->SetVisibility(true);
 	simpleAlSVisAtt->SetForceSolid(true);
 	PCLTubeLV->SetVisAttributes(simpleAlSVisAtt);

	// Region for cuts
	PCLTubeLV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(PCLTubeLV);

	bCreated = true;
	return bCreated;
}
bool CML2Acc2::vacuumWindow()
{
	bool bCreated = false;

	G4Material *Be=G4NistManager::Instance()->FindOrBuildMaterial("G4_Be");

	// Region for cuts
	G4Region *regVol = new G4Region("BeWindow");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*cm);
	regVol->SetProductionCuts(cuts);


	G4Tubs* BeWTube = new G4Tubs("BeWindowTube", 0., 50.*mm, 0.12*mm, 0.*deg, 360.*deg);
	G4LogicalVolume *BeWTubeLV = new G4LogicalVolume(BeWTube, Be, "BeWTubeLV", 0, 0, 0);
	new G4PVPlacement(0, G4ThreeVector(0.,0.,90.*mm), "BeWTubePV", BeWTubeLV, PVWorld, false, 0);

	G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::Yellow());
	simpleAlSVisAtt->SetVisibility(true);
	BeWTubeLV->SetVisAttributes(simpleAlSVisAtt);
	BeWTubeLV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(BeWTubeLV);

	bCreated = true;
	return bCreated;
}
bool CML2Acc2::flatteningFilter()
{
	G4String iName;
	char a[10];
	bool bCreated = false;
	switch (idEnergy)
	{
	case 6:
	{
		G4double distanceLast = 116.3386 - 0.932*25.4/2.*mm;

		G4Material *Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
		const int PointsNumber=22;
		// Coordinates of flattening filter in inches
		G4double xRadius[PointsNumber]={0.0000000001, .025, .050, .075, .1, .15, .2, .25, .3, .35, .4, .5, .6, .7, .8, .9, 1., 1.1, 1.205, 1.3, 1.325, 1.5};
		G4double yHeight[PointsNumber]={.932, .92, .907, .892, .874, .83, .782, .736, .69, .645, .601, .516, .437, .361, .294, .22, .173, .118, .08, .08, .125, .125};

		int i;
		// From inches to millimeters
		for (i = 0; i<PointsNumber; i++)
		{
			xRadius[i]*=25.4*mm;
			yHeight[i]*=25.4*mm;
		}
		G4double angleStart=0.;
		G4double angleStop=360.0;

		G4double halfHeigth, rMaxInf, rMaxSup, rMinInf, rMinSup;
		G4String name;
		G4ThreeVector centre;
		G4LogicalVolume *logVol;
		G4Cons *cone;
//		G4VPhysicalVolume *phVol;
		G4VisAttributes* simpleAlSVisAtt;
		G4Region *regVol= new G4Region("flatfilterR");
		G4ProductionCuts* cuts = new G4ProductionCuts;
		cuts->SetProductionCut(0.2*cm);
		regVol->SetProductionCuts(cuts);

		for (i = 0 ; i < 19; i++)
		{
			sprintf(a,"%d", i);
			iName = (G4String)a;
			halfHeigth = (yHeight[i]-yHeight[i+1])/2. ;
			if (i == 18)
			{
				halfHeigth = yHeight[i]/2.;
			}

			rMaxInf = xRadius[i];
			rMaxSup = xRadius[i+1];
			rMinInf = 0.;
			rMinSup = 0.;
			centre.set(0.,0.,distanceLast+halfHeigth); distanceLast+=halfHeigth*2.;
			name = "ffConeG"+iName;
			cone = new G4Cons(name, rMinInf, rMaxInf, rMinSup, rMaxSup, halfHeigth, angleStart, angleStop);
			name = "ffConeLV"+iName;
			logVol = new G4LogicalVolume(cone, Cu, name, 0, 0, 0);
			name = "ffConePV"+iName;
			new G4PVPlacement(0, centre, name, logVol, PVWorld, false, 0);

			// Region for cuts
			logVol->SetRegion(regVol);
			regVol->AddRootLogicalVolume(logVol);
			// Visualization
			simpleAlSVisAtt= new G4VisAttributes(G4Colour::Cyan());
			simpleAlSVisAtt->SetVisibility(true);
			logVol->SetVisAttributes(simpleAlSVisAtt);
		}

		halfHeigth = yHeight[19]/2. ;
		rMaxInf = xRadius[21];
		rMaxSup = xRadius[21];
		rMinInf = 0.;
		rMinSup = 0.;
		centre.set(0.,0., distanceLast+halfHeigth);
		sprintf(a,"%d", i);
		iName = (G4String)a;
		name = "ffConeG"+iName;
		cone = new G4Cons(name, rMinInf, rMaxInf, rMinSup, rMaxSup, halfHeigth, angleStart, angleStop);
		name = "ffConeLV"+iName;
		logVol = new G4LogicalVolume(cone, Cu, name, 0, 0, 0);
		name = "ffConePV"+iName;
		new G4PVPlacement(0, centre, name, logVol, PVWorld, false, 0);

		// Region for cuts
		logVol->SetRegion(regVol);
		regVol->AddRootLogicalVolume(logVol);
		// Visualization
		simpleAlSVisAtt = new G4VisAttributes(G4Colour::Cyan());
		simpleAlSVisAtt -> SetVisibility(true);
		logVol -> SetVisAttributes(simpleAlSVisAtt);

		halfHeigth = (yHeight[20]-yHeight[19])/2. ;
		rMaxInf = xRadius[21];
		rMaxSup = xRadius[21];
		rMinInf = xRadius[20];
		rMinSup = xRadius[19];
		centre.set(0.,0., distanceLast-halfHeigth);
		sprintf(a,"%d", ++i);
		iName = (G4String)a;
		name = "ffConeG"+iName;
		cone = new G4Cons(name, rMinInf, rMaxInf, rMinSup, rMaxSup, halfHeigth, angleStart, angleStop);
		name = "ffConeLV"+iName;
		logVol = new G4LogicalVolume(cone, Cu, name, 0, 0, 0);
		name = "ffConePV"+iName;
		new G4PVPlacement(0, centre, name, logVol,PVWorld, false, 0);

		// Region for cuts
		logVol->SetRegion(regVol);
		regVol->AddRootLogicalVolume(logVol);
		// Visualization
		simpleAlSVisAtt= new G4VisAttributes(G4Colour::Cyan());
		simpleAlSVisAtt->SetVisibility(true);
		logVol->SetVisAttributes(simpleAlSVisAtt);

		bCreated = true;
		break;
	}
		case 15:
		{
			G4double distanceLast = 125.+.125*25.4 - 0.744*25.4/2.- 0.744*25.4/2.*mm;

			G4Material *Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
			const int PointsNumber = 21;
			// Coordinate in pollici del flattening filter come da disegno
			G4double xRadius[PointsNumber] = {0.0000000001, 0.040, 0.078,
					0.105, 0.131, 0.160, 0.205, 0.248, 0.343, 0.436, 0.531,
					0.628, 0.727, 0.829, 0.880, 0.932, 0.983, 1.037, 1.250,
					1.350, 1.5};
			G4double yHeight[PointsNumber] = {0.744, 0.718, 0.682, 0.653,
					0.622, 0.593, 0.547, 0.506, 0.427, 0.354, 0.287, 0.225,
					0.168, 0.119, 0.090, 0.067, 0.047, 0.030, 0.030, 0.125,
					0.125};

			// From inches to millimeters
			int i;
			for (i = 0; i < PointsNumber; i++)
			{
				xRadius[i]*=25.4*mm;
				yHeight[i]*=25.4*mm;
			}
			G4double angleStart = 0.;
			G4double angleStop = 360.0;

			G4double halfHeigth, rMaxInf, rMaxSup, rMinInf, rMinSup;
			G4String name;
			G4ThreeVector centre;
			G4Cons *cone;
			G4LogicalVolume *logVol;
//			G4VPhysicalVolume *phVol;
			G4VisAttributes* simpleAlSVisAtt;
			G4Region *regVol= new G4Region("flatfilterR");
			G4ProductionCuts* cuts = new G4ProductionCuts;
			cuts->SetProductionCut(0.2*cm);
			regVol->SetProductionCuts(cuts);

			for (i=0;i<18;i++)
			{
				sprintf(a,"%d", i);
				iName = (G4String)a;
				halfHeigth = (yHeight[i]-yHeight[i+1])/2. ;
				if (i == 17)
				{
					halfHeigth = yHeight[i]/2.;
				}
				rMaxInf = xRadius[i];
				rMaxSup = xRadius[i+1];
				rMinInf = 0.;
				rMinSup = 0.;
				centre.set(0.,0.,distanceLast+halfHeigth); distanceLast+=halfHeigth*2.;
				name = "ffConeG"+iName;
				cone = new G4Cons(name, rMinInf, rMaxInf, rMinSup, rMaxSup, halfHeigth, angleStart, angleStop);
				name = "ffConeLV"+iName;
				logVol = new G4LogicalVolume(cone, Cu, name, 0, 0, 0);
				name = "ffConePV"+iName;
				new G4PVPlacement(0, centre, name, logVol, PVWorld, false, 0);

				// Region for cuts
				logVol -> SetRegion(regVol);
				regVol -> AddRootLogicalVolume(logVol);
				// Visualization
				simpleAlSVisAtt = new G4VisAttributes(G4Colour::Cyan());
				simpleAlSVisAtt -> SetVisibility(true);
				logVol -> SetVisAttributes(simpleAlSVisAtt);
			}

			halfHeigth = yHeight[18]/2. ;
			rMaxInf = xRadius[20];
			rMaxSup = xRadius[20];
			rMinInf = 0.;
			rMinSup = 0.;
			centre.set(0.,0., distanceLast+halfHeigth);
			sprintf(a,"%d", i);
			iName = (G4String)a;
			name = "ffConeG"+iName;
			cone = new G4Cons(name, rMinInf, rMaxInf, rMinSup, rMaxSup, halfHeigth, angleStart, angleStop);
			name = "ffConeLV"+iName;
			logVol = new G4LogicalVolume(cone, Cu, name, 0, 0, 0);
			name = "ffConePV"+iName;
			new G4PVPlacement(0, centre, name, logVol, PVWorld, false, 0);

			// Region for cuts
			logVol -> SetRegion(regVol);
			regVol -> AddRootLogicalVolume(logVol);
			// Visualization
			simpleAlSVisAtt = new G4VisAttributes(G4Colour::Cyan());
			simpleAlSVisAtt -> SetVisibility(true);
			logVol -> SetVisAttributes(simpleAlSVisAtt);

			halfHeigth = (yHeight[19]-yHeight[18])/2. ;
			rMaxInf = xRadius[20];
			rMaxSup = xRadius[20];
			rMinInf = xRadius[19];
			rMinSup = xRadius[18];
			centre.set(0.,0., distanceLast-halfHeigth);
			sprintf(a,"%d", ++i);
			iName = (G4String)a;
			name = "ffConeG"+iName;
			cone = new G4Cons(name, rMinInf, rMaxInf, rMinSup, rMaxSup, halfHeigth, angleStart, angleStop);
			name = "ffConeLV"+iName;
			logVol = new G4LogicalVolume(cone, Cu, name, 0, 0, 0);
			name = "ffConePV"+iName;
			new G4PVPlacement(0, centre, name, logVol,PVWorld, false, 0);

			// Region for cuts
			logVol -> SetRegion(regVol);
			regVol -> AddRootLogicalVolume(logVol);
			// Visualization
			simpleAlSVisAtt = new G4VisAttributes(G4Colour::Cyan());
			simpleAlSVisAtt -> SetVisibility(true);
			logVol -> SetVisAttributes(simpleAlSVisAtt);

			bCreated = true;
			break;
		}
	}
	return bCreated;
}
bool CML2Acc2::ionizationChamber()
{
	bool bCreated = false;


	G4Material *KAPTON = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
	G4Region *regVol = new G4Region("ionizationChamber");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts -> SetProductionCut(0.1*cm);
	regVol -> SetProductionCuts(cuts);

	G4VisAttributes* simpleAlSVisAtt;
	// Region for cuts

	G4Tubs* ICTubeW = new G4Tubs("ionizationChamberTube", 0., 3.75*2.54*10.*mm, 0.005*25.4*mm, 0.*deg, 360.*deg);
	G4Tubs* ICTubeP = new G4Tubs("ionizationChamberTube", 0., 3.75*2.54*10.*mm, 0.002*25.4*mm, 0.*deg, 360.*deg);

	G4ThreeVector centre;
	// W1
	centre.set(0.,0.,148.35*mm);
	G4LogicalVolume *PCUTubeW1LV = new G4LogicalVolume(ICTubeW, KAPTON, "ionizationChamberTubeW1LV", 0, 0, 0);
	new G4PVPlacement(0, centre, "ionizationChamberTubeW1PV", PCUTubeW1LV, PVWorld, false, 0);
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAtt -> SetVisibility(true);
	PCUTubeW1LV -> SetVisAttributes(simpleAlSVisAtt);
	PCUTubeW1LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(PCUTubeW1LV);

	// P1
	centre.set(0.,0.,150.73*mm);
	G4LogicalVolume *PCUTubeP1LV = new G4LogicalVolume(ICTubeP, KAPTON, "ionizationChamberTubeP1LV", 0, 0, 0);
	new G4PVPlacement(0, centre, "ionizationChamberTubeP1PV", PCUTubeP1LV, PVWorld, false, 0);
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Yellow());
	simpleAlSVisAtt -> SetVisibility(true);
	PCUTubeP1LV -> SetVisAttributes(simpleAlSVisAtt);
	PCUTubeP1LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(PCUTubeP1LV);

	// W2
	centre.set(0.,0.,155.5*mm);
	G4LogicalVolume *PCUTubeW2LV = new G4LogicalVolume(ICTubeW, KAPTON, "ionizationChamberTubeW2LV", 0, 0, 0);
	new G4PVPlacement(0, centre, "ionizationChamberTubeW2PV", PCUTubeW2LV, PVWorld, false, 0);
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAtt -> SetVisibility(true);
	PCUTubeW2LV -> SetVisAttributes(simpleAlSVisAtt);
	PCUTubeW2LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(PCUTubeW2LV);

	// P2
	centre.set(0.,0.,153.12*mm);
	G4LogicalVolume *PCUTubeP2LV = new G4LogicalVolume(ICTubeP, KAPTON, "ionizationChamberTubeP2LV", 0, 0, 0);
	new G4PVPlacement(0, centre, "ionizationChamberTubeP2PV", PCUTubeP2LV, PVWorld, false, 0);
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Yellow());
	simpleAlSVisAtt->SetVisibility(true);
	PCUTubeP2LV->SetVisAttributes(simpleAlSVisAtt);
	PCUTubeP2LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(PCUTubeP2LV);

	// W3
	centre.set(0.,0.,162.65*mm);
	G4LogicalVolume *PCUTubeW3LV = new G4LogicalVolume(ICTubeW, KAPTON, "ionizationChamberTubeW3LV", 0, 0, 0);
	new G4PVPlacement(0, centre, "ionizationChamberTubeW3PV", PCUTubeW3LV, PVWorld, false, 0);
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAtt -> SetVisibility(true);
	PCUTubeW3LV -> SetVisAttributes(simpleAlSVisAtt);
	PCUTubeW3LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(PCUTubeW3LV);

	// P3
	centre.set(0.,0.,157.88*mm);
	G4LogicalVolume *PCUTubeP3LV = new G4LogicalVolume(ICTubeP, KAPTON, "ionizationChamberTubeP3LV", 0, 0, 0);
	new G4PVPlacement(0, centre, "ionizationChamberTubeP3PV", PCUTubeP3LV, PVWorld, false, 0);
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Yellow());
	simpleAlSVisAtt -> SetVisibility(true);
	PCUTubeP3LV -> SetVisAttributes(simpleAlSVisAtt);
	PCUTubeP3LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(PCUTubeP3LV);

	// P4
	centre.set(0.,0.,160.27*mm);
	G4LogicalVolume *PCUTubeP4LV = new G4LogicalVolume(ICTubeP, KAPTON, "ionizationChamberTubeP4LV", 0, 0, 0);
	new G4PVPlacement(0, centre, "ionizationChamberTubeP4PV", PCUTubeP4LV, PVWorld, false, 0);
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Yellow());
	simpleAlSVisAtt->SetVisibility(true);
	PCUTubeP4LV->SetVisAttributes(simpleAlSVisAtt);
	PCUTubeP4LV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(PCUTubeP4LV);

	bCreated = true;
	return bCreated;
}
bool CML2Acc2::mirror()
{
	bool bCreated = false;

	G4Material *MYLAR = G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR");

	// Region for cuts
	G4Region *regVol = new G4Region("Mirror");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts -> SetProductionCut(0.1*cm);
	regVol -> SetProductionCuts(cuts);

	G4Tubs* MirrorTube = new G4Tubs("MirrorTube", 0., 86.*mm, 1.*mm, 0.*deg, 360.*deg);
	G4LogicalVolume *MirrorTubeLV = new G4LogicalVolume(MirrorTube, MYLAR, "MirrorTubeLV", 0, 0, 0);
	G4RotationMatrix *cRotation = new G4RotationMatrix();
	cRotation -> rotateY(35.0*deg);
	new G4PVPlacement(cRotation, G4ThreeVector(0., 0., 220.*mm), "MirrorTubePV", MirrorTubeLV,PVWorld, false, 0);

	G4VisAttributes* simpleAlSVisAtt = new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAtt -> SetVisibility(true);
	MirrorTubeLV -> SetVisAttributes(simpleAlSVisAtt);
	MirrorTubeLV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(MirrorTubeLV);

	bCreated = true;
	return bCreated;
}
bool CML2Acc2::Jaw1X()
{
	bool bCreated = false;

	G4Material *W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
	G4String name = "Jaws1X";

	G4Region *regVol= new G4Region(name+"R");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(1.*cm);
	regVol->SetProductionCuts(cuts);

	G4VisAttributes* simpleAlSVisAtt;
	G4ThreeVector centre, halfSize;
	G4RotationMatrix *cRotation = new G4RotationMatrix();
	centre.set(0.,0.,(367.+78./2.)*mm);
	halfSize.set(50.*mm, 90.*mm, 78./2.*mm);
	G4Box *box = new G4Box(name+"Box", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	G4LogicalVolume *logVol = new G4LogicalVolume(box, W, name+"LV", 0, 0, 0);
	SetJawAperture(1, centre, halfSize,cRotation);
	new G4PVPlacement(cRotation, centre, name+"PV", logVol, PVWorld, false, 0);

	// Region for cuts
	logVol->SetRegion(regVol);
	regVol->AddRootLogicalVolume(logVol);

	// Visibility
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAtt->SetVisibility(true);
	logVol->SetVisAttributes(simpleAlSVisAtt);

	bCreated = true;
	return bCreated;
}
bool CML2Acc2::Jaw2X()
{
	bool bCreated=false;
	G4Material *W=G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
	G4String name="Jaws2X";
	G4Region *regVol= new G4Region(name+"R");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(1.*cm);
	regVol->SetProductionCuts(cuts);

	G4VisAttributes* simpleAlSVisAtt;
	G4ThreeVector centre, halfSize;
	G4RotationMatrix *cRotation = new G4RotationMatrix();
	centre.set(0.,0.,(367.+78./2.)*mm);
	halfSize.set(50.*mm, 90.*mm, 78./2.*mm);
	G4Box *box = new G4Box(name+"Box", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	G4LogicalVolume *logVol = new G4LogicalVolume(box, W, name+"LV", 0, 0, 0);
	SetJawAperture(2, centre, halfSize, cRotation);
	new G4PVPlacement(cRotation, centre, name+"PV", logVol, PVWorld, false, 0);

	// Region for cuts
	logVol->SetRegion(regVol);
	regVol->AddRootLogicalVolume(logVol);

	// Visibility
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAtt->SetVisibility(true);
	logVol->SetVisAttributes(simpleAlSVisAtt);

	bCreated = true;
	return bCreated;
}
bool CML2Acc2::Jaw1Y()
{
	bool bCreated = false;
	G4Material *W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
	G4String name="Jaws1Y";

	G4Region *regVol= new G4Region(name+"R");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(1.*cm);
	regVol->SetProductionCuts(cuts);

	G4VisAttributes* simpleAlSVisAtt;
	G4ThreeVector centre, halfSize;
	G4RotationMatrix *cRotation=new G4RotationMatrix();
	centre.set(0.,0.,(280.+78./2.)*mm);
	halfSize.set(90.*mm, 50.*mm, 78./2.*mm);
	G4Box *box = new G4Box(name+"Box", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	G4LogicalVolume *logVol = new G4LogicalVolume(box, W, name+"LV", 0, 0, 0);
	SetJawAperture(3, centre, halfSize, cRotation);
	new G4PVPlacement(cRotation, centre, name+"PV", logVol, PVWorld, false, 0);

	// Region for cuts
	logVol->SetRegion(regVol);
	regVol->AddRootLogicalVolume(logVol);

	// Visibility
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Red());
	simpleAlSVisAtt->SetVisibility(true);
	logVol->SetVisAttributes(simpleAlSVisAtt);

	bCreated = true;
	return bCreated;
}
bool CML2Acc2::Jaw2Y()
{
	bool bCreated = false;

	G4Material *W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
	G4String name = "Jaws2Y";
	G4Region *regVol = new G4Region(name+"R");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(1.*cm);
	regVol->SetProductionCuts(cuts);

	G4ThreeVector centre, halfSize;
	G4RotationMatrix *cRotation=new G4RotationMatrix();
	centre.set(0.,0.,(280.+78./2.)*mm);
	halfSize.set(90.*mm, 50.*mm, 78./2.*mm);
	G4Box *box = new G4Box(name+"Box", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	G4LogicalVolume *logVol = new G4LogicalVolume(box, W, name+"LV", 0, 0, 0);
	SetJawAperture(4, centre, halfSize, cRotation);
	new G4PVPlacement(cRotation, centre, name+"PV", logVol, PVWorld, false, 0);

	// Region for cuts
	logVol->SetRegion(regVol);
	regVol->AddRootLogicalVolume(logVol);

	// Visibility

	G4VisAttributes* simpleAlSVisAtt = new G4VisAttributes(G4Colour::Red());
	simpleAlSVisAtt -> SetVisibility(true);
	logVol -> SetVisAttributes(simpleAlSVisAtt);

	bCreated = true;
	return bCreated;
}
bool CML2Acc2::MLC()
{
	G4String iName;
	char a[12]; 
	bool bCreated = false;
	//    material
	G4Material *W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
	G4VisAttributes* simpleAlSVisAtt;
	// Region for cuts
	G4Region *regVol = new G4Region("MLCR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.5*cm);
	regVol->SetProductionCuts(cuts);

	G4ThreeVector box1Size, box2Size;
	
	G4ThreeVector centreStart;
	centreStart.set(0.,0.,(482.5+535.5)/2.*mm);

	box1Size.set(5.2883/2.*mm, 167.488/2.*mm, 61.2953/2.*mm);
	G4double circleRadiuos=80.0*mm;
	
	G4double alfa=std::asin(box1Size.getZ()/circleRadiuos);
	G4double externalYPartSize=circleRadiuos*(1-std::cos(alfa));
	
	box2Size.set(box1Size.getX(), externalYPartSize/2., box1Size.getZ());

	// single leaf
	G4Box *box1Leaf = new G4Box("LeafBox1", box1Size.getX(), box1Size.getY(), box1Size.getZ());
	G4Box *box2leaf = new G4Box("LeafBox2", box2Size.getX(), box2Size.getY(), box2Size.getZ());
    G4Tubs *cyLeaf  = new G4Tubs("LeafCylinder", 0, circleRadiuos, box1Size.getX(),0,CLHEP::twopi);
	
	G4RotationMatrix *rm=new G4RotationMatrix();
	rm->rotateY(90.*deg);
	G4DisplacedSolid *cyLeaf90Y=new G4DisplacedSolid("LeafCylinder90Y",cyLeaf,rm,G4ThreeVector(0,0,0));
	G4DisplacedSolid *cyLeafTr=new G4DisplacedSolid("LeafCylinderTr",cyLeaf90Y,0,G4ThreeVector(0., -(circleRadiuos-box2Size.getY()), 0.));
	G4IntersectionSolid* bx2CyleafTr =new G4IntersectionSolid("bx2CyleafTr", box2leaf, cyLeafTr);
	G4DisplacedSolid *bx2CyleafTrTr=new G4DisplacedSolid("bx2CyleafTrTr",bx2CyleafTr,0,G4ThreeVector(0., +(box1Size.getY()+box2Size.getY()-.1), 0.)); // -.1 to guarantee the volumes touch each other
	
	G4UnionSolid *leafSolidA = new G4UnionSolid("SingleLeafA", bx2CyleafTrTr, box1Leaf);
	rm = new G4RotationMatrix();
	rm -> rotateZ(180.*deg);
	G4DisplacedSolid *leafSolidB = new G4DisplacedSolid("SingleLeafB", leafSolidA, rm, G4ThreeVector(0,0,0));
	
	
	G4ThreeVector halfSize;
	halfSize.set(0., box2Size.getY()*2.-.1,0.);
	halfSize+=box1Size;

	G4LogicalVolume *leafLVA = new G4LogicalVolume(leafSolidA, W, "leafSolidALV", 0, 0, 0);
	G4LogicalVolume *leafLVB = new G4LogicalVolume(leafSolidB, W, "leafSolidBLV", 0, 0, 0);

	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAtt->SetVisibility(true);
	leafLVA->SetVisAttributes(simpleAlSVisAtt);
	leafLVA->SetRegion(regVol);
	regVol->AddRootLogicalVolume(leafLVA);
	
	simpleAlSVisAtt= new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAtt->SetVisibility(true);
	leafLVB->SetVisAttributes(simpleAlSVisAtt);
	leafLVB->SetRegion(regVol);
	regVol->AddRootLogicalVolume(leafLVB);


	G4String PVname;
//	G4VPhysicalVolume *leafPhys;
	int i;
	int j = 0;

	G4ThreeVector centre;
	centre= centreStart + G4ThreeVector(-38.*halfSize.getX(), 0.,0.);

	for (i = 1; i < (int)vec_leavesA.size(); i++)
	{
		sprintf(a,"%d", i);
		iName = (G4String)a;
		PVname = "leafA"+iName;
		centre.setX(centre.getX()+halfSize.getX()*2.);
		centre.setY(-halfSize.getY()-vec_leavesA[i]);
		new G4PVPlacement(0, centre, PVname, leafLVA, PVWorld, false, i);
		j++;
	}
	centre=centreStart+G4ThreeVector(-38.*halfSize.getX(), 0.,0.);
	for (i = 1; i < (int)vec_leavesB.size(); i++)
	{
		sprintf(a,"%d", i);
		iName = (G4String)a;
		PVname = "leafB"+iName;
		centre.setX(centre.getX()+halfSize.getX()*2.);
		centre.setY(+halfSize.getY()+vec_leavesB[i]);
		new G4PVPlacement(0, centre, PVname, leafLVB, PVWorld, false, i);
		j++;
	}
	bCreated = true;
	return bCreated;
}

