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


#include "ML2WorldConstruction.hh"
#include "G4SystemOfUnits.hh"

CML2WorldConstruction::CML2WorldConstruction():acceleratorEnv(0),phantomEnv(0),PVWorld(0),phaseSpace(0),backScatteredPlane(0)
{
    phantomEnv = CML2PhantomConstruction::GetInstance();
    acceleratorEnv = CML2AcceleratorConstruction::GetInstance();
    bWorldCreated = false;
    bOnlyVisio = 0;
}

CML2WorldConstruction::~CML2WorldConstruction(void)
{
    delete phaseSpace;
    delete backScatteredPlane;
}

CML2WorldConstruction* CML2WorldConstruction::instance = 0;

CML2WorldConstruction* CML2WorldConstruction::GetInstance()
{
    if (instance == 0)
    {
        instance = new CML2WorldConstruction();
    }
    return instance;
}

G4VPhysicalVolume* CML2WorldConstruction::Construct()
{
    return PVWorld;
}

bool CML2WorldConstruction::create(SInputData *inputData, bool bOV)
{
    // create the world box
    bOnlyVisio = bOV;
    G4double halfSize = 3000.*mm;
    G4Material *Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4Box *worldB = new G4Box("worldG", halfSize, halfSize, halfSize);
    G4LogicalVolume *worldLV = new G4LogicalVolume(worldB, Vacuum, "worldL", 0, 0, 0);
    G4VisAttributes* simpleWorldVisAtt = new G4VisAttributes(G4Colour::Black());
    simpleWorldVisAtt -> SetVisibility(true);
    worldLV -> SetVisAttributes(simpleWorldVisAtt);
    PVWorld = new G4PVPlacement(0,  G4ThreeVector(0.,0.,0.), "worldPV", worldLV, 0, false, 0);

    // create the accelerator-world box
    if (!acceleratorEnv -> Construct(PVWorld, bOV))
    {
        G4cout << "\n\n The macro file '" << inputData->generalData.StartFileInputData <<
        		"' refers to a not defined accelerator.\n" << acceleratorEnv->getAcceleratorName() <<
				"\n\nSTOP\n\n" << G4endl;
        return false;
    }

    // create the phantom-world box
    if ( !phantomEnv->Construct(PVWorld,
                        inputData->voxelSegmentation.nX,
                        inputData->voxelSegmentation.nY,
                        inputData->voxelSegmentation.nZ,
			bOV) )
    {
        G4cout << "\n\n The macro file '" << inputData->generalData.StartFileInputData <<
        		"' refers to a not defined phantom.\n" << phantomEnv->getPhantomName() <<
				"\n\nSTOP\n\n" << G4endl;
        return false;
    }

    // if the bSavePhaseSpace flag is true create a phase plane
/*   
 if (inputData -> generalData.bSavePhaseSpace)
    {
        phaseSpace = new CML2PhaseSpaces();
        if (inputData -> generalData.bForcePhaseSpaceBeforeJaws)
        {
        	inputData -> generalData.centrePhaseSpace.setZ(acceleratorEnv->getZ_Value_PhaseSpaceBeforeJaws());
        }

        phaseSpace -> createPlane(idSD_PhaseSpace,
        		inputData->generalData.max_N_particles_in_PhSp_File,
        		inputData->generalData.seed,
				inputData->generalData.nMaxParticlesInRamPlanePhaseSpace,
				acceleratorEnv->getPhysicalVolume(), "PhSp",
				inputData->generalData.PhaseSpaceOutFile,
				inputData->generalData.bSavePhaseSpace,
				inputData->generalData.bStopAtPhaseSpace,
				inputData->generalData.centrePhaseSpace,
				inputData->generalData.halfSizePhaseSpace,
				&inputData->primaryParticleData,
				acceleratorEnv->getAcceleratorIsoCentre()); // phase space plane, yellow
    }

    // create a killer plane to destroy the particles back scattered from the target
    backScatteredPlane = new CML2PhaseSpaces();
    backScatteredPlane -> createPlane(acceleratorEnv->getPhysicalVolume(),
    		"killerPlane", G4ThreeVector(0, 0, -50*mm), G4ThreeVector(200*mm, 200*mm, 1*mm)); // killer plane, cyan
*/
    bWorldCreated = true;
    return bWorldCreated;
}
void CML2WorldConstruction::checkVolumeOverlap()
{
    // loop inside all the daughters volumes
	G4cout<< G4endl;
    //        bool bCheckOverlap;
    //        bCheckOverlap=false;

    int nSubWorlds, nSubWorlds2;
    for (int i=0; i<(int) PVWorld->GetLogicalVolume()->GetNoDaughters(); i++)
    {
        PVWorld->GetLogicalVolume()->GetDaughter(i)->CheckOverlaps();
        nSubWorlds=(int) PVWorld->GetLogicalVolume()->GetDaughter(i)->GetLogicalVolume()->GetNoDaughters();
        for (int j=0; j<nSubWorlds; j++)
        {
            PVWorld->GetLogicalVolume()->GetDaughter(i)->GetLogicalVolume()->GetDaughter(j)->CheckOverlaps();
            nSubWorlds2=(int) PVWorld->GetLogicalVolume()->GetDaughter(i)->GetLogicalVolume()->GetDaughter(j)->GetLogicalVolume()->GetNoDaughters();
            for (int k=0; k<nSubWorlds2; k++)
            {
                PVWorld->GetLogicalVolume()->GetDaughter(i)->GetLogicalVolume()->GetDaughter(j)->GetLogicalVolume()->GetDaughter(k)->CheckOverlaps();
            }
        }
    }
    G4cout<< G4endl;
}
bool CML2WorldConstruction::newGeometry()
{

    G4bool bNewRotation = false;
    G4bool bNewCentre = false;
    G4bool bNewGeometry = false;
    bNewCentre = phantomEnv -> applyNewCentre();
    G4RotationMatrix *rmInv = acceleratorEnv -> rotateAccelerator();
    if (rmInv!=0)
    {
        CML2PrimaryGenerationAction::GetInstance()->setRotation(rmInv);
        bNewRotation = true;
    }
    if (bNewRotation || bNewCentre)
    {
    	bNewGeometry = true;
    }

    return bNewGeometry;
}

