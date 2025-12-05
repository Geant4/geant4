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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4RegionStore.hh"
#include "G4VisAttributes.hh"
#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
    //instantiate the messenger
    fMessenger = new DetectorConstructionMessenger(this);

    //Crystal size
    fCrystalSize.setX(20*CLHEP::mm);
    fCrystalSize.setY(20*CLHEP::mm);
    fCrystalSize.setZ(0.0305*CLHEP::mm);

    //Crystal planes or axes considered
    fLattice = "(111)";

    //Detector size
    fDetectorSize.setX(10*CLHEP::cm);
    fDetectorSize.setY(10*CLHEP::cm);
    fDetectorSize.setZ(0.1*CLHEP::mm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    //Check overlap option
    G4bool checkOverlaps = true;

    //Materials
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
    G4Material* silicon = nist->FindOrBuildMaterial("G4_Si");

    //to use a diamond crystal
    G4Element* elC = nist->FindOrBuildElement("C");
    G4Material* diamond = new G4Material("G4_Diamond", 3.520*CLHEP::g/CLHEP::cm3, 1);
    diamond->AddElement(elC, 1);

    //World
    G4Box* solidWorld = new G4Box("World", 0.2*CLHEP::m, 0.2*CLHEP::m, 30.*CLHEP::m);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");  
    G4VPhysicalVolume* physWorld = new G4PVPlacement
                                    (0,                     // no rotation
                                     G4ThreeVector(),       // centre position
                                     logicWorld,            // its logical volume
                                     "World",               // its name
                                     0,                     // its mother volume
                                     false,                 // no boolean operation
                                     0,                     // copy number
                                     checkOverlaps);        // overlaps checking
    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());


    // --------------- Crystal ------------------------------------

    //Select crystal material
    fCrystalMaterial = nist->FindOrBuildMaterial(fCrystalMaterialStr);

    //Crystal rotation angle (also the angle of crystal planes vs the beam)
    G4RotationMatrix* crystalRotationMatrix = new G4RotationMatrix;
    crystalRotationMatrix->rotateY(-fAngleX);
    crystalRotationMatrix->rotateX(-fAngleY);

    //Setting crystal position
    G4ThreeVector posCrystal = G4ThreeVector(0., 0., fCrystalSize.z()/2.);

    //crystal volume
    G4Box* solidCrystal = new G4Box("Crystal",
                                    fCrystalSize.x()/2,
                                    fCrystalSize.y()/2,
                                    fCrystalSize.z()/2.);
    
    fLogicCrystal = new G4LogicalVolume(solidCrystal,
                                        fCrystalMaterial,
                                        "Crystal");
    new G4PVPlacement(crystalRotationMatrix,
                      posCrystal,
                      fLogicCrystal,
                      "Crystal",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    
    if (fActivateChannelingModel)
    {
      //crystal region (necessary for the FastSim model)
      G4Region* regionCh = new G4Region("Crystal");
      regionCh->AddRootLogicalVolume(fLogicCrystal);
    }

    //visualization attributes
    G4VisAttributes* crystalVisAttribute =
        new G4VisAttributes(G4Colour(1., 0., 0.));
    crystalVisAttribute->SetForceSolid(true);
    fLogicCrystal->SetVisAttributes(crystalVisAttribute);

    //print Crystal info
    G4cout << "Crystal material: " << fCrystalMaterial->GetName() << G4endl;
    G4cout << "Crystal size: " << fCrystalSize.x()/CLHEP::mm
           << "x" << fCrystalSize.y()/CLHEP::mm
           << "x" << fCrystalSize.z()/CLHEP::mm << " mm3" << G4endl;

    if (fActivateChannelingModel)
    {
        G4cout << "G4ChannelingFastSimModel activated" << G4endl;
        G4cout << "Crystal bending angle: " << fBendingAngle << " rad" << G4endl;
        G4cout << "Crystal Lattice: " << fLattice << G4endl;
        G4cout << "Crystal angleX: " << fAngleX << " rad" << G4endl;
        G4cout << "Crystal angleY: " << fAngleY << " rad" << G4endl;
        G4cout << "ActivateRadiationModel: " << fActivateRadiationModel << G4endl;

        if(fVirtualCollimatorHalfSize<CLHEP::halfpi)
        {
            G4cout << "Setting virtual collimator angular radius: "
                   << fVirtualCollimatorHalfSize/CLHEP::mrad << " mrad" << G4endl;
        }
        else
        {
            G4cout << "No virtual collimator set." << G4endl;
        }

        if(fCrystalInternalGeometryPath != "")
        {
            G4cout << "Reading crystal internal geometry activated: " << G4endl;
            G4cout << "reading from the file: " << fCrystalInternalGeometryPath << G4endl;
        }
        else
        {
            if (fCrystallineUndulatorAmplitude > DBL_EPSILON &&
                fCrystallineUndulatorPeriod > DBL_EPSILON) {
                G4cout << "Crystalline undulator activated: " << G4endl;
                G4cout << "undulator amplitude: "
                       << fCrystallineUndulatorAmplitude/CLHEP::nm
                       << " nm" << G4endl;
                G4cout << "undulator period: "
                       << fCrystallineUndulatorPeriod/CLHEP::mm
                       << " mm" << G4endl;
                G4cout << "undulator phase: "
                       << fCrystallineUndulatorPhase
                       << " rad" << G4endl;
            }
        }

        G4cout << G4endl;
    }
    else
    {
        G4cout << "G4ChannelingFastSimModel is not activated" << G4endl << G4endl;
    }

    // --------------- Detector -----------------------------------
    //Setting detector position
    G4ThreeVector posDetector =
            G4ThreeVector(0, 0, fDetectorFrontPosZ+fDetectorSize.z()/2.);

    //particle detector volume
    G4Box* detector = new G4Box("Detector",
                                fDetectorSize.x()/2,
                                fDetectorSize.y()/2,
                                fDetectorSize.z()/2.);
    
    G4LogicalVolume* logicDetector = new G4LogicalVolume(detector,
                                                         silicon,
                                                         "Detector");
    new G4PVPlacement(0, 
                      posDetector, 
                      logicDetector,
                      "Detector", 
                      logicWorld, 
                      false, 
                      0, 
                      checkOverlaps);

    //visualization attributes
    G4VisAttributes* detectorVisAttribute =
        new G4VisAttributes(G4Colour(1., 0., 1., 0.3));
    detectorVisAttribute->SetForceSolid(true);
    logicDetector->SetVisAttributes(detectorVisAttribute);
    
    //always return the physical World
    return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  if (fActivateChannelingModel)
  {
    // --------------- fast simulation ----------------------------
    //extract the region of the crystal from the store
    G4RegionStore* regionStore = G4RegionStore::GetInstance();
    G4Region* regionCh = regionStore->GetRegion("Crystal");

    //create the channeling model for this region
    G4ChannelingFastSimModel* channelingModel =
        new G4ChannelingFastSimModel("ChannelingModel", regionCh);

    //activate the channeling model
    channelingModel->Input(fCrystalMaterial, fLattice, fPotentialPath);
    //setting bending angle of the crystal planes (default is 0)
    channelingModel->GetCrystalData()->SetBendingAngle(fBendingAngle, fLogicCrystal);

    //reading internal crystal geometry from file has a priority vs its setup from parameters;
    //it is used to setup a realistic undulator, but may be used also for a bent crystal geometry
    if(fCrystalInternalGeometryPath != "")
    {
        channelingModel->GetCrystalData()->
            SetCrystallineUndulatorParameters(fLogicCrystal,fCrystalInternalGeometryPath);
    }
    else
    {
        //setting crystalline undulator parameters
        //NOTE: they are incompatible with a bent crystal
        if (fCrystallineUndulatorAmplitude > DBL_EPSILON &&
            fCrystallineUndulatorPeriod > DBL_EPSILON)
        {
            channelingModel->GetCrystalData()->SetCrystallineUndulatorParameters(
                fCrystallineUndulatorAmplitude,
                fCrystallineUndulatorPeriod,
                fCrystallineUndulatorPhase,
                fLogicCrystal);
        }
    }

    /*
    Set the multiple of critical channeling angles (Lindhard angles) which defines
    the angular cut of the model (otherwise standard Geant4 is active):
    a number too low reduces the accuracy of the model;
    a number too high sometimes drastically reduces the simulation speed.
    The default value is 100, while for many problems 10-20 is ok.
    CAUTION: If you set this value to 1 meaning the angular cut = 1*Lindhard angle,
    this will cut off the physics of overbarrier motion, still important even if
    the particle is not in channeling. This will provide incorrect results.
    */
    channelingModel->SetDefaultLindhardAngleNumberHighLimit(fLindhardAngles);
    //you may set a particular limit for a certain particle type
    //(has a priority vs default):
    channelingModel->SetLindhardAngleNumberHighLimit(fLindhardAnglesProton,"proton");
    channelingModel->SetLindhardAngleNumberHighLimit(fLindhardAnglesAntiproton,
                                                     "anti_proton");
    channelingModel->SetLindhardAngleNumberHighLimit(fLindhardAnglesPiPlus,"pi+");
    channelingModel->SetLindhardAngleNumberHighLimit(fLindhardAnglesPiMinus,"pi-");
    channelingModel->SetLindhardAngleNumberHighLimit(fLindhardAnglesPositron,"e+");
    channelingModel->SetLindhardAngleNumberHighLimit(fLindhardAnglesElectron,"e-");
    channelingModel->SetLindhardAngleNumberHighLimit(fLindhardAnglesMuPlus,"mu+");
    channelingModel->SetLindhardAngleNumberHighLimit(fLindhardAnglesMuMinus,"mu-");
    //channelingModel->SetLindhardAngleNumberHighLimit(your_value,"your_particle");

    /*
    Set the low kinetic energy cut for the model:
    too low energy may reduce the simulation speed (sometimes drastically),
    too high energy can cut off a useful physics for radiation losses.
    A recommended value depends a lot on the case. Usually it should be above 100 MeV,
    since below quantum channeling effects may become important, though lower energies
    are not forbidden. For energies considerably above 1 GeV and a crystal thick enough
    for multiphoton radiation emission,a lower energy limit of 300-500 MeV is recommended.
    */
    channelingModel->SetDefaultLowKineticEnergyLimit(fParticleMinKinEnergy);
    //you may set a particular limit for a certain particle type
    //(has a priority vs default):
    channelingModel->SetLowKineticEnergyLimit(fProtonMinKinEnergy,"proton");
    channelingModel->SetLowKineticEnergyLimit(fAntiprotonMinKinEnergy,"anti_proton");
    channelingModel->SetLowKineticEnergyLimit(fPiPlusMinKinEnergy,"pi+");
    channelingModel->SetLowKineticEnergyLimit(fPiMinusMinKinEnergy,"pi-");
    channelingModel->SetLowKineticEnergyLimit(fPositronMinKinEnergy,"e+");
    channelingModel->SetLowKineticEnergyLimit(fElectronMinKinEnergy,"e-");
    channelingModel->SetLowKineticEnergyLimit(fMuPlusMinKinEnergy,"mu+");
    channelingModel->SetLowKineticEnergyLimit(fMuMinusMinKinEnergy,"mu-");
    //channelingModel->SetLowKineticEnergyLimit(your_value,"your_particle");

    /*
    activate the radiation model (do it only when you want to take into account the
    radiation production in an oriented crystal; it reduces simulation speed.)
    */
    if (fActivateRadiationModel)
    {
        channelingModel->RadiationModelActivate();
        G4cout << "Radiation model activated" << G4endl;

        /*
        Set the number of the photons used in Monte Carlo integral in Baier-Katkov:
        too low number reduces the accuracy of the model;
        too high number reduces the calculation speed.
        In most of the cases 150 is a minimal secure number.
        */
        channelingModel->GetRadiationModel()->
                SetSamplingPhotonsNumber(fSamplingPhotonsNumber);

        /*
        Increase the statistics of sampling photons in a certain energy range.
        By default it is not active! In many cases it is not necessary.
        It is very useful for a soft spectrum part, when it is considerably below
        the charged particle energy, in order to increase the accuracy. It is possible
        to apply as many ranges as you want.
        NOTE: usually important for crystalline undulator
        CAUTION: insert only an integer number as a multiple of a photon statistics
        and ONLY > 1 .
        CAUTION: this energy range must not be beyond the range of radiation energies
        (i.e. below minimum photon energy (see below)) and must not intersect another
        range with an increased statistics if any.
        CAUTION: this is a multiple of the statistics of sampling photons randomly get in
        this energy range => make sure the total number of sampling photons is high enough
        to regularly get in this energy range, otherwise the statistics will not increase.
        */
        if(fTimesPhotonStatistics>1)
        {channelingModel->GetRadiationModel()->
                    AddStatisticsInPhotonEnergyRegion(fMinPhotonEnergyAddStat,
                                                      fMaxPhotonEnergyAddStat,
                                                      fTimesPhotonStatistics);}

        /*
        Adjust the angular distribution of the sampling photons by
        changing the multiple of the opening radiation angle 1/gamma:
        the model should work correctly in the range 2-5;
        too small multiple reduces the accuracy;
        too high value requires more sampling photons.
        */
        channelingModel->GetRadiationModel()->
                SetRadiationAngleFactor(fRadiationAngleFactor);

        /*
        Set the minimal energy of radiated photon: generally it depends on which part
        of the spectrum you are interested in:
        too low number vs the charged particle energy may require more sampling photons
        in a total or in a particular energy range;
        too high number may reduce the accuracy in radiation energy loss.
        Generally 1 MeV is a recommended value.
        */
        channelingModel->GetRadiationModel()->SetMinPhotonEnergy(fMinPhotonEnergy);

        /*
        Set the maximal energy in the spectrum to be written into the output file.
        Note: the minimal energy written is equal to fMinPhotonEnergy.
        Note: unlike the minimal energy, the maximal one is just a scoring parameter,
        it does not modify the simulations.
        */
        if(fMaxPhotonEnergySpectrum > fMinPhotonEnergy + DBL_EPSILON)
        {
            channelingModel->GetRadiationModel()->SetMaxPhotonEnergy(fMaxPhotonEnergySpectrum);
        }
        else
        {
            G4cout << "Warning: the maximal energy in BK spectrum <= the minimal energy." << G4endl;
            G4cout << "The maximal energy is default now." << G4endl;
            G4cout << " "<< G4endl;
        }

        /*
        Set the maximal energy in the spectrum to be written into the output file
        */
        channelingModel->GetRadiationModel()->SetNBinsSpectrum(fNBinsSpectrum);

        /*
        Set the angular size of virtual round collimator to accumulate the spectrum the output file
        to be written into the output file. Default is infinite => no collimator.
        */
        channelingModel->GetRadiationModel()->
            SetRoundVirtualCollimator(fVirtualCollimatorHalfSize,-fAngleX,fAngleY);

        /*
        Set the number of trajectory steps after which the radiation probability
        check (whether the probability is below or above of the threshold) is performed;
        at the first iteration also the sampling photons are generated by using
        the angles of this first part of the trajectory as an argument:
        too higher number of steps reduces the accuracy due to considerable excess
        of the single radiation probability threshold;
        too low number may reduce the accuracy of the angular distribution of
        sampling photons.
        Generally the range between 1000-10000 steps is recommended.
        */
        channelingModel->GetRadiationModel()->
                SetNSmallTrajectorySteps(fNSmallTrajectorySteps);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

