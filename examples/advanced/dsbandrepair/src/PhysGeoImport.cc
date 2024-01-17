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
/// \file PhysGeoImport.cc
/// \brief Implementation of the PhysGeoImport class

#include "PhysGeoImport.hh"
#include "Randomize.hh"
#include "G4Sphere.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysGeoImport::PhysGeoImport() :
    fFactor(1.),
    fGeoName("")
{
    DefineMaterial();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysGeoImport::PhysGeoImport(G4bool isVisu) :
    fIsVisu(isVisu),
    fFactor(1.),
    fGeoName("")
{
    if (fIsVisu) G4cout<<" For Checking Visualization!"<<G4endl;
    DefineMaterial();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* PhysGeoImport::CreateLogicVolume(const G4String& fileName,G4String& voxelName)
{
    // The idea is to return a pointer to the mother volume of the input file as output.
    // Create a temporary material map to assign materials

    std::map<G4String, G4Material*> materials;

    // Loop on each material to build the map
    for(G4int i=0, ie=fMaterialVect.size(); i<ie; ++i)
    {
        G4String matName("");

        if(fMaterialVect[i]->GetName()=="adenine_PU") matName = "adenine";
        else if(fMaterialVect[i]->GetName()=="thymine_PY") matName = "thymine";
        else if(fMaterialVect[i]->GetName()=="guanine_PU") matName = "guanine";
        else if(fMaterialVect[i]->GetName()=="cytosine_PY") matName = "cytosine";
        else if(fMaterialVect[i]->GetName()=="backbone_THF") matName = "deoxyribose";
        else if(fMaterialVect[i]->GetName()=="backbone_TMP") matName = "phosphate";
        else if(fMaterialVect[i]->GetName()=="G4_WATER") matName = "water";
        else if(fMaterialVect[i]->GetName()=="homogeneous_dna") matName = "homogeneous_dna";
        else if(fMaterialVect[i]->GetName()=="G4_Galactic") matName = "vacuum";
        else
        {
            matName = fMaterialVect[i]->GetName();
        }

        materials[matName] = fMaterialVect[i];
    }


    // Parse the input file

    voxelName = ParseFile(fileName);

    // Create the volumes

    G4String boxNameSolid = fGeoName+"_solid";
    G4Box* box_solid = new G4Box(boxNameSolid, fSize/2, fSize/2, fSize/2);

    G4String boxNameLogic = fGeoName+"_logic";
    G4LogicalVolume* box_logic = new G4LogicalVolume(box_solid, fpWater, boxNameLogic);

    // sort the molecules
    std::sort(fMolecules.begin(), fMolecules.end() );
   
    // Loop on all the parsed molecules
    for(G4int i=0, ie=fMolecules.size(); i<ie; ++i)
    {
        // Retrieve general molecule informations
        //
        G4String name = fMolecules[i].fName;
        G4String materialName = "water";// fMolecules[i].fMaterial;
        G4double radius = fMolecules[i].fRadius;
        G4double waterRadius = fMolecules[i].fRadiusWater;
        G4ThreeVector moleculePosition = fMolecules[i].fPosition;
        G4int copyNum = fMolecules[i].fCopyNumber;

        // Water hydration shell volume part

        G4Orb* moleculeWater_solid = 0;
        G4VSolid* moleculeWaterCut_solid = 0;
        G4LogicalVolume* moleculeWater_logic = 0;

        // If water radius != 0 then we have a water hydration shell
        G4double tol = 0.0001;
        if(waterRadius > (0 + tol)*nm)
        {
            G4Material* mat = 0;
            
            mat=materials[materialName];

            G4String nameWaterSolid = name+"_water_solid";
            moleculeWater_solid = new G4Orb(nameWaterSolid, waterRadius);
            moleculeWaterCut_solid = CreateCutSolid(moleculeWater_solid, fMolecules[i], fMolecules, false);

            G4String nameWaterLogic = name+"_water_logic";
            moleculeWater_logic = new G4LogicalVolume(moleculeWaterCut_solid, mat, nameWaterLogic);

            G4String nameWaterPhys = name+"_water_phys";
            new G4PVPlacement(0, moleculePosition, moleculeWater_logic, nameWaterPhys, box_logic, false, copyNum);
        }

        // Dna volume part

        G4Orb* molecule_solid = 0;
        G4VSolid* moleculeCut_solid = 0;
        G4LogicalVolume* molecule_logic = 0;

        G4String nameSolid = fMolecules[i].fName+"_solid";
        molecule_solid = new G4Orb(nameSolid, radius);
        moleculeCut_solid = CreateCutSolid(molecule_solid, fMolecules[i], fMolecules, true);

        G4String nameLogic = name+"_logic";
        molecule_logic = new G4LogicalVolume(moleculeCut_solid, materials[materialName], nameLogic);

        // If there was a water hydration shell volume then the current dna volume is
        // placed within it and, thus, its relative coordinates are 0,0,0.
        if(waterRadius > (0 + tol)*nm)
        {
            G4ThreeVector position(0.,0.,0.);
            G4String namePhys = name+"_phys";
            new G4PVPlacement(0, position, molecule_logic, namePhys, moleculeWater_logic, false, copyNum);
        }
        // If not, coordinates are those of the molecule
        else
        {
            G4ThreeVector position = moleculePosition;
            G4String namePhys = name+"_phys";
            new G4PVPlacement(0, position, molecule_logic, namePhys, box_logic, false, copyNum);
        }
    }

    // Clear the containers
    fMolecules.clear();
    fRadiusMap.clear();
    fWaterRadiusMap.clear();

    return box_logic;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String PhysGeoImport::ParseFile(const G4String& fileName)
{
    G4cout<<"Start parsing of "<<fileName<<G4endl;

    // Clear the containers
    fMolecules.clear();
    fRadiusMap.clear();
    fWaterRadiusMap.clear();

    // Setup the input stream
    std::ifstream file(fileName.c_str());

    // Check if the file was correctly opened
    if(!file.is_open())
    {
        // Geant4 exception
        G4String msg = fileName+" could not be opened";
        G4Exception("PhysGeoImport::ParseFile()", "Geo_InputFileNotOpened", FatalException, msg);
    }

    // Define the line string variable
    G4String line;
    G4int preBpCpNo = -1;
    G4int preHistoneCpNo =-1;
    G4int NoBp = 0;
    G4int NoHistone = 0;
    // Read the file line per line
    while(std::getline(file, line) )
    {
        // Check the line to determine if it is empty
        if(line.empty()) continue; // skip the line if it is empty

        // Data string stream
        std::istringstream issLine(line);

        // String to determine the first letter/word
        G4String firstItem;

        // Put the first letter/word within the string
        issLine >> firstItem;

        // Check first letter to determine if the line is data or comment
        if(firstItem=="#") continue; // skip the line if it is comment

        // Use the file
        else if(firstItem=="_Name")
        {
            G4String name;
            issLine >> name;

            fGeoName = name;
        }
        else if(firstItem=="_Size")
        {
            G4double size;
            issLine >> size;
            size *= fFactor*nm;

            fSize = size;
        }
        else if(firstItem == "_Version")
        {

        }
        else if(firstItem=="_Number")
        {

        }
        else if(firstItem=="_Radius")
        {
            G4String name;
            issLine >> name;

            G4double radius;
            issLine >> radius;
            radius *= fFactor*nm;

            G4double waterRadius;
            issLine >> waterRadius;
            waterRadius *= fFactor*nm;

            fRadiusMap[name] = radius;
            fWaterRadiusMap[name] = waterRadius;
        }
        else if(firstItem=="_pl")
        {
            G4String name;
            issLine >> name;

            G4String material;
            issLine >> material;

            G4int strand;
            issLine >> strand;

            G4int copyNumber;
            issLine >> copyNumber;

            if (name == "histone") {
                if (copyNumber != preHistoneCpNo ) {
                    NoHistone ++;
                    preHistoneCpNo = copyNumber;
                }
            } else {
                if (copyNumber != preBpCpNo) {
                    NoBp++;
                    preBpCpNo = copyNumber;
                }
            }

            G4double x;
            issLine >> x;
            x *= fFactor*nm;

            G4double y;
            issLine >> y;
            y *= fFactor*nm;

            G4double z;
            issLine >> z;
            z *= fFactor*nm;

            Molecule molecule(name, copyNumber, G4ThreeVector(x, y, z), fRadiusMap[name], fWaterRadiusMap[name], material, strand);

            fMolecules.push_back(molecule);
        }
        else
        {
            // Geant4 exception
            G4String msg = firstItem+" is not defined in the parser. Check the input file: "+fileName+".";
            G4Exception("PhysGeoImport::ParseFile()", "Geo_WrongParse", FatalException, msg);
        }
    }

    // Close the file once the reading is done
    file.close();

    G4cout<<"End parsing of "<<fileName<<G4endl;
    fVoxelNbHistoneMap.insert({fGeoName,NoHistone});
    fVoxelNbBpMap.insert({fGeoName,NoBp});
    return fGeoName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VSolid* PhysGeoImport::CreateCutSolid(G4Orb *solidOrbRef,
                                               Molecule &molRef,
                                               std::vector<Molecule> &molList,
                                               G4bool in)
{
    // The idea behing this method is to cut overlap volumes by selecting one of them (the reference) and checking all the other volumes (the targets).
    // If a reference and a target volumes are close enough to overlap they will be cut.
    // The reference is already selected when we enter this method.

    // Use the tiny space to differentiate the frontiers (may not be necessary)
    G4double tinySpace = 0.001*fFactor*nm;

    // Cutted solid to be returned
    G4SubtractionSolid* solidCut(NULL);

    // Some flags
    G4bool isCutted = false;
    G4bool isOurVol = false;

    // Radius of the molecule to cut
    G4double radiusRef;
    if(molRef.fRadiusWater==0) radiusRef = molRef.fRadius;
    else radiusRef = molRef.fRadiusWater;

    // Reference volume position
    G4ThreeVector posRef = molRef.fPosition;
    // Before looping on all the volumes we check if the current 
    //reference volume overlaps with its container voxel boundaries.
    if(std::abs(posRef.x() ) + radiusRef > fSize/2 // along x
            || std::abs(posRef.y() ) + radiusRef > fSize/2 // along y
            || std::abs(posRef.z() ) + radiusRef > fSize/2) // along z
    {
        // If we enter here, then the reference volume overlaps with the boundaries of the container voxel
        // Box used to cut
        G4Box* solidBox = new G4Box("solid_box_for_cut", fSize/2, fSize/2, fSize/2);
        G4ThreeVector posBox;

        // Create a dummy rotation matrix
        G4RotationMatrix *rotMat = new G4RotationMatrix;
        rotMat->rotate(0, G4ThreeVector(0,0,1) );

        // To choose the cut direction

        // Up
        if(std::abs( posRef.y() + radiusRef ) > fSize/2 )
        {
           posBox = -posRef +  G4ThreeVector(0,fSize,0);

            // If the volume is cutted for the first time
            if(!isCutted)
            {
                solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, solidBox, rotMat, posBox);
                isCutted = true;
            }
            // For the other times
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, posBox);
        }

        // Down
        if(std::abs( posRef.y() - radiusRef ) > fSize/2 )
        {
            posBox = -posRef + G4ThreeVector(0,-fSize,0);

            // If the volume is cutted for the first time
            if(!isCutted)
            {
                solidCut = new G4SubtractionSolid("solidCut", solidOrbRef,solidBox,rotMat, posBox);
            isCutted = true;
        }
            // For the other times
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, posBox);
        }

        // Left
        if(std::abs( posRef.x() + radiusRef ) > fSize/2 )
        {
            posBox = -posRef + G4ThreeVector(fSize,0,0);

            // If the volume is cutted for the first time
            if(!isCutted)
            {
                solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, solidBox, rotMat, posBox);
                isCutted = true;
            }
            // For the other times
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, posBox);
        }

        // Right
        if(std::abs( posRef.x() - radiusRef ) > fSize/2 )
        {
            posBox = -posRef + G4ThreeVector(-fSize,0,0);

            // If the volume is cutted for the first time
            if(!isCutted)
            {
                solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, solidBox, rotMat, posBox);
                isCutted = true;
            }
            // For the other times
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, posBox);
        }

        // Forward
        if(std::abs( posRef.z() + radiusRef ) > fSize/2 )
        {
            posBox = -posRef + G4ThreeVector(0,0,fSize);

            // If the volume is cutted for the first time
            if(!isCutted)
            {
                solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, solidBox, rotMat, posBox);
                isCutted = true;
            }
            // For the other times
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, posBox);
        }

        // Backward
        if(std::abs( posRef.z() - radiusRef ) > fSize/2 )
        {
            posBox = -posRef + G4ThreeVector(0,0,-fSize);

            // If the volume is cutted for the first time
            if(!isCutted)
            {
                solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, solidBox, rotMat, posBox);
                isCutted = true;
            }
            // For the other times
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, posBox);
        }
    }

    // Look the other volumes of the voxel
    // Loop on all the target volumes (other volumes with potential overlaps)
    for(G4int i=0, ie=molList.size(); i<ie; ++i)
    {
        G4ThreeVector posTar = molList[i].fPosition;

        G4double rTar = posRef.z();
        G4double zTar = posTar.z();

        if(zTar>rTar+20*fFactor*nm)
        {
            break;
        }
        else if(zTar<rTar-20*fFactor*nm)
        {
            continue;
        }

        // Retrieve current target sphere informations
        G4double radiusTar;
        if(molList[i].fRadiusWater==0) radiusTar = molList[i].fRadius;
        else radiusTar = molList[i].fRadiusWater;

        // Compute the distance reference-target
        G4double distance = std::abs( (posRef - posTar).getR() );

        // Use the distance to check if the current target is also the reference.
        // This can only happen once per loop.
        if(distance==0 && !isOurVol)
        {
            // Target volume is also reference volume.

            // Set the flag
            isOurVol = true;

            // Next iteration
            continue;
        }
        // If the condition is correct more than one time then there is a mistake somewhere.
        else if(distance == 0 && isOurVol)
        {
            G4ExceptionDescription desmsg;
            desmsg << "DetectorConstruction::CreateCutSolid: Two volumes are placed at the same position.";
            G4Exception("ChemGeoImport::BuildGeometry", "", FatalException, desmsg);
        }

        // If the volumes are differents then we want to know if they are
        // close enough to overlap and, thus, to intiate a cut.
        else if(distance <= radiusRef+radiusTar)
        {
            // Volumes are close enough, there will be a cut

            // Box used to cut
            G4Box* solidBox = new G4Box("solid_box_for_cut", 2*radiusTar, 2*radiusTar, 2*radiusTar);

            // This part is tricky.
            // The goal is to calculate the position of the intersection center

            // diff vector to from ref to tar
            G4ThreeVector diff = posTar - posRef;

            // Find the intersection point and add to it half the length of the box used to cut
            G4double d = (std::pow(radiusRef,2)-std::pow(radiusTar,2)+std::pow(distance,2) ) / 
                        (2*distance) + solidBox->GetZHalfLength() - tinySpace;

            // If we are in another volume we choose to double the tinySpace to 
            //differentiate without ambiguities the inner and outer volume frontiers.
            // (may not be necessary)
            if(in) d -= 2*tinySpace;

            // Position of the box in order to achieve the cut.
            // "* ( diff/diff.getR() )" is necessary to get a vector in the right direction as output
            G4ThreeVector pos = d *( diff/diff.getR() );

            // Get the rotation angles because the box used to cut needs to be rotated
            // to give the right "cut angle".
            G4double phi = std::acos(pos.getZ()/pos.getR());
            G4double theta = std::acos( pos.getX() / ( pos.getR()*std::cos(pi/2.-phi) ) );

            if(pos.getY()<0) theta = -theta;

            G4ThreeVector rotAxisForPhi(1*fFactor*nm,0.,0.);
            rotAxisForPhi.rotateZ(theta+pi/2);

            // Create the rotation matrix
            G4RotationMatrix *rotMat = new G4RotationMatrix;
            rotMat->rotate(-phi, rotAxisForPhi);

            // Rotate it again
            G4ThreeVector rotZAxis(0.,0.,1*fFactor*nm);
            rotMat->rotate(theta, rotZAxis);

            // If the volume is cutted for the first time
            if(!isCutted) solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, 
            solidBox, rotMat, pos);

            // For the other times
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, pos);

            // Set the cut flag
            isCutted = true;
        }
    }

    // If there was at least one cut then we return the cutted volume
    if(isCutted) return solidCut;

    // Otherwise, we return the original volume
    else return solidOrbRef;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysGeoImport::DefineMaterial()
{
    G4NistManager * man = G4NistManager::Instance();
    fpWater = man->FindOrBuildMaterial("G4_WATER");
    fMaterialVect.push_back(fpWater);
    fVacuum = man->FindOrBuildMaterial("G4_Galactic");
    fMaterialVect.push_back(fVacuum);

    G4double z, a, density;
    G4String name, symbol;
    G4int nComponents, nAtoms;

    a = 12.0107*g/mole;
    G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a);
    a = 1.00794*g/mole;
    G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);
    a = 15.9994*g/mole;
    G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);
    a = 14.00674*g/mole;
    G4Element* elN  = new G4Element(name="Nitrogen"  ,symbol="N" , z= 7., a);
    a = 30.973762*g/mole;
    G4Element* elP  = new G4Element(name="Phosphorus"  ,symbol="P" , z= 15., a);

    density = 1.*g/cm3;

    fTHF = new G4Material("THF", density, nComponents=3);
    fTHF->AddElement(elC, nAtoms=4);
    fTHF->AddElement(elH, nAtoms=8);
    fTHF->AddElement(elO, nAtoms=1);

    fPY = new G4Material("PY", density, nComponents=3);
    fPY->AddElement(elC, nAtoms=4);
    fPY->AddElement(elH, nAtoms=4);
    fPY->AddElement(elN, nAtoms=2);

    fPU = new G4Material("PU", density, nComponents=3);
    fPU->AddElement(elC, nAtoms=5);
    fPU->AddElement(elH, nAtoms=4);
    fPU->AddElement(elN, nAtoms=4);

    fTMP = new G4Material("TMP", density,4);
    fTMP->AddElement(elC, nAtoms=3);
    fTMP->AddElement(elH, nAtoms=6);
    fTMP->AddElement(elP, nAtoms=1);
    fTMP->AddElement(elO, nAtoms=4);

    fMaterialVect.push_back(fTHF);
    fMaterialVect.push_back(fPY);
    fMaterialVect.push_back(fPU);
    fMaterialVect.push_back(fTMP);

    // Mixture creation
    fSugarMixt = new G4Material("sugarMixt", density, nComponents=2);
    fSugarMixt->AddMaterial(fTHF,0.5);
    fSugarMixt->AddMaterial(fTMP,0.5);
    fMaterialVect.push_back(fSugarMixt);

    // Real DNA materials creation
    //
    // Density calculation:
    //
    // d=m/V and M=m/n <=> m = M*n
    // <=> d = M*n/V
    // <=> d = (M*nb)/(Na*V)

    G4double MBackbone = (5.*12.0104 + 10.*1.00794 + 4.*15.9994)*g/mole;
    G4double VBackbone = 1.823912/20. * std::pow(1.e-9, 3.) *m3;
    G4double densityBaTHF = (MBackbone/(g/mole) * 1. )/(6.022e23*VBackbone/cm3) *g/cm3; // g/cm3

    fDeoxyribose = new G4Material("backbone_THF", densityBaTHF, nComponents=3);
    fDeoxyribose->AddElement(elC, nAtoms=5);
    fDeoxyribose->AddElement(elH, nAtoms=10);
    fDeoxyribose->AddElement(elO, nAtoms=4);
    fMaterialVect.push_back(fDeoxyribose);

    MBackbone = (3.*1.00794 + 1.*30.973762 + 4.*15.9994)*g/mole;
    VBackbone = 1.19114/20. * std::pow(1.e-9, 3.) *m3;
    G4double densityBaTMP = (MBackbone/(g/mole) * 1. )/(6.022e23*VBackbone/cm3) *g/cm3; // g/cm3

    fPhosphate = new G4Material("backbone_TMP", densityBaTMP,3);
    fPhosphate->AddElement(elH, nAtoms=3);
    fPhosphate->AddElement(elP, nAtoms=1);
    fPhosphate->AddElement(elO, nAtoms=4);
    fMaterialVect.push_back(fPhosphate);

    G4double MCytosine = (4.*12.0104 + 5.*1.00794 + 3.*14.0067 + 1.*15.9994)*g/mole;
    G4double VBase = 1.8527205/20. * std::pow(1.e-9, 3.) *m3;
    G4double densityCy = (MCytosine/(g/mole) * 1. )/(6.022e23*VBase/cm3) *g/cm3; // g/cm3

    fCytosine_PY = new G4Material("cytosine_PY", densityCy, nComponents=4);
    fCytosine_PY->AddElement(elC, nAtoms=4);
    fCytosine_PY->AddElement(elH, nAtoms=5);
    fCytosine_PY->AddElement(elN, nAtoms=3);
    fCytosine_PY->AddElement(elO, nAtoms=1);
    fMaterialVect.push_back(fCytosine_PY);

    G4double MThy = (5.*12.0104 + 6.*1.00794 + 2.*14.0067 + 2.*15.9994)*g/mole;
    VBase = 1.8527205/20. * std::pow(1.e-9, 3.) *m3;
    densityCy = (MThy/(g/mole) * 1. )/(6.022e23*VBase/cm3) *g/cm3; // g/cm3

    fThymine_PY = new G4Material("thymine_PY", densityCy, nComponents=4);
    fThymine_PY->AddElement(elC, nAtoms=5);
    fThymine_PY->AddElement(elH, nAtoms=6);
    fThymine_PY->AddElement(elN, nAtoms=2);
    fThymine_PY->AddElement(elO, nAtoms=2);
    fMaterialVect.push_back(fThymine_PY);

    G4double MAdenine = (5.*12.0104 + 5.*1.00794 + 5.*14.0067)*g/mole;
    VBase = 1.8527205/20. * std::pow(1.e-9, 3.) *m3;
    G4double densityAd = (MAdenine/(g/mole) * 1. )/(6.022e23*VBase/cm3) *g/cm3; // g/cm3

    fAdenine_PU = new G4Material("adenine_PU", densityAd, nComponents=3);
    fAdenine_PU->AddElement(elC, nAtoms=5);
    fAdenine_PU->AddElement(elH, nAtoms=5);
    fAdenine_PU->AddElement(elN, nAtoms=5);
    fMaterialVect.push_back(fAdenine_PU);

    MAdenine = (5.*12.0104 + 5.*1.00794 + 5.*14.0067 + 1.*15.999)*g/mole;
    VBase = 1.8527205/20. * std::pow(1.e-9, 3.) *m3;
    densityAd = (MAdenine/(g/mole) * 1. )/(6.022e23*VBase/cm3) *g/cm3; // g/cm3

    fGuanine_PU = new G4Material("guanine_PU", densityAd, nComponents=4);
    fGuanine_PU->AddElement(elC, nAtoms=5);
    fGuanine_PU->AddElement(elH, nAtoms=5);
    fGuanine_PU->AddElement(elN, nAtoms=5);
    fGuanine_PU->AddElement(elO, nAtoms=1);
    fMaterialVect.push_back(fGuanine_PU);

    fHomogeneous_dna = new G4Material("homogeneous_dna", 19.23e-21/(12.30781*1.e-21)  *g/cm3, nComponents=7);//21
    fHomogeneous_dna->AddMaterial(fpWater, 0.37);
    fHomogeneous_dna->AddMaterial(fDeoxyribose, 0.23);
    fHomogeneous_dna->AddMaterial(fPhosphate, 0.17);
    fHomogeneous_dna->AddMaterial(fAdenine_PU, 0.06);
    fHomogeneous_dna->AddMaterial(fGuanine_PU, 0.07);
    fHomogeneous_dna->AddMaterial(fCytosine_PY, 0.05);
    fHomogeneous_dna->AddMaterial(fThymine_PY, 0.05);

    fMaterialVect.push_back(fHomogeneous_dna);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* PhysGeoImport::CreateNucleusLogicVolume(const G4String& fileName)
{
    // Here, we parse the nucleus file and build a logical Geant4 volume from it.

    G4cout<<"Start the nucleus creation..."<<G4endl;

    G4VSolid* nucleusSolid = 0;
    G4LogicalVolume* nucleusLogic = 0;

    G4bool isNucleusBuild (false);

    // Open the world file
    std::ifstream nucleusFile;
    nucleusFile.open(fileName.c_str());

    // Check if the file was correctly opened
    if(!nucleusFile.is_open())
    {
        // Geant4 exception
        G4String msg = fileName+" could not be opened";
        G4Exception("PhysGeoImport::ParseFile()", "Geo_InputFileNotOpened", FatalException, msg);
    }

    // Define the line string variable
    G4String line;

    // Read the file line per line
    while(std::getline(nucleusFile, line) && !isNucleusBuild)
    {
        // Check the line to determine if it is empty
        if(line.empty()) continue; // skip the line if it is empty

        // Data string stream
        std::istringstream issLine(line);

        // String to determine the first letter/word
        G4String firstItem;

        // Put the first letter/word within the string
        issLine >> firstItem;

        // Check first letter to determine if the line is data or comment
        if(firstItem=="#") continue; // skip the line if it is comment

        // Use the file

        else if(firstItem == "_Name")
        {
            issLine >> fNucleusName;
        }

        else if(firstItem=="_Type")
        {
            issLine >> fNucleusType;

            if (fNucleusType == "Ellipsoid")
            {
                G4float semiX, semiY, semiZ;

                issLine >> semiX >> semiY >> semiZ;

                semiX *= fFactor*nm;
                semiY *= fFactor*nm;
                semiZ *= fFactor*nm;

                fNucleusData["SemiX"] = semiX;
                fNucleusData["SemiY"] = semiY;
                fNucleusData["SemiZ"] = semiZ;

                nucleusSolid = new G4Ellipsoid(fNucleusName, semiX, semiY, semiZ);
                nucleusLogic = new G4LogicalVolume(nucleusSolid, fpWater, "nucleus_logic");
                G4double r3 = semiX*semiY*semiZ/m3;
                fNucleusVolume = 4.0*pi*r3/3.0;
            } else if (fNucleusType == "EllipticCylinder") {
            	G4float semiX, semiY, semiZ;
		
                issLine >> semiX >> semiY >> semiZ;
		
                semiX *= fFactor*nm;
                semiY *= fFactor*nm;
                semiZ *= fFactor*nm;
		
                fNucleusData["SemiX"] = semiX;
                fNucleusData["SemiY"] = semiY;
                fNucleusData["SemiZ"] = semiZ;

                nucleusSolid = new G4EllipticalTube(fNucleusName, semiX, semiY, semiZ);
                nucleusLogic = new G4LogicalVolume(nucleusSolid, fpWater, "nucleus_logic");
                G4double r3 = 2*semiX*semiY*semiZ/m3; // multiplied by 2 cuz semiZ is hafl of height
                fNucleusVolume = pi*r3;
            } else if (fNucleusType == "Spherical") {
                G4float radius;
                issLine >> radius;
                radius *= fFactor*nm;
                fNucleusData["SemiX"] = radius;
                fNucleusData["SemiY"] = radius;
                fNucleusData["SemiZ"] = radius;
                nucleusSolid = new G4Orb(fNucleusName,radius);
                nucleusLogic = new G4LogicalVolume(nucleusSolid, fpWater, "nucleus_logic");
                G4double r3 = radius*radius*radius/m3;
                fNucleusVolume = 4.0*pi*r3/3.0;
            } else {
                G4String msg =fNucleusType+" is not registered.";
                G4Exception("PhysGeoImport::CreateNucleusLogicVolume()", 
                "Geo_InputFile_Unknown_Nucleus", FatalException, msg);
            }

            isNucleusBuild = true;
        }
    }

    nucleusFile.close();

    if(!isNucleusBuild)
    {
        G4String msg = "Nucleus data were not found for "+fNucleusType;
        G4Exception("PhysGeoImport::CreateNucleusLogicVolume()", 
        "Geo_InputFile_Unknown_Nucleus", FatalException, msg);
    }

    return nucleusLogic;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<Voxel>* PhysGeoImport::CreateVoxelsData(const G4String& fileName)
{
    fTotalNbBpPlacedInGeo = 0;
    fTotalNbHistonePlacedInGeo = 0;
    G4cout<<"Start the voxel data generation..."<<G4endl;

    std::vector<Voxel>* voxels = new std::vector<Voxel>;

    // Clear any previously loaded world
    voxels->clear();

    voxels->reserve(1000000);

    // Open the world file
    std::ifstream nucleusFile;
    nucleusFile.open(fileName.c_str());

    // Check if the file was correctly opened
    if(!nucleusFile.is_open())
    {
        // Geant4 exception
        G4String msg = fileName+" could not be opened";
        G4Exception("PhysGeoImport::ParseFile()", "Geo_InputFileNotOpened", FatalException, msg);
    }

    // Initialise the copy number to 0
    G4int voxelCopyNumber = 0;

    // Define the line string variable
    G4String line;

    auto whatchromatintype = [] (Voxel::VoxelType voxt) -> ChromatinType {
        if (voxt == Voxel::Straight || voxt == Voxel::Left || voxt == Voxel::Right ||
            voxt == Voxel::Up || voxt == Voxel::Down) return ChromatinType::fHeterochromatin;
        else if (voxt == Voxel::Straight2 || voxt == Voxel::Left2 || voxt == Voxel::Right2 ||
            voxt == Voxel::Up2 || voxt == Voxel::Down2) return ChromatinType::fEuchromatin;
        else return ChromatinType::fUnspecified;   
    };
    // Read the file line per line
    while(std::getline(nucleusFile, line) )
    {
        // Check the line to determine if it is empty
        if(line.empty()) continue; // skip the line if it is empty

        // Data string stream
        std::istringstream issLine(line);

        // String to determine the first letter/word
        G4String firstItem;

        // Put the first letter/word within the string
        issLine >> firstItem;

        // Check first letter to determine if the line is data or comment
        if(firstItem=="#") continue; // skip the line if it is comment

        // Use the file

        else if(firstItem == "_pl")
        {
            // Set the flags

            G4String name;
            issLine >> name;

            G4int chromoCpN;
            issLine >> chromoCpN;

            G4int domainCpN;
            issLine >> domainCpN;

            G4double x;
            issLine >> x;
            x *= fFactor*nm;

            G4double y;
            issLine >> y;
            y *= fFactor*nm;

            G4double z;
            issLine >> z;
            z *= fFactor*nm;
            
            G4double xx;
            issLine >> xx;
            
            G4double yx;
            issLine >> yx;
            
            G4double zx;
            issLine >> zx;
            
            G4double xy;
            issLine >> xy;
            
            G4double yy;
            issLine >> yy;
            
            G4double zy;
            issLine >> zy;
            
            G4double xz;
            issLine >> xz;
            
            G4double yz;
            issLine >> yz;
            
            G4double zz;
            issLine >> zz;
            

            G4RotationMatrix* rot = new G4RotationMatrix(G4ThreeVector(xx,xy,xz),
            G4ThreeVector(yx,yy,yz),G4ThreeVector(zx,zy,zz));

            Voxel::VoxelType type = Voxel::Other;

            if(name=="voxelStraight" || name=="VoxelStraight")
            {
                type = Voxel::Straight;
            }
            else if(name=="voxelUp" || name=="VoxelUp")
            {
                type = Voxel::Up;
            }
            else if(name=="voxelDown" || name=="VoxelDown")
            {
                type = Voxel::Down;
            }
            else if(name=="voxelRight" || name=="VoxelRight")
            {
                type = Voxel::Right;
            }
            else if(name=="voxelLeft" || name=="VoxelLeft")
            {
                type = Voxel::Left;
            }
            else if(name=="voxelStraight2" || name=="VoxelStraight2")
            {
                type = Voxel::Straight2;
            }
            else if(name=="voxelUp2" || name=="VoxelUp2")
            {
                type = Voxel::Up2;
            }
            else if(name=="voxelDown2" || name=="VoxelDown2")
            {
                type = Voxel::Down2;
            }
            else if(name=="voxelRight2" || name=="VoxelRight2")
            {
                type = Voxel::Right2;
            }
            else if(name=="voxelLeft2" || name=="VoxelLeft2")
            {
                type = Voxel::Left2;
            }
            else
            {
                G4ExceptionDescription msg ;
                msg << "The voxel named "<<name<<" is not registered in the parser";
                G4Exception("PhysGeoImport::CreateVoxelsData", "", FatalException, msg, "");
            }

            voxels->push_back(Voxel(voxelCopyNumber, chromoCpN, domainCpN, type, 
            G4ThreeVector(x,y,z), rot) );
            fTotalNbBpPlacedInGeo += fVoxelNbBpMap[name];
            fTotalNbHistonePlacedInGeo += fVoxelNbHistoneMap[name];
            voxelCopyNumber++;
            auto chromatintype =  whatchromatintype(type);
            if (fChromatinTypeCount.find(chromatintype) == fChromatinTypeCount.end()) 
                fChromatinTypeCount.insert({chromatintype,1});
            else fChromatinTypeCount[chromatintype]++;
        }
    }

    nucleusFile.close();

    G4cout<<"End the voxel data generation"<<G4endl;

    return voxels;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......