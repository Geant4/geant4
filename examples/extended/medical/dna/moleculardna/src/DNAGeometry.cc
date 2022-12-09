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
/// file: DNAGeometry.cc
/// brief:
/*
 * This class builds the DNA geometry. It interacts in a non-standard way
 * with the DNAWorld class, to create a physical DNA geometry
 * in a parallel world.
 *
 */
#include "DNAGeometry.hh"
#include "DNAGeometryMessenger.hh"
#include "OctreeNode.hh"
#include "PlacementVolumeInfo.hh"
#include "VirtualChromosome.hh"
#include "G4Box.hh"
#include "G4Ellipsoid.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PhysicalConstants.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4NistManager.hh"
#include <fstream>
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool compareLVByName::operator()(const G4LogicalVolume* lhs,
                                 const G4LogicalVolume* rhs) const
{
  return lhs->GetName() < rhs->GetName();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAGeometry::DNAGeometry()
  : G4VDNAMolecularGeometry()
{
  fpMessenger        = new DNAGeometryMessenger(this);
  fpChromosomeMapper = new ChromosomeMapper();

  // Create Damage Model and add some default values
  fpDamageModel = new DamageModel();
  G4NistManager* man  = G4NistManager::Instance();
  fpSugarMaterial     = man->FindOrBuildMaterial("G4_DNA_DEOXYRIBOSE");
  fpMediumMaterial    = man->FindOrBuildMaterial("G4_WATER");
  fpPhosphateMaterial = man->FindOrBuildMaterial("G4_DNA_PHOSPHATE");
  fpGuanineMaterial   = man->FindOrBuildMaterial("G4_DNA_GUANINE");
  fpAdenineMaterial   = man->FindOrBuildMaterial("G4_DNA_ADENINE");
  fpThymineMaterial   = man->FindOrBuildMaterial("G4_DNA_THYMINE");
  fpCytosineMaterial  = man->FindOrBuildMaterial("G4_DNA_CYTOSINE");
  fpHistoneMaterial   = man->FindOrBuildMaterial("G4_WATER");
  fpDNAPhysicsWorld   = new DNAWorld();
  fpDrawCellAttrs     = new G4VisAttributes();
  fpDrawCellAttrs->SetColor(G4Color(0, 0, 1, 0.2));
  fpDrawCellAttrs->SetForceSolid(true);
  // set default molecule sizes
  fMoleculeSizes["phosphate"] =
    G4ThreeVector(2.282354, 2.282354, 2.282354) * angstrom;
  fMoleculeSizes["sugar"] =
    G4ThreeVector(2.632140, 2.632140, 2.632140) * angstrom;
  fMoleculeSizes["guanine"] =
    G4ThreeVector(3.631503, 3.799953, 1.887288) * angstrom;
  fMoleculeSizes["cytosine"] =
    G4ThreeVector(3.597341, 3.066331, 1.779361) * angstrom;
  fMoleculeSizes["adenine"] =
    G4ThreeVector(3.430711, 3.743504, 1.931958) * angstrom;
  fMoleculeSizes["thymine"] =
    G4ThreeVector(4.205943, 3.040448, 2.003359) * angstrom;
  fMoleculeSizes["histone"] = G4ThreeVector(25, 25, 25) * angstrom;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAGeometry::~DNAGeometry() { delete fpMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAGeometry::BuildDNA(G4LogicalVolume* vol)
{
  // TODO: Add Assertion tests to make sure files have been loaded
  fpDNAVolumeChem = vol;

  // TODO: Make a complete copy of vol and it's physical volume
  auto* volAsEllipsoid = dynamic_cast<G4Ellipsoid*>(vol->GetSolid());
  auto* volAsBox             = dynamic_cast<G4Box*>(vol->GetSolid());

  G4VSolid* cloneSolid;
  if(volAsEllipsoid != nullptr)
  {
    cloneSolid = new G4Ellipsoid(*((G4Ellipsoid*) vol->GetSolid()));
  }
  else if(volAsBox != nullptr)
  {
    cloneSolid = new G4Box(*((G4Box*) vol->GetSolid()));
  }
  else
  {
    G4ExceptionDescription errmsg;
    errmsg
      << "An invalid physical volume shape has been used to hold the DNA volume"
      << G4endl;
    G4Exception("MolecularDnaGeometry", "ERR_INVALID_DNA_SHAPE", FatalException,
                errmsg);
    return;
  }

  fpDNAVolumePhys = new G4LogicalVolume(cloneSolid, vol->GetMaterial(),
                                        "DNAPhysLV");  // dousatsu
  FillParameterisedSpace();
  fpDNAPhysicsWorld->SetDNAVolumePointer(fpDNAVolumePhys);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAGeometry::FillParameterisedSpace()
{
  // Create map for each logical volume with its parametrisation
  // LV, placement_index, copy number, position, rotation
  std::vector<placement> placementVector;
  std::map<G4LogicalVolume*, G4int> currentCopyNumber;
  G4int numberOfChains = 0;

  // read in and create the LV's specified in the macro
  std::map<G4String, G4LogicalVolume*> namesToLVs;
  for(const auto& it : fVoxelNames)
  {
    G4String thisVoxelName = it.first;
    G4String thisVoxelFile = it.second;

    auto voxel                = LoadVoxelVolume(thisVoxelName, thisVoxelFile);
    G4LogicalVolume* lv       = voxel.first;
    fInfoMap[lv]              = voxel.second;
    namesToLVs[thisVoxelName] = lv;
    currentCopyNumber[lv]     = 0;
    numberOfChains            = this->GetPVInfo(lv)->GetNumberOfChains();
    if(!(numberOfChains == 1 || numberOfChains == 4 || numberOfChains == 8))
    {
      G4cout << "Geometry contains volumes with none of 1/4/8"
             << " chains. This can cause errors" << G4endl;
    }
  }

  // Go through the voxel types to build a placement vector
  for(unsigned ii = 0; ii < fVoxelTypes.size(); ii++)
  {
    G4LogicalVolume* lv = namesToLVs[fVoxelTypes[ii]];

    placement thisPlacement =
      std::make_tuple(lv, fVoxelIndices[ii], currentCopyNumber[lv]++,
                      fVoxelPositions[ii], fVoxelRotations[ii]);
    placementVector.push_back(thisPlacement);
  }

  // Do each placement.
  std::map<G4String, int64_t> histonesPerChromosome;   // dousatsu
  std::map<G4String, int64_t> basePairsPerChromosome;  // dousatsu
  // std::map<G4String, long long> basePairsPerChromosome;//ORG
  for(auto it = placementVector.begin(); it != placementVector.end(); it++)
  {
    G4LogicalVolume* thisLogical = std::get<0>(*it);
    if(thisLogical == nullptr)
    {
      G4Exception("DNAGeometry::FillParameterisedSpace",
                  "ERR_LOG_VOL_EXISTS_FALSE", FatalException,
                  "At least one of the placement names"
                  " specified doesn't exist");
    }
    G4int thisIndex  = std::get<1>(*it);
    G4int thisCopyNo = std::get<2>(*it);

    G4ThreeVector thisPosition = std::get<3>(*it);
    G4ThreeVector thisRotation = std::get<4>(*it);
    std::stringstream ss;
    ss << thisLogical->GetName() << "-" << thisIndex;

    auto inv = new G4RotationMatrix();
    inv->rotateX(thisRotation.getX());
    inv->rotateY(thisRotation.getY());
    inv->rotateZ(thisRotation.getZ());
    auto pRot = new G4RotationMatrix(inv->inverse());
    delete inv;

    G4bool isPhysicallyPlaced = false;
    VirtualChromosome* thisChromosome =
      this->GetChromosomeMapper()->GetChromosome(thisPosition);
    if(thisChromosome != nullptr)
    {
      isPhysicallyPlaced = true;
      // Place for Physics
      new G4PVPlacement(pRot, thisPosition, thisLogical, ss.str(),
                        fpDNAVolumePhys, false, thisCopyNo,
                        this->GetOverlaps());
      // Place for Chemistry
      new G4PVPlacement(
        pRot, thisPosition, this->GetMatchingChemVolume(thisLogical), ss.str(),
        fpDNAVolumeChem, false, thisCopyNo, this->GetOverlaps());

      // dousatsu==============
      if(histonesPerChromosome.find(thisChromosome->GetName()) ==
         histonesPerChromosome.end())
      {
        histonesPerChromosome[thisChromosome->GetName()] = 0;
      }
      histonesPerChromosome[thisChromosome->GetName()] +=
        this->GetPVInfo(thisLogical)->GetTotalHistones();
      // dousatsu==============

      if(basePairsPerChromosome.find(thisChromosome->GetName()) ==
         basePairsPerChromosome.end())
      {
        basePairsPerChromosome[thisChromosome->GetName()] = 0;
      }
      basePairsPerChromosome[thisChromosome->GetName()] +=
        this->GetPVInfo(thisLogical)->GetTotalBasePairs();
    }

    // add the placement to the local register, even if it isn't in a
    // chromosome, as this tracks base pair index number.
    if(numberOfChains == 4 || numberOfChains == 8)
    {
      AddFourChainPlacement(it, placementVector.begin(), isPhysicallyPlaced);
    }
    else if(numberOfChains == 1)
    {
      AddSingleChainPlacement(it, placementVector.begin(), isPhysicallyPlaced);
    }
    else
    {
      G4ExceptionDescription errmsg;
      errmsg << "Having none of 1/4/8 chains in a placement is not "
             << "supported, the application will crash" << G4endl;
      G4Exception("DNAGeometry::FillParameterisedSpace",
                  "Number of chains not supported", FatalException, errmsg);
    }
  }

  // dousatsu------------------------------
  G4cout << "Histones placed per chromosome:" << G4endl;
  G4cout << "Chromosome               Histones" << G4endl;
  G4cout << "-----------------------------------" << G4endl;
  for(auto& it : histonesPerChromosome)
  {
    G4cout <<it.first.c_str()<<" : "<< it.second << G4endl;
  }
  // dousatsu------------------------------

  G4cout << "Base Pairs placed per chromosome:" << G4endl;
  G4cout << "Chromosome               Base Pairs" << G4endl;
  G4cout << "-----------------------------------" << G4endl;
  for(auto& it : basePairsPerChromosome)
  {
    G4cout <<it.first.c_str()<<" : "<< it.second << G4endl;
  }
  G4cout << "-------------------------------" << G4endl;

  if(this->GetVerbosity() > 0)
  {
    this->PrintOctreeStats();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAGeometry::AddSingleChainPlacement(
  const std::vector<placement>::iterator it,
  const std::vector<placement>::iterator, G4bool isPhysicallyPlaced)
{
  G4LogicalVolume* thisLogical = std::get<0>(*it);
  G4bool twist                 = fVoxelTwist[thisLogical->GetName()];
  AddNewPlacement(thisLogical, { { 0, 1, 2, 3, 4, 5, 6, 7 } }, twist,
                  isPhysicallyPlaced);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAGeometry::AddFourChainPlacement(
  const std::vector<placement>::iterator it,
  const std::vector<placement>::iterator begin, G4bool isPhysicallyPlaced)
{
  G4LogicalVolume* thisLogical = std::get<0>(*it);
  G4int thisIndex              = std::get<1>(*it);
  // Use info to build the placement transforms.
  std::array<int, 8> global_chain{};
  G4bool twist = fVoxelTwist[thisLogical->GetName()];
  if(it == begin)
  {
    global_chain = { { 0, 1, 2, 3, 4, 5, 6, 7 } };
  }
  else
  {
    // work out global chain
    placement previous         = *(it - 1);
    G4ThreeVector oldRotation  = std::get<4>(previous);
    G4ThreeVector thisRotation = std::get<4>(*it);

    G4RotationMatrix rotnew;
    rotnew.rotateX(thisRotation.getX());
    rotnew.rotateY(thisRotation.getY());
    rotnew.rotateZ(thisRotation.getZ());

    rotnew.rectify();

    G4RotationMatrix rotold;
    rotold.rotateX(oldRotation.getX());
    rotold.rotateY(oldRotation.getY());
    rotold.rotateZ(oldRotation.getZ());

    G4ThreeVector zdiff = rotold.colZ() - rotnew.colZ();
    if(zdiff.mag2() > 0.1)
    {
      rotold = rotold.rotate(halfpi, rotold.colY());
    }
    rotold.rectify();
    // rotnew.colz == rotold.colz now.
    // There exists now a rotation around rotnew.colZ()
    // that transforms rotold to rotnew.
    G4RotationMatrix relative = rotnew * rotold.inverse();
    relative.rectify();
    G4ThreeVector axis = relative.getAxis();
    G4double angle     = relative.getDelta();
    // get the sign of the rotation, and correct for quadrant.
    G4int sign = (axis.dot(rotnew.colZ()) >= 0) ? 1 : -1;
    if(sign < 0) {
      angle = 2 * pi - angle;}
    G4int zrots = ((G4int)(2.05 * angle / pi)) % 4;  // ORG

    global_chain = { { (this->GetGlobalChain(thisIndex - 1, 0) + zrots) % 4,
                       (this->GetGlobalChain(thisIndex - 1, 1) + zrots) % 4,
                       (this->GetGlobalChain(thisIndex - 1, 2) + zrots) % 4,
                       (this->GetGlobalChain(thisIndex - 1, 3) + zrots) % 4,
                       4 + (this->GetGlobalChain(thisIndex - 1, 4) + zrots) % 4,
                       4 + (this->GetGlobalChain(thisIndex - 1, 5) + zrots) % 4,
                       4 + (this->GetGlobalChain(thisIndex - 1, 6) + zrots) % 4,
                       4 + (this->GetGlobalChain(thisIndex - 1, 7) + zrots) %
                             4 } };
    twist        = false;
  }
  this->AddNewPlacement(thisLogical, global_chain, twist, isPhysicallyPlaced);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAGeometry::AddNewPlacement(const G4LogicalVolume* lv,
                                  std::array<int, 8> global_chain, G4bool twist,
                                  G4bool isPhysicallyPlaced)
{
  placement_transform pt;
  std::array<int64_t, 8> start{};  // dousatsu
  std::array<int64_t, 8> end{};    // dousatsu
  // std::array<long long, 8> start;//ORG
  // std::array<long long, 8> end;//ORG
  if(fPlacementTransformations.empty())
  {
    start = { { 0, 0, 0, 0, 0, 0, 0, 0 } };
    end   = { { this->GetPVInfo(lv)->GetPairsOnChain(global_chain[0]),
              this->GetPVInfo(lv)->GetPairsOnChain(global_chain[1]),
              this->GetPVInfo(lv)->GetPairsOnChain(global_chain[2]),
              this->GetPVInfo(lv)->GetPairsOnChain(global_chain[3]) } };
  }
  else
  {
    placement_transform previous       = fPlacementTransformations.back();
    std::array<int64_t, 8> previousEnd = std::get<2>(previous);  // dousatsu
    // std::array<long long, 8> previousEnd = std::get<2>(previous);//ORG
    start = { { previousEnd[0], previousEnd[1], previousEnd[2], previousEnd[3],
                previousEnd[4], previousEnd[5], previousEnd[6],
                previousEnd[7] } };
    end   = {
      { previousEnd[0] + this->GetPVInfo(lv)->GetPairsOnChain(global_chain[0]),
        previousEnd[1] + this->GetPVInfo(lv)->GetPairsOnChain(global_chain[1]),
        previousEnd[2] + this->GetPVInfo(lv)->GetPairsOnChain(global_chain[2]),
        previousEnd[3] + this->GetPVInfo(lv)->GetPairsOnChain(global_chain[3]),
        previousEnd[4] + this->GetPVInfo(lv)->GetPairsOnChain(global_chain[4]),
        previousEnd[5] + this->GetPVInfo(lv)->GetPairsOnChain(global_chain[5]),
        previousEnd[6] + this->GetPVInfo(lv)->GetPairsOnChain(global_chain[6]),
        previousEnd[7] + this->GetPVInfo(lv)->GetPairsOnChain(global_chain[7]) }
    };
  }
  pt = std::make_tuple(global_chain, start, end, twist, isPhysicallyPlaced);
  fPlacementTransformations.push_back(pt);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::pair<G4LogicalVolume*, PlacementVolumeInfo*> DNAGeometry::LoadVoxelVolume(
  const G4String& volumeName, const G4String& filename)
{
  if(this->GetVerbosity() > 0)
  {
    G4cout << "Loading Voxel: " << volumeName << G4endl;
  }

  auto thisInfo = new PlacementVolumeInfo();

  G4double vxdim = 0.5 * fVoxelSideLength.getX();
  G4double vydim = 0.5 * fVoxelSideLength.getY();
  G4double vzdim = 0.5 * fVoxelSideLength.getZ();

  // Make Volumes, and put them into the maps for phys/chem comparison
  auto physicsPhysical = new G4Box(volumeName, vxdim, vydim, vzdim);
  auto physicsLogical =
    new G4LogicalVolume(physicsPhysical, fpMediumMaterial, volumeName);
  auto chemPhysical = new G4Box(volumeName, vxdim, vydim, vzdim);
  auto chemLogical =
    new G4LogicalVolume(chemPhysical, fpMediumMaterial, volumeName);

  if(this->GetDrawCellVolumes())
  {
    chemLogical->SetVisAttributes(fpDrawCellAttrs);
  }

  fChemToPhys[chemLogical]    = physicsLogical;
  fPhysToChem[physicsLogical] = chemLogical;

  physicsLogical->SetSmartless(this->GetSmartless());

  auto thisOctree        = new OctreeNode(G4ThreeVector(0, 0, 0),
                                   G4ThreeVector(vxdim, vydim, vzdim), 1);
  auto thisHistoneOctree = new OctreeNode(
    G4ThreeVector(0, 0, 0), G4ThreeVector(vxdim, vydim, vzdim), 1);

  // open and load file
  if(!utility::Path_exists(filename))
  {
    G4ExceptionDescription errmsg;
    errmsg << "The file: " << filename << " could not be found." << G4endl;
    G4Exception("DNAGeometry::LoadVoxelVolume", "File Not Found",
                FatalException, errmsg);
  }
  std::ifstream infile(filename, std::ifstream::in);
  G4String currentline;
  G4String lower_name;
  G4ThreeVector pos;
  G4ThreeVector rot;
  G4ThreeVector mol_size;
  std::vector<molecule_t> molecules;
  std::vector<molecule_t> histones;
  G4bool file_contains_size     = false;
  G4int line_fields             = 0;
  G4int line_number             = 0;
  G4int uncommented_line_number = 0;

  if(infile.is_open())
  {
    while(getline(infile, currentline))
    {
      if(currentline[0] != '#')
      {
        std::vector<G4String> line = utility::Split(currentline, ' ');
        // validation
        if(uncommented_line_number == 0)
        {
          line_fields        = line.size();
          file_contains_size = (line_fields == 14);
        }
        if((line_fields != 10) && (line_fields != 14))
        {
          G4ExceptionDescription errmsg;
          errmsg << filename << " should have 10 or 14 fields, only got "
                 << line_fields << G4endl;
          G4Exception("DNAGeometry::LoadVoxelVolume", "BadFileFormat",
                      FatalErrorInArgument, errmsg);
        }
        if(line_fields != (G4int) line.size())
        {
          G4ExceptionDescription errmsg;
          errmsg << "Unexpected number of fields on line " << line_number
                 << " of file " << filename << G4endl << "Expected "
                 << line_fields << " fields, but got " << line.size() << G4endl;
          G4Exception("DNAGeometry::LoadVoxelVolume", "BadFileFormat",
                      FatalErrorInArgument, errmsg);
        }
        molecule_t thisMolecule;

        thisMolecule.fname      = (G4String) line.at(0);
        thisMolecule.fshape     = (G4String) line.at(1);
        thisMolecule.fchain_id  = (G4int) std::stoi(line.at(2));
        thisMolecule.fstrand_id = (G4int) std::stoi(line.at(3));
        thisMolecule.fbase_idx  = (int64_t) std::stoll(line.at(4));

        lower_name = (G4String) line.at(0);
        G4StrUtil::to_lower(lower_name);
        // lower_name.toLower();

        G4int size_offset = (file_contains_size) ? 3 : 0;
        if(file_contains_size)
        {
          mol_size = G4ThreeVector(std::stod(line.at(5)) * angstrom,
                                   std::stod(line.at(6)) * angstrom,
                                   std::stod(line.at(7)) * angstrom);
        }
        else
        {
          if(fMoleculeSizes.count(lower_name) == 0)
          {
            G4ExceptionDescription errmsg;
            errmsg << "A molecule size has not been defined for " << lower_name
                   << G4endl;
            G4Exception("DNAGeometry::LoadVoxelVolume", "MoleculeNotDefined",
                        FatalErrorInArgument, errmsg);
          }
          mol_size = fMoleculeSizes.at(lower_name);
        }
        pos = G4ThreeVector(std::stod(line.at(5 + size_offset)) * angstrom,
                            std::stod(line.at(6 + size_offset)) * angstrom,
                            std::stod(line.at(7 + size_offset)) * angstrom);

        rot = G4ThreeVector(std::stod(line.at(8 + size_offset)),
                            std::stod(line.at(9 + size_offset)),
                            std::stod(line.at(10 + size_offset)));

        thisMolecule.fsize     = mol_size;
        thisMolecule.fposition = pos;
        thisMolecule.frotation = rot;

        if(lower_name == "histone")
        {
          histones.push_back(thisMolecule);
        }
        else
        {
          molecules.push_back(thisMolecule);
        }
        uncommented_line_number++;
      }
      line_number++;
    }
    infile.close();
  }
  else
  {
    G4ExceptionDescription errmsg;
    errmsg << "Voxel Position file read error" << G4endl
           << "Name: " << volumeName << G4endl << "Filename: " << filename
           << G4endl << "Format should be: "
           << "# NAME SHAPE CHAIN_ID STRAND_ID BP_INDEX SIZE_X SIZE_Y "
           << "SIZE_Z POS_X POS_Y POS_Z ROT_X ROT_Y ROT_Z" << G4endl << "or"
           << G4endl
           << "# NAME SHAPE CHAIN_ID STRAND_ID BP_INDEX POS_X POS_Y POS_Z "
              "ROT_X ROT_Y ROT_Z"
           << G4endl
           << "Separator is space, lines with '#' as the first  character are "
              "commented"
           << G4endl;

    G4Exception("DNAGeometry::LoadVoxelVolume", "BadFile", FatalErrorInArgument,
                errmsg);
  }

  // We can place the elements now, specify an order that they'll be placed
  // for the union solid cuts.
  G4VPhysicalVolume* pv_placement = nullptr;

  // Place Histones
  G4int his_arr_size = histones.size();
  for(G4int ii = 0; ii != his_arr_size; ++ii)
  {
    molecule_t thisMolecule = histones[ii];
    if(this->GetVerbosity() > 4)
    {
      G4cout << "Placing molecule " << ii << ": " << thisMolecule.fname
             << G4endl;
    }
    // Identify number of histones in voxel volume
    pv_placement = PlaceHistone(physicsLogical, thisMolecule);
    thisHistoneOctree->AddPhysicalVolume(pv_placement);
  }
  thisInfo->SetNHistones(his_arr_size);
  thisInfo->SetHistoneOctree(thisHistoneOctree);

  // Place Molecules
  G4int mol_arr_size = molecules.size();
  for(G4int ii = 0; ii != mol_arr_size; ii++)
  {
    molecule_t thisMolecule = molecules[ii];
    if(this->GetVerbosity() > 4)
    {
      G4cout << "Placing molecule " << ii << ": " << thisMolecule.fname
             << G4endl;
    }
    // Identify number of molecules on chain
    if(thisInfo->GetPairsOnChain(thisMolecule.fchain_id) <
       (thisMolecule.fbase_idx + 1))
    {
      thisInfo->SetPairsOnChain(thisMolecule.fchain_id,
                                thisMolecule.fbase_idx + 1);
    }

    // Place Volume
    if(thisMolecule.fname == "Sugar")
    {
      molecule_t phosphate = molecules[ii - 1];
      // G4int phosph_idx = (thisMolecule.strand_id == 0) ? ii - 1 : ii - 1;
      G4int phosph_idx = ii - 1;
      if((phosph_idx >= 0) && (phosph_idx < (G4int) mol_arr_size))
      {
        phosphate = molecules[phosph_idx];
        if(phosphate.fchain_id != thisMolecule.fchain_id)
        {
          // There is no next sugar, just the voxel wall
          phosphate      = molecule_t();
          phosphate.fname = "wall";
        }
      }
      else
      {
        // There is no next sugar, just the voxel wall
        phosphate      = molecule_t();
        phosphate.fname = "wall";
      }
      pv_placement = PlaceSugar(physicsLogical, thisMolecule, phosphate);
    }
    else if(thisMolecule.fname == "Phosphate")
    {
      molecule_t sugar;
      G4int sugar_idx = (thisMolecule.fstrand_id == 0) ? ii + 7 : ii - 5;
      if((sugar_idx >= 0) && (sugar_idx < mol_arr_size))
      {
        sugar = molecules[sugar_idx];
        if(sugar.fchain_id != thisMolecule.fchain_id)
        {
          // There is no next sugar, just the voxel wall
          sugar      = molecule_t();
          sugar.fname = "wall";
        }
      }
      else
      {
        // There is no next sugar, just the voxel wall
        sugar      = molecule_t();
        sugar.fname = "wall";
      }
      pv_placement = PlacePhosphate(physicsLogical, thisMolecule, sugar);
    }
    else if((thisMolecule.fname == "Guanine") ||
            (thisMolecule.fname == "Adenine") ||
            (thisMolecule.fname == "Cytosine") ||
            (thisMolecule.fname == "Thymine"))
    {
      molecule_t sugar = molecules[ii - 1];
      molecule_t oppBasePair;
      molecule_t nextBasePair;

      G4int nextBPidx;
      G4bool endOfChain = false;
      G4bool nextBPexists;
      if(thisMolecule.fstrand_id == 0)
      {
        oppBasePair  = molecules[ii + 3];
        nextBPexists = (ii + 6 < mol_arr_size);
        if(nextBPexists)
        {
          endOfChain = (thisMolecule.fchain_id != molecules[ii + 6].fchain_id);
        }
        nextBPidx = (nextBPexists && !endOfChain) ? ii + 6 : ii - 6;
      }
      else
      {
        oppBasePair  = molecules[ii - 3];
        nextBPexists = (ii - 6 >= 0);
        if(nextBPexists)
        {
          endOfChain = (thisMolecule.fchain_id != molecules[ii - 6].fchain_id);
        }
        nextBPidx = (nextBPexists && !endOfChain) ? ii - 6 : ii + 6;
      }
      nextBasePair = molecules[nextBPidx];
      pv_placement = PlaceBase(physicsLogical, thisMolecule, oppBasePair,
                               nextBasePair, sugar);
    }
    else
    {
      G4ExceptionDescription errmsg;
      errmsg << "Molecule name: " << thisMolecule.fname << " is invalid."
             << G4endl << "Aborting" << G4endl;
      G4Exception("DNAGeometry::LoadVoxelVolume", "Invalid molecule name",
                  FatalException, errmsg);
    }
    if(this->GetOverlaps())
    {
      if(pv_placement->CheckOverlaps())
      {
        // NOTE: if there are overlaps, the simulation will probably
        // crash. This is the desired behaviour (Fail Loud)
        G4cout << "Overlap, translation: "
               << pv_placement->GetFrameTranslation() << G4endl;
        G4cout << "Overlap, rotationX: "
               << pv_placement->GetFrameRotation()->colX() << G4endl;
        G4cout << "Overlap, rotationY: "
               << pv_placement->GetFrameRotation()->colY() << G4endl;
        G4cout << "Overlap, rotationZ: "
               << pv_placement->GetFrameRotation()->colZ() << G4endl;
      }
    }
    if(this->GetDrawCellVolumes())
    {
      // do not draw DNA, as we are drawing the cell volume
      pv_placement->GetLogicalVolume()->SetVisAttributes(
        G4VisAttributes::GetInvisible());
    }
    thisOctree->AddPhysicalVolume(pv_placement);
  }
  thisInfo->SetOctree(thisOctree);
  return std::make_pair(physicsLogical, thisInfo);
}

void DNAGeometry::FillVoxelVectors()
{
  G4String filename = fFractalCurveFile;
  if(!utility::Path_exists(filename))
  {
    G4ExceptionDescription errmsg;
    errmsg << "The file: " << filename << " could not be found." << G4endl;
    G4Exception("DNAGeometry::FillVoxelVectors", "File Not Found",
                FatalException, errmsg);
  }
  std::ifstream infile(filename, std::ifstream::in);
  G4String currentline;
  G4ThreeVector pos;
  G4ThreeVector rot;

  fVoxelTypes.clear();
  fVoxelIndices.clear();
  fVoxelPositions.clear();
  fVoxelRotations.clear();

  if(infile.is_open())
  {
    while(getline(infile, currentline))
    {
      if(currentline[0] != '#')
      {
        try
        {
          G4int pi_or_one            = fAnglesAsPi ? pi : 1;
          std::vector<G4String> line = utility::Split(currentline, ' ');
          fVoxelIndices.push_back(std::stoi(line.at(0)));
          fVoxelTypes.push_back((G4String) line.at(1));
          pos = G4ThreeVector(std::stod(line.at(2)) * fFractalScaling.getX(),
                              std::stod(line.at(3)) * fFractalScaling.getY(),
                              std::stod(line.at(4)) * fFractalScaling.getZ());
          rot = G4ThreeVector(std::stod(line.at(5)), std::stod(line.at(6)),
                              std::stod(line.at(7)));
          fVoxelPositions.push_back(pos);
          fVoxelRotations.push_back(rot * pi_or_one);
        } catch(int exception)
        {
          G4cout << "Voxel Position file read error" << G4endl;
          G4cout << "Got Line: " << currentline << G4endl;
          G4cout << "Format should be: "
                 << "IDX KIND POS_X POS_Y POS_Z EUL_PSI EUL_THETA"
                 << " EUL_PHI" << G4endl;
          G4cout << "Separator is space, lines with '#' as the first "
                 << "character are commented" << G4endl;
          std::exit(1);
        }
      }
    }
    infile.close();
  }
  else
  {
    G4cout << "Voxel Position file read error on open" << G4endl;
    G4cout << "Format should be: "
           << "IDX KIND POS_X POS_Y POS_Z EUL_PSI EUL_THETA EUL_PHI" << G4endl;
    G4cout << "Separator is space, lines with '#' as the first "
           << "character are commented" << G4endl;
  }
}

OctreeNode* DNAGeometry::GetTopOctreeNode(G4LogicalVolume* lv) const
{
  const PlacementVolumeInfo* info = GetPVInfo(lv);
  return (info == nullptr) ? nullptr : info->GetOctree();
}

const PlacementVolumeInfo* DNAGeometry::GetPVInfo(
  const G4LogicalVolume* lv) const
{
  if(fInfoMap.find(lv) != fInfoMap.end())
  {
    return fInfoMap.at(lv);
  }
  else if(fChemToPhys.find(lv) != fChemToPhys.end())
  {
    return fInfoMap.at(fChemToPhys.at(lv));
  }
  else
  {
    return nullptr;
  }
}

void DNAGeometry::PrintOctreeStats()
{
  G4cout << "Name         Vols     EndNodes  Max(vols/node)" << G4endl;
  G4cout << "----------------------------------------------" << G4endl;
  OctreeNode* currentOctree;
  for(auto& it : fInfoMap)
  {
    char buffer[80];
    currentOctree      = it.second->GetOctree();
    G4int vols         = currentOctree->GetContents().size();
    const char* lvname = it.first->GetName().c_str();
    G4int mvols        = currentOctree->GetMaxContents();
    G4int nodes        = currentOctree->GetNumberOfTerminalNodes();
    std::snprintf(buffer, 80, "%-12s %-7i  %-7i   %-9i", lvname, vols, nodes,
                  mvols);
    G4cout << buffer << G4endl;
  }
  G4cout << "-------------------------------------------" << G4endl;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

G4VPhysicalVolume* DNAGeometry::PlacePhosphate(G4LogicalVolume* physicsLogical,
                                               const molecule_t& thisMolecule,
                                               const molecule_t& sugar)
{
  // Define Vis Attributes
  auto vattrs = new G4VisAttributes(G4Color::Yellow());
  vattrs->SetForceSolid(true);

  std::stringstream ss;
  ss << thisMolecule.fname << "-" << thisMolecule.fchain_id << "-"
     << thisMolecule.fstrand_id << "-" << thisMolecule.fbase_idx;
  G4String name = ss.str();
  // Initial Rotation of the sphere
  G4RotationMatrix rot;
  rot.rotateX(thisMolecule.frotation.getX());
  rot.rotateY(thisMolecule.frotation.getY());
  rot.rotateZ(thisMolecule.frotation.getZ());

  G4VPhysicalVolume* pv_placement;
  G4RotationMatrix* rot_pointer;
  if(sugar.fname == "wall")
  {
    if(this->GetVerbosity() >= 2)
    {
      G4cout << "Wall placement" << G4endl;
    }
    G4ThreeVector abspos(std::abs(thisMolecule.fposition.getX()),
                         std::abs(thisMolecule.fposition.getY()),
                         std::abs(thisMolecule.fposition.getZ()));
    G4ThreeVector overlaps = abspos + thisMolecule.fsize - fVoxelSideLength;

    G4Ellipsoid* mol;
    // NOTE: will work to cut molecules at boundary only when the chain is
    // leaving the box.
    G4double z_cut = -1 * utility::Min(overlaps);
    if(z_cut < 0)  // all molecules fit
    {
      mol = new G4Ellipsoid(name, thisMolecule.fsize.getX(),
                            thisMolecule.fsize.getY(), thisMolecule.fsize.getZ
                                                       ());
    }
    else  // need to cut molecule to fit in box
    {
      mol = new G4Ellipsoid(name, thisMolecule.fsize.getX(),
                            thisMolecule.fsize.getY(), thisMolecule.fsize
                                                           .getZ(),
                            -thisMolecule.fsize.getZ(), z_cut);
    }

    rot_pointer   = new G4RotationMatrix(rot.inverse());
    auto molLogic = new G4LogicalVolume(mol, fpPhosphateMaterial, name);

    molLogic->SetVisAttributes(vattrs);
    pv_placement =
      new G4PVPlacement(rot_pointer, thisMolecule.fposition, molLogic, name,
                        physicsLogical, false, 0, this->GetOverlaps());
  }
  else
  {
    // Case where there is a sugar molecule that we would intersect with
    G4ThreeVector z_direction = sugar.fposition - thisMolecule.fposition;
    G4double separation       = z_direction.mag();
    G4double z_cut            = separation - sugar.fsize.getX();

    G4ThreeVector local_z = G4ThreeVector(0, 0, 1);
    G4ThreeVector z_new   = z_direction / z_direction.mag();

    G4ThreeVector ortho = local_z.cross(z_new);
    G4double dp         = local_z.dot(z_new);
    G4double angle      = -std::atan2(ortho.mag(), dp);

    G4RotationMatrix second_rot(ortho, angle);
    rot_pointer = new G4RotationMatrix(second_rot);

    auto mol = new G4Ellipsoid(
      name, thisMolecule.fsize.getX(), thisMolecule.fsize.getY(),
      thisMolecule.fsize.getZ(), -thisMolecule.fsize.getZ(), z_cut);

    auto molLogic = new G4LogicalVolume(mol, fpPhosphateMaterial, name);

    molLogic->SetVisAttributes(vattrs);
    pv_placement =
      new G4PVPlacement(rot_pointer, thisMolecule.fposition, molLogic, name,
                        physicsLogical, false, 0, this->GetOverlaps());
  }

  return pv_placement;
}

G4VPhysicalVolume* DNAGeometry::PlaceSugar(G4LogicalVolume* physicsLogical,
                                           const molecule_t& thisMolecule,
                                           const molecule_t& phosphate)
{
  // Define Vis Attributes
  auto vattrs = new G4VisAttributes(G4Color::Red());
  vattrs->SetForceSolid(true);

  // Name:
  std::stringstream ss;
  ss << thisMolecule.fname << "-" << thisMolecule.fchain_id << "-"
     << thisMolecule.fstrand_id << "-" << thisMolecule.fbase_idx;
  G4String name = ss.str();
  // Initial Rotation of the sphere
  G4RotationMatrix rot;
  rot.rotateX(thisMolecule.frotation.getX());
  rot.rotateY(thisMolecule.frotation.getY());
  rot.rotateZ(thisMolecule.frotation.getZ());

  G4VPhysicalVolume* pv_placement;
  G4RotationMatrix* rot_pointer;
  if(phosphate.fname == "wall")
  {
    if(this->GetVerbosity() >= 2)
    {
      G4cout << "Wall placement" << G4endl;
    }
    G4ThreeVector abspos(std::abs(thisMolecule.fposition.getX()),
                         std::abs(thisMolecule.fposition.getY()),
                         std::abs(thisMolecule.fposition.getZ()));
    G4ThreeVector overlaps = abspos + thisMolecule.fsize - fVoxelSideLength;

    G4Ellipsoid* mol;
    // NOTE: will work to cut molecules at boundary only when the chain is
    // leaving the box.
    G4double z_cut = -1 * utility::Min(overlaps);
    if(z_cut < 0)  // all molecules fit
    {
      mol = new G4Ellipsoid(name, thisMolecule.fsize.getX(),
                            thisMolecule.fsize.getY(), thisMolecule.fsize.getZ
                                                       ());
    }
    else  // need to cut molecule to fit in box
    {
      mol = new G4Ellipsoid(name, thisMolecule.fsize.getX(),
                            thisMolecule.fsize.getY(), thisMolecule.fsize
                                                           .getZ(),
                            -thisMolecule.fsize.getZ(), z_cut);
    }

    rot_pointer   = new G4RotationMatrix(rot.inverse());
    auto molLogic = new G4LogicalVolume(mol, fpSugarMaterial, name);

    molLogic->SetVisAttributes(vattrs);
    pv_placement =
      new G4PVPlacement(rot_pointer, thisMolecule.fposition, molLogic, name,
                        physicsLogical, false, 0, this->GetOverlaps());
  }
  else
  {
    // Get information about the phsophate
    G4ThreeVector z_direction = phosphate.fposition - thisMolecule.fposition;
    G4double separation       = z_direction.mag();
    G4double z_cut            = separation - phosphate.fsize.getX();

    G4ThreeVector local_z = G4ThreeVector(0, 0, 1);
    G4ThreeVector z_new   = z_direction / z_direction.mag();

    G4ThreeVector ortho = local_z.cross(z_new);
    G4double dp         = local_z.dot(z_new);
    G4double angle      = -std::atan2(ortho.mag(), dp);

    G4RotationMatrix first_rot(ortho, angle);
    rot_pointer = new G4RotationMatrix(first_rot);

    auto mol = new G4Ellipsoid(
      name, thisMolecule.fsize.getX(), thisMolecule.fsize.getY(),
      thisMolecule.fsize.getZ(), -thisMolecule.fsize.getZ(), z_cut);

    auto molLogic = new G4LogicalVolume(mol, fpSugarMaterial, name);

    molLogic->SetVisAttributes(vattrs);

    pv_placement =
      new G4PVPlacement(rot_pointer, thisMolecule.fposition, molLogic, name,
                        physicsLogical, false, 0, this->GetOverlaps());
  }
  return pv_placement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DNAGeometry::PlaceBase(G4LogicalVolume* physicsLogical,
                                          const molecule_t& thisMolecule,
                                          const molecule_t& oppBasePair,
                                          const molecule_t& nextBasePair,
                                          const molecule_t& sugar)
{
  std::stringstream ss;
  ss << thisMolecule.fname << "-" << thisMolecule.fchain_id << "-"
     << thisMolecule.fstrand_id << "-" << thisMolecule.fbase_idx;
  G4String name = ss.str();
  // Initial Rotation of the sphere
  G4RotationMatrix rot;
  rot.rotateX(thisMolecule.frotation.getX());
  rot.rotateY(thisMolecule.frotation.getY());
  rot.rotateZ(thisMolecule.frotation.getZ());

  G4ThreeVector sugar_dirn = sugar.fposition - thisMolecule.fposition;
  G4ThreeVector oppBPdirn  = oppBasePair.fposition - thisMolecule.fposition;
  G4ThreeVector nextBPdirn = nextBasePair.fposition - thisMolecule.fposition;
  // Do a cut for the sugar, do a scaling for the base pair
  G4double sugar_cut = sugar_dirn.mag() - sugar.fsize.getX();

  G4double base_scaling = oppBPdirn.mag() - oppBasePair.fsize.getY();
  base_scaling          = std::min(1., base_scaling / thisMolecule.fsize.getY
                                             ());

  // Find the new constrained values
  G4double thisShapeX = base_scaling * thisMolecule.fsize.getX();
  G4double thisShapeY = base_scaling * thisMolecule.fsize.getY();
  G4double thisShapeZ = base_scaling * thisMolecule.fsize.getZ();
  thisShapeZ          = std::min(thisShapeZ, 1.7 * angstrom);

  // first rotation to turn the ellipse to run along sugar_dirn
  G4ThreeVector z_old = G4ThreeVector(0, 0, 1);
  G4ThreeVector z_new = sugar_dirn / sugar_dirn.mag();

  G4ThreeVector ortho = z_old.cross(z_new);
  G4double angle      = std::acos(z_new.dot(z_old));

  // Now find out quadrant of angle
  G4RotationMatrix first_rot;
  G4RotationMatrix first_rot_1(ortho, angle);
  G4RotationMatrix first_rot_2(ortho, angle + halfpi);
  G4RotationMatrix first_rot_3(ortho, angle + pi);
  G4RotationMatrix first_rot_4(ortho, angle + 3 * halfpi);

  if(first_rot_1(z_old).dot(z_new) > 0.99)
  {
    first_rot = first_rot_1;
  }
  else if(first_rot_2(z_old).dot(z_new) > 0.99)
  {
    first_rot = first_rot_2;
  }
  else if(first_rot_3(z_old).dot(z_new) > 0.99)
  {
    first_rot = first_rot_3;
  }
  else if(first_rot_4(z_old).dot(z_new) > 0.99)
  {
    first_rot = first_rot_4;
  }
  else
  {
    G4cout << "Could not identify first rotation" << G4endl;
  }

  // second rotation
  // Get y vector aligned with next_base_pair_direction through a rotation
  // around sugar_dirn
  G4ThreeVector y  = first_rot(G4ThreeVector(0, 1, 0));
  G4ThreeVector up = rot(G4ThreeVector(0, 0, 1));

  G4double interval = 0.1;  // precision for finding theta
  G4double theta    = -2 * interval;
  G4double prev1    = -DBL_MAX;
  G4double prev2    = -DBL_MAX;
  G4double val      = -DBL_MAX;
  G4RotationMatrix z_rotation(z_new, theta);
  while(theta <= twopi + interval)
  {
    prev2 = prev1;
    prev1 = val;
    theta += interval;
    z_rotation.rotate(interval, z_new);
    val = z_rotation(y).dot(up);
    if((prev1 < val) && (prev1 < prev2))
    {
      // prev1 contains a local minimum
      theta -= interval;
      break;
    }
  }
  G4RotationMatrix second_rot(first_rot.inverse()(z_new), theta);

  auto rot_pointer =
    new G4RotationMatrix(second_rot.inverse() * first_rot.inverse());
  auto mol = new G4Ellipsoid(name, thisShapeX, thisShapeZ, thisShapeY,
                             -thisShapeY, sugar_cut);

  G4LogicalVolume* molLogic = nullptr;
  G4Color color;
  if(thisMolecule.fname == "Guanine")
  {
    color    = G4Color::Green();
    molLogic = new G4LogicalVolume(mol, fpGuanineMaterial, name);
  }
  else if(thisMolecule.fname == "Adenine")
  {
    color    = G4Color::Blue();
    molLogic = new G4LogicalVolume(mol, fpAdenineMaterial, name);
  }
  else if(thisMolecule.fname == "Cytosine")
  {
    color    = G4Color::Cyan();
    molLogic = new G4LogicalVolume(mol, fpCytosineMaterial, name);
  }
  else if(thisMolecule.fname == "Thymine")
  {
    color    = G4Color::Magenta();
    molLogic = new G4LogicalVolume(mol, fpThymineMaterial, name);
  }
  else
  {
    G4ExceptionDescription errmsg;
    errmsg << "Molecule name: " << thisMolecule.fname << " is invalid." <<
        G4endl
           << "Aborting" << G4endl;
    G4Exception("DNAGeometry::PlaceBase", "Invalid molecule name",
                FatalException, errmsg);
  }

  // Define Vis Attributes
  auto vattrs = new G4VisAttributes(color);
  vattrs->SetForceSolid(true);
  molLogic->SetVisAttributes(vattrs);

  G4VPhysicalVolume* pv_placement =
    new G4PVPlacement(rot_pointer, thisMolecule.fposition, molLogic, name,
                      physicsLogical, false, 0, this->GetOverlaps());

  return pv_placement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DNAGeometry::PlaceHistone(G4LogicalVolume* physicsLogical,
                                             const molecule_t& thisMolecule)
{
  G4VPhysicalVolume* pv_placement;

  // Define Vis Attributes
  auto vattrs = new G4VisAttributes(G4Color::Gray());
  vattrs->SetColor(G4Colour(0, 0, 1, 0.2));
  vattrs->SetForceSolid(true);

  // Name:
  std::stringstream ss;
  ss << thisMolecule.fname << "-" << thisMolecule.fchain_id << "-"
     << thisMolecule.fstrand_id << "-" << thisMolecule.fbase_idx;
  G4String name = ss.str();

  //// Does the Histone fit in the voxel (assumes square voxels)
  // G4double max_coord = std::max(std::abs(thisMolecule.position.getX()),
  //                              std::abs(thisMolecule.position.getY()))
  // max_coord = std::max(max_coord, std::abs(thisMolecule.position.getZ()));
  //// if overlap
  // G4double radius = utility::max(thisMolecule.size);
  // G4double radiusx=0,radiusx=0,radiusx=0;
  // if ((max_coord + radius) > fVoxelSideLength)
  //{
  //    radius = 0.9*(fVoxelSideLength - max_coord)
  //    G4cout<<" Warning : Histone overlap with mother voxel volume "<<G4endl;
  //}

  // Initial Rotation of the sphere
  G4RotationMatrix rot;
  rot.rotateX(thisMolecule.frotation.getX());
  rot.rotateY(thisMolecule.frotation.getY());
  rot.rotateZ(thisMolecule.frotation.getZ());

  G4RotationMatrix* rot_pointer;
  auto mol =
    new G4Ellipsoid(name, thisMolecule.fsize.getX(), thisMolecule.fsize.getY(),
                    thisMolecule.fsize.getZ());

  fHistoneInteractionRadius = thisMolecule.fsize.mag();

  rot_pointer = new G4RotationMatrix(rot.inverse());

  auto molLogic = new G4LogicalVolume(mol, fpHistoneMaterial, name);
  molLogic->SetVisAttributes(vattrs);

  pv_placement =
    new G4PVPlacement(rot_pointer, thisMolecule.fposition, molLogic, name,
                      physicsLogical, false, 0, this->GetOverlaps());

  return pv_placement;
}
// dousatsu --------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAGeometry::FindNearbyMolecules(
  const G4LogicalVolume* lv, const G4ThreeVector& localPosition,
  std::vector<G4VPhysicalVolume*>& daughters_pv_out, double searchRange)
{
  const PlacementVolumeInfo* pvInfo = GetPVInfo(lv);
  if(pvInfo == nullptr)
  {
    return;
  }
  // never search more than the radical kill distance
  searchRange            = std::min(searchRange, fRadicalKillDistance);
  OctreeNode* octreeNode = pvInfo->GetOctree();
  octreeNode->SearchOctree(localPosition, daughters_pv_out, searchRange);

  //     G4double r = std::pow(std::pow(localPosition.getX(), 2) +
  //     std::pow(localPosition.getY(), 2), 0.5)/nm; if (r > 6)
  //     {
  //       G4cout << "candidates found in radius " << searchRange/nm << " nm : "
  //              << daughters_pv_out.size() << G4endl;
  //       G4cout << "Radius from center of point: "
  //              << r << " nm" << G4endl;
  //       G4cout << "----------------------------------------------------" <<
  //       G4endl;
  //     }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool DNAGeometry::IsInsideHistone(const G4LogicalVolume* lv,
                                    const G4ThreeVector& localPosition) const
{
  const PlacementVolumeInfo* pvInfo = GetPVInfo(lv);
  if(pvInfo == nullptr)
  {
    G4cout << " Warning : no PV for IsInsideHistone " << G4endl;
    return false;
  }

  if(!fUseHistoneScav) {
    return false;}

  OctreeNode* octreeNode = pvInfo->GetHistoneOctree();
  const std::vector<G4VPhysicalVolume*> result =
    octreeNode->SearchOctree(localPosition, fHistoneInteractionRadius);
  return !result.empty();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAGeometry::BasePairIndexTest()
{
  G4cout << "------------------------" << G4endl;
  G4cout << "Begin BasePairIndexTest" << G4endl;
  G4bool pass = true;
  for(G4int ii = 0; ii != GetNumberOfChains(); ii++)
  {
    G4cout << "Testing Chain " << ii << G4endl;
    std::vector<int64_t> test_vector;  // dousatsu
    // std::vector<long long> test_vector;//ORG
    int64_t start_idx;
    int64_t end_idx;  // dousatsu
    // long long start_idx, end_idx;//ORG
    for(G4int jj = 0; jj != (G4int) fPlacementTransformations.size(); jj++)
    {
      start_idx = GetStartIdx(jj, ii);
      end_idx   = GetEndIdx(jj, ii);

      test_vector.push_back(start_idx);
      test_vector.push_back(end_idx);
    }
    for(auto it = test_vector.begin(); it != test_vector.end(); it++)
    {
      if(it == test_vector.begin())  // exception because start of vector
      {
        if(it + 1 != test_vector.end())
        {
          if(*it > *(it + 1))
          {
            G4ExceptionDescription errmsg;
            errmsg << "BasePairIndexTest fails at start of test. "
                   << "The index " << *(it + 1) << " occurs after " << (*it)
                   << "." << G4endl;
            G4Exception("DNAGeometry::BasePairIndexTest", "ERR_TEST_FAIL",
                        FatalException, errmsg);
          }
        }
      }
      else if(it == test_vector.end() - 1)  // exception because vec end
      {
        if(*it < *(it - 1))
        {
          G4ExceptionDescription errmsg;
          errmsg << "BasePairIndexTest fails. "
                 << "The index " << *(it - 1) << " occurs before " << (*it)
                 << "." << G4endl;
          G4Exception("DNAGeometry::BasePairIndexTest", "ERR_TEST_FAIL",
                      FatalException, errmsg);
        }
      }
      else  // default case
      {
        if(*it < *(it - 1))
        {
          G4ExceptionDescription errmsg;
          errmsg << "BasePairIndexTest fails. "
                 << "The index " << *(it - 1) << " occurs before " << (*it)
                 << "." << G4endl;
          G4Exception("DNAGeometry::BasePairIndexTest", "ERR_TEST_FAIL",
                      FatalException, errmsg);
        }
        if(*it > *(it + 1))
        {
          G4ExceptionDescription errmsg;
          errmsg << "BasePairIndexTest fails. "
                 << "The index " << *(it + 1) << " occurs after " << (*it)
                 << "." << G4endl;
          G4Exception("DNAGeometry::BasePairIndexTest", "ERR_TEST_FAIL",
                      FatalException, errmsg);
        }
      }
    }
  }
  if(pass)
  {
    G4cout << "BasePairIndexTest passed" << G4endl;
    G4cout << "------------------------" << G4endl;
  }
}

/**
 * Global unique id is an integer with the following bases
 * moleculeID: base ::molecule::Count
 * PlacementIdx: base this->GetNumberOfPlacements()
 * StrandID: base 2
 * ChainID: base GetNumberOfChains()
 * BasePairIdx: base GetMaxBPIdx()
 */
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int64_t DNAGeometry::GetGlobalUniqueID(G4VPhysicalVolume* dnavol,
                               const G4VTouchable* touch) const
{
  const G4String& dnaName        = dnavol->GetName();
  std::array<G4String, 4> pv_arr = utility::Get_four_elements(dnaName, '-');
  molecule mol                   = utility::GetMoleculeEnum(pv_arr.at(0));
  G4int chainIdx                 = std::stoi(pv_arr.at(1));
  G4int strandIdx                = std::stoi(pv_arr.at(2));
  int64_t baseIdx = std::stoll(pv_arr.at(3));  // dousatsu
  // G4int baseIdx   = std::stoi(pv_arr.at(3));//ORG

  const G4String& placementName = touch->GetVolume()->GetName();
  G4int placeIdx =
    std::stoi(utility::Get_seperated_element(placementName, '-', 1));

  G4int chains     = this->GetNumberOfChains();      // dousatsu
  G4int placements = this->GetNumberOfPlacements();  // dousatsu
  G4int molecules  = ::molecule::Count;              // dousatsu
  // long long chains = this->GetNumberOfChains();//ORG
  // long long placements = this->GetNumberOfPlacements();//ORG
  // long long molecules  = ::molecule::Count;//ORG

  chainIdx = this->GetGlobalChain(placeIdx, chainIdx);
  // long long globalPair = baseIdx + this->GetStartIdx(placeIdx, chainIdx);
  if(this->GetStrandsFlipped(placeIdx)) {
    strandIdx = (strandIdx + 1) % 2;}

  int64_t val1 = mol;
  int64_t val2 = molecules * placeIdx;
  int64_t val3 = molecules * placements * strandIdx;
  int64_t val4 = molecules * placements * 2 * chainIdx;
  int64_t val5 = molecules * placements * 2 * chains * baseIdx;
  int64_t sum  = val1 + val2 + val3 + val4 + val5;
  if(sum < 0)
  {
    G4ExceptionDescription errmsg;
    errmsg << "v1: " << val1 << G4endl << "v2: " << val2 << G4endl
           << "v3: " << val3 << G4endl << "v4: " << val4 << G4endl
           << "v5: " << val5 << G4endl << "sum: " << sum << G4endl
           << "placeIdx: " << placeIdx << G4endl << "placements: " << placements
           << G4endl << "strandIdx: " << strandIdx << G4endl
           << "chainIdx: " << chainIdx << G4endl << "chains: " << chains
           << G4endl << "baseIdx: " << baseIdx << G4endl << "mol: " << mol
           << G4endl;
    G4Exception("DNAGeometry::GetGlobalUniqueID", "ERR_BAD_ID", FatalException,
                errmsg);
  }
  return sum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// copy of the GetGlobalUniqueID except with integer inputs rather than strings
int64_t
DNAGeometry::GetGlobalUniqueIDTest(int64_t mol, int64_t placeIdx, int64_t chainIdx,
                                   int64_t strandIdx, int64_t baseIdx) const
{
  int64_t chains     = this->GetNumberOfChains();      // dousatsu
  int64_t placements = this->GetNumberOfPlacements();  // dousatsu
  int64_t molecules  = ::molecule::Count;              // dousatsu
  // long long chains = this->GetNumberOfChains();//ORG
  // long long placements = this->GetNumberOfPlacements();//ORG
  // long long molecules  = ::molecule::Count;//ORG

  chainIdx = this->GetGlobalChain(placeIdx, chainIdx);
  if(this->GetStrandsFlipped(placeIdx)) {
    strandIdx = (strandIdx + 1) % 2;}

  int64_t val1 = mol;
  int64_t val2 = molecules * placeIdx;
  int64_t val3 = molecules * placements * strandIdx;
  int64_t val4 = molecules * placements * 2 * chainIdx;
  int64_t val5 = molecules * placements * 2 * chains * baseIdx;
  int64_t sum  = val1 + val2 + val3 + val4 + val5;
  if(sum < 0)
  {
    G4ExceptionDescription errmsg;
    errmsg << "v1: " << val1 << G4endl << "v2: " << val2 << G4endl
           << "v3: " << val3 << G4endl << "v4: " << val4 << G4endl
           << "v5: " << val5 << G4endl << "sum: " << sum << G4endl
           <<" molecules : "<<molecules<<G4endl
           << "placeIdx: " << placeIdx << G4endl
           << "placements: " << placements<< G4endl
           << "strandIdx: " << strandIdx << G4endl
           << "chainIdx: " << chainIdx << G4endl << "chains: " << chains
           << G4endl << "baseIdx: " << baseIdx << G4endl << "mol: " << mol
           << G4endl;
    G4Exception("DNAGeometry::GetGlobalUniqueIDTest",  // dousatsu
                "ERR_BAD_ID", FatalException, errmsg);
  }
  return sum;
}

void DNAGeometry::UniqueIDTest()
{
  G4cout << "Unique ID Test starting." << G4endl;
  G4int molecules = ::molecule::Count;

  G4int chains     = this->GetNumberOfChains();
  G4int placements = this->GetNumberOfPlacements();
  G4int basepairs =
    100000;  // Maximum number of basepairs to consider per block

  G4bool pass   = true;
  G4int counter = 0;
  for(G4int pment = 0; pment != (placements); ++pment)
  {
    for(G4int chain = 0; chain != (chains); ++chain)
    {
      for(G4int mol = 0; mol != (molecules); ++mol)
      {
        // repeat 10 times along the strand/chain/molecule
        for(int64_t ii = 0; ii != 10; ++ii)
        {
          counter++;
          G4int strand = G4UniformRand() > 0.5 ? 0 : 1;
          G4int basepair = (G4int) basepairs * G4UniformRand();  // dousatsu
          int64_t unique_id =                                    // dousatsu
                               // long long unique_id = //ORG
            this->GetGlobalUniqueIDTest(mol, pment, chain, strand, basepair);

          G4int real_strand = (this->GetStrandsFlipped(pment)) ?  // ORG
                                (strand + 1) % 2
                                                               : strand;
          G4int real_chain = this->GetGlobalChain(pment, chain);  // ORG
          int64_t real_bp =
            basepair + this->GetStartIdx(pment, real_chain);  // dousatsu
          // long long real_bp = basepair + this->GetStartIdx(pment,
          // real_chain);//ORG

          G4int obs_mol       = this->GetMoleculeFromUniqueID(unique_id);
          G4int obs_placement = this->GetPlacementIndexFromUniqueID(unique_id);
          G4int obs_strand    = this->GetStrandIDFromUniqueID(unique_id);
          G4int obs_chain     = this->GetChainIDFromUniqueID(unique_id);
          int64_t obs_basepair =
            this->GetBasePairFromUniqueID(unique_id);  // dousatsu
          // long long obs_basepair = this->GetBasePairFromUniqueID
          // (unique_id);//ORG

          if((pment != obs_placement) || (real_strand != obs_strand) ||
             (real_chain != obs_chain) || (obs_mol != mol) ||
             (real_bp != obs_basepair))
          {
            pass = false;
            G4cerr << "Unique ID test failed at ID: " << unique_id << G4endl;
            G4cerr << "real_placement: " << pment
                   << " observed: " << obs_placement << G4endl;
            G4cerr << "real_strand: " << real_strand
                   << " observed: " << obs_strand << G4endl;
            G4cerr << "real_chain: " << real_chain << " observed: " << obs_chain
                   << G4endl;
            G4cerr << "real_mol: " << mol << " observed: " << obs_mol << G4endl;
            G4cerr << "real_bp: " << real_bp << " observed: " << obs_basepair
                   << G4endl;
          }
        }
      }
    }
  }
  if(!pass)
  {
    G4Exception("DNAGeometry::UniqueIDTest", "ERR_TEST_FAIL", FatalException,
                "Algorithm to generate unique ID's failed");
  }
  G4cout << "Unique ID Test completed, after running " << counter << " tests"
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DNAGeometry::GetMaterialFromUniqueID(int64_t idx) const
{
  molecule mol = GetMoleculeFromUniqueID(idx);
  switch(mol)
  {
    case ::molecule::SUGAR:
      return fpSugarMaterial;
    case ::molecule::PHOSPHATE:
      return fpPhosphateMaterial;
    case ::molecule::GUANINE:
      return fpGuanineMaterial;
    case ::molecule::ADENINE:
      return fpAdenineMaterial;
    case ::molecule::THYMINE:
      return fpThymineMaterial;
    case ::molecule::CYTOSINE:
      return fpCytosineMaterial;
    case ::molecule::UNSPECIFIED:
      G4cout << "DNAGeometry::GetMaterialFromUniqueID: "
             << "Encountered an unspecified material." << G4endl;
      return fpMediumMaterial;
    default:
      G4cout << "DNAGeometry::GetMaterialFromUniqueID: "
             << "Encountered a bad material enumerator." << G4endl;
      return fpMediumMaterial;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DNAGeometry::GetPlacementIndexFromUniqueID(
  int64_t idx) const  // dousatsu
{
  // remove molecule
  idx -= idx % (int64_t)::molecule::Count;
  // remove everything above molecule
  idx = idx % (int64_t)(::molecule::Count * this->GetNumberOfPlacements());
  return (G4int) idx / (::molecule::Count);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DNAGeometry::GetStrandIDFromUniqueID(int64_t idx) const  // dousatsu
{
  // remove molecule
  idx -= (int64_t)(idx % ::molecule::Count);
  // remove placement
  idx -= (int64_t)(idx % (::molecule::Count * this->GetNumberOfPlacements()));
  // remove everything above strand
  idx =
    (int64_t)(idx % (::molecule::Count * this->GetNumberOfPlacements() * 2));
  // recover strand
  return (G4int) idx / (::molecule::Count * this->GetNumberOfPlacements());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DNAGeometry::GetChainIDFromUniqueID(int64_t idx) const  // dousatsu
{
  // remove molecule
  idx -= idx % (int64_t)::molecule::Count;
  // remove placement
  idx -= idx % (int64_t)(::molecule::Count * this->GetNumberOfPlacements());
  // remove strand
  idx -= idx % (int64_t)(::molecule::Count * this->GetNumberOfPlacements() * 2);
  // remove everything above chain
  idx = idx % (int64_t)(::molecule::Count * this->GetNumberOfPlacements() * 2 *
                        this->GetNumberOfChains());
  // recover chain
  G4int chain = (G4int) idx / (int64_t)(::molecule::Count *
                                        this->GetNumberOfPlacements() * 2);
  return chain;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int64_t DNAGeometry::GetBasePairFromUniqueID(int64_t idx) const
{
  // remove molecule
  idx -= idx % ::molecule::Count;

  // remove placement
  int64_t placement_ =
    idx % (int64_t)(::molecule::Count * this->GetNumberOfPlacements());
  idx -= placement_;

  // find placement for later
  placement_ = placement_ / (int64_t)(::molecule::Count);

  // remove strand
  idx -= idx % (int64_t)(::molecule::Count * this->GetNumberOfPlacements() * 2);

  // remove chain
  int64_t chain =
    idx % (int64_t)(::molecule::Count * this->GetNumberOfPlacements() * 2 *
                    this->GetNumberOfChains());
  idx -= chain;

  // find chain for later
  chain =
    chain / (int64_t)(::molecule::Count * this->GetNumberOfPlacements() * 2);

  // only base pairs left
  int64_t bp =
    idx / (int64_t)(::molecule::Count * this->GetNumberOfPlacements() * 2 *
                    this->GetNumberOfChains());

  // use chain and placement to get BP offset
  bp = bp + this->GetStartIdx(placement_, chain);
  return bp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int64_t DNAGeometry::GetMaxBPIdx() const  // dousatsu
{
  int64_t mx = 0;  // dousatsu
  for(auto val : std::get<2>(*fPlacementTransformations.end()))
  {
    if(val > mx) {
      mx = val;}
  }
  return mx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DNAGeometry::GetLocalChain(G4int vol_idx, G4int global_chain) const
{
  auto chains = std::get<0>(fPlacementTransformations[vol_idx]);
  G4int local = 0;
  for(auto chain : chains)
  {
    if(global_chain == chain) {
      return local;}
    local++;
  }
  return -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DNAGeometry::GetMatchingPhysVolume(
  G4LogicalVolume* chem) const
{
  auto it = fChemToPhys.find(chem);
  if(it != fChemToPhys.end())
  {
    return it->second;
  }
  else if(fPhysToChem.find(chem) != fPhysToChem.end())
  {
    G4cout << "A phyics volume was used to find a chem volume" << G4endl;
    return chem;
  }
  else
  {
    return nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DNAGeometry::GetMatchingChemVolume(
  G4LogicalVolume* physics) const
{
  auto it = fPhysToChem.find(physics);
  if(it != fPhysToChem.end())
  {
    return it->second;
  }
  else if(fChemToPhys.find(physics) != fChemToPhys.end())
  {
    G4cout << "A chem volume was used to find a physics volume" << G4endl;
    return physics;
  }
  else
  {
    return nullptr;
  }
}

void DNAGeometry::ChromosomeTest() { fpChromosomeMapper->Test(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DNAGeometry::AddChangeMoleculeSize(const G4String& name, const
                                         G4ThreeVector& size)
{
  G4String lower_name = name;
  // lower_name.toLower();
  G4StrUtil::to_lower(lower_name);
  if(!(fEnableCustomMoleculeSizes))
  {
    G4cerr << "Custom Molecule Sizes need to be enabled first" << G4endl;
  }
  else
  {
    if(fMoleculeSizes.count(lower_name) == 0)
    {
      // add molecule
      fMoleculeSizes[lower_name] = size;
    }
    else
    {
      // error message for replacement
      G4cout << lower_name
             << " is already defined in molecule size, I will update"
             << G4endl << " this molecule's size from"
             << fMoleculeSizes[lower_name] << G4endl << " to " << size
             << G4endl;
      fMoleculeSizes[lower_name] = size;
    }
  }
}