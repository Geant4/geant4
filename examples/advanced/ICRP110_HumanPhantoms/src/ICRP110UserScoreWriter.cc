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
// Code developed by:
// S.Guatelli, M. Large and A. Malaroda, University of Wollongong
//
//Original code from geant4/examples/extended/runAndEvent/RE03
//
#include <vector>
#include <map>
#include "ICRP110UserScoreWriter.hh"
#include "ICRP110ScoreWriterMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDParticleFilter.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VScoringMesh.hh"

ICRP110UserScoreWriter::ICRP110UserScoreWriter():
G4VScoreWriter() 
{
 fMessenger = new ICRP110ScoreWriterMessenger(this);
 fSex = "female"; //Default phantom sex is female
 fSection = "head"; // Default phantom section is head
}

ICRP110UserScoreWriter::~ICRP110UserScoreWriter() 
{
  delete fMessenger;
}

void ICRP110UserScoreWriter::DumpQuantityToFile(const G4String & psName, const G4String & fileName, const G4String & option) 
{
    using MeshScoreMap = G4VScoringMesh::MeshScoreMap;

    if(verboseLevel > 0) 
      {
      G4cout << "ICRP110UserScorer-defined DumpQuantityToFile() method is invoked." << G4endl; 
      }

    // change the option string into lowercase to the case-insensitive.
    G4String opt = option;
    std::transform(opt.begin(), opt.end(), opt.begin(), (int (*)(int))(tolower));
    
    // confirm the option
    if(opt.size() == 0) opt = "csv";

//--------------------------------------------------------------------//
//----------------Create Scoring Mesh Output Text File----------------//
//--------------------------------------------------------------------//
// First we create use the scoring mesh to create a default output text
// file containing 4 columns: voxel number along x, y, z, and edep deposited
// in that voxel (in J). This file is to be called "PhantomMesh_Edep.txt". 

std::ofstream ofile(fileName);
  
if(!ofile) 
{
   G4cerr << "ERROR : DumpToFile : File open error -> " << fileName << G4endl;
   return;
}
  ofile << "# mesh name: " << fScoringMesh -> GetWorldName() << G4endl;

// retrieve the map
MeshScoreMap fSMap = fScoringMesh -> GetScoreMap();
  
auto msMapItr = fSMap.find(psName);
  
if(msMapItr == fSMap.end()) 
  {
   G4cerr << "ERROR : DumpToFile : Unknown quantity, \""<< psName 
   << "\"." << G4endl;
   return;
  }

std::map<G4int, G4StatDouble*> * score = msMapItr -> second-> GetMap();
  
ofile << "# primitive scorer name: " << msMapItr -> first << G4endl;

  // declare dose array and initialize to zero.
  std::vector<double> ScoringMeshEdep;
  for(G4int y = 0; y < fNMeshSegments[0]*fNMeshSegments[1]*fNMeshSegments[2]; y++) ScoringMeshEdep.push_back(0.);

ofile << std::setprecision(16); // for double value with 8 bytes
  
for(G4int x = 0; x < fNMeshSegments[0]; x++) {
   for(G4int y = 0; y < fNMeshSegments[1]; y++) {
     for(G4int z = 0; z < fNMeshSegments[2]; z++){
        // Retrieve dose in each scoring mesh bin/voxel
        G4int idx = GetIndex(x, y, z);
        std::map<G4int, G4StatDouble*>::iterator value = score -> find(idx);
        if (value != score -> end()) ScoringMeshEdep[idx] += (value->second->sum_wx())/(joule);
       }
      }
     }

ofile << std::setprecision(6);

ofile << std::setprecision(16); // for double value with 8 bytes
  
for(G4int x = 0; x < fNMeshSegments[0]; x++) {
   for(G4int y = 0; y < fNMeshSegments[1]; y++) {
     for(G4int z = 0; z < fNMeshSegments[2]; z++){
         
         G4int idx = GetIndex(x, y, z);
         
         if (ScoringMeshEdep[idx] != 0){
         ofile << x << '\t' << y << '\t' << z << '\t' << ScoringMeshEdep[idx] << G4endl;
         }
         //Store x,y,z and dose for each voxel in output text file. 
     
       }
      }
     }

// Close the output ASCII file
ofile.close();

//----------------------------------------------------------------------------------------//
//-----Read Data.dat File to determine the name and number of slice files to open---------//
//----------------------------------------------------------------------------------------//
// Using the macro commands from the .in files, the UserScoreWriter identifies which
// phantom sex and section has been constructed. We then store the names of the individual z-slices 
// which have been called upon in the detector construction when creating the phantom.

G4int NSlices = 0;
G4int NXVoxels = 0;
G4int NYVoxels = 0;
std::ifstream DataFile;

  G4cout << "Phantom Sex: " << fSex << G4endl;
  G4cout << "Phantom Section: " << fSection << G4endl;

G4String male = "male";
G4String female = "female";
//G4String 

  //Determine Phantom Sex and Section which was Simulated
  if (fSex == "male")
  {
    if (fSection == "head")
    {
      DataFile.open("ICRPdata/MaleHead.dat");
      G4cout << "Selecting data file ICRPdata/MaleHead.dat..." << G4endl;
    }
    if (fSection == "trunk")
    {
      DataFile.open("ICRPdata/MaleTrunk.dat");
      G4cout << "Selecting data file ICRPdata/MaleTrunk.dat..." << G4endl;
    }
    if (fSection == "full")
    {
      DataFile.open("ICRPdata/MaleData.dat");
      G4cout << "Selecting data file ICRPdata/MaleData.dat..." << G4endl;
    }
  }
  else if(fSex == "female")
  {
    if (fSection == "head")
    {
      DataFile.open("ICRPdata/FemaleHead.dat");
      G4cout << "Selecting data file ICRPdata/FemaleHead.dat..." << G4endl;
    }
    if (fSection == "trunk")
    {
      DataFile.open("ICRPdata/FemaleTrunk.dat");
      G4cout << "Selecting data file ICRPdata/FemaleTrunk.dat..." << G4endl;
    }
    if (fSection == "full")
    {
      DataFile.open("ICRPdata/FemaleData.dat");
      G4cout << "Selecting data file ICRPdata/FemaleData.dat..." << G4endl;
    }
  }
  else
  {
    G4cout << "Phantom Sex or section not correctly specified to ICRP110UserScoreWriter" << G4endl;
  }

		//Check if file opens
		if(DataFile.good() != 1 )
		{
			G4cout << "Problem Reading Data File" << G4endl;
		}
		else 
		{
			G4cout << "Opening Data.dat File..." << G4endl; 
		}
   
DataFile >> NSlices;
G4cout << "Number of Phantom Slices Simulated = " << NSlices << G4endl;

DataFile >> NXVoxels >> NYVoxels;
G4cout << "Number of X Voxels per slice = " << NXVoxels << G4endl;
G4cout << "Number of Y Voxels per slice = " << NYVoxels << G4endl;


G4int VoxelsPerSlice = 0;
VoxelsPerSlice = NXVoxels * NYVoxels;

//Skip lines 4-60 of Data.dat as they hold no useful information for this code
	for (G4int i = 0; i < 58; i++){
		DataFile.ignore(256, '\n'); 
	}

// Read file names to open from Data.dat (i.e. those that were used in simulation)
 
std::vector<G4String> SliceName; //char SliceName[NSlices][20]; //Stores name and number of phantom slices used in simulation

  for(G4int i = 0; i < NSlices; i++){
		SliceName.push_back("empty");
	}

	for (G4int i = 0; i < NSlices; i++){
		DataFile >> SliceName[i];
	}


//--------------------------------------------------------------------//
//----------------------Read Phantom Slice Files----------------------//
//--------------------------------------------------------------------//
// Reads each of the phantom z-slice files identified by the above code
// and stores the organIDs within each voxel sequentially. 
// Later, we will compare the dose in each voxel to the organ ID of each 
// voxel to calculate total dose in each organ. 

G4int ARRAY_SIZE = VoxelsPerSlice * NSlices;

std::vector<G4int> OrganIDs; 

std::ifstream PhantomFile;
 
for (G4int ii = 0; ii < NSlices; ii++){

G4String sliceVar = SliceName[ii];
G4String slice;

  if (fSex == "male")
  {
    slice = "ICRPdata/ICRP110_g4dat/AM/"+sliceVar;
  }
  else if(fSex == "female")
  {
    slice = "ICRPdata/ICRP110_g4dat/AF/"+sliceVar;
  }

  PhantomFile.open(slice.c_str());

		//Check if file opens
		if(PhantomFile.good() != 1 )
		{
			G4cout << "Problem Reading Phantom Slice File:" << SliceName[ii] << G4endl;
		}
		else 
		{
			G4cout << "Opening Phantom Slice File: " << SliceName[ii] << G4endl; 
		}

	G4int FNVoxelX;
	G4int FNVoxelY;
	G4int FNVoxelZ;

		//Input first 3 numbers of file - They are not organ IDs	
		PhantomFile >> FNVoxelX >> FNVoxelY >> FNVoxelZ;
  
  G4int aa = 0;
  
 for (G4int i = 0; i < VoxelsPerSlice; i++)
	{
		PhantomFile >> aa;
    OrganIDs.push_back(aa);
	}

  PhantomFile.close();
}

//---------------------------------------------------------------//
//--------------Read Data from Scoring Mesh Text File------------//
//---------------------------------------------------------------//
// Opens and reads initial/default scoring mesh output file which
// was created at the beginning of the code. We now store the edep in each
// voxel to cross reference against the organ ID of each voxel for 
// calculations of total dose in each organ.

std::ifstream MeshFile(fileName);

	//Check if file opens
	if(MeshFile.good() != 1 )
	{
		G4cout << "Problem Reading Data File: " << fileName << G4endl;
	}
	else {
		G4cout << "Opening File: " << fileName << G4endl; 
	}


//-----Reads Phantom Mesh Text File and Stores Data in 6 different Vectors-------//

G4int col = 4;
G4int lines = VoxelsPerSlice*NSlices;

//Ignore first 2 lines of PhantomMesh.txt as they are text headers
MeshFile.ignore(256, '\n');
MeshFile.ignore(256, '\n');

std::vector<G4int> X_MeshID; //Stores X-position of all scoring mesh voxels
std::vector<G4int> Y_MeshID; //Stores Y-position
std::vector<G4int> Z_MeshID; //Stores Z-position
std::vector<G4double> Edep; //Stores edep in voxels

G4int nX = 0; //Number along X of scoring mesh voxel
G4int nY = 0; //Number along Y
G4int nZ = 0; //Number along Z
G4double EDep = 0.0; //edep deposited in individual voxels

for (G4int i=0; i< lines; i++){
	for (G4int j=0; j < col;){
 
		MeshFile >> nX; //Reads number along X of current scoring mesh voxel
		X_MeshID.push_back(nX); //Stores it sequentially in vector X_MeshID
      j++;
     
		MeshFile >> nY; // Reads number along Y
		Y_MeshID.push_back(nY); // Stores in vector
      j++;
  
	  MeshFile >> nZ; // Reads number along Z
		Z_MeshID.push_back(nZ); // Stores in vector
      j++;
      
		MeshFile >> EDep; // Reads edep in each voxel
		Edep.push_back(EDep); // Stores in vector
      j++;
 
 }
}

MeshFile.close();

//--------------------------------------------------------------------//
//---------Reads AM/AF_organs.dat file and stores info about----------//
//------------------------the phantom organs--------------------------//
//--------------------------------------------------------------------//

std::ifstream PhantomOrganNames;

  if (strcmp(fSex.c_str(), male.c_str()) == 0)
  {
      PhantomOrganNames.open ("ICRPdata/ICRP110_g4dat/P110_data_V1.2/AM/AM_organs.dat");
    	//Check if file opens
      	if(PhantomOrganNames.good() != 1 )
      	{
      		G4cout << "Problem reading AM_organs.dat" << G4endl;
      	}
      	else
        {
      		G4cout << "Reading AM_organs.dat" << G4endl; 
      	}
  }
  else if(strcmp(fSex.c_str(), female.c_str()) == 0)
  {
      PhantomOrganNames.open ("ICRPdata/ICRP110_g4dat/P110_data_V1.2/AF/AF_organs.dat");
    	//Check if file opens
      	if(PhantomOrganNames.good() != 1 )
      	{
      		G4cout << "Problem reading AF_organs.dat" << G4endl;
      	}
      	else
        {
      		G4cout << "Reading AF_organs.dat" << G4endl; 
      	}
  }

	
PhantomOrganNames.ignore(256, '\n');
PhantomOrganNames.ignore(256, '\n');
PhantomOrganNames.ignore(256, '\n');
PhantomOrganNames.ignore(256, '\n');

G4String str;
std::vector<G4String> OrganNames;

  OrganNames.push_back("0     Air"); // Register air surrounding phantom
    //Fill with organ IDs of the phantom from ICRP data files (located in /ICRPdata/ICRP110_g4dat/P110_data_V1.2/)
    while(getline(PhantomOrganNames, str))
    {
      OrganNames.push_back(str);
    }
  OrganNames.push_back("141    Phantom Top/Bottom Skin Layer"); //Registers top and bottom slices of phantom made entirely of
  // skin. The skin in these layers has organ ID 141 to differentiate it from other skin, and is given its
  // own organ ID so that the user can choose whether to include it or not. 
  
//-------------------------------------------------------------------//
//---------Reads OrganMasses.dat file and stores info about----------//
//--------------------the phantom organ massess----------------------//
//-------------------------------------------------------------------//

G4int NOrganIDs = OrganNames.size();

std::ifstream OrganMasses;

      OrganMasses.open ("ICRPdata/OrganMasses.dat");
    	//Check if file opens
      	if(OrganMasses.good() != 1 )
      	{
      		G4cout << "Problem reading OrganMasses.dat" << G4endl;
      	}
      	else
        {
      		G4cout << "Reading OrganMasses.dat" << G4endl; 
      	}


OrganMasses.ignore(256, '\n'); //Igonore first line as it is a header

std::vector<G4int> iteratorID;
std::vector<G4double> MaleOrganMasses;
std::vector<G4double> FemaleOrganMasses;

G4int itID = 0;
G4double massM = 0.0;
G4double massF = 0.0;  
  
  for (G4int i = 0; i < NOrganIDs; i++)
    {
        OrganMasses >> itID;
        iteratorID.push_back(itID);
        
        OrganMasses >> massM;
        MaleOrganMasses.push_back(massM);
        
        OrganMasses >> massF;
        FemaleOrganMasses.push_back(massF);
    }

//---------------------------------------------------------------------------//
//----------------------Writes Outputs of code to File-----------------------//
//-------------------------------OrganDeps.out------------------------------//
//---------------------------------------------------------------------------//
// As the final step, we compare the edep in each voxel with the voxels organID 
// and sum the edep in voxels with identical organIDs to obtain total edep in 
// each organ. We then divide total edep in each organ by their respective organ
// mass to give total dose received in each organ (in Gy). 
// All this information is then output to the file "ICRP110.out".

std::ofstream OutputFile2;

G4int VoxelNumber = 0; 
G4int OrganIndex = 0;
G4cout << "NOrganIDs: " << NOrganIDs << G4endl;
std::vector <G4double> OrganDep;
std::vector <G4double> OrganDose;

G4double a = 0.0;
G4double b = 0.0;

for (G4int i = 0; i < NOrganIDs; i++)
{
  OrganDep.push_back(a);
  OrganDose.push_back(b);
}


for (G4int i = 0; i < ARRAY_SIZE; i++){
	VoxelNumber = X_MeshID[i] + NXVoxels * Y_MeshID[i] + VoxelsPerSlice * Z_MeshID[i]; 
  OrganIndex = OrganIDs[VoxelNumber];
  OrganDep[OrganIndex] += Edep[i];
}

//Calculate dose in each organ by dividing edep in each organ by organ masses (in kg)
  if (strcmp(fSex.c_str(), male.c_str()) == 0)
  {
    for (G4int i = 0; i < NOrganIDs; i++)
      {
        OrganDose[i] = (MaleOrganMasses[i] == 0 ) ? 0 : OrganDep[i]/(MaleOrganMasses[i] * 1e-3); 
      }
  }
  else if(strcmp(fSex.c_str(), female.c_str()) == 0)
  {
    for (G4int i = 0; i < NOrganIDs; i++)
      {
        OrganDose[i] = (FemaleOrganMasses[i] == 0 ) ? 0 : OrganDep[i]/(FemaleOrganMasses[i] * 1e-3); 
      }  
  }

OutputFile2.open ("ICRP110.out");

	//Check if file opens
	if(OutputFile2.good() != 1 )
	{
		G4cout << "Problem writing output to ICRP110.out" << G4endl;
	}
	else {
		G4cout << "Writing output to ICRP110.out" << G4endl; 
	}


G4double TotalDep = 0.0;
G4double TotalDose = 0.0;

OutputFile2 << G4endl; 
OutputFile2 << '\t' << "-------------------------------- " << G4endl;
OutputFile2 << '\t' << "OrganID" << '\t' << "Edep (J)" << '\t' << "Dose (Gy)" << G4endl;
OutputFile2 << '\t' << "-------------------------------- " << G4endl;


for (G4int i = 1; i < NOrganIDs; i++)
{
  if (OrganDep[i] != 0)
  {
    if (i != 140) //Skip dose deposited in air inside body
    {
      OutputFile2 << '\t' << i << " |" << '\t' << '\t' << OrganDep[i] << '\t' << OrganDose[i] << G4endl;
    }
  }
}

OutputFile2 << "----------------------------------------------------------------------------" << G4endl;
OutputFile2 << "-------------------------------ORGAN INFO-----------------------------------" << G4endl;
OutputFile2 << "-----------------(of organs where edep/dose was recorded)-------------------" << G4endl;
OutputFile2 << "----------------------------------------------------------------------------" << G4endl;
OutputFile2 << "ID" << '\t' << '\t' << "Organ Name" << '\t' << "    " << '\t' << "    " << '\t' << "    " << '\t' << "Material ID" << '\t'  << '\t' << "Density (g/cm^3) " << G4endl;

for(G4int i = 1; i < NOrganIDs; i++)
{
  if (OrganDep[i] != 0)
  {
    if (i != 140) //Skip dose deposited in air inside body
    {
      OutputFile2 << OrganNames[i] << G4endl;
    }
  }
}

//Sum total dose over all organs
for (G4int i = 1; i < NOrganIDs; i++)
{
    if (i != 140) //Skip dose deposited in air inside body
    {
      TotalDep += OrganDep[i];
      TotalDose += OrganDose[i];
    }
}

OutputFile2 << G4endl;
OutputFile2 << "Total Edep over all organs = " << TotalDep << " J" << G4endl;
OutputFile2 << "Total dose absorbed over all organs = " << TotalDose << " Gy" << G4endl;


OutputFile2 << G4endl;
OutputFile2 << "----------------------------------------------------------------------------" << G4endl;
OutputFile2 << "----------------ORGAN ENERGY DEPOSITIONS AND ABSORBED DOSE------------------" << G4endl;
OutputFile2 << "-----------(for all organs [includes air - OrganIDs = 0, 140])--------------" << G4endl;
OutputFile2 << "--------------([and top/bottom skin layer - OrganIDs = 141])----------------" << G4endl;
OutputFile2 << "----------------------------------------------------------------------------" << G4endl;
OutputFile2 << "OrganID" << '\t' << "Edep (J) " << '\t' << "Dose (Gy) " << G4endl;
OutputFile2 << "-------------------------------" << G4endl;

for (G4int i = 0; i < NOrganIDs; i++)
{
    OutputFile2 << i << "  | " << '\t' << '\t' << OrganDep[i] << '\t' << OrganDose[i] << G4endl;
}

OutputFile2 << "Total energy depositied over all organs = " << TotalDep << " J" << G4endl;
OutputFile2 << "Total absorbed dose over all organs = " << TotalDose << " Gy " << G4endl;
G4cout << "Total energy deposited over all Organs within the Phantom is " << TotalDep << " J" << G4endl;
G4cout << "Total absorbed dose over all phantom organs is " << TotalDose << " Gy " << G4endl;

OutputFile2.close();
}


// Sets the sex of the phantom as defined through the messenger class
void ICRP110UserScoreWriter::SetPhantomSex(G4String newSex)
{
  fSex = newSex;
  
  if (fSex == "male")
    {
      G4cout << ">> Male Phantom identified by UserScoreWriter." << G4endl;
    }
  if (fSex == "female")
    {
      G4cout << ">> Female Phantom identified by UserScoreWriter." << G4endl;
    }
  if ((fSex != "female") && (fSex != "male"))
    G4cout << fSex << " can not be defined!" << G4endl;
}

// Sets the section of the phantom as defined through the messenger class
void ICRP110UserScoreWriter::SetPhantomSection(G4String newSection)
{
  fSection = newSection;
  
  if (fSection == "head")
    {
      G4cout << ">> Partial Head Phantom identified by UserScoreWriter." << G4endl;
    }
  if (fSection == "trunk")
    {
      G4cout << ">> Partial Trunk Phantom identified by UserScoreWriter." << G4endl;
    }
  if (fSection == "full")
    {
      G4cout << ">> Custom/Full Phantom identified by UserScoreWriter." << G4endl;
    }
  if ((fSection != "head") && (fSection != "trunk") && (fSection != "full"))
    G4cout << fSection << " can not be defined!" << G4endl;
}
