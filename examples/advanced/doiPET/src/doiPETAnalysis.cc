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

//GEANT4 - Depth-of-Interaction enabled Positron emission tomography (PET) advanced example 

//Authors and contributors

// Author list to be updated, with names of co-authors and contributors from National Institute of Radiological Sciences (NIRS)

// Abdella M. Ahmed (1, 2), Andrew Chacon (1, 2), Harley Rutherford (1, 2),
// Hideaki Tashima (3), Go Akamatsu (3), Akram Mohammadi (3), Eiji Yoshida (3), Taiga Yamaya (3)
// Susanna Guatelli (2), and Mitra Safavi-Naeini (1, 2)

// (1) Australian Nuclear Science and Technology Organisation, Australia
// (2) University of Wollongong, Australia
// (3) National Institute of Radiological Sciences, Japan

//Implemetation of the doiPETAnalysis.cc class
//This implementation mimics readout (or digitizer) of the PET scanner. To mimic realistic PET detector, the signals are blurred. Blurring 
//parameters are given in inputParameters.txt file. Crystal based energy resolution and quantum efficiency has been applied. Deadtime (on 
//each detector block and axially multiplexed detector) is also applied before the event is rejected by the energy window. The units for 
//blurring parameters are in keV (for energy) and ns (nano sec) for dead time. If the units are different, exception will be thrown and the
//program quits. In this class, an ideal PMT is assumed to be placed at the corners of the crystal block. First, the ideal interaction position
//(obtained by G4) is used to determine the distance of the PMT from the interaction point. The signal on each PMT depends on the lateral (2D)
//distance from the PMTs. Light sharing method (reflector based) DOI identification method has been used. If the crystal ID is out of bound, 
//error message will be displayed and the event will be rejected. The output file is single based list-mode ASCII file and can be then be 
//processed into coinsidence list-mode data. As, an option, binary output method is also given.  
//Explanation is given for the methods provided. 


#include "doiPETAnalysis.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <iomanip>
#include "Randomize.hh"
#include "G4SPSRandomGenerator.hh"
#include "doiPETAnalysisMessenger.hh"


doiPETAnalysis* doiPETAnalysis::instance=0;

/////////// Constructor /////////////////////////////////////////////
doiPETAnalysis::doiPETAnalysis()
{
	fAnalysisMessenger = new doiPETAnalysisMessenger(this);

	//Set energy window
	lowerThreshold = 400*keV;
	upperThreshold = 600*keV;


	//give default initial activity. Activity strength is changed in the .mac file
	InitialActivity = 1000000*becquerel;

	//In NEMA NU2, all test is done with F-18
	halfLife = 109.771*60 * s;//The halfLife of a given isotope can be changed via the run.mac file

	//
	totalTime = 0 * s;
	prev_totalTime = 0 * s;
	prev_eventID = 0;

	//
	//Initialize crystal ID
	crystalIDNew_tan = -1;
	crystalIDNew_axial = -1;
	crystalIDNew_DOI = -1;

	//
	scatterIndex = 0;

	//
	numberOfPixel_tan = 32;
	numberOfPixel_axial = 32;

	//Default value for deadtime.
	block_DeadTime = 256*ns;
	module_DeadTime = 0*ns;
	//

	//Crystal blurring parameters. One detector has 1024 crystals. All the crystals have different energy resolution. 
	//So, a range of energy resolution is applied between minumun and maximum values. 
	//The energy resolution can be set in the inputParameter.txt file
	crystalResolutionMin = 0.13;//13%
	crystalResolutionMax = 0.17;//17%

	crystalEnergyRef = 511 * keV;//Energy of reference in which the energy resolution of the crystal is computed

	//The quantum efficiency models the probability for the event to be detected by the photo-detector.
	//The quantum efficiency can be set inputParameter.txt file
	crystalQuantumEfficiency = 1;//100% 
	//

	//intialize deadtime for blocks and modules
	numberOfBlocks_total = numberOfRings * numberOfDetector_perRing; 
	blockTime = new double[numberOfBlocks_total];//for each individual block.
	moduleTime = new double[numberOfBlocks_total];//for axially multiplexed detectors.

	//Initialize the deadtime for each detector and axially multiplexed detector (also called modules)
	for(G4int i = 0; i<numberOfBlocks_total; i++){
		blockTime [i] = 0.0;
		moduleTime [i] = 0.0;
	}

	//Initialize type of output. The default output is single events
	getSinglesData  = false; //default value 
	getCoincidenceData = false;
	numberOfHit = 0;

	//This value is based on the assumption that the shift due to the reflector is half distance from the interaction position to the air gap.
	shiftCoeff = 0.5;
}
////////// Destructor ///////////////////////////////////////////////
doiPETAnalysis::~doiPETAnalysis()
{
	delete fAnalysisMessenger;
	delete [] blockTime;
	delete [] moduleTime;
}

////////// GetInstance /////////////////////////////////////////////
doiPETAnalysis* doiPETAnalysis::GetInstance()
{
	if(instance==0) instance = new doiPETAnalysis();
	return instance;
}
void doiPETAnalysis::Delete()
{
	delete instance;
}

//If there is energy deposition in the phantom by the photon, the scatter index is 1, otherwise it is 0
//Use this for checking
void doiPETAnalysis::GetScatterIndexInPhantom(G4double edepInPhantom){
	if(edepInPhantom>0)scatterIndex = 1;
	else scatterIndex = 0;
}

//Get the source position if the process is annihilation.
//Use this for checking
void doiPETAnalysis::SetSourcePosition(G4ThreeVector spos){
	spositionX = spos.x();
	spositionY = spos.y();
	spositionZ = spos.z();
}


//Set the event ID
//Use this for checking. eventID should not be used to sort coincidence events if realistic simulation is done
void doiPETAnalysis::SetEventID(G4int evID){
	eventID = evID;
}

//
void doiPETAnalysis::GetSizeOfDetector(G4double detSizeDoi, G4double detSizeTan, G4double detSizeAxial){
	sizeOfDetector_DOI = detSizeDoi;
	sizeOfDetector_axial = detSizeTan;
	sizeOfDetector_tangential = detSizeAxial;
}

//
void doiPETAnalysis::SetActivity(G4double newActivity){
	InitialActivity = newActivity;
	G4cout<<"Initial activity: "<<InitialActivity/becquerel<<" Bq."<<G4endl;
}
void doiPETAnalysis::SetIsotopeHalfLife(G4double newHalfLife){
	halfLife = newHalfLife;
	G4cout<<"Half life of the isotope "<<halfLife/s<<" sec."<<G4endl;
}

//The time is based on random time intervals between events. See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3267383/
//This time mimics acquisition time of a PET scanner for a given number of particles. 
void doiPETAnalysis::CalulateAcquisitionTime(){
	//Calculate the strength of activity at t=totaltime using decay equation 
	activityNow = InitialActivity * std::exp(-((0.693147/halfLife)*totalTime)); //ln(2) = 0.693147181

	//Activity based time interval. 
	timeInterval = -std::log(G4UniformRand())*(1./activityNow);
	totalTime = timeInterval+prev_totalTime;
	prev_totalTime = totalTime; 
}

//Apply energy blurring on the crystals. The value of the energy blurring with respect to a reference energy is given in the inputParameter.txt file
G4double doiPETAnalysis::QuantumEffifciency(G4double edep, G4int blkID, G4int cysID)
{
	if(fixedResolution){
		crystalResolution = energyResolution_fixed;
	}
	else{
		crystalResolution = energyResolution_cryDependent[blkID][cysID];
	}
	crystalCoeff = crystalResolution * std::sqrt(crystalEnergyRef);

	G4double QE = G4UniformRand();

	//The quantum efficiency models the probability for the event to be detected by the photo-detector. It can be changed in the inputParameter.txt file
	if(QE <= crystalQuantumEfficiency)
	{
		edep_AfterCrystalBlurring = G4RandGauss::shoot(edep,crystalCoeff*std::sqrt(edep)/2.35);
	}
	else {
		//not detected by the photodetector, eventhough there was an interaction
		edep_AfterCrystalBlurring = 0 *keV;
	}
	return edep_AfterCrystalBlurring;
}

///////// ReadOut ///////////////////////////////////

void doiPETAnalysis::ReadOut(G4int blkID, G4int cryID, G4double interTime, G4double timeAnnih, G4ThreeVector interPos, G4double edep)
{
	blockID = blkID;
	crystalID = cryID;
	interactionTime = interTime;
	time_annihil = timeAnnih;
	interactionPos = interPos;
	totalEdep = edep;

	//Get the time of flight. This is the duration from the annihilation process to the detection of the photons by the scintillator. 
	time_tof = interactionTime - time_annihil;

	//time of the event when detected (timerTag)
	timeStamp = totalTime + time_tof;


	//************************************** Apply dead-time ********************************************//
	//Apply paralizable dead-time in the block beofore events are rejected by the energy window
	if(std::fabs(timeStamp - blockTime[blockID]) >=  block_DeadTime){ //If true, the event is accepted
		blockTime[blockID] = timeStamp;
	}
	else {
		//If the time difference is less than the processing time of the detector (dead time), then the dead time (blockTime) of the block is extended.
		blockTime[blockID] = timeStamp;

		//the event is then rejected
		//continue;
		return;
	}

	//Apply Non-paralyzable dead-time on axially multiplexed detectors (4 detectors are arranged axailly)
	//If the time difference is less than the processing time of the module,  the event is rejected without extending the dead time of the module
	if(std::fabs(timeStamp - moduleTime[blockID]) > module_DeadTime){

		//The following finds the block id's of four blocks which are arranged axially
		for (G4int r_ring = 0; r_ring < numberOfRings; r_ring++){
			if (blockID >= r_ring*numberOfDetector_perRing && blockID <(r_ring + 1)*numberOfDetector_perRing){
				for (G4int m_module = 0; m_module < numberOfRings; m_module++){

					//Set the time of the module (four blocks) the same
					moduleTime[blockID + (m_module - r_ring)*numberOfDetector_perRing] = timeStamp;
				}
			}
		}
	}
	else return;


	/////////////////////////   Write qualified single events based the energy deposition in the detector   ///////////

	if(totalEdep>lowerThreshold && totalEdep<upperThreshold ){

		//identifiy the layer
		DOI_ID =  G4int(crystalID/(numberOfCrystal_tangential * numberOfCrystal_axial));

		//identify the crystal id for each Layer. Now, crystalID_2D can take  0,1, ... numberOfCrystal_tangential x numberOfCrystal_axial 
		crystalID_2D = crystalID - (DOI_ID*numberOfCrystal_tangential * numberOfCrystal_axial);

		//identify the crystal ID in the tangential and axial direction
		crystalID_axial = crystalID_2D/numberOfCrystal_axial;
		crystalID_tangential = crystalID_2D%numberOfCrystal_tangential;

		//Calculate local position of the crystal with respect to the detector. Only the lateral distances (tangential (y) and axial (z) are needed.)
		//G4double posCrystalX = (DOI_ID-((G4double)numberOfCrystal_DOI)/2 + 0.5)*(sizeOfCrystal_DOI + crystalGap_DOI) + interactionPos.x();
		G4double posCrystalY = (crystalID_tangential-((G4double)numberOfCrystal_tangential)/2 + 0.5)*(sizeOfCrystal_tangential + crystalGap_tangential) + interactionPos.y();
		G4double posCrystalZ = (crystalID_axial-((G4double)numberOfCrystal_axial)/2 + 0.5)*(sizeOfCrystal_axial + crystalGap_axial) + interactionPos.z();


		//shiftCoeff = 0.5 is used. This value is based on the assumption that the shift due to the reflector is half distance from the interaction position to the air gap. 
		AngerLogic(DOI_ID, crystalID_tangential, crystalID_axial, posCrystalY, posCrystalZ, totalEdep, shiftCoeff);//
		

		//Single event output. Coincidence events can then be made using the single events.
		if(getSinglesData) WriteOutput();

		//Coincidence output
		if(getCoincidenceData){
			eventID_coin.push_back(eventID);
			blockID_coin.push_back(blockID);
			cryID_axial_coin.push_back(crystalID_axial);
			cryID_tan_coin.push_back(crystalID_tangential);
			edep_coin.push_back(totalEdep);
			cryDOI_coin.push_back(DOI_ID);
			time_coin.push_back(timeStamp);

			numberOfHit++;

			if(numberOfHit == 2){ //two events within the energy window are qualified.
				WriteOutput();
				ResetNumberOfHits();
			}
		}
	}
}



////////// Clear ///////////////////////////////////////////////////
void doiPETAnalysis::ResetNumberOfHits()
{
	numberOfHit = 0;
	eventID_coin.clear();
	blockID_coin.clear();
	cryID_axial_coin.clear();
	cryID_tan_coin.clear();
	edep_coin.clear();
	cryDOI_coin.clear();
	time_coin.clear();

}

//
void doiPETAnalysis::Open(G4String fileName)
{
	if(getSinglesData){
		asciiFileName = fileName + "Singles.data";
	}
	if(getCoincidenceData){
		asciiFileName = fileName + "Coincidence.data";
	}

	ofs.open(asciiFileName.c_str());
	if(!ofs.is_open()){
		G4cerr<<"=== \n File opening Error to write the output ===="<<G4endl;
		exit(0);
	}
	//
#ifdef USEROOT
	if(getSinglesData){ 
		rootFileName = fileName+"Singles.root";
		tSingles = new TTree("tSingles","SinglesTree");
		tSingles->Branch("eventID",&eventID,"eventID/I");
		tSingles->Branch("blockID",&blockID,"blockID/I");
		tSingles->Branch("crystalID_axial",&crystalID_axial,"crystalID_axial/I");
		tSingles->Branch("crystalID_tangential",&crystalID_axial,"crystalID_tangential/I");
		tSingles->Branch("DOI_ID",&DOI_ID0,"DOI_ID/I");
		tSingles->Branch("timeStamp",&timeStamp,"timeStamp/D");
		tSingles->Branch("totalEdep",&totalEdep,"totalEdep/D");
	}

	if(getCoincidenceData){
		rootFileName = fileName+"Coincidence.root";
		tCoincidence = new TTree("tCoincidence","CoincidenceTree");
		//First Single
		tCoincidence->Branch("eventID0",&eventID0,"eventID0/I");
		tCoincidence->Branch("blockID0",&blockID0,"blockID0/I");
		tCoincidence->Branch("crystalID_axial0",&crystalID_axial0,"crystalID_axial0/I");
		tCoincidence->Branch("crystalID_tangential0",&crystalID_axial0,"crystalID_tangential0/I");
		tCoincidence->Branch("DOI_ID0",&DOI_ID0,"DOI_ID0/I");
		tCoincidence->Branch("timeStamp0",&timeStamp0,"timeStamp0/D");
		tCoincidence->Branch("totalEdep0",&totalEdep0,"totalEdep0/D");

		//Second Single
		tCoincidence->Branch("eventID1",&eventID1,"eventID1/I");
		tCoincidence->Branch("blockID1",&blockID1,"blockID1/I");
		tCoincidence->Branch("crystalID_axial1",&crystalID_axial1,"crystalID_axial1/I");
		tCoincidence->Branch("crystalID_tangential1",&crystalID_axial1,"crystalID_tangential1/I");
		tCoincidence->Branch("DOI_ID1",&DOI_ID1,"DOI_ID1/I");
		tCoincidence->Branch("timeStamp1",&timeStamp1,"timeStamp1/D");
		tCoincidence->Branch("totalEdep1",&totalEdep1,"totalEdep1/D");
	}
#endif
	//
}

void doiPETAnalysis::WriteOutput(){
	if(getSinglesData){
		ofs<<eventID<<" "<<blockID<<" "<<crystalID_axial<<" "<<crystalID_tangential<<" "<<DOI_ID<<" "<<std::setprecision(17)<<timeStamp/s<<" "<<std::setprecision(7)<<totalEdep/keV<<G4endl;
#ifdef USEROOT
		tSingles->Fill();
#endif
	}
	if(getCoincidenceData){
		//2 singles will qualify to be in coincidence within the energy window.
		for(G4int i=0; i<2; i++){

			//First Single
			if(i==0){
				eventID0				= eventID_coin[0];
				blockID0				= blockID_coin[0];
				crystalID_axial0		= cryID_axial_coin[0];
				crystalID_tangential0	= cryID_tan_coin[0];
				DOI_ID0					= cryDOI_coin[0];
				timeStamp0				= time_coin[0];
				totalEdep0				= edep_coin[0];
			}
			if(i==1){
				//Second Single
				eventID1				= eventID_coin[1];
				blockID1				= blockID_coin[1];
				crystalID_axial1		= cryID_axial_coin[1];
				crystalID_tangential1	= cryID_tan_coin[1];
				DOI_ID1					= cryDOI_coin[1];
				timeStamp1				= time_coin[1];
				totalEdep1				= edep_coin[1];
			}
		}

		ofs<<eventID0<<" "<<blockID0<<" "<<crystalID_axial0<<" "<<crystalID_tangential0<<" "<<DOI_ID0<<" "<<std::setprecision(17)<<timeStamp0/s<<" "<<std::setprecision(7)<<totalEdep0/keV<<" "
			<<eventID1<<" "<<blockID1<<" "<<crystalID_axial1<<" "<<crystalID_tangential1<<" "<<DOI_ID1<<" "<<std::setprecision(17)<<timeStamp1/s<<" "<<std::setprecision(7)<<totalEdep1/keV<<G4endl;

#ifdef USEROOT
		tCoincidence->Fill();
#endif
	}

}

//
///////// Close /////////////////////////////////////////////////////
void doiPETAnalysis::Close()
{
	//close ascii file
	ofs.close();

	//
#ifdef USEROOT
	TFile f(rootFileName.c_str(),"RECREATE");
	if(getSinglesData){
		tSingles->Write();
		delete tSingles;
	}
	if(getCoincidenceData){
		tCoincidence->Write();
		delete tCoincidence;
	}

	f.Close();


#endif

}

//Place the photomultiplier tube (PMT) at each corner of the detector. 
//The positions of the PMT is with respect to the axis of the detector block
//All the PMTs are placed at the same doi (x) position 
//(at +sizeOfDetector_DOI/2 which is at the top of the detector). 

//The PMT is placed at each corner of the crystal block and is assumed to be an ideal PMT.
//The signal (energy deposition) of each PMT depends on  the distance of the respective 
// PMT from the interaction point
void doiPETAnalysis::PMTPosition(){

	sizeOfDetector_DOI = (numberOfCrystal_DOI * sizeOfCrystal_DOI) + (numberOfCrystal_DOI - 1)*crystalGap_DOI;
	sizeOfDetector_axial = (numberOfCrystal_axial * sizeOfCrystal_axial) + (numberOfCrystal_axial - 1)*crystalGap_axial;
	sizeOfDetector_tangential = (numberOfCrystal_tangential * sizeOfCrystal_tangential) + (numberOfCrystal_tangential - 1)*crystalGap_tangential;

	//Position of PMT1. 
	posPMT1x = sizeOfDetector_DOI/2;//mm
	posPMT1y = -sizeOfDetector_tangential/2;
	posPMT1z = -sizeOfDetector_axial/2;

	//Position of PMT2
	posPMT2x = sizeOfDetector_DOI/2;
	posPMT2y = sizeOfDetector_tangential/2;
	posPMT2z = -sizeOfDetector_axial/2;

	//Position of PMT3 
	posPMT3x = sizeOfDetector_DOI/2;
	posPMT3y = -sizeOfDetector_tangential/2;
	posPMT3z = sizeOfDetector_axial/2;

	//Position of PMT4 
	posPMT4x = sizeOfDetector_DOI/2;
	posPMT4y = sizeOfDetector_tangential/2;
	posPMT4z = sizeOfDetector_axial/2;

	G4cout<<"PMT positions: "<<G4endl;
	G4cout<<"PMT1 (mm) ("<<posPMT1x<<", "<<posPMT1y<<", "<<posPMT1z<<")"<<G4endl;
	G4cout<<"PMT2 (mm) ("<<posPMT2x<<", "<<posPMT2y<<", "<<posPMT2z<<")"<<G4endl;
	G4cout<<"PMT3 (mm) ("<<posPMT3x<<", "<<posPMT3y<<", "<<posPMT3z<<")"<<G4endl;
	G4cout<<"PMT4 (mm) ("<<posPMT4x<<", "<<posPMT4y<<", "<<posPMT4z<<")"<<G4endl;

}

//The blurring parameters are given and can be changed in the inputParameter.txt file
void doiPETAnalysis::BlurringParameters(){
	char inputChar[256];
	std::string inputLine;
	G4String value[7];
	std::string  filename = "inputParameter.txt";
	ifs.open(filename.c_str());
	if(!ifs.good()){
		G4cerr<<"File opening Error: Could not open "<<filename<<G4endl;
		exit(0);
	}
	while(!ifs.eof()){
		ifs.getline(inputChar,256);
		inputLine = inputChar;
		if(inputChar[0]!='#' && inputLine.length()!=0 ){
			if( (std::string::size_type)inputLine.find("block_DeadTime:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				tmpStream >> value[0] >> value[1] >> value[2];
				block_DeadTime = atof(value[1].c_str());
				if(value[2] != "ns"){
					G4cerr<<" Dead time unit is not in nano seconds (ns), Make it in 'ns' "<<G4endl; 
					exit(0);
				}
				block_DeadTime = block_DeadTime*ns;
				G4cout<<"Dead time of the detector: "<<block_DeadTime <<" ns."<<G4endl;
			}
			if( (std::string::size_type)inputLine.find("module_DeadTime:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				tmpStream >> value[0] >> value[1] >> value[2];
				module_DeadTime = atof(value[1].c_str());
				if(value[2] != "ns"){
					G4cerr<<" Dead time unit is not in nano seconds (ns), Make it in 'ns' "<<G4endl; 
					exit(0);
				}
				module_DeadTime = module_DeadTime*ns;
				G4cout<<"Dead time of the module (axially multiplexed detectors): "<<module_DeadTime <<" ns."<<G4endl;
			}
			//
			if( (std::string::size_type)inputLine.find("crystalResolutionMin:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				tmpStream >> value[0] >> value[1];
				crystalResolutionMin = atof(value[1].c_str());
				G4cout<<"crystal Resolution (Min.): "<<crystalResolutionMin*100<< " %." <<G4endl;
			}
			if( (std::string::size_type)inputLine.find("crystalResolutionMax:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				tmpStream >> value[0] >> value[1];
				crystalResolutionMax = atof(value[1].c_str());
				G4cout<<"crystal Resolution (Max.): "<<crystalResolutionMax*100<<" %"<<G4endl;
			}

			//
			if( (std::string::size_type)inputLine.find("fixedResolution:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				tmpStream >> value[0] >> value[1];
				if(value[1]=="true"){
					fixedResolution = true;
					energyResolution_fixed = (crystalResolutionMin + crystalResolutionMax)*0.5;
					G4cout<<"Fixed crystal resolution is used. "<<G4endl;
				}
				else {
					fixedResolution = false;
					//Store into a file if needed.
					//std::string fname = "crystalDependentResolution.txt";
					//std::ofstream outFname(fname.c_str());

					G4cout<<" \n Crystal dependent resolution is used. preparing look-up table .... "<<G4endl;
					energyResolution_cryDependent.resize(numberOfBlocks_total,std::vector<G4double>(numberOfCrystal_tangential*numberOfCrystal_axial*numberOfCrystal_DOI,0));	
					for(G4int i_blk = 0; i_blk < numberOfBlocks_total; i_blk++){
						for(G4int i_cry = 0; i_cry < numberOfCrystal_tangential*numberOfCrystal_axial*numberOfCrystal_DOI; i_cry++){
							energyResolution_cryDependent[i_blk][i_cry] = crystalResolutionMin + (crystalResolutionMax - crystalResolutionMin)*G4UniformRand();
							//store into a file
							//outFname<<i_blk<<" "<<i_cry<<" "<<energyResolution_cryDependent[i_blk][i_cry]<<G4endl;
						}
					}
					G4cout<<"Done. \n"<<G4endl;
					//outFname.close();
				}

			}
			//

			if( (std::string::size_type)inputLine.find("crystalEnergyRef:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				tmpStream >> value[0] >> value[1] >> value[2];
				crystalEnergyRef = atof(value[1].c_str());
				if(value[2] != "keV"){
					G4cerr<<" The unit of reference energy is not in keV, Make it in 'keV' "<<G4endl; 
					exit(0);
				}
				crystalEnergyRef = crystalEnergyRef*keV;
				G4cout<<"Energy of refernce: "<<crystalEnergyRef/keV<<" keV."<<G4endl;
			}
			if( (std::string::size_type)inputLine.find("crystalQuantumEfficiency:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				tmpStream >> value[0] >> value[1];
				crystalQuantumEfficiency = atof(value[1].c_str());
				G4cout<<"Quantum Efficiency "<<crystalQuantumEfficiency*100<< " % "<<G4endl;
			}
			if( (std::string::size_type)inputLine.find("lowerThreshold:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				tmpStream >> value[0] >> value[1] >> value[2];
				lowerThreshold = atof(value[1].c_str());
				if(value[2] != "keV"){
					G4cerr<<" The unit of Lower energy threshold is not in keV, Make it in 'keV' "<<G4endl; 
					exit(0);
				}
				lowerThreshold = lowerThreshold*keV;
				G4cout<<"Lower energy threshold: "<<lowerThreshold/keV<<" keV."<<G4endl;

			}
			if( (std::string::size_type)inputLine.find("upperThreshold:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				tmpStream >> value[0] >> value[1] >> value[2];
				upperThreshold = atof(value[1].c_str());
				if(value[2] != "keV"){
					G4cerr<<" The unit of Upper energy threshold is not in keV, Make it in 'keV' "<<G4endl; 
					exit(0);
				}
				upperThreshold = upperThreshold*keV;
				G4cout<<"Upper energy threshold: "<<upperThreshold/keV<<" keV."<<G4endl;

			}

			if( (std::string::size_type)inputLine.find("numberOfPixel_2D_Pixel:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				tmpStream >> value[0] >> value[1] >> value[2];
				numberOfPixel_axial = atof(value[1].c_str());
				numberOfPixel_tan = atof(value[2].c_str());
				G4cout<<"Number of pixels for a 2D position histogram of the response: "<<numberOfPixel_tan<<" x "<< numberOfPixel_axial <<G4endl;
			}

			//
			if( (std::string::size_type)inputLine.find("TypeOfOutput:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				tmpStream >> value[0] >> value[1];
				if(value[1]=="singlesOutput"){
					getSinglesData = true;
					G4cout<<"Single mode output enabled. "<<G4endl;
				}
				else if(value[1]=="coincidenceOutput") {
					getCoincidenceData = true;
					G4cout<<"Coicidence mode output enabled. "<<G4endl;
				}

			}

		}
	}
	ifs.close();
}

//The following function reads the reflector pattern for each layer. 
//Each layer has different patterns along the tangetial and axial positions.
//For defualt reflector pattern, see https://link.springer.com/article/10.1007/s12194-013-0231-4
//The patter of the reflectors can be changed in the inputParameter.txt file
//The pattern is given as 0 and 1. If there is reflector the value is 1 and if there is no reflector, the value is 0.

void doiPETAnalysis::ReadReflectorPattern(){
	G4cout<<" Reflector pattern is being read "<<G4endl;
	//
	std::vector<std::string> stringReflectorValue;
	//
	char inputChar[256];
	std::string inputLine;

	//open inputParameter.txt to read reflector pattern. 
	std::string  filename = "inputParameter.txt";

	G4String refValue;

	ifs.open(filename.c_str());
	if(!ifs.good()){
		G4cerr<<"File opening Error: Could not open "<<filename<<G4endl;
		exit(0);
	}
	while(!ifs.eof()){
		ifs.getline(inputChar,256);
		inputLine = inputChar;

		//The reflector patter in read from the inputparamter.txt file			
		if(inputChar[0]!='#' && inputLine.length()!=0 ){

			//Reflector patter for Layer1 in the tangential direction
			if( (std::string::size_type)inputLine.find("reflectorLayer1_Tangential:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				while(tmpStream >> refValue){
					stringReflectorValue.push_back(refValue);
					if(stringReflectorValue.size()>1){
						G4int tmp_value = atoi(stringReflectorValue[stringReflectorValue.size()-1].c_str());
						ireflectorLayer1_Tangential.push_back(tmp_value);
					}
				}
			}
			stringReflectorValue.clear();

			//Reflector patter for Layer1 in the axial direction
			if( (std::string::size_type)inputLine.find("reflectorLayer1_Axial:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				while(tmpStream >> refValue){
					stringReflectorValue.push_back(refValue);
					if(stringReflectorValue.size()>1){
						G4int tmp_value = atoi(stringReflectorValue[stringReflectorValue.size()-1].c_str());
						ireflectorLayer1_Axial.push_back(tmp_value);
					}
				}
			}
			stringReflectorValue.clear();

			//Reflector patter for Layer2 in the tangential direction
			if( (std::string::size_type)inputLine.find("reflectorLayer2_Tangential:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				while(tmpStream >> refValue){
					stringReflectorValue.push_back(refValue);
					if(stringReflectorValue.size()>1){
						G4int tmp_value = atoi(stringReflectorValue[stringReflectorValue.size()-1].c_str());
						ireflectorLayer2_Tangential.push_back(tmp_value);
					}
				}
			}
			stringReflectorValue.clear();

			//Reflector patter for Layer2 in the axial direction
			if( (std::string::size_type)inputLine.find("reflectorLayer2_Axial:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				while(tmpStream >> refValue){
					stringReflectorValue.push_back(refValue);
					if(stringReflectorValue.size()>1){
						G4int tmp_value = atoi(stringReflectorValue[stringReflectorValue.size()-1].c_str());
						ireflectorLayer2_Axial.push_back(tmp_value);
					}
				}
			}
			stringReflectorValue.clear();

			//Reflector patter for Layer3 in the tangential direction
			if( (std::string::size_type)inputLine.find("reflectorLayer3_Tangential:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				while(tmpStream >> refValue){
					stringReflectorValue.push_back(refValue);
					if(stringReflectorValue.size()>1){
						G4int tmp_value = atoi(stringReflectorValue[stringReflectorValue.size()-1].c_str());
						ireflectorLayer3_Tangential.push_back(tmp_value);
					}
				}
			}
			stringReflectorValue.clear();

			//Reflector patter for Layer3 in the axial direction
			if( (std::string::size_type)inputLine.find("reflectorLayer3_Axial:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				while(tmpStream >> refValue){
					stringReflectorValue.push_back(refValue);
					if(stringReflectorValue.size()>1){
						G4int tmp_value = atoi(stringReflectorValue[stringReflectorValue.size()-1].c_str());
						ireflectorLayer3_Axial.push_back(tmp_value);
					}
				}
			}
			stringReflectorValue.clear();

			//Reflector patter for Layer4 in the tangential direction
			if( (std::string::size_type)inputLine.find("reflectorLayer4_Tangential:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				while(tmpStream >> refValue){
					stringReflectorValue.push_back(refValue);
					if(stringReflectorValue.size()>1){
						G4int tmp_value = atoi(stringReflectorValue[stringReflectorValue.size()-1].c_str());
						ireflectorLayer4_Tangential.push_back(tmp_value);
					}
				}
			}
			stringReflectorValue.clear();

			//Reflector patter for Layer4 in the axial direction
			if( (std::string::size_type)inputLine.find("reflectorLayer4_Axial:")!=std::string::npos){
				std::istringstream tmpStream(inputLine);
				while(tmpStream >> refValue){
					stringReflectorValue.push_back(refValue);
					if(stringReflectorValue.size()>1){
						G4int tmp_value = atoi(stringReflectorValue[stringReflectorValue.size()-1].c_str());
						ireflectorLayer4_Axial.push_back(tmp_value);
					}
				}
			}
			stringReflectorValue.clear();
		}//#
	}//while(eof)


	//prepare Look up table for crystal identification. 
	G4cout<<"DOI look-up table is being prepared. "<<G4endl;
	std::string outputFileName = "check_2Dposition.txt";// excuted only once
	std::ofstream outFile(outputFileName.c_str());

	G4double crystalPositionY;
	G4double crystalPositionZ;
	doi_table.resize(numberOfCrystal_tangential*numberOfCrystal_axial*numberOfCrystal_DOI,0);

	for(G4int i_DOI = 0; i_DOI<numberOfCrystal_DOI; i_DOI++){
		//crystalPositionX=(i_DOI-((float)numberOfCrystal_DOI)/2 + 0.5)*(sizeOfCrystal_DOI + crystalGap_DOI); //Becuase only lateral distances are used
		for(G4int i_axial=0; i_axial< numberOfCrystal_axial;i_axial++){
			crystalPositionZ = (i_axial-((float)numberOfCrystal_axial)/2 + 0.5)*(sizeOfCrystal_axial + crystalGap_axial);
			for(G4int i_tan=0; i_tan<numberOfCrystal_tangential;i_tan++){
				crystalPositionY=(i_tan-((float)numberOfCrystal_tangential)/2 + 0.5)*(sizeOfCrystal_tangential + crystalGap_tangential);
				AngerLogic(i_DOI, i_tan, i_axial, crystalPositionY, crystalPositionZ, 1, 0.5);
				outFile<<i_DOI<<" "<<i_axial<<" "<<i_tan<<" "<<crystalID_in2D_posHist<<" "<<PositionAngerZ<<" "<<PositionAngerY<<G4endl;
				doi_table[crystalID_in2D_posHist]=i_DOI;
			}
		}
	}
	ifs.close();
}


//Based on ideal photomultiplier tube (PMT) placement, the interaction position of the photon with the detector is calculated using Anger Logic method. 
//The reflectors shifts the response by some distance so that the response can be projected into 2D position histogram. 
//From this 2D position histogram, the new crystal ID (in 3D along the tangential (y), axial (z) and DOI (x)) (after Anger Logic method is applied) can be obtained.
//If the crystal ID after Anger method apllied is out of the give number of crystals (in 3D), then an error message is displayed and the event will be rejected.
void doiPETAnalysis::AngerLogic(G4int i_doi, G4int i_tan, G4int i_axial, G4double posCrystalY, G4double posCrystalZ, G4double Edep, G4double shiftDis)
{
	//1z and 2z are at the same z distance; 3z and 4z are at the same z distance
	//The signal (the energy deposition) is devided into the four PMTs depending on their lateral distances (in the axial and tangential directions) from the interaction position

	//The following is based on symetrical placment of the PMTs 

	//Calculate the axial (z) distance from the position of interaction to each PMT
	dist1z = std::fabs(posPMT1z - posCrystalZ);
	dist2z = std::fabs(posPMT2z - posCrystalZ);
	dist3z = std::fabs(posPMT3z - posCrystalZ);
	dist4z = std::fabs(posPMT4z - posCrystalZ);

	//Calculate the resultant distance
	//dist1z = dist2z, and dist3z = dist4z, so only take two of them or take the average
	//distz = ((dist1z + dist2z) + (dist3z + dist4z))/2; 
	distz = dist1z + dist3z; 

	//1y and 3y are at the same y distance; and 2y and 4y are at the same y distance
	//Calculate the tangential (y)  distance from the position of interaction to each PMT
	dist1y = std::fabs(posPMT1y - posCrystalY);
	dist2y = std::fabs(posPMT2y - posCrystalY);
	dist3y = std::fabs(posPMT3y - posCrystalY);
	dist4y = std::fabs(posPMT4y - posCrystalY);

	//Calculate the resultant distance
	//dist1y = dist3y, and dist2y = dist4y, so only take two of them or take the average
	//disty = ((dist1y + dist3y) + (dist2y+dist4y))/2;
	disty = dist1y + dist2y;

	//signalPMT1z = signalPMT2z, and signalPMT3z = signalPMT4z
	signalPMT1z = Edep * dist3z/(dist1z + dist3z);
	signalPMT3z = Edep * dist1z/(dist1z + dist3z);

	signalPMT2z = Edep * dist4z/(dist2z + dist4z);
	signalPMT4z = Edep * dist2z/(dist2z + dist4z);


	//signalPMT1y = signalPMT3y, and signalPMT2y = signalPMT4y
	signalPMT1y = Edep * dist2y/(dist1y + dist2y);
	signalPMT2y = Edep * dist1y/(dist1y + dist2y);

	signalPMT3y = Edep * dist4y/(dist3y + dist4y);
	signalPMT4y = Edep * dist3y/(dist3y + dist4y);

	//Calculate the signal on each PMT from the 'component' signal
	signalPMT1 = (signalPMT1z +  signalPMT1y)*0.5; 
	signalPMT2 = (signalPMT2z +  signalPMT2y)*0.5;
	signalPMT3 = (signalPMT3z +  signalPMT3y)*0.5;
	signalPMT4 = (signalPMT4z +  signalPMT4y)*0.5;


	signalZplus = (signalPMT3 + signalPMT4);
	signalZminus = (signalPMT1 + signalPMT2);
	signalYplus = (signalPMT2 + signalPMT4);
	signalYminus = (signalPMT1 + signalPMT3);


	//Position of interaction is calculated based on Anger logic method. 
	//To get the position by Anger calculation, the result should be multiplied by the dimenion of the total distance.
	PositionAngerZ = (signalZplus - signalZminus)/(signalZplus + signalZminus)*distz; 
	PositionAngerY = (signalYplus - signalYminus)/(signalYplus + signalYminus)*disty;

	//Find local position in the interacting cystal to estimate the shift due to reflector
	//double localPosX = posCrystalX - (i_doi-((G4double)numberOfCrystal_DOI)/2 + 0.5)*(sizeOfCrystal_DOI + crystalGap_DOI);
	G4double localPosY = posCrystalY - (G4double)(i_tan-((G4double)numberOfCrystal_tangential)/2 + 0.5)*(sizeOfCrystal_tangential + crystalGap_tangential);
	G4double localPosZ = posCrystalZ - (G4double)(i_axial-((G4double)numberOfCrystal_axial)/2 + 0.5)*(sizeOfCrystal_axial + crystalGap_axial);


	G4double crystalPitch_tan = sizeOfCrystal_tangential + crystalGap_tangential;
	G4double crystalPitch_axial = sizeOfCrystal_axial + crystalGap_axial;

	//For detectors with reflector insertion (light sharing), the response is shifted depending on the reflector patter.
	//Here, it is assumed that the shift of the response is equal to half of the distance from the interaction position to the airgap in the lateral (transversal direction)

	//If reflector is only in the left side of the crystal, then response shift is to the right side (away from the reflector).
	//If reflector is only in the right side of the crystal, then response shift is to the left side (away from the reflector).

	//Response shift for 1st Layer
	if(i_doi == 0){
		//If reflector is only in one (left) side of the crystal, then response shifts to the right side (away from the reflector)
		if(ireflectorLayer1_Tangential[i_tan] == 1 && ireflectorLayer1_Tangential[i_tan + 1] == 0) PositionAngerY += (crystalPitch_tan/2 - localPosY)*shiftDis;

		//If reflector is only in one (right) side of the crystal, then response shifts to the left side (away from the reflector)
		if(ireflectorLayer1_Tangential[i_tan] == 0 && ireflectorLayer1_Tangential[i_tan + 1] == 1) PositionAngerY -= (crystalPitch_tan/2 - localPosY)*shiftDis;

		if(ireflectorLayer1_Axial[i_axial] == 1 && ireflectorLayer1_Axial [i_axial + 1] == 0) PositionAngerZ += (crystalPitch_axial/2 - localPosZ)*shiftDis;
		if(ireflectorLayer1_Axial[i_axial] == 0 && ireflectorLayer1_Axial [i_axial + 1] == 1) PositionAngerZ -= (crystalPitch_axial/2 - localPosZ)*shiftDis;
	}
	if(i_doi == 1){ //Response shift for 2nd Layer
		if(ireflectorLayer2_Tangential[i_tan] == 1 && ireflectorLayer2_Tangential[i_tan + 1] == 0) PositionAngerY += (crystalPitch_tan/2 - localPosY)*shiftDis;
		if(ireflectorLayer2_Tangential[i_tan] == 0 && ireflectorLayer2_Tangential[i_tan + 1] == 1) PositionAngerY -= (crystalPitch_tan/2 - localPosY)*shiftDis;

		if(ireflectorLayer2_Axial[i_axial] == 1 && ireflectorLayer2_Axial [i_axial + 1] == 0) PositionAngerZ += (crystalPitch_axial/2 - localPosZ)*shiftDis;
		if(ireflectorLayer2_Axial[i_axial] == 0 && ireflectorLayer2_Axial [i_axial + 1] == 1) PositionAngerZ -= (crystalPitch_axial/2 - localPosZ)*shiftDis;
	}
	if(i_doi == 2){ //Response shift for 3rd Layer
		if(ireflectorLayer3_Tangential[i_tan] == 1 && ireflectorLayer3_Tangential[i_tan + 1] == 0) PositionAngerY += (crystalPitch_tan/2 - localPosY)*shiftDis;
		if(ireflectorLayer3_Tangential[i_tan] == 0 && ireflectorLayer3_Tangential[i_tan + 1] == 1) PositionAngerY -= (crystalPitch_tan/2 - localPosY)*shiftDis;

		if(ireflectorLayer3_Axial[i_axial] == 1 && ireflectorLayer3_Axial [i_axial + 1] == 0) PositionAngerZ += (crystalPitch_axial/2 - localPosZ)*shiftDis;
		if(ireflectorLayer3_Axial[i_axial] == 0 && ireflectorLayer3_Axial [i_axial + 1] == 1) PositionAngerZ -= (crystalPitch_axial/2 - localPosZ)*shiftDis;
	}
	if(i_doi == 3){ //Response shift for 4th Layer
		if(ireflectorLayer4_Tangential[i_tan] == 1 && ireflectorLayer4_Tangential[i_tan + 1] == 0) PositionAngerY += (crystalPitch_tan/2 - localPosY)*shiftDis;
		if(ireflectorLayer4_Tangential[i_tan] == 0 && ireflectorLayer4_Tangential[i_tan + 1] == 1) PositionAngerY -= (crystalPitch_tan/2 - localPosY)*shiftDis;

		if(ireflectorLayer4_Axial[i_axial] == 1 && ireflectorLayer4_Axial [i_axial + 1] == 0) PositionAngerZ += (crystalPitch_axial/2 - localPosZ)*shiftDis;
		if(ireflectorLayer4_Axial[i_axial] == 0 && ireflectorLayer4_Axial [i_axial + 1] == 1) PositionAngerZ -= (crystalPitch_axial/2 - localPosZ)*shiftDis;
	}

	//The main purpose of shifting the response is to be able to project the response of all the crytal elements into a 2D position histogram so that we can identify the DOI layer 
	//by comparing with a look-up-table which is prepared based on the reflector insertion. 

	//The crystal ID in 2D position histogram along the axial (z) direction. It can have values of: 0, 1, .. , 31, in 32x32 pixel position histogram
	crystalID_in2D_posHist_axial = (G4int)((G4double)PositionAngerZ/(G4double)(crystalPitch_axial*0.5) + (G4double)(numberOfPixel_axial - 1.0)*0.5 + 0.5);//0.5 is added for round off

	//The crystal ID in 2D position histogram along the tangential (y) direction. It can have values of: 0, 1, .. , 31, in 32x32 pixel position histogram
	crystalID_in2D_posHist_tan =   (G4int)((G4double)PositionAngerY/(G4double)(crystalPitch_tan*0.5) + (G4double)(numberOfPixel_tan - 1.0)*0.5 + 0.5);//y_ID

	//continuous crystal ID in the 2D position histogram. It will be from 0 to 1023 (in the case of 16x16x4 crystal array). 
	crystalID_in2D_posHist = crystalID_in2D_posHist_axial + crystalID_in2D_posHist_tan * numberOfPixel_tan;//32;


	//Now, lets find the crystal ID in 3D after applying Anger Logic calculation. NOTE that its value can be the same as the original crystal ID or not.

	//Crystal ID along the tangential diraction after Anger Logic calculation
	crystalIDNew_tan = (G4int)(crystalID_in2D_posHist_tan/2);

	//Crystal ID along the axial diraction after Anger Logic calculation
	crystalIDNew_axial = (G4int)(crystalID_in2D_posHist_axial/2);

	////Crystal ID along the DOI diraction after Anger Logic calculation
	if(crystalID_in2D_posHist>numberOfCrystal_tangential*numberOfCrystal_axial*numberOfCrystal_DOI) return;
	crystalIDNew_DOI = doi_table[crystalID_in2D_posHist];

	//If the crsytal ID is beyond the given the number of crystal in the detector, the following is excecuted and the event will be rejected
	if(crystalIDNew_tan < 0 || crystalIDNew_axial < 0 || crystalIDNew_DOI < 0 ||
		crystalIDNew_tan >= numberOfCrystal_tangential || crystalIDNew_axial >= numberOfCrystal_axial || crystalIDNew_DOI >= numberOfCrystal_DOI){
			return;
	}

	CrystalIDAfterAngerLogic(crystalIDNew_tan,crystalIDNew_axial,crystalIDNew_DOI);	
}

void doiPETAnalysis::CrystalIDAfterAngerLogic(G4int i_tan, G4int i_axial, G4int i_doi){
	crystalID_tangential = i_tan;
	crystalID_axial = i_axial;
	DOI_ID = i_doi;
}
