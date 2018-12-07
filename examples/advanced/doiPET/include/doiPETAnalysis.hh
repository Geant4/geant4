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


#ifndef doiPETAnalysis_h 
#define doiPETAnalysis_h  1

#include "doiPETGlobalParameters.hh"
#include "globals.hh"
#include <vector>
#include <time.h>
#include <map>
#include <set>
#include "G4ThreeVector.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <algorithm> 

#ifdef USEROOT
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "TRandom3.h"
#include "TFile.h"
#include "TNtuple.h"
#include "g4root.hh"
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TSystem.h"
#pragma GCC diagnostic pop
#endif

class doiPETAnalysisMessenger;

//class InteractionInformation; 

class doiPETAnalysis
{
private:
	doiPETAnalysis();

public:
	~doiPETAnalysis();
	static doiPETAnalysis* GetInstance();
	void FindInteractingCrystal();
	void Open(G4String);
	void Close();
	void Delete();
	void ResetNumberOfHits();
	void Write(/*G4int, G4int, G4int, G4double*/);
	void WriteOutput();

	//void GetIntractionInfomation(InteractionInformation*);

	void GetParentParticleName(G4String);
	void GetSizeOfDetector (G4double, G4double, G4double);
	void GetScatterIndexInPhantom(G4double);

	void SetSourcePosition(G4ThreeVector);//
	void SetEventID(G4int);

	void BlurringParameters();
	void GetTimeOfAnnihilation(G4double);

	void PMTPosition();
	void AngerLogic(G4int, G4int, G4int, G4double, G4double, G4double, G4double);
	void ReadReflectorPattern();

	void SetActivity(G4double);
	void SetIsotopeHalfLife(G4double);
	void CrystalIDAfterAngerLogic(G4int, G4int, G4int);
	void TypeOfOutput(G4String);//Single or coincidence list-mode data
	void CalulateAcquisitionTime();
	//G4double QuantumEffifciency(G4double);
	G4double QuantumEffifciency(G4double, G4int, G4int);
	void ReadOut(G4int, G4int, G4double, G4double, G4ThreeVector, G4double);

private:
	static doiPETAnalysis* instance;
	doiPETAnalysisMessenger* fAnalysisMessenger;
	//std::multimap< G4int, InteractionInformation* > mapBlockInteraction;
	std::set<G4int> setBlockInteraction;

	G4double upperThreshold, lowerThreshold;	



	//G4ThreeVector sourcePosition;



	//
	G4int scatterIndex;

	G4String parentParticleName;//

	//
	G4int numberofInteractions;
	G4int countCoincidence;

	G4int numberOfBlocks_total;

	G4double sizeOfDetector_DOI,sizeOfDetector_axial,sizeOfDetector_tangential;

	//Virtual position of the PMT
	G4double signalPMT1, signalPMT2, signalPMT3, signalPMT4;

	G4double posPMT1x, posPMT2x, posPMT3x, posPMT4x;
	G4double posPMT1y, posPMT2y, posPMT3y, posPMT4y;
	G4double posPMT1z, posPMT2z, posPMT3z, posPMT4z;

	//
	G4double signalPMT1z, signalPMT2z, signalPMT3z, signalPMT4z;
	G4double signalPMT1y, signalPMT2y, signalPMT3y, signalPMT4y;

	//
	G4double signalZplus, signalZminus; 
	G4double signalYplus, signalYminus;
	//

	G4double dist1z, dist2z, dist3z, dist4z, distz;
	G4double dist1y, dist2y, dist3y, dist4y, disty;

	G4double shiftCoeff;

	G4double PositionAngerZ, PositionAngerY;

	//reflector pattern
	std::vector<G4int> ireflectorLayer1_Tangential;
	std::vector<G4int> ireflectorLayer1_Axial;
	std::vector<G4int> ireflectorLayer2_Tangential;
	std::vector<G4int> ireflectorLayer2_Axial;
	std::vector<G4int> ireflectorLayer3_Tangential;
	std::vector<G4int> ireflectorLayer3_Axial;
	std::vector<G4int> ireflectorLayer4_Tangential;
	std::vector<G4int> ireflectorLayer4_Axial;
	std::vector<G4int> doi_table;
	//

	//The number of pixes for the 2D position histogram after Anger Logic calculation 
	G4int numberOfPixel_axial;
	G4int numberOfPixel_tan;

	//source position
	G4double spositionX;
	G4double spositionY;
	G4double spositionZ;

	//interaction position with respect to the crystal axis
	G4ThreeVector interactionPos;


	G4double interactionTime;

	G4int crystalID;//contineous crystal ID in 3D
	G4int crystalID_2D;

	G4int prev_eventID;

	//Single output
	G4int eventID;
	G4int blockID;
	G4int crystalID_axial;
	G4int crystalID_tangential;
	G4int DOI_ID;
	G4double timeStamp;
	G4double totalEdep;

	//coincidence output
	G4int eventID0,					eventID1;
	G4int blockID0,					blockID1;
	G4int crystalID_axial0,			crystalID_axial1;
	G4int crystalID_tangential0,	crystalID_tangential1;
	G4int DOI_ID0,					DOI_ID1;
	G4double timeStamp0,			timeStamp1;
	G4double totalEdep0,			totalEdep1;

	//choice for the user
	G4bool getSinglesData;
	G4bool getCoincidenceData;

	G4String outputData;
	G4int numberOfHit;
	std::vector<G4int> eventID_coin;
	std::vector<G4double> edep_coin;
	std::vector<G4int>blockID_coin;
	std::vector<G4int> cryID_axial_coin;
	std::vector<G4int> cryID_tan_coin;
	std::vector<G4int> cryDOI_coin;
	std::vector<G4double> time_coin;

	//Crystal IDs after Anger Logic calculation
	G4int crystalIDNew_DOI, crystalIDNew_tan, crystalIDNew_axial;

	//Crystal ID in the 2D position histogram along the axial and tangetial direction 
	G4int crystalID_in2D_posHist_axial, crystalID_in2D_posHist_tan;

	//continous crystal ID after after Anger Logic. 
	G4int crystalID_in2D_posHist;


	//Crystal blurring
	G4double crystalResolution;
	G4double crystalResolutionMin;//
	G4double crystalResolutionMax;//

	//G4bool variableResolution;
	G4bool fixedResolution;

	G4double energyResolution_fixed;
	std::vector<std::vector<G4double>> energyResolution_cryDependent;

	G4double crystalEnergyRef;//This 511 keV
	G4double crystalQuantumEfficiency;//
	G4double edep_AfterCrystalBlurring;
	G4double crystalCoeff;
	G4double sigma_energyResolution;

	G4double totalTime;
	G4double prev_totalTime;
	G4double timeInterval;
	G4double time_annihil;
	G4double time_tof;
	G4double block_DeadTime;
	G4double module_DeadTime;

	G4double *blockTime;
	G4double *moduleTime;

	//
	G4double activityNow;
	G4double InitialActivity;
	G4double halfLife;

	G4String simulationType;

	//
	//for output file to write results
	std::ofstream ofs;
	G4String asciiFileName;
	G4String rootFileName;

#ifdef USEROOT
	TTree* tSingles;
	TTree* tCoincidence;
	//TH1F*hb;
#endif

	//input file to read reflector pattern
	std::ifstream ifs;
};

#endif
