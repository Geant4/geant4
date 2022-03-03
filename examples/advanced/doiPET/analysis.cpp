//The following code is for anylizing axial sensitivity from coincidence list-mode data
//It takes source position of the event (of those which are detected by the scanner) and analises axial sensitivty.
//by Abdella M. Ahmed, 2020

#define _USE_MATH_DEFINES
#include <iostream>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <string>
#include <limits>
#include <stdio.h>
#include <random>
#include <string>

#define UseROOT //ASCII or UseROOT

#ifdef UseROOT
//for root
#include "TRandom3.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TSystem.h"
#endif




#define AxialLength 216 //(mm), Axial length of the scanner

//Change this number based of the number of primary particles simulated in the "run.mac" file
#define NumberOfPrimaries 3000000 //Number of particles simulated

using namespace std;
int main(){

	cout<<"============Axial Sensitivity Analysis =========================="<<endl;

	int eventID0, blockID0, crystalID_axial0, crystalID_tangential0, DOI_ID0; 
	double timeStamp0, totalEdep0;
    int eventID1, blockID1, crystalID_axial1, crystalID_tangential1, DOI_ID1;
	double timeStamp1, totalEdep1;
	double  spositionX, spositionY, spositionZ; //source position

	double z_offset = 0.0;//Axial offset position where the plane is located
	double planeWidth = 3;// (mm)
	int planeNumber;

	float total_sensitivity = 0.0;

	int numberOfSlices = int (AxialLength/planeWidth);
	cout<<"Number of axial planes (slices) are: " <<numberOfSlices<<endl;
	double Counts_per_plane[numberOfSlices];

	ofstream OutFile("axial_sensitivity.csv");
	string filename = "resultCoincidence.data";
	string filepath = "";//provide the file path if it is stored in a different location.
	ifstream InFile(filepath+filename);

	if (!InFile.is_open())
	{
		cout << "Unable to open input file to read .... " << endl;
		return 1;
	}
	if (!OutFile.is_open())
	{
		cout << "Unable to open out file to write ....";
		return 1;
	}

	OutFile << "PlaneNmber" << "," << "Z(mm)" << "," << "Sensitivity(%)" << endl;

	for (int i_plane = 0; i_plane < numberOfSlices; i_plane++){
		Counts_per_plane[i_plane] = 0;
	}
#ifdef ASCII
    cout<<"\nASCII coincidence list-mode data is being analised..."<<endl;
	while (InFile >> eventID0 >> blockID0 >> crystalID_axial0 >> crystalID_tangential0 >> DOI_ID0 >> timeStamp0 >> totalEdep0 >> 
		             eventID1 >> blockID1 >> crystalID_axial1 >> crystalID_tangential1 >> DOI_ID1 >> timeStamp1 >> totalEdep1 >> 
					 spositionX >> spositionY >> spositionZ){
		
		planeNumber = int(spositionZ/planeWidth + numberOfSlices/2 - 0.5 + 0.5);
		Counts_per_plane[planeNumber]++;
		total_sensitivity++;
	}
#endif
	
	
#ifdef UseROOT
	cout<<"\nROOT coincidence list-mode data is being analised..."<<endl;
	TFile *f = new TFile("resultCoincidence.root","READ");
	TTree *Singles = (TTree*)gDirectory->Get("Coincidence");

	Singles->SetBranchAddress("eventID0",&eventID0);
	Singles->SetBranchAddress("blockID0",&blockID0);
	Singles->SetBranchAddress("crystalID_axial0",&crystalID_axial0);
	Singles->SetBranchAddress("crystalID_tangential0",&crystalID_tangential0);
	Singles->SetBranchAddress("DOI_ID0",&DOI_ID0);
	Singles->SetBranchAddress("timeStamp0",&timeStamp0);
	Singles->SetBranchAddress("totalEdep0",&totalEdep0);

	Singles->SetBranchAddress("eventID1",&eventID1);
	Singles->SetBranchAddress("blockID1",&blockID1);
	Singles->SetBranchAddress("crystalID_axial1",&crystalID_axial1);
	Singles->SetBranchAddress("crystalID_tangential1",&crystalID_tangential1);
	Singles->SetBranchAddress("DOI_ID1",&DOI_ID1);
	Singles->SetBranchAddress("timeStamp1",&timeStamp1);
	Singles->SetBranchAddress("totalEdep1",&totalEdep1);

	Singles->SetBranchAddress("spositionX",&spositionX);
	Singles->SetBranchAddress("spositionY",&spositionY);
	Singles->SetBranchAddress("spositionZ",&spositionZ);

	//
	int nentries = 0;
	nentries = Singles->GetEntries();
	for(int entry=0; entry<nentries; entry++){
		Singles->GetEntry(entry);      
		planeNumber = int(spositionZ/planeWidth + numberOfSlices/2 - 0.5 + 0.5);
		Counts_per_plane[planeNumber]++;
		total_sensitivity++;
	}
#endif 

	
	//Save it into CSV file
	for (int i_plane = 0; i_plane < numberOfSlices; i_plane++){
		//Axial mid-position of the plane or slice
		z_offset = (i_plane - numberOfSlices/2 + 0.5)*planeWidth;

		OutFile << i_plane << "," << z_offset << "," << (Counts_per_plane[i_plane]/NumberOfPrimaries)*100 << endl;
	}
	cout<<"\nSensitivity evaluation has completed."<<endl;
	cout<<"\nTotal Sensitivity (Number of coinsidence per total number of pair of photons): "<<(total_sensitivity/NumberOfPrimaries)*100 << "%"<<endl;
	InFile.close();
	OutFile.close();

	return 0;
}
