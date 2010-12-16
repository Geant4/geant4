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
// --------------------------------------------------------------------
// Class implementation
// --------------------------------------------------------------------

#include "G4TSimHelper.h"

G4TSimHelper *gSimHelper = new G4TSimHelper();


ClassImp(G4TSimHelper)

using namespace std;
using namespace ROOT;
using namespace TMath;

//______________________________________________________________________________
void G4TSimHelper::LoadLibraries()
{
	//cout << "Loading Libraries..." << endl;

	gSystem->Load("libRIO.so");
	gSystem->Load("libGeom.so");
	gSystem->Load("libMatrix.so");
	gSystem->Load("libHist.so");
	gSystem->Load("libGraf.so");
	gSystem->Load("libGpad.so");
	gSystem->Load("libGX11.so");
	gSystem->Load("libGui.so");
	gSystem->Load("libGX11TTF.so");
	gSystem->Load("libThread.so");
	gSystem->Load("libNet.so");
	gSystem->Load("libTree.so");
	gSystem->Load("libGraf3d.so");
	gSystem->Load("libTreePlayer.so");
	gSystem->Load("libPostscript.so");
	gSystem->Load("libEG.so");
	gSystem->Load("libGuiBld.so");
	gSystem->Load("libASImage.so");

}

//______________________________________________________________________________
TFile* G4TSimHelper::ExecuteTest(Int_t np , Int_t nn ,  Int_t e , Int_t runNumber , Int_t nbEvents, const TString& dir)
{
	cout << "Executing Chips..." << endl;

	TString pf = dir + "chipstest_default.in";
	cout << "Using parameters: " + pf << endl;

	G4TModelParams ModelParams;
	ModelParams.Load(pf);

	TFile* lastFile = NULL;

	if(runNumber > 0){
		Float_t eConst = 938.272;
		Int_t pdg  = np * 1000 + nn;
		Int_t pdgt = 90000000 + pdg;
		Float_t energy = e + eConst;
		Float_t moment = sqrt((energy * energy) - (eConst * eConst));

		cout << "Number of events to simulate: " << nbEvents << endl;

		// preparing parameters
		ModelParams.GetData().pdgtg = pdgt;
		ModelParams.GetData().nevnt = nbEvents;
		ModelParams.GetData().pdgpr = 2212;

		// Make automatic energy calculation
		ModelParams.GetData().enb = 0;
		ModelParams.GetData().momb = moment;

		TString pfW = TString::Format("%schipstest.in", dir.Data());
		ModelParams.Save(pfW);

		cout << "Starting Geant4 simulations..." << endl;
		gSystem->Exec("$G4INSTALL/bin/$G4SYSTEM/test49");
		gSystem->Sleep(2000);

	}

	lastFile = MakeTree();


	return lastFile;
}

//______________________________________________________________________________
TFile* G4TSimHelper::MakeTree()
{
	Int_t nevts = GetEventsNumber();
	cout << "Total number of events = " << nevts << endl;

	TFile* file = new TFile("chips.root", "recreate");
	TTree* TuplInclTree = new TTree("TuplInclTree", "TuplIncl Tree");

	TuplInclTree->ReadFile("tuplincl.out", "Nevt:MtotD:MtotR:Mprot:Mneut:Mdeut:Mtrit:MHe3:MHe4:Mgam:Mpim:Mpip:Mpi0:MKp:MK0:MKm:MaK0:Meta:Metap:Mrhom:Mrhop:Mrho0:Momega:Mphi:MKS0:MKSC:MaKS0:MaKSC:Mf2:Ma2m:Ma2p:Ma20:Mf2p:ND:PDG:NS:NZ:NN:m:Px:Py:Pz:E");
	cout << "Number of entries = " << TuplInclTree->GetEntries() << endl;

	// save the data
	file->Write();
	return file;
}

//______________________________________________________________________________
Int_t G4TSimHelper::GetEventsNumber()
{
	ifstream 	asciiFile;
	string   	asciiLine;
	Int_t		result = 0;

	// open the file
	asciiFile.open(fEventsNumberFileName.Data());
	if (!asciiFile) {
		cout << "Error opening file"<< fEventsNumberFileName.Data() << endl;
		return 0;
	}

	// process
	getline(asciiFile, asciiLine);
	result = atoi(asciiLine.data());

	asciiFile.close();
	return result;

}
