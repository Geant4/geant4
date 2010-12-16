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

#include "G4TDataItem.h"


ClassImp(G4TDataItem)

using namespace std;
using namespace ROOT;
using namespace TMath;



//______________________________________________________________________________
void G4TDataItem::LoadFromASCII(const TString& filename, Bool_t load3)
{

	gROOT->cd();

	TString tname = GetHeader();
	TString vname = "t:s";

	TTree* DataTree = new TTree(tname.Data(), tname.Data());
	fData = DataTree;

	if(load3)
	{
		vname = "t:s:d";
		DataTree->ReadFile(filename.Data(),vname.Data());
	}
	if(fHeader.fSecondaryParticlePDG == 1000000010 || fHeader.fSecondaryParticlePDG == 2112) // neutrons
	{

		vname = "t:s:r";
		DataTree->ReadFile(filename.Data(),vname.Data());

		// wait for tree
		while(DataTree == NULL)
			gSystem->Sleep(10);

		// TODO: OPTIMIZE (need a better way to do this)
		// ve/cop $SIGMA(r1[a]*s1[a]/100.) d1[a]
		Float_t dval = 0;
		Float_t rval = 0;
		Float_t sval = 0;
		TBranch* branchR = DataTree->GetBranch("r");
		TBranch* branchS = DataTree->GetBranch("s");
		TBranch* branchD = DataTree->Branch("d", &dval, "d");

		branchR->SetAddress(&rval);
		branchS->SetAddress(&sval);
		Int_t rEntries = DataTree->GetEntries();
		for(Int_t x = 0; x< rEntries; x++)
		{
			branchR->GetEntry(x);
			branchS->GetEntry(x);
			dval = rval * sval / 100;
			branchD->Fill();
		}

	}
	else // all other secondaries
	{

		DataTree->ReadFile(filename.Data(),vname.Data());

		// wait for tree
		while(DataTree == NULL)
			gSystem->Sleep(10);

		// TODO: OPTIMIZE (need a better way to do this)
		// e/cop $SIGMA(s[j][a]/20.) d[j][a]
		Float_t dval = 0;
		Float_t sval = 0;
		TBranch* branchS = DataTree->GetBranch("s");
		TBranch* branchD = DataTree->Branch("d", &dval, "d");

		branchS->SetAddress(&sval);
		Int_t rEntries = DataTree->GetEntries();
		for(Int_t x = 0; x< rEntries; x++)
		{
			branchS->GetEntry(x);
			dval = sval / 20;
			branchD->Fill();
		}
	}


}

//______________________________________________________________________________
TString G4TDataItem::HeaderToString(DataItemObjectHeader_t header) const
{
	TString result = TString::Format("%s_%d_%d_%g_%g_%d_%d_%d_%d",
		  gParticlesDAL->GetParticleName(fHeader.fSecondaryParticlePDG).Data(),
		  fHeader.fCutVar,
		  fHeader.fCutUnits,
		  fHeader.fCutValue,
		  fHeader.fCutDelta,
		  fHeader.fFunctionVar,
		  fHeader.fFunctionUnits,
		  fHeader.fArgumentVar,
		  fHeader.fArgumentUnits);
	return result;
}

//______________________________________________________________________________
DataItemObjectHeader_t G4TDataItem::StringToHeader(TString headerStr) const
{
	DataItemObjectHeader_t header;

	TObjArray* tokens = headerStr.Tokenize("_");
	for(int i = 0; i< tokens->GetEntriesFast(); ++i)
	{
		TObjString *T = (TObjString*) (*tokens)[i];
		TString value = T->GetString();


		if(i == 0) header.fSecondaryParticlePDG = gParticlesDAL->GetPDG(value.Data());
		if(i == 1) header.fCutVar = (CutEnum)atoi(value.Data());
		if(i == 2) header.fCutUnits = (UnitsEnum)atoi(value.Data());
		if(i == 3) header.fCutValue = atof(value.Data());
		if(i == 4) header.fCutDelta = atof(value.Data());
		if(i == 5) header.fFunctionVar = (FuncEnum)atoi(value.Data());
		if(i == 6) header.fFunctionUnits = (UnitsEnum)atoi(value.Data());
		if(i == 7) header.fArgumentVar = (ArgEnum)atoi(value.Data());
		if(i == 8) header.fArgumentUnits = (UnitsEnum)atoi(value.Data());
	}
	return header;
}

//______________________________________________________________________________
TString	G4TDataItem::GetFormula()
{
	Double_t mass = gParticlesDAL->GetParticleMass(fHeader.fSecondaryParticlePDG);
	TString P = TString::Format("(sqrt(t * (2. * %g + t)))", mass);
	TString DFormula = TString::Format("(s/%s) : t : (d/%s)", P.Data(), P.Data());
	return DFormula;
}

//______________________________________________________________________________
DataItemLimit_t G4TDataItem::GetLimits()
{
	DataItemLimit_t result;
	Bool_t			initialized = false;
	if(fData != 0){
		fData->Draw(this->GetFormula().Data(),"","goff");

		Long64_t nb = fData->GetSelectedRows();
		Double_t* p = fData->GetV1();
		for(Int_t i = 0; i < nb; ++i)
		{
			 if(!initialized)
			 {
				 result.fMax = *p;
				 result.fMin = *p;
				 initialized = true;
			 }

			 if(*p > result.fMax) result.fMax = *p;
			 if(*p < result.fMin) result.fMin = *p;



			 p++;
		}

	}
	//cout << "limits: min = " << result.fMin << " max = " << result.fMax << endl;

	return result;
}

//______________________________________________________________________________
Double_t G4TDataItem::GetAngleInRadians() const
{
	// Set Angle
	Double_t a = GetCutValue();
	Double_t arad = 0;
	if (GetCutUnits() != Radians){
		arad = a;
	}if(GetCutUnits() == Degrees ){
		arad = 3.141593 * a / 180;
	}else{
		cout << "ERROR: Angle must be in degrees or radians, aborting..." << endl;
	}
	return arad;
}

//______________________________________________________________________________
TString	G4TDataItem::GetHistogramName(Int_t additionalIndex ) const
{
	if(additionalIndex == -1)
		return TString::Format("Histogram_%s", HeaderToString(fHeader).Data());
	return TString::Format("Histogram_%s_%d", HeaderToString(fHeader).Data(), additionalIndex);
}

//______________________________________________________________________________
TString G4TDataItem::GetHeader() const
{
	  return HeaderToString(fHeader);
}


//______________________________________________________________________________
Int_t G4TDataItem::GetSecondaryParticlePDG() const
{
	return fHeader.fSecondaryParticlePDG;
}

//______________________________________________________________________________
CutEnum G4TDataItem::GetCutVar() const
{
	return fHeader.fCutVar;
}

//______________________________________________________________________________
UnitsEnum G4TDataItem::GetCutUnits() const
{
	return fHeader.fCutUnits;
}

//______________________________________________________________________________
Double_t G4TDataItem::GetCutDelta() const
{
	return fHeader.fCutDelta;
}

//______________________________________________________________________________
Double_t G4TDataItem::GetCutValue() const
{
	return fHeader.fCutValue;
}

//______________________________________________________________________________
FuncEnum G4TDataItem::GetFunctionVar() const
{
	return fHeader.fFunctionVar;
}

//______________________________________________________________________________
UnitsEnum G4TDataItem::GetFunctionUnits() const
{
	return fHeader.fFunctionUnits;
}

//______________________________________________________________________________
ArgEnum G4TDataItem::GetArgumentVar() const
{
	return fHeader.fArgumentVar;
}

//______________________________________________________________________________
UnitsEnum G4TDataItem::GetArgumentUnits() const
{
	return fHeader.fArgumentUnits;
}

//______________________________________________________________________________
TH1F* G4TDataItem::GetHistogram() const
{
	return fHistogram;
}

//______________________________________________________________________________
TTree* G4TDataItem::GetData() const
{
	return fData;
}

//______________________________________________________________________________
void G4TDataItem::SetHistogram(TH1F *fHistogram)
{
	this->fHistogram = fHistogram;
}

//______________________________________________________________________________
void G4TDataItem::SetData(TTree *fData)
{
	  this->fData = fData;
}
