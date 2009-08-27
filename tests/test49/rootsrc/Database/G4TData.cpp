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

#include "G4TData.h"


ClassImp(G4TData)

using namespace std;
using namespace ROOT;
using namespace TMath;

//______________________________________________________________________________
void G4TData::Save()
{
	TString hfname = fDirectory + GetHeader() + ".root";
	//cout << "Saving the DataObject to file " << hfname << "..." << endl;

	// prepare the file
	TFile hfile(hfname, "recreate");
	if(hfile.IsZombie()){
		return;
	}

	TObjArray hlist(0); // array to hold the list of objects

	for(UInt_t i = 0; i < fItems.size(); ++i)
	{
		TString subHeader = fItems[i]->GetHeader();
		TTree* data = fItems[i]->GetData();
		hlist.Add(data);
	}

	TList* paramsList = new TList();
	TNamed* CrossSection = new TNamed("CrossSection", TString::Format("%g", fCrossSection).Data());
	TNamed* NumberOfEvents = new TNamed("NumberOfEvents", TString::Format("%d", fNumberOfEvents).Data());

	paramsList->Add(CrossSection);
	paramsList->Add(NumberOfEvents);

	hfile.WriteObject(paramsList,"Parameters","");



	hlist.Write();
	hfile.Close();

	TString name = this->GetHeader() + ".root";
	gCatalog->Insert(name);
}

//______________________________________________________________________________
void G4TData::Load(Int_t secondaryPDG)
{
	TString hfname = fDirectory + HeaderToString(fHeader) + ".root";
	TFile* hfile = new TFile(hfname.Data());
	if(hfile->IsZombie()){
		return;
	}

	TString secName;
	if(secondaryPDG != 0)
	{
		secName = gParticlesDAL->GetParticleName(secondaryPDG);
	}

	TIter next(hfile->GetListOfKeys());
	TKey* key;

	while ((key = (TKey*)next()))
	{
		TString name(key->GetName());
		if(name == "Parameters")
		{
			TList* paramsList = (TList*)hfile->Get(name.Data());
			fCrossSection 	=  atof(((TNamed*)paramsList->At(0))->GetTitle());
			fNumberOfEvents	=  atoi(((TNamed*)paramsList->At(1))->GetTitle());
			continue;
		}

		if(secondaryPDG != 0 && !name.BeginsWith(secName)){
			continue;
		}

		G4TDataItem* item = this->AddItem(name);
		item->SetData( (TTree*) hfile->Get(name.Data()) );

	}

	fIsLoaded = true;
}



//______________________________________________________________________________
G4TDataItem* G4TData::AddItem(
	Int_t 		SecondaryParticlePDG,
	CutEnum 	CutVar,
	UnitsEnum 	CutUnits,
	Double_t	CutValue,
	Double_t	CutDelta,
	FuncEnum	FunctionVar,
	UnitsEnum	FunctionUnits,
	ArgEnum		ArgumentVar,
	UnitsEnum	ArgumentUnits)
{
	G4TDataItem* item = new G4TDataItem(SecondaryParticlePDG, CutVar, CutUnits, CutValue, CutDelta, FunctionVar,
			FunctionUnits, ArgumentVar,ArgumentUnits);
	fItems.push_back(item);
	return item;
}


//______________________________________________________________________________
G4TDataItem* G4TData::AddItem(TString const& headerStr)
{
	G4TDataItem* item = new G4TDataItem(headerStr);
	fItems.push_back(item);
	return item;
}


//______________________________________________________________________________
TString G4TData::HeaderToString(DataObjectHeader_t header) const
{
	if(header.fTypeVar != E_Kin || header.fTypeUnits != MeV)
	{
		cout << "Error: HeaderToString() typeVar and/or typeUnits conversion is not implemented! Use E_Kin with MeV" << endl;
		return "";
	}
	TString result = TString::Format("%s_%g_%s_%s",
		  gParticlesDAL->GetParticleName(header.fProjectilePDG).Data() ,
		  header.fTypeValue,
		  gParticlesDAL->GetParticleName(header.fTargetPDG).Data(),
		  header.fModelName.Data() );
	return result;
}

//______________________________________________________________________________
DataObjectHeader_t G4TData::StringToHeader(TString headerStr) const
{
	DataObjectHeader_t header;
	TString workStr(headerStr.ReplaceAll(".root",""));

	header.fTypeUnits = MeV;
	header.fTypeVar	  = E_Kin;

	TObjArray* tokens = workStr.Tokenize("_");
	for(int i = 0; i< tokens->GetEntriesFast(); ++i)
	{
		TObjString *T = (TObjString*) (*tokens)[i];
		TString value = T->GetString();

		if(i == 0) header.fProjectilePDG 	= gParticlesDAL->GetPDG(value.Data());
		if(i == 1) header.fTypeValue 		= atof(value.Data());
		if(i == 2) header.fTargetPDG	 	= gParticlesDAL->GetPDG(value.Data());
		if(i == 3) header.fModelName		= value.Data();
	}
	return header;
}

//______________________________________________________________________________
TString G4TData::GetHeader() const
{
	  return HeaderToString(fHeader);
}


//______________________________________________________________________________
void G4TData::PrepareHistograms(Double_t hnbin, Double_t hlxmin, Double_t hlxmax, Int_t particleIdx, Int_t additionalIndex)
{
	// Keep them in memory right now
	gROOT->cd();


	for(UInt_t i = 0; i < fItems.size(); ++i)
	{
		G4TDataItem* item = fItems[i];
		Double_t a 		= item->GetCutValue();
		TString hname 	= item->GetHistogramName(additionalIndex);
		TString htitle	= TString::Format("T_{n} at %g deg (MeV)", a );

		//cout << "Creating histogram " << hname.Data() << " at angle " << a  << endl;
		TObject* existing = gROOT->FindObject(hname.Data());
		if(existing != 0)
			delete[] existing;
		TH1F* hist = new TH1F(hname.Data(), htitle.Data(), hnbin, hlxmin, hlxmax);
		item->SetHistogram(hist);
	}

}

//______________________________________________________________________________
Color_t G4TData::GetRenderColor() const
{
	return fRenderColor;
}

//______________________________________________________________________________
void G4TData::SetRenderColor(Color_t fRenderColor)
{
	this->fRenderColor = fRenderColor;
}

//______________________________________________________________________________
TString G4TData::GetModelName() const
{
	return fHeader.fModelName;
}

//______________________________________________________________________________
vector<G4TDataItem*>	G4TData::GetItems()
{
	return fItems;
}

//______________________________________________________________________________
vector<Double_t> G4TData::GetCutValues()
{
	vector<Double_t> resultVector;
	for(UInt_t i = 0; i < fItems.size(); ++i)
	{
		G4TDataItem* item = fItems[i];
		Double_t cutValue = item->GetCutValue();

	   vector<Double_t>::iterator result;
	   result = find( resultVector.begin(), resultVector.end(), cutValue );

	   if( result == resultVector.end() ) {
		   // not found, insert
		   resultVector.push_back(cutValue);
	   }
	}
	return resultVector;
}


//______________________________________________________________________________
vector<Int_t> G4TData::GetSecondaryPDGs()
{
	vector<Int_t> resultVector;
	for(UInt_t i = 0; i < fItems.size(); ++i)
	{
		G4TDataItem* item = fItems[i];
		Int_t pdg = item->GetSecondaryParticlePDG();

	   vector<Int_t>::iterator result;
	   result = find( resultVector.begin(), resultVector.end(), pdg );

	   if( result == resultVector.end() ) {
		   // not found, insert
		   resultVector.push_back(pdg);
	   }
	}
	return resultVector;
}

//______________________________________________________________________________
vector<Double_t> G4TData::GetCutValuesForSecondary(Int_t secondaryPDG)
{
	vector<Double_t> resultVector;
	for(UInt_t i = 0; i < fItems.size(); ++i)
	{
		G4TDataItem* item = fItems[i];
		Int_t pdg = item->GetSecondaryParticlePDG();

		if(pdg == secondaryPDG)
		{
			Double_t cutValue = item->GetCutValue();
			vector<Double_t>::iterator result;
			result = find( resultVector.begin(), resultVector.end(), cutValue );

			if( result == resultVector.end() ) {
			   // not found, insert
			   resultVector.push_back(cutValue);
			}
		}
	}
	return resultVector;
}

//______________________________________________________________________________
vector<G4TDataItem*>	G4TData::GetItemsForSecondary(Int_t secondaryPDG)
{
	vector<G4TDataItem*> resultVector;
	for(UInt_t i = 0; i < fItems.size(); ++i)
	{
		G4TDataItem* item = fItems[i];
		Int_t pdg = item->GetSecondaryParticlePDG();

		if(pdg == secondaryPDG)
		{
			vector<G4TDataItem*>::iterator result;
			result = find( resultVector.begin(), resultVector.end(), item );

			if( result == resultVector.end() ) {
			   // not found, insert
			   resultVector.push_back(item);
			}
		}
	}
	return resultVector;
}

//______________________________________________________________________________
vector<G4TDataItem*>	G4TData::GetItemsForSecondary(vector<Int_t> secondaries)
{
	vector<G4TDataItem*> resultVector;
	for(UInt_t i = 0; i< secondaries.size(); ++i)
	{
		vector<G4TDataItem*> items = GetItemsForSecondary(secondaries[i]);
		for(UInt_t j = 0; j < items.size(); ++j)
		{
			resultVector.push_back(items[j]);
		}
	}
	return resultVector;
}

//______________________________________________________________________________
G4TDataItem* G4TData::GetItem(Int_t secondaryPDG, Double_t cutValue)
{
	G4TDataItem* result = 0;
	for(UInt_t i = 0; i < fItems.size(); ++i)
	{
		G4TDataItem* item = fItems[i];
		Int_t pdg = item->GetSecondaryParticlePDG();
		Double_t cut = item->GetCutValue();

		if(pdg == secondaryPDG && cut == cutValue)
		{
			result = item;
			break;
		}
	}
	return result;
}


//______________________________________________________________________________
DataItemLimit_t	G4TData::GetLimits(Int_t secondaryIdx, Int_t padsPerRow)
{
	DataItemLimit_t rowLimits;
	vector<Int_t> fragments = this->GetSecondaryPDGs();
	Bool_t initialized = false;

	// Auto scale feature
	if(secondaryIdx % padsPerRow == 0)
	{
		UInt_t max = secondaryIdx + padsPerRow;
		if(max > fragments.size()) max = fragments.size();

		for(UInt_t k = secondaryIdx; k < max; ++k)
		{
			Int_t pdgCode	= fragments[k];
			vector<G4TDataItem*> citems = this->GetItemsForSecondary(pdgCode);
			for(UInt_t l = 0; l < citems.size(); ++l)
			{
				DataItemLimit_t limits = citems[l]->GetLimits();
				if(!initialized)
				{
					 rowLimits.fMax = limits.fMax;
					 rowLimits.fMin = limits.fMin;
					 initialized = true;
				}

				if(limits.fMax > rowLimits.fMax ) rowLimits.fMax = limits.fMax;
				if(limits.fMin < rowLimits.fMin ) rowLimits.fMin = limits.fMin;
			}
		}
	}

	//cout << " >>> limits: min = " << rowLimits.fMin << " max = " << rowLimits.fMax << endl;

	return rowLimits;
}

//______________________________________________________________________________
DataItemLimit_t	G4TData::GetLimits()
{
	DataItemLimit_t rowLimits;
	vector<Int_t> fragments = this->GetSecondaryPDGs();
	Bool_t initialized = false;

	// Auto scale feature
	for(UInt_t k = 0; k < fragments.size(); ++k)
	{
		Int_t pdgCode	= fragments[k];
		vector<G4TDataItem*> citems = this->GetItemsForSecondary(pdgCode);
		for(UInt_t l = 0; l < citems.size(); ++l)
		{
			DataItemLimit_t limits = citems[l]->GetLimits();
			if(!initialized)
			{
				 rowLimits.fMax = limits.fMax;
				 rowLimits.fMin = limits.fMin;
				 initialized = true;
			}

			if(limits.fMax > rowLimits.fMax ) rowLimits.fMax = limits.fMax;
			if(limits.fMin < rowLimits.fMin ) rowLimits.fMin = limits.fMin;
		}
	}


	return rowLimits;
}

//______________________________________________________________________________
Bool_t G4TData::IsLoaded() const
{
	  return fIsLoaded;
}

//______________________________________________________________________________
Double_t G4TData::GetCrossSection() const
{
	  return fCrossSection;
}

//______________________________________________________________________________
void G4TData::SetCrossSection(Double_t fCrossSection)
{
	  this->fCrossSection = fCrossSection;
}

//______________________________________________________________________________
Int_t G4TData::GetNumberOfEvents() const
{
	  return fNumberOfEvents;
}

//______________________________________________________________________________
void G4TData::SetNumberOfEvents(Int_t fNumberOfEvents)
{
	  this->fNumberOfEvents = fNumberOfEvents;
}

//______________________________________________________________________________
TString G4TData::GetDirectory() const
{
	return fDirectory;
}

//______________________________________________________________________________
void G4TData::SetDirectory(TString fDirectory)
{
	this->fDirectory = fDirectory;
}

