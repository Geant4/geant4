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
// Class G4TData
//
// Class description:
//
// A class to represent a publication or simulation. Contains a header
// and a vector of G4TDataItem.
//
// History:
// Roman Atachiants, 18/08/2009 - initial version
//
// --------------------------------------------------------------------

#ifndef G4TData_H_
#define G4TData_H_

#include "../CommonHeaders.h"
#include "G4TDataItem.h"
#include "G4TParticlesDAL.h"
#include "G4TCatalog.h"

struct DataObjectHeader_t
{
	  Int_t							fProjectilePDG;		// PDG code for a projectile
	  Int_t							fTargetPDG;			// PDG code for a target
	  Bool_t						fIsPublication;		// if false, then ModelName is required
	  TString 						fModelName;			// required only if isPublication is true
	  ArgEnum						fTypeVar;			// Energy in p90*.kumac
	  Double_t						fTypeValue; 		// 90 MeV in p90
	  UnitsEnum						fTypeUnits;			// MeV in p90
};

class G4TData : public TObject {

  protected:

	  // Other fields
	  Color_t						fRenderColor;
	  vector<G4TDataItem*>	fItems;
	  Double_t						fCrossSection;
	  Int_t							fNumberOfEvents;
	  Bool_t						fIsLoaded;
	  TString						fDirectory;

	  // Methods
	  TString						HeaderToString(DataObjectHeader_t header) const;
	  DataObjectHeader_t			StringToHeader(TString headerStr) const;

  public:
	  // Header
	  DataObjectHeader_t			fHeader;			// The header of the object

	  G4TData(TString const& headerString) {
		  fHeader = StringToHeader(headerString);
		  fDirectory = "./";
		  fIsLoaded = false;
		  if(fHeader.fIsPublication) fHeader.fModelName = "data";
	  }
	  G4TData(
			  Int_t projectilePDG,
			  Int_t targetPDG,
			  Bool_t isPublication,
			  TString const& modelName,
			  ArgEnum typeVar,
			  Double_t typeValue,
			  UnitsEnum typeUnits,
			  Color_t color
		  )  {
				  fHeader.fProjectilePDG = projectilePDG;
				  fHeader.fTargetPDG = targetPDG;
				  fHeader.fIsPublication = isPublication;
				  fHeader.fModelName = modelName;
				  fHeader.fTypeVar = typeVar;
				  fHeader.fTypeValue = typeValue;
				  fHeader.fTypeUnits = typeUnits;

				  fDirectory = "./";
				  fIsLoaded = false;
				  if(fHeader.fIsPublication) fHeader.fModelName = "data";

				  fRenderColor = color;
			  }
	  virtual ~G4TData () {}


	  void 		Save();
	  void 		Load(Int_t secondaryPDG = 0);
	  void 		PrepareHistograms(Double_t hnbin, Double_t hlxmin, Double_t hlxmax, Int_t particleIdx = 0/* 0 for ALL */, Int_t additionalIndex = -1 );



	  vector<G4TDataItem*>			GetItems();
	  vector<Double_t>				GetCutValues();
	  vector<Int_t>					GetSecondaryPDGs();
	  vector<Double_t>				GetCutValuesForSecondary(Int_t secondaryPDG);
	  vector<G4TDataItem*>			GetItemsForSecondary(Int_t secondaryPDG);
	  vector<G4TDataItem*>			GetItemsForSecondary(vector<Int_t> secondaries);
	  G4TDataItem*					GetItem(Int_t secondaryPDG, Double_t cutValue);
	  TString 						GetModelName() const;
	  TString						GetHeader() const;
	  DataItemLimit_t				GetLimits(Int_t secondaryIdx, Int_t padsPerRow);
	  DataItemLimit_t				GetLimits();
	  Bool_t 						IsLoaded() const;


	  // Getters/Setters
	  Color_t 						GetRenderColor() const;
	  void 							SetRenderColor(Color_t);
	  TString						GetDirectory() const;
	  void							SetDirectory(TString fDirectory);
	  Double_t 						GetCrossSection() const;
	  void 							SetCrossSection(Double_t fCrossSection);
	  Int_t							GetNumberOfEvents() const;
	  void 							SetNumberOfEvents(Int_t fNumberOfEvents);


	  G4TDataItem*			AddItem( Int_t 	SecondaryParticlePDG,
											 CutEnum 	CutVar,
											 UnitsEnum CutUnits,
											 Double_t	CutValue,
											 Double_t	CutDelta,
											 FuncEnum	FunctionVar,
											 UnitsEnum	FunctionUnits,
											 ArgEnum	ArgumentVar,
											 UnitsEnum	ArgumentUnits
									);
	  G4TDataItem*			AddItem( TString const& headerStr);


	  ClassDef(G4TData, 1)  //The class for Geant4 Model Data handling
};



#endif




