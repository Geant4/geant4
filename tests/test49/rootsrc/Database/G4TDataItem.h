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
// Class G4TDataItem
//
// Class description:
//
// A class to represent a sub-item of a publication or simulation.
// Contains and saves a TTree (data). Also used to hold the histogram
// (for plotting), but the histogram itself is not saved.
//
// History:
// Roman Atachiants, 18/08/2009 - initial version
//
// --------------------------------------------------------------------

#ifndef G4TDataItem_H_
#define G4TDataItem_H_

#include "../CommonHeaders.h"
#include "G4TParticlesDAL.h"


enum FuncEnum 	{ dS_over_dE, dN_over_dE, dO_over_dE, dN_over_dEdO, dO_over_dEdO};
enum CutEnum	{ Theta, CosTheta, LogTgThetaOverTwo, PL, PT };
enum UnitsEnum	{ Degrees, Radians, MeV, MeVoverC };
enum ArgEnum 	{ E_Kin, Momentum, Energy };


struct DataItemObjectHeader_t
{
	Int_t		fSecondaryParticlePDG; // the PDG code for the secondary particle
	CutEnum 	fCutVar;			   // the cut (angle) variable
	UnitsEnum	fCutUnits;			   // the units for the cut
	Double_t	fCutValue;			   // the value for the cut variable (using determined units)
	Double_t	fCutDelta;			   // the delta value for the cut
	FuncEnum	fFunctionVar;		   // the function used for the scaling/transformation
	UnitsEnum	fFunctionUnits;		   // the units of transformation function
	ArgEnum		fArgumentVar;		   // the argument type
	UnitsEnum	fArgumentUnits;		   // the argument units
};

struct DataItemLimit_t
{
	Double_t fMin;
	Double_t fMax;

	DataItemLimit_t() : fMin(0), fMax(0) {}
};


class G4TDataItem : public TObject {

  protected:
	  // Data
	  TH1F*		fHistogram;
	  TTree*	fData;

	  // Methods
	  TString						HeaderToString(DataItemObjectHeader_t header) const;
	  DataItemObjectHeader_t		StringToHeader(TString headerStr) const;

  public:
	  // Header
	  DataItemObjectHeader_t fHeader;

	  G4TDataItem(TString const& headerStr) {
		  fHeader = StringToHeader(headerStr);
	  }
	  G4TDataItem(
			  Int_t 	SecondaryParticlePDG,
			  CutEnum 	CutVar,
			  UnitsEnum CutUnits,
			  Double_t	CutValue,
			  Double_t	CutDelta,
			  FuncEnum	FunctionVar,
			  UnitsEnum	FunctionUnits,
			  ArgEnum	ArgumentVar,
			  UnitsEnum	ArgumentUnits
		  ) {
				  fHeader.fSecondaryParticlePDG = SecondaryParticlePDG;
				  fHeader.fCutVar = CutVar;
				  fHeader.fCutUnits = CutUnits;
				  fHeader.fCutValue = CutValue;
				  fHeader.fCutDelta = CutDelta;
				  fHeader.fFunctionVar = FunctionVar;
				  fHeader.fFunctionUnits = FunctionUnits;
				  fHeader.fArgumentVar = ArgumentVar;
				  fHeader.fArgumentUnits = ArgumentUnits;
				  fHistogram = 0;
				  fData = 0;

	  }
	  virtual ~G4TDataItem () {}

	  void 				LoadFromASCII(const TString& file, Bool_t load3 = false);

	  Int_t 			GetSecondaryParticlePDG() const;
	  CutEnum 			GetCutVar() const;
	  UnitsEnum 		GetCutUnits() const;
	  Double_t 			GetCutValue() const;
	  Double_t 			GetAngleInRadians() const;
	  Double_t 			GetCutDelta() const;
	  FuncEnum 			GetFunctionVar() const;
	  UnitsEnum 		GetFunctionUnits() const;
	  ArgEnum 			GetArgumentVar() const;
	  UnitsEnum 		GetArgumentUnits() const;
	  TH1F*				GetHistogram() const;
	  TTree*			GetData() const;
	  TString			GetFormula();
	  DataItemLimit_t	GetLimits();

	  TString			GetHeader() const;
	  TString			GetHistogramName(Int_t additionalIndex = -1) const;

	  void 		SetHistogram(TH1F *fHistogram);
	  void		SetData(TTree *fData);

	  ClassDef(G4TDataItem, 1)
};



#endif




