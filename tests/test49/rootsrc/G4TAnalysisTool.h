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
// Class G4TAnalysisTool
//
// Class description:
//
// The analysis tool, takes a publication, searches for the simulation
// files and renders the model comparison plot.
//
// History:
// Roman Atachiants, 18/08/2009 - initial version
//
// --------------------------------------------------------------------

#ifndef G4TAnalysisTool_H_
#define G4TAnalysisTool_H_

#include "CommonHeaders.h"
#include "G4TTool.h"
#include "G4TModelParams.h"
#include "Helpers/G4TPlotHelper.h"
#include "Helpers/G4TSimHelper.h"
#include "Database/G4TData.h"
#include "Database/G4TDataBase.h"

class G4TAnalysisTool : public G4TTool {

  protected:
	  Int_t 					fNP; // Number of Protons
	  Int_t 					fNN; // Number of Neutrons
	  Double_t 					fKineticEnergy;

	  Int_t						fTargetPDG;
	  Int_t						fProjectilePDG;
	  Int_t						fSecondaryPDG;
	  TString					fModelName;
	  Int_t						fCrossSection;

	  vector<G4TData*>			fSimulations;



	  void			SetSecondaryToAnalyze(Int_t idx, Int_t np);
	  void 			InternalExecute(Int_t secondaryPDGorIdx = 0, Int_t np = 13, Int_t nn = 14, Int_t e = 90, const TString& pq = "32-rb20-hp",
			  Int_t nzone = 2, Int_t nvex = 26, const TString& dir = "./"  );

  public:

	  virtual ~G4TAnalysisTool() {}
	  G4TAnalysisTool() { }

	  int Run(TString const& publicationFile, Int_t secondaryPDGorIdx = 0, TString const& printQue = "32-rb20-hp" );

	  ClassDef(G4TAnalysisTool, 1)  //The class for Geant4 Analyzing (p90.kumac base)
};

R__EXTERN G4TAnalysisTool *gAnalysisTool;

#endif




