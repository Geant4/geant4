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
// Class G4TSimulationTool
//
// Class description:
//
// The simulation tool, performs a simulation (prepares chipstest.in
// and launches test19 executable), does the cuts and also plots the
// output histogram.
//
// History:
// Roman Atachiants, 18/08/2009 - initial version
//
// --------------------------------------------------------------------

#ifndef G4TSimulationTool_H_
#define G4TSimulationTool_H_

#include "CommonHeaders.h"

#include "G4TTool.h"
#include "G4TModelParams.h"
#include "Helpers/G4TSimHelper.h"
#include "Helpers/G4TPlotHelper.h"
#include "Database/G4TData.h"
#include "Database/G4TDataItem.h"
#include "Database/G4TDataBase.h"


using namespace std;
using namespace ROOT;
using namespace TMath;


class G4TSimulationTool : public G4TTool {

  protected:
	  Int_t					fTargetPDG;
	  Int_t					fProjectilePDG;
	  TString				fModelName;
	  Int_t					fRunsNumber;
	  Int_t					fEventsNumber;
	  Int_t					fCrossSection;
	  Bool_t				fUseExistingData;

	  G4TSimHelper			fSimHelper;
	  G4TData*				fSimulation;

	  void Plot(Int_t np, TTree* inclTree, Int_t e, Int_t cs, const TString& file);

	  void InternalExecute(Int_t np = 6, Int_t nn = 6, Int_t e = 90, const TString& pq = "32-rb20-hp", Int_t nzone = 2, Int_t nvex = 26,
	  		const TString& dir = "./"  );

  public:

	  G4TSimulationTool() {}
	  virtual ~G4TSimulationTool () {}

	  int Run(TString const& publicationFile, Int_t crossSection = 450, const TString& modelName = "preco", Int_t runsNumber = 25 ,
			  Bool_t useExistingData = false, TString const& dbPath = "./database/", TString const& printQue = "32-rb20-hp");

	  ClassDef(G4TSimulationTool, 1)  //The class for Geant4 Models Testing
};

R__EXTERN G4TSimulationTool *gSimulationTool;

#endif



