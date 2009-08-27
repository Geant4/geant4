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
// Class G4TTool
//
// Class description:
//
// An abstract (base) class for the tools (simulation and analysis).
//
// History:
// Roman Atachiants, 18/08/2009 - initial version
//
// --------------------------------------------------------------------


#ifndef G4TTool_H_
#define G4TTool_H_

#include "CommonHeaders.h"

#include "G4TModelParams.h"
#include "Helpers/G4TSimHelper.h"
#include "Helpers/G4TPlotHelper.h"
#include "Database/G4TData.h"
#include "Database/G4TDataItem.h"
#include "Database/G4TDataBase.h"

using namespace std;
using namespace ROOT;
using namespace TMath;


class G4TTool : public TObject {

  protected:

		Double_t fHxmin;
		Double_t fHxmax;
		Double_t fHymin;
		Double_t fHymax;
		Double_t fHfbin;
		Double_t fHnbin;
		Double_t fDangle;
		Double_t fPi;
		Double_t fDanrad;
		Double_t fHlxmin;
		Double_t fHlxmax;

		TH1F* fHZL;
		TH1F* fHDT;

		G4TData*		fPublication;

		virtual void Initialize();
		virtual void PrepareHistograms(Double_t hnbin, Double_t hlxmin, Double_t hlxmax);
		virtual void RenderHSolid(TH1F* hist, Int_t hf, Int_t hn, Double_t m, Color_t color = 2, Bool_t noDots = false);


  public:

		G4TTool() {

		}
		virtual ~G4TTool () {}


		ClassDef(G4TTool, 1)  //Base class for tools
};

#endif



