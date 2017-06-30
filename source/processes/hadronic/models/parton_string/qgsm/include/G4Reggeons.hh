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
#ifndef G4Reggeons_h
#define G4Reggeons_h 1
//
// $Id: G4Reggeons.hh 96940 2016-05-18 09:44:20Z ribon $
//
#include "G4Proton.hh"
#include "G4Neutron.hh"

#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"

#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZero.hh"
#include "G4KaonZeroShort.hh"
#include "G4KaonZeroLong.hh"
#include "G4Gamma.hh"

#include "Randomize.hh"

class G4Reggeons
{

  public:
  	G4Reggeons(const G4ParticleDefinition * );

        G4double Get_Cprojectile();

        G4double Get_Ctarget();

	~G4Reggeons();

        void SetS(G4double S);
	void CalculateXs();

        G4double Chi_pomeron(G4double Mult, G4double B);
        G4double Chi_reggeon(G4double Mult, G4double B);

        G4double GetTotalX();
        G4double GetTotalXp();        
        G4double GetTotalXr();

        G4double GetElasticX();
        G4double GetPrDiffX();
        G4double GetTrDiffX();
        G4double GetDDiffX();

        G4double GetInelX();
        G4double GetND_X();
        G4double GetNDp_X();
        G4double GetNDr_X();

        void GetProbabilities(G4double B, G4int Mode,
			      G4double & Pint,
			      G4double & Pprd, G4double & Ptrd, G4double & Pdd, 
			      G4double & Pnd, G4double & Pnvr);

	G4int ncPomerons();

  private:         
	enum  { ALL, WITHOUT_R, NON_DIFF };

        G4ParticleDefinition * Target=G4Proton::Proton();

	G4double Alpha_pomeron, Alphaprime_pomeron, Gamma_pomeron,    Rsquare_pomeron, S0_pomeron;
        G4double Alpha_pomeronHard,                 Gamma_pomeronHard;

        G4double Freggeon_Alpha, Freggeon_Alphaprime, Freggeon_Gamma, Freggeon_Rsquare, Freggeon_C, FParity;
        G4double Wreggeon_Alpha, Wreggeon_Alphaprime, Wreggeon_Gamma, Wreggeon_Rsquare, Wreggeon_C, WParity;

	G4double C_pomeron;                              // = pomeron_Cpr * pomeron_Ctr             // Uzhi Oct. 2016
	G4double Cpr_pomeron;                            // shower enhancement in projectile vertex // Uzhi Oct. 2016
	G4double Ctr_pomeron;                            // shower enhancement in target     vertex // Uzhi Oct. 2016

	G4double Sint=0.;

        G4double chiPin;                                 // Pomeron inelastic phase needed for NcutPomeron
                 					 // calculations.
	G4double Xtotal  , XtotalP, XtotalR;
	G4double Xelastic, Xpr_Diff, Xtr_Diff, XDDiff; 
        G4double Xinel   , Xnd, XndP, XndR;
};
#endif
