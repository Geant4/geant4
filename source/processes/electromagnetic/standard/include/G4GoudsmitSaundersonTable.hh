//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software Liscense,  included in the file *
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
// $Id: G4GoudsmitSaundersonTable.hh 103884 2017-05-03 08:04:50Z gcosmo $
//
// -----------------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4GoudsmitSaundersonTable
//
// Author:        Mihaly Novak / (Omrane Kadri)
//
// Creation date: 20.02.2009
//
// Class description:
//   Class to handle multiple scattering angular distributions precomputed by
//   using Kawrakow-Bielajew Goudsmit-Saunderson MSC model based on the screened
//   Rutherford DCS for elastic scattering of electrons/positrons [1,2]. This
//   class is used by G4GoudsmitSaundersonMscModel to sample the angular
//   deflection of electrons/positrons after travelling a given path.
//
// Modifications:
// 04.03.2009 V.Ivanchenko cleanup and format according to Geant4 EM style
// 18.05.2015 M. Novak This class has been completely replaced (only the original
//            class name was kept; class description was also inserted):
//            A new version of Kawrakow-Bielajew Goudsmit-Saunderson MSC model
//            based on the screened Rutherford DCS for elastic scattering of
//            electrons/positrons has been introduced[1,2]. The corresponding MSC
//            angular distributions over a 2D parameter grid have been recomputed
//            and the CDFs are now stored in a variable transformed (smooth) form
//            together with the corresponding rational interpolation parameters.
//            The new version is several times faster, more robust and accurate
//            compared to the earlier version (G4GoudsmitSaundersonMscModel class
//            that use these data has been also completely replaced)
// 28.04.2017 M. Novak: the GS angular distributions has been recomputed, the 
//            data size has been reduced from 16 MB down to 5 MB by using a new 
//            representation, the class has been modified significantly due to 
//            this new data representation.
//
// References:
//   [1] A.F.Bielajew, NIMB, 111 (1996) 195-208
//   [2] I.Kawrakow, A.F.Bielajew, NIMB 134(1998) 325-336
//
// -----------------------------------------------------------------------------


#ifndef G4GoudsmitSaundersonTable_h
#define G4GoudsmitSaundersonTable_h 1

#include <vector>

#include "G4Types.hh"


class G4GoudsmitSaundersonTable {

public:
  G4GoudsmitSaundersonTable();
 ~G4GoudsmitSaundersonTable();

  void Initialise();

  // structure to store one GS transformed angular distribution (for a given s/lambda_el,s/lambda_elG1)
  struct GSMSCAngularDtr {
    G4int     fNumData;    // # of data points
    G4double  fQScale;
    G4double *fUValues;    // array of transformed variables
    G4double *fParamA;     // array of interpolation parameters a
    G4double *fParamB;     // array of interpolation parameters b
  };

  void LoadMSCData();
  GSMSCAngularDtr* GetOne(G4int indx) {return fGSMSCAngularDistributions1[indx];}

  void   Sampling(G4double lambdaval, G4double qval, G4double scra, 
                  G4double &cost, G4double &sint);
  G4double SampleCosTheta(G4double lambdaval, G4double qval, G4double scra, 
                          G4double rndm1, G4double rndm2, G4double rndm);
  G4double SampleCosTheta1(G4double lambdaval, G4double qval, G4double scra, 
                           G4double rndm1, G4double rndm2, G4double rndm);
  G4double SampleCosTheta2(G4double lambdaval, G4double qval, G4double scra,
                           G4double rndm1, G4double rndm2, G4double rndm);
  G4double GetScreeningParam(G4double G1);

  // material dependent MSC parameters (computed at initialisation) regarding
  // Moliere's screening parameter
  G4double GetMoliereBc(G4int matindx){return (*fgMoliereBc)[matindx];}
  G4double GetMoliereXc2(G4int matindx){return (*fgMoliereXc2)[matindx];}

private:
  // initialisation of material dependent Moliere's MSC parameters
  void InitMoliereMSCParams();


 private:
   static bool             gIsInitialised;         // are the precomputed angular distributions already loaded in?
   static constexpr G4int    gLAMBNUM = 64;        // # L=s/lambda_el in [fLAMBMIN,fLAMBMAX]
   static constexpr G4int    gQNUM1   = 15;        // # Q=s/lambda_el G1 in [fQMIN1,fQMAX1] in the 1-st Q grid
   static constexpr G4int    gQNUM2   = 32;        // # Q=s/lambda_el G1 in [fQMIN2,fQMAX2] in the 2-st Q grid
   static constexpr G4int    gNUMSCR1 = 201;       // # of screening parameters in the A(G1) function
   static constexpr G4int    gNUMSCR2 = 51;        // # of screening parameters in the A(G1) function
   static constexpr G4double gLAMBMIN = 1.0;       // minimum s/lambda_el
   static constexpr G4double gLAMBMAX = 100000.0;  // maximum s/lambda_el
   static constexpr G4double gQMIN1   = 0.001;     // minimum s/lambda_el G1 in the 1-st Q grid
   static constexpr G4double gQMAX1   = 0.99;      // maximum s/lambda_el G1 in the 1-st Q grid
   static constexpr G4double gQMIN2   = 0.99;      // minimum s/lambda_el G1 in the 1-st Q grid
   static constexpr G4double gQMAX2   = 7.99;      // maximum s/lambda_el G1 in the 1-st Q grid
   // precomputed A(G1) function with its interpolation parameters
   static constexpr G4double gSCRMIN1  = 1.93214991408357e-12;
   static constexpr G4double gSCRMAX1  = 2.42974344203683e-01;
   static constexpr G4double gSCRMAX2  = 5.50564555556202e+01;
   //
   static const G4double gG1Values1[];
   static const G4double gScrAValues1[];
   static const G4double gScrBValues1[];
   static const G4double gG1Values2[];
   static const G4double gScrAValues2[];
   static const G4double gScrBValues2[];


   G4double fLogLambda0;          // ln(gLAMBMIN)
   G4double fLogDeltaLambda;      // ln(gLAMBMAX/gLAMBMIN)/(gLAMBNUM-1)
   G4double fInvLogDeltaLambda;   // 1/[ln(gLAMBMAX/gLAMBMIN)/(gLAMBNUM-1)]
   G4double fInvDeltaQ1;          // 1/[(gQMAX1-gQMIN1)/(gQNUM1-1)]
   G4double fDeltaQ2;             // [(gQMAX2-gQMIN2)/(gQNUM2-1)]
   G4double fInvDeltaQ2;          // 1/[(gQMAX2-gQMIN2)/(gQNUM2-1)]
   // for the precumputed A(G1) function
   G4double fLogG1FuncMin1;
   G4double fInvLogDeltaG1Func1;
   G4double fLogG1FuncMin2;
   G4double fInvLogDeltaG1Func2;

   // vector to store all GS transformed angular distributions
   std::vector<GSMSCAngularDtr*> fGSMSCAngularDistributions1;
   std::vector<GSMSCAngularDtr*> fGSMSCAngularDistributions2;

   //@{
   /** Precomputed \f$ b_lambda_{c} $\f and \f$ \chi_c^{2} $\f material dependent
   *   Moliere parameters that can be used to compute the screening parameter,
   *   the elastic scattering cross section (or \f$ \lambda_{e} $\f) under the
   *   screened Rutherford cross section approximation. (These are used in
   *   G4GoudsmitSaundersonMscModel if fgIsUsePWATotalXsecData is FALSE.)
   */
   static std::vector<G4double> *fgMoliereBc;
   static std::vector<G4double> *fgMoliereXc2;

};

#endif

