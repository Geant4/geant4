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
// Roman Atachiants, 18/08/2009: version for 1 hardwired data set
// Modified:
// Mikhail Kosov, 11/05/2010: generalized for universal DB
//
// --------------------------------------------------------------------

#ifndef G4TDataItem_H_
#define G4TDataItem_H_

#include "../CommonHeaders.h"
#include "G4TParticlesDAL.h"

// New items in the end   -1/9-   -2/10-   -3/11-   -4/12-    -5/13-     -6/14-   -7-
enum FunctEnum{dN_dEdO, dS_dEdO, dN_dpdO, dS_dpdO, dN_pdEdO, dS_pdEdO, EdN_d3p, EdS_d3p,
               EdS_Ad3p, dN_dydpT, dS_dydpT, dN_dxdpT, dS_dxdpT, dN_dxdET, dS_dxdET,
               dN_dydET, dS_dydET, dN_detadET, dS_detadET, dS_dt, dS_du};            // 21
//              ^15^       ^16^      ^17^        ^18^       ^19^   ^20^
//              -0-   -1-   -2-   -3- -4-  -5-  -6-   -7-    -8-        -9-
enum ArgEnum  {T_kin,p_mom,E_tot, p_L,p_T, x_F, y_R, Theta,CosTheta,LogTgHalfTheta}; // 10
//               -0-       -1-      -2-        -3-      -4-    -5-     -6-     -7-
enum FunUnEnum{_MeV_sr, mb_MeV_sr,_MeV2_sr, mb_MeV2_sr, _MeV, mb_MeV, _MeV2, mb_MeV2,
               _GeV_sr, mb_GeV_sr,_GeV2_sr, mb_GeV2_sr, _GeV, mb_GeV, _GeV2, mb_GeV2,
               _1, mb_1};   //      -10-      -11-      -12-  -13-     -14-   -15-   // 18
//             -16/8--17/9-
// !!! Be consistent with G4TData::HeaderToString() and ::StringToHeader() !!!
//             -0-  -1-  -2-  -3-   -4-      -5-     -6-
enum UnitsEnum{MeV, GeV, TeV, keV, Degrees, NoDim, Radians};                         //  7
//              -0-      -1-        -2-       -3-
enum ErrorType{NoError, Absolute, Relative, InPerCent};                              //  4
//              -0-    -1-   -2-
enum SigmaEnum{mBarn, Barn, mkBarn};                                                 //  3


struct DataItemObjectHeader_t
{
  Int_t     fSecondaryParticlePDG; // the PDG code for the secondary particle
  ArgEnum   fCutVar;               // the cut (angle) variable
  UnitsEnum fCutUnits;             // the units for the cut
  Double_t  fCutValue;             // the value for the cut variable (in defined units)
  Double_t  fCutDelta;             // the delta value for the cut (in defined units)
  FunctEnum fFunctionVar;          // the function (CS) used in the data table
  FunUnEnum fFunctionUnits;        // the units of data function
  ErrorType fFunctErrorType;       // The type of errors of the function
  ArgEnum   fArgumentVar;          // the argument type of the data
  UnitsEnum fArgumentUnits;        // the argument units of the data
};

class G4TDataItem : public TObject
{
  protected:
  // Body
  TH1F*  fHistogram;    // Histograms
  TTree* fData;         // Data from ASCII file
  //Information from the top level Header for kinematics calculations
  Double_t fPrMomentum; // projectile momentum (in MeV/c)
  Double_t fPrEnergy;   // projectile Energy (in MeV)
  Double_t fTgMass;     // target mass (in MeV/c2)

  // Methods
  TString      HeaderToString(DataItemObjectHeader_t header) const;
  DataItemObjectHeader_t  StringToHeader(TString headerStr) const;

  public:
  // Header
  DataItemObjectHeader_t fHeader;

  G4TDataItem(TString const& headerStr, Double_t PrMomentum, Double_t PrEnergy,
                      Double_t TgMass);

  G4TDataItem(Int_t     SecondaryParticlePDG,
              ArgEnum   CutVar,
              UnitsEnum CutUnits,
              Double_t  CutValue,
              Double_t  CutDelta,
              FunctEnum FunctionVar,
              FunUnEnum FunctionUnits,
              ErrorType FunctErrorType,
              ArgEnum   ArgumentVar,
              UnitsEnum ArgumentUnits,
              Double_t  PrMomentum,
              Double_t  PrEnergy,
              Double_t  TgMass);

  virtual ~G4TDataItem () {}

  void     LoadFromASCII(const TString& file, Bool_t load3 = false);

  Int_t           GetSecondaryParticlePDG() const {return fHeader.fSecondaryParticlePDG;}
  TString         GetHeader() const               {return HeaderToString(fHeader);}
  ArgEnum         GetCutVar() const               {return fHeader.fCutVar;}
  UnitsEnum       GetCutUnits() const             {return fHeader.fCutUnits;}
  Double_t        GetCutValue() const             {return fHeader.fCutValue;}
  Double_t        GetCutDelta() const             {return fHeader.fCutDelta;}
  FunctEnum       GetFunctionVar() const          {return fHeader.fFunctionVar;}
  FunUnEnum       GetFunctionUnits() const        {return fHeader.fFunctionUnits;}
  ErrorType       GetFunctErrorType() const       {return fHeader.fFunctErrorType;}
  ArgEnum         GetArgumentVar() const          {return fHeader.fArgumentVar;}
  UnitsEnum       GetArgumentUnits() const        {return fHeader.fArgumentUnits;}
  Double_t        GetProjMomentMeV() const        {return fPrMomentum;}
  Double_t        GetProjEnergyMeV() const        {return fPrEnergy;}
  Double_t        GetTargMassMeV() const          {return fTgMass;}
  TH1F*           GetHistogram() const            {return fHistogram;}
  TTree*          GetData() const                 {return fData;}

  void            SetHistogram(TH1F *fHist)       {fHistogram = fHist;}
  void            SetData(TTree *fDat)            {fData = fDat;}

  Double_t                      GetAngleInRadians() const; // Convert whatever to radians
  TString                       GetInvariantFunction();
  Double_t                      GetMaxT();
  std::pair<Double_t, Double_t> GetLimits();
  TString                       GetHistogramName(Int_t anIndex = -1) const;

  ClassDef(G4TDataItem, 1)
};

#endif
