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
// Class G4TModelParams
//
// Class description:
//
// TClass that takes care of parsing (loading, saving) the chipstest.in
// file, uses also chipstest_default.in in order to take the defaults and transform them.
//
// History:
// Created by Roman Atachiants, 18/08/2009
// Modified:
// Mikhail Kosov, 25/05/2010: include mode name in the string
//
// --------------------------------------------------------------------

#ifndef G4TMODELPARAMS_H_
#define G4TMODELPARAMS_H_

#include "CommonHeaders.h"

struct ParamsData_t
{
  Float_t temp;
  Float_t ssse;
  Float_t eepr;
  Float_t nop;
  Float_t momb;
  Float_t enb;
  Int_t pdgpr;
  Int_t pdgtg;
  Int_t nevnt;
  Float_t freeN;
  Float_t freeD;
  Float_t clustP;
  Float_t rMed;
  Float_t solA;
  TString mdName;
};

class G4TModelParams : public TObject
{
  private:
  ParamsData_t fData;

  inline void Tokenize(const string& str, vector<string>& tokens,
                       const string& delimiters = " ");
  public:
   // Selectors and Modifiers
   ParamsData_t& GetData(){ return fData; }
   void SetData(ParamsData_t const& data){ fData = data; }

   void Load(const TString& pf = "chipstest.in");
   void Save(const TString& pf = "chipstest.in");

   G4TModelParams()
   {
     fData.temp = 180.00;
     fData.ssse = 0.30;
     fData.eepr = 0.30;
     fData.nop = 152;
     fData.momb = 0;
     fData.enb = 0;
     fData.pdgpr = 2212;
     fData.pdgtg = 90013014;
     fData.nevnt = 2000;
     fData.freeN = 0.40;
     fData.freeD = 0.20;
     fData.clustP = 5.00;
     fData.rMed = 1.00;
     fData.solA = 0.50;
     fData.mdName = "chips";
   }
   virtual ~G4TModelParams () {}

   ClassDef(G4TModelParams, 1)  //The class for Geant4 parameters handling
};

#endif
