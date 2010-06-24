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
// Class G4TSimHelper
//
// Class description:
//
// This class makes use of G4TModelParams class in order to prepare the
// parameters and launch a simulation.
// The testing program test19 should be located in:
//        $G4INSTALL/bin/$G4SYSTEM/test19
//
// History:
// Created by Roman Atachiants, 18/08/2009
// Modified:
// Mikhail Kosov, 25/05/2010: transfer the model & pPDG name for G4 simulation
//
// ---------------------------------------------------------------------------
#ifndef G4TSIMHELPER_H_
#define G4TSIMHELPER_H_

#include "../CommonHeaders.h"
#include "../G4TModelParams.h"

class G4TSimHelper : public TObject
{
  private:

  TString fEventsNumberFileName;
  TFile* MakeTree();

  public:

  G4TSimHelper() {fEventsNumberFileName = "histnevt.out";}
  virtual ~G4TSimHelper () {}

  //Static Members
  static void LoadLibraries();

  //Public Members
  TString const& GetEventsNumberFileName(){ return fEventsNumberFileName; }
  void SetEventsNumberFileName(const TString& filename) {fEventsNumberFileName = filename;}

  Int_t GetEventsNumber();
  //--------------------------------------- Tkin=90 MeV (mom is in MeV/c) ------------
  TFile* ExecuteTest(Int_t pPDG=2212, Int_t tPDG=90013014,  Double_t mom=421.,
                     Int_t runNumber=25, Int_t nbEvents=50000, const TString& dir="" ,
                     const TString& model="chips");

  ClassDef(G4TSimHelper, 1)  //The class for Geant4 Simulation
};

R__EXTERN G4TSimHelper *gSimHelper;

#endif




