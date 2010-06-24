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
// and a vector of G4TDataItem's (one per secondary/"angle").
//
// History:
// Created by Roman Atachiants, 18/08/2009
// Modified:
// Mikhail Kosov, 11.05.2010: Transfer the kin data to each data item
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
  Int_t     fProjectilePDG; // PDG code for a projectile
  Int_t     fTargetPDG;     // PDG code for a target
  Bool_t    fIsPublication; // if false, then ModelName is required
  TString   fModelName;     // required only if isPublication is false
  TString   fComment;       // Comment for Publication and Postfix for Simulation
  ArgEnum   fTypeVar;       // e.g. Energy in p90
  Double_t  fTypeValue;     // e.g. 90 in p90
  UnitsEnum fTypeUnits;     // e.g. MeV in p90
  Double_t  fSigValue;      // total interaction cross-section
  SigmaEnum fSigUnits;      // units of the interaction cross-section
};

class G4TData : public TObject
{
  protected:
  // Body
  Color_t              fRenderColor;
  vector<G4TDataItem*> fItems;
  Double_t             fCrossSection;
  Int_t                fNumberOfEvents;
  Bool_t               fIsLoaded;
  TString              fDirectory;
  Double_t             fPrMomentum; // projectile momentum (in MeV/c)
  Double_t             fPrEnergy;   // projectile Energy (in MeV)
  Double_t             fTgMass;     // target mass (in MeV/c2)

  // Methods
  TString              HeaderToString(DataObjectHeader_t header) const;
  DataObjectHeader_t   StringToHeader(TString headerStr) const;
  void                 CalcKinValues();

  public:
  // Header
  DataObjectHeader_t   fHeader;   // The header of the object

  G4TData(TString const& headerString)
  {
    fHeader = StringToHeader(headerString);
    fDirectory = "./";
    fIsLoaded = false;
    if(fHeader.fIsPublication) fHeader.fModelName = "data";
    CalcKinValues();
  }
  G4TData( Int_t projectilePDG,
           Int_t targetPDG,
           Bool_t isPublication,
           TString const& modelName,
           TString const& postComment,
           ArgEnum typeVar,
           Double_t typeValue,
           UnitsEnum typeUnits,
           Double_t sigValue,
           SigmaEnum sigUnits,
           Color_t color)
  {
    fHeader.fProjectilePDG = projectilePDG;
    fHeader.fTargetPDG     = targetPDG;
    fHeader.fIsPublication = isPublication;
    fHeader.fModelName     = modelName;
    fHeader.fComment       = postComment;
    fHeader.fTypeVar       = typeVar;
    fHeader.fTypeValue     = typeValue;
    fHeader.fTypeUnits     = typeUnits;
    fHeader.fSigValue      = sigValue;
    fHeader.fSigUnits      = sigUnits;

    fDirectory = "./";
    fIsLoaded = false;
    if(fHeader.fIsPublication) fHeader.fModelName = "data";

    fRenderColor = color;
    CalcKinValues();
  }

  virtual ~G4TData () {}

  void   Save();
  void   Load(Int_t secondaryPDG = 0);
  void   BookHistograms( Int_t    hnbin,
                         Double_t hlxmin,
                         Double_t hlxmax,
                         Int_t particleIdx = 0/* 0 for ALL */,
                         Int_t additionalIndex = -1);
  vector<G4TDataItem*>          GetItems() {return fItems;}
  vector<Double_t>              GetCutValues();
  vector<Int_t>                 GetSecondaryPDGs();
  vector<Double_t>              GetCutValuesForSecondary(Int_t secondaryPDG);
  vector<G4TDataItem*>          GetItemsForSecondary(Int_t secondaryPDG);
  vector<G4TDataItem*>          GetItemsForSecondary(vector<Int_t> secondaries);
  G4TDataItem*                  GetItem(Int_t secondaryPDG, Double_t cutValue);
  TString                       GetModelName() const {return fHeader.fModelName;}
  TString                       GetComment() const {return fHeader.fComment;}
  TString                       GetHeader() const {return HeaderToString(fHeader);}
  Double_t                      GetMaxT(Int_t secondaryIdx, Int_t padsPerRow);
  Double_t                      GetMaxT();
  std::pair<Double_t, Double_t> GetLimits(Int_t secondaryIdx, Int_t padsPerRow);
  std::pair<Double_t, Double_t> GetLimits();
  Bool_t                        IsLoaded() const {return fIsLoaded;}

  // Selectors & Modifiers
  Color_t       GetRenderColor() const              {return fRenderColor;}
  void          SetRenderColor(Color_t fRColor)     {fRenderColor = fRColor;}
  TString       GetDirectory() const                {return fDirectory;}
  void          SetDirectory(TString fDir)          {fDirectory = fDir;}
  Double_t      GetCrossSection() const             {return fCrossSection;}
  void          SetCrossSection(Double_t fCS)       {fCrossSection = fCS;}
  Double_t      GetProjEnergy()                     {return fPrEnergy;}
  Double_t      GetProjMomentum()                   {return fPrMomentum;}
  Double_t      GetTargMass()                       {return fTgMass;}
  Int_t         GetNumberOfEvents() const           {return fNumberOfEvents;}
  void          SetNumberOfEvents(Int_t fNOfEvents) {fNumberOfEvents = fNOfEvents;}

  G4TDataItem*  AddItem( Int_t     SecondaryParticlePDG,
                         ArgEnum   CutVar,
                         UnitsEnum CutUnits,
                         Double_t  CutValue,
                         Double_t  CutDelta,
                         FunctEnum FunctionVar,
                         FunUnEnum FunctionUnits,
				 ErrorType FunctErType,
                         ArgEnum   ArgumentVar,
                         UnitsEnum ArgumentUnits );
  G4TDataItem*  AddItem( TString const& headerStr);

  ClassDef(G4TData, 1)  //The class for Geant4 Model Data handling
};



#endif




