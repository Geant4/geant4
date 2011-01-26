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
// Class G4TParticlesDAL
//
// Class description:
//
// Data access layer for the particles pdg. This class handles nuclear
// mass calculation, requests to get the particles by name or A and Z.
//
// History:
// Created by Roman Atachiants, 18/08/2009
// Modified:
// Mikhail Kosov, 10/05/2010: warnings & mesons are added for p400
//
// --------------------------------------------------------------------

#ifndef G4TParticlesDAL_H_
#define G4TParticlesDAL_H_

#include "../CommonHeaders.h"

struct MTableEntry_t
{
  Int_t fZ;
  TString fName;

  MTableEntry_t() {}
  MTableEntry_t(Int_t z, TString const& name) : fZ(z), fName(name) {}
};


class G4TParticlesDAL : public TObject
{
  private:

    vector<MTableEntry_t> fMendeleevTable; // Mendeleev table to get name-> or Z->name
    Bool_t    fNuclearMassesRead; // A flag for the nuclear masses database has been read
    TGeoManager*   fGeoManager;   // Pointer to the geometry manager
    TGeoElementTable*  fTable;    // Pointer to the database of elements

    Double_t   ComputeNuclearMass(Double_t AtomicMass, Double_t A, Int_t Z);

  public:

    G4TParticlesDAL() : fNuclearMassesRead(false), fGeoManager(0), fTable(0) { }
    virtual ~G4TParticlesDAL () {}

    void     ReadNuclearMasses();
   TGeoElementRN* GetParticle(Int_t PDG);
   Double_t   GetParticleMass(Int_t PDG);
   TString   GetParticleName(Int_t PDG, Bool_t inLatex = false);
   TString   GetElementName(Int_t Z);
   TString   GetFileName(Int_t PDG);
   Int_t    GetPDG(TString const& particleName);
   Int_t    GetPDG(Int_t Z, Int_t A, Int_t S=0);
   Int_t    GetA(Int_t PDG);
   Int_t    GetA(TString const& elementName);
   Int_t    GetZ(Int_t PDG);
   Int_t    GetZ(TString const& elementName);
   Int_t    GetN(Int_t PDG);


   TString   GetCut(Int_t PDG);


   ClassDef(G4TParticlesDAL, 1)  //The class for Geant4 Testing Database DAL
};


R__EXTERN G4TParticlesDAL *gParticlesDAL;

#endif




