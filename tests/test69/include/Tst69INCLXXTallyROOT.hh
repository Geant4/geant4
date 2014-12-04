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
#ifndef Tst69INCLXXTallyROOT_hh
#define Tst69INCLXXTallyROOT_hh

#include "G4INCLXXVInterfaceTally.hh"
#include "Rtypes.h"

#include "TString.h"
#include "TTree.h"
#include "TFile.h"

class Tst69INCLXXTallyROOT : public G4INCLXXVInterfaceTally {
  public:
    Tst69INCLXXTallyROOT(TString const &physList);
    virtual ~Tst69INCLXXTallyROOT();
    virtual void Open();
    virtual void Close();
    virtual void Tally(G4HadProjectile const &aTrack, G4Nucleus const &theNucleus, G4HadFinalState const &result);

  protected:
    TString filename;
    TFile *file;
    TTree *reac;

    static const size_t maxNParticles = 1000;

    Short_t Ap;
    Short_t Zp;
    Float_t Ep;
    Short_t At;
    Short_t Zt;
    Short_t nParticles;
    Short_t A[maxNParticles];
    Short_t Z[maxNParticles];
    Float_t EKin[maxNParticles];
    Float_t px[maxNParticles];
    Float_t py[maxNParticles];
    Float_t pz[maxNParticles];
    Float_t theta[maxNParticles];
    Float_t phi[maxNParticles];
};

#endif
