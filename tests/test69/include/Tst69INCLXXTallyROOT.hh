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
