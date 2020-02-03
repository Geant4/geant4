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
//
// $Id: G4NeutrinoNucleusModel.hh 90228 2015-05-21 08:49:57Z gcosmo $
//
// Geant4 Header : G4NeutrinoNucleusModel
//
// Author : V.Grichine 12.2.19
//  
// Modified:
//
// Class Description
// Default model for muon neutrino-nucleus charge current scattering; 
// Class Description - End

#ifndef G4NeutrinoNucleusModel_h
#define G4NeutrinoNucleusModel_h 1
 
#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4NucleiProperties.hh"
#include "G4LorentzVector.hh"

class G4ParticleDefinition;
class G4PreCompoundModel;
// class G4CascadeInterface;
// class G4BinaryCascade;
// class G4TheoFSGenerator;
// class G4LundStringFragmentation;
// class G4ExcitedStringDecay;
// class G4INCLXXInterface;
class G4Nucleus;
class G4Fragment;
class G4GeneratorPrecompoundInterface;
class G4ExcitationHandler;

class G4NeutrinoNucleusModel : public G4HadronicInteraction
{
public:

  G4NeutrinoNucleusModel(const G4String& name = "neutrino-nucleus");

  virtual ~G4NeutrinoNucleusModel();

  virtual G4bool IsApplicable(const G4HadProjectile & aTrack, 
  			      G4Nucleus & targetNucleus);

  virtual G4HadFinalState * ApplyYourself(const G4HadProjectile & aTrack, 
					  G4Nucleus & targetNucleus)=0;

 //////// fragmentation functions /////////////////////////

  void ClusterDecay( G4LorentzVector & lvX, G4int qX);

  void MesonDecay( G4LorentzVector & lvX, G4int qX);

  void FinalBarion( G4LorentzVector & lvB, G4int qB, G4int pdgB);

  void RecoilDeexcitation( G4Fragment& fragment);

  void FinalMeson( G4LorentzVector & lvM, G4int qM, G4int pdgM);

  void CoherentPion( G4LorentzVector & lvP, G4int pdgP, G4Nucleus & targetNucleus);


  // set/get class fields

  void SetCutEnergy(G4double ec){fCutEnergy=ec;};
  G4double GetCutEnergy(){return fCutEnergy;};

  G4double GetNuEnergy(){return fNuEnergy;};
  G4double GetQtransfer(){return fQtransfer;};
  G4double GetQ2(){return fQ2;};
  G4double GetXsample(){return fXsample;};

  G4int    GetPDGencoding(){return fPDGencoding;};
  G4bool   GetCascade(){return fCascade;};
  G4bool   GetString(){return fString;};

  G4double GetCosTheta(){return fCosTheta;};
  G4double GetEmu(){return fEmu;};
  G4double GetEx(){return fEx;};
  G4double GetMuMass(){return fMu;};
  G4double GetW2(){return fW2;};
  G4double GetM1(){return fM1;};
  G4double GetMr(){return fMr;};
  G4double GetTr(){return fTr;};
  G4double GetDp(){return fDp;};

  G4LorentzVector GetLVl(){return fLVl;};
  G4LorentzVector GetLVh(){return fLVh;};
  G4LorentzVector GetLVt(){return fLVt;};
  G4LorentzVector GetLVcpi(){return fLVcpi;};

  G4double GetMinNuMuEnergy(){ return fMu + 0.5*fMu*fMu/fM1 + 4.*CLHEP::MeV; }; // kinematics + accuracy for sqrts

  G4double ThresholdEnergy(G4double mI, G4double mF, G4double mP) // for cluster decay
  { 
    G4double w = std::sqrt(fW2);
    return w + 0.5*( (mP+mF)*(mP+mF)-(w+mI)*(w+mI) )/mI;
  };
  G4double FinalMomentum(G4double mI, G4double mF, G4double mP, G4LorentzVector lvX); // for cluster decay

  // nucleon binding

  G4double FermiMomentum( G4Nucleus & targetNucleus);
  G4double NucleonMomentum( G4Nucleus & targetNucleus);
  
  G4double GetEx( G4int A, G4bool fP );
  G4double GgSampleNM(G4Nucleus & nucl);
  
  G4int    GetEnergyIndex(G4double energy);
  G4double GetNuMuQeTotRat(G4int index, G4double energy);

  G4int    GetOnePionIndex(G4double energy);
  G4double GetNuMuOnePionProb(G4int index, G4double energy);
  
  virtual void ModelDescription(std::ostream&) const;

protected:

  G4ParticleDefinition* theMuonMinus;
  G4ParticleDefinition* theMuonPlus;
 
  G4double fSin2tW;    // sin^2theta_Weinberg
  G4double fCutEnergy; // minimal recoil electron energy detected

  G4int fNbin, fIndex, fEindex, fXindex, fQindex, fOnePionIndex, fPDGencoding;
  G4bool fCascade, fString, fProton, f2p2h, fBreak;

  G4double fNuEnergy, fQ2, fQtransfer, fXsample;

  G4double fM1, fM2, fMt, fMu, fW2,  fMpi, fW2pi, fMinNuEnergy, fDp, fTr;

  G4double fEmu, fEmuPi, fEx, fMr, fCosTheta, fCosThetaPi; // final lepton

  G4LorentzVector fLVh, fLVl, fLVt, fLVcpi;

  G4GeneratorPrecompoundInterface* fPrecoInterface;
  G4PreCompoundModel*              fPreCompound;
  G4ExcitationHandler*             fDeExcitation;


  G4Nucleus* fRecoil;

  static const G4int fResNumber;
  static const G4double fResMass[6]; // [fResNumber];

  static const G4int fClustNumber;

  static const G4double fMesMass[4];
  static const G4int    fMesPDG[4];

  static const G4double fBarMass[4];
  static const G4int    fBarPDG[4];

  static const G4double fNuMuResQ[50][50];
  

  static const G4double fNuMuEnergy[50];
  static const G4double fNuMuQeTotRat[50];
  static const G4double fOnePionEnergy[58];
  static const G4double fOnePionProb[58];
  
};



#endif
