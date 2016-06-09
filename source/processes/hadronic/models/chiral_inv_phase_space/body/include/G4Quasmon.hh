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
// $Id$
//
//      ---------------- G4Quasmon ----------------
//             by Mikhail Kossov, July 1999.
//      class for a Quasmon used by the CHIPS Model
// ------------------------------------------------------------
// Short description: If partons from the G4QPartonPair are close in
// rapidity, they create Quasmons, but if they are far in the rapidity
// space, they can not interact directly. Say the bottom parton (quark)
// has rapidity 0, and the top parton (antiquark) has rapidity 8, then
// the top quark splits in two by radiating gluon, and each part has
// rapidity 4, then the gluon splits in quark-antiquark pair (rapidity
// 2 each), and then the quark gadiates anothe gluon and reachs rapidity
// 1. Now it can interact with the bottom antiquark, creating a Quasmon
// or a hadron. The intermediate partons is the string ladder.
// ---------------------------------------------------------------------

#ifndef G4Quasmon_h
#define G4Quasmon_h 1

// Standard G4-headers
#include "G4ios.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
//CHIPS-headers
#include "G4QCHIPSWorld.hh"
#include "G4QChipolino.hh"
#include "G4QHadronVector.hh"
#include "G4QNucleus.hh"

class G4Quasmon 
{
public:
  G4Quasmon(G4QContent qQCont    = G4QContent(0,0,0,0,0,0),
            G4LorentzVector q4M  = G4LorentzVector(0.,0.,0.,0.),
            G4LorentzVector ph4M = G4LorentzVector(0.,0.,0.,0.));// Direct Constructor
  G4Quasmon(const G4Quasmon& right);                             // Copy Quasmon by object
  G4Quasmon(G4Quasmon* right);                                   // Copy Quasmon by pointer

  ~G4Quasmon();                                                  // Public Destructor

  // Overloaded Operators
  const G4Quasmon& operator=(const G4Quasmon& right);
  G4bool operator==(const G4Quasmon &right) const;
  G4bool operator!=(const G4Quasmon &right) const;

  // Static functions
  static void SetParameters(G4double temper=180., G4double ssin2g=.3, G4double etaetap=.3);
  static void SetTemper(G4double temperature);
  static void SetSOverU(G4double ssin2g);
  static void SetEtaSup(G4double etaetap);
  static void OpenElectromagneticDecays();
  static void CloseElectromagneticDecays();

  //Selectors
  G4double          GetTemper()       const;
  G4double          GetSOverU()       const;
  G4double          GetEtaSup()       const;
  G4LorentzVector   Get4Momentum()    const;
  G4QContent        GetQC()           const;
  G4QPDGCode        GetQPDG()         const;
  G4int             GetStatus()       const;
  G4int             GetCharge()       const;
  G4int             GetBaryonNumber() const;
  G4int             GetStrangeness()  const;

  //Modifiers
  void Set4Momentum(G4LorentzVector Q4M) {q4Mom=Q4M;} // Set new value for the Quasmon 4mom
  void SetQC(G4QContent QQC)             {valQ=QQC;}  // Set new Quark Cont for the Quasmon
  void Boost(const G4LorentzVector& theBoost);        // Boosts hadron's 4Momentum using 4M
  void Boost(const G4ThreeVector& B){q4Mom.boost(B);} // Boosts 4-Momentum using v/c
  // Public wrapper for HadronizeQuasmon(,)
  G4QHadronVector*  Fragment(G4QNucleus& nucEnviron, G4int nQ = 1);
  G4QHadronVector*  DecayQuasmon();                   // Decay Quasmon if it's Res or Chipo
  G4QHadronVector*  DecayQHadron(G4QHadron* hadron);  // Decay QHadron if it's Res or Chipo
  void              ClearOutput();                    // Clear but not destroy the output
  void              InitQuasmon(const G4QContent& qQCont, const G4LorentzVector& q4M);
  void              IncreaseBy(const G4Quasmon* pQuasm); // as operator+= but by pointer
  void              IncreaseBy(G4QContent& qQCont, const G4LorentzVector& q4M);
  void              ClearQuasmon();                   // Clear Quasmon (status=0)
  void              KillQuasmon();                    // Kill Quasmon (status=0)
  G4int             CalculateNumberOfQPartons(G4double qMass);

private:  
  G4QHadronVector   HadronizeQuasmon(G4QNucleus& qEnv, G4int nQ=1); // + new Neuclear Envir
  G4double          GetRandomMass(G4int PDGCode, G4double maxM);
  void              ModifyInMatterCandidates();
  void              CalculateHadronizationProbabilities(G4double excE, G4double kQ,
                                                        G4LorentzVector k4M, G4bool piF,
                                                        G4bool gaF, G4bool first=false);
  void              FillHadronVector(G4QHadron* qHadron);
  G4int             RandomPoisson(G4double meanValue);
  G4double          GetQPartonMomentum(G4double mMinResidual2, G4double mCandidate2);
  G4bool            CheckGroundState(G4bool corFlag=false); // Forbid correction by default
  void              KillEnvironment(); // Kill Environment (Z,N,S=0,LV=0)

// Body
private:
  // Static Parameters
  static G4double    Temperature;   // Quasmon Temperature
  static G4double    SSin2Gluons;   // Percent of ssbar sea in a constituen gluon
  static G4double    EtaEtaprime;   // Part of eta-prime in all etas
  static G4bool      WeakDecays;    // Flag for opening WeakDecays (notUsed: allAreClosed)
  static G4bool      ElMaDecays;    // Flag for opening ElectroMagDecays (true by default)
  // Hadronic input
  G4LorentzVector    q4Mom;         // 4-momentum of the Quasmon +++++
  G4QContent         valQ;          // Quark Content of the Quasmon  +++++
  G4QNucleus         theEnvironment;// Nuclear (or Vacuum) Environment around the Quasmon
  // Output
  G4int status; // -1-Panic,0-Done,1-FilledSomething,2-DidNothing,3-StopedByCB,4-JustBorn
  G4QHadronVector    theQHadrons;   // Vector of generated secondary hadrons +++++
  // Internal working objects (@@ it is possible to sacrifice some of them in future)
  G4QCHIPSWorld*     theWorld;      // Pointer to the CHIPS World
  G4LorentzVector    phot4M;        // Gamma 4-momentum for interaction with a quark-parton
  G4int              nBarClust;     // Maximum barion number of clusters in the Environment
  G4int              nOfQ;          // a#of quark-partons in theQuasmon (to accelerate)
  G4QCandidateVector theQCandidates;// Vector of possible secondary hadrons/clusters(*del*)
  G4double           f2all;         // Ratio of free nucleons to free+freeInDense nucleons
  G4double           rEP;           // E+p for the Residual Coloured Quasmon im LS +++++
  G4double           rMo;           // p for the Residual Coloured Quasmon im LS +++++
  G4double           totMass;       // Mass of totalCompoundSys: curQuasmon+curEnvironment
  G4int              bEn;           // BaryoNumber of the Limit Active Nuclear Environment
  G4double           mbEn;          // Mass of the LimActNucEnv
  G4LorentzVector    bEn4M;         // 4-momentum of the LimitActiveNuclearEnviron
  G4QContent         bEnQC;         // QuarkContent of the LimitActiveNuclEnv
};

inline G4bool G4Quasmon::operator==(const G4Quasmon &rhs) const {return this == &rhs;}
inline G4bool G4Quasmon::operator!=(const G4Quasmon &rhs) const {return this != &rhs;}
inline G4double        G4Quasmon::GetTemper()    const {return Temperature;}
inline G4double        G4Quasmon::GetSOverU()    const {return SSin2Gluons;}
inline G4double        G4Quasmon::GetEtaSup()    const {return EtaEtaprime;}
inline G4LorentzVector G4Quasmon::Get4Momentum() const {return q4Mom;}
inline G4QContent      G4Quasmon::GetQC()        const {return valQ;}
inline G4int           G4Quasmon::GetCharge()       const{return valQ.GetCharge();}
inline G4int           G4Quasmon::GetBaryonNumber() const{return valQ.GetBaryonNumber();}
inline G4int           G4Quasmon::GetStrangeness()  const{return valQ.GetStrangeness();}
inline G4QPDGCode      G4Quasmon::GetQPDG()      const {return G4QPDGCode(valQ);}
inline G4int           G4Quasmon::GetStatus()    const {return status;}
inline void            G4Quasmon::ClearOutput()        
      {std::for_each(theQHadrons.begin(), theQHadrons.end(), DeleteQHadron());
       theQHadrons.clear();
      }
inline void            G4Quasmon::KillEnvironment()
                           {theEnvironment=G4QNucleus(0,0,0,G4LorentzVector(0.,0.,0.,0.));}
inline G4double        G4Quasmon::GetRandomMass(G4int PDG, G4double maxM)
{
  G4QParticle* part = theWorld->GetQParticle(PDG);
  return G4QHadron(part, maxM).GetMass();
}

inline void G4Quasmon::IncreaseBy(const G4Quasmon* pQuasm)
{
  valQ  += pQuasm->GetQC();
  q4Mom += pQuasm->Get4Momentum();
  status= 3;
}
inline void G4Quasmon::IncreaseBy(G4QContent& qQCont, const G4LorentzVector& q4M)
{
  valQ  += qQCont;
  q4Mom += q4M;
  status= 3;
}

inline void G4Quasmon::InitQuasmon(const G4QContent& qQCont, const G4LorentzVector& q4M)
{
  valQ  = qQCont;
  q4Mom = q4M;
  status= 3;
}

inline void G4Quasmon::ClearQuasmon()
{
  static const G4QContent zeroQC(0,0,0,0,0,0);
  static const G4LorentzVector nothing(0.,0.,0.,0.);
  phot4M= nothing;
  valQ  = zeroQC;
  q4Mom = nothing;
  status= 0;
  std::for_each(theQCandidates.begin(), theQCandidates.end(), DeleteQCandidate());
  theQCandidates.clear();

}

inline void G4Quasmon::KillQuasmon()
{
  ClearQuasmon();
  ClearOutput();
}

#endif





