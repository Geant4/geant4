// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Quasmon.hh,v 1.11 2000-09-25 07:26:21 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4Quasmon_h
#define G4Quasmon_h 1


// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4Quasmon ----------------
//             by Mikhail Kossov, July 1999.
//      class for a Quasmon used by the CHIPS Model
// ------------------------------------------------------------

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
#include "G4QCandidateVector.hh"
//#include "G4QParentClusterVector.hh"

class G4Quasmon 
{
public:
  G4Quasmon(G4QContent qQCont    = G4QContent(0,0,0,0,0,0),
            G4LorentzVector q4M  = G4LorentzVector(0.,0.,0.,0.),
            G4LorentzVector ph4M = G4LorentzVector(0.,0.,0.,0.));// Direct Constructor
  G4Quasmon(const G4Quasmon& right);                             // Copy Quasmon by object
  G4Quasmon(G4Quasmon* right);                                   // Copy Quasmon by pointer

  ~G4Quasmon();

  // Overloaded Operators
  const G4Quasmon& operator=(const G4Quasmon& right);
  int operator==(const G4Quasmon &right) const;
  int operator!=(const G4Quasmon &right) const;

  // Static functions
  static void SetParameters(G4double temperature, G4double ssin2g, G4double etaetap);
  static void SetTemper(G4double temperature);
  static void SetSOverU(G4double ssin2g);
  static void SetEtaSup(G4double etaetap);

  //Selectors
  G4double          GetTemper()    const;
  G4double          GetSOverU()    const;
  G4double          GetEtaSup()    const;
  G4LorentzVector   Get4Momentum() const;
  G4QContent        GetQC()        const;
  G4QPDGCode        GetQPDG()      const;
  G4int             GetStatus()    const;

  //Modifiers
  G4QHadronVector*  Fragment(G4QNucleus& nucEnviron); // Pub-wrapper for HadronizeQuasmon()
  void              ClearOutput();                    // Clear but not destroy the output
  void              InitQuasmon(const G4QContent& qQCont, const G4LorentzVector& q4M);
  void              IncreaseBy(const G4Quasmon* pQuasm); // as operator+= but by pointer
  void              KillQuasmon();                    // Kill Quasmon (status=0)

private:  
  G4QHadronVector    HadronizeQuasmon(G4QNucleus& qEnv); // Return new neuclear environment (!)
  G4double           GetRandomMass(G4int PDGCode, G4double maxM);
  G4double           CoulombBarrier(const G4double& tZ, const G4double& tA, const G4double& cZ,
                                    const G4double& cA);
  void               ModifyInMatterCandidates();
  void               InitCandidateVector(G4int maxMes, G4int maxBar, G4int maxClust);
  void               CalculateNumberOfQPartons(G4double qMass);
  void               PrepareClusters();
  void               PrepareCandidates(G4int j);
  void               CalculateHadronizationProbabilities(G4double kQ, G4double kLS, G4int j);
  void               FillHadronVector(G4QHadron* qHadron);
  G4int              RandomPoisson(G4double meanValue);
  G4double           GetQPartonMomentum(G4double mMinResidual2, G4double mCandidate2);
  G4bool             DecayOutHadron(G4QHadron* qHadron);

// Body
private:
  // Static Parameters
  static G4double    Temperature;     // Quasmon Temperature
  static G4double    SSin2Gluons;     // Percent of ssbar sea in a constituen gluon
  static G4double    EtaEtaprime;     // Part of eta-prime in all etas
  // Hadronic input
  G4LorentzVector    q4Mom;           // 4-momentum of the Quasmon +++++
  G4QContent         valQ;            // Quark Content of Quasmon  +++++
  G4QNucleus         theEnvironment;  // Nuclear (or Vacuum) Environment for the Quasmon
  // Output
  G4int              status;          // Status of Quasmon (0-Done,1-FilledSomething,2-DidNothing)
  G4QHadronVector    theQHadrons;     // Vector of generated secondary hadrons +++++
  // Internal working parameters
  G4QCHIPSWorld      theWorld;        // The CHIPS World
  G4LorentzVector    phot4M;          // Gamma 4-momentum for interaction with a quark-parton
  G4int              nBarClust;       // Maximum barion number of clusters @@ ??
  G4int              nOfQ;            // number of quark-partons just to accelerate +++++
  G4QCandidateVector theQCandidates;  // Vector of possible secondary hadrons (***delete***) @@ ??
  G4double           f2all;           // Ratio of free nucleons to free+dense nucleons @@ ??
  G4double           rEP;             // E+p for the Residual Coloured Quasmon im LS +++++
  G4double           rMo;             // p for the Residual Coloured Quasmon im LS +++++
};

inline int G4Quasmon::operator==(const G4Quasmon &right) const {return this == &right;}
inline int G4Quasmon::operator!=(const G4Quasmon &right) const {return this != &right;}
inline G4double        G4Quasmon::GetTemper()    const {return Temperature;}
inline G4double        G4Quasmon::GetSOverU()    const {return SSin2Gluons;}
inline G4double        G4Quasmon::GetEtaSup()    const {return EtaEtaprime;}
inline G4LorentzVector G4Quasmon::Get4Momentum() const {return q4Mom;}
inline G4QContent      G4Quasmon::GetQC()        const {return valQ;}
inline G4QPDGCode      G4Quasmon::GetQPDG()      const {return G4QPDGCode(valQ);}
inline G4int           G4Quasmon::GetStatus()    const {return status;}
inline void            G4Quasmon::ClearOutput()        {theQHadrons.clearAndDestroy();}
inline G4double        G4Quasmon::GetRandomMass(G4int PDG, G4double maxM)
{
  G4QParticle* part = theWorld.GetQParticle(PDG);
  return G4QHadron(part, maxM).GetMass();
}

inline void G4Quasmon::IncreaseBy(const G4Quasmon* pQuasm)
{
  valQ  += pQuasm->GetQC();
  q4Mom += pQuasm->Get4Momentum();
  status= 3;
}

inline void G4Quasmon::InitQuasmon(const G4QContent& qQCont, const G4LorentzVector& q4M)
{
  valQ  = qQCont;
  q4Mom = q4M;
  status= 3;
}

inline void G4Quasmon::KillQuasmon()
{
  static const G4QContent zeroQC(0,0,0,0,0,0);
  static const G4LorentzVector nothing(0.,0.,0.,0.);
  valQ  = zeroQC;
  q4Mom = nothing;
  status= 0;
}

#endif





