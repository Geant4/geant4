// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Quasmon.hh,v 1.5 2000-08-22 09:05:13 hpw Exp $
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
#include "G4QParticleVector.hh"
#include "G4QChipolino.hh"
#include "G4QHadronVector.hh"
#include "G4QCandidateVector.hh"
#include "G4QParentClusterVector.hh"

class G4Quasmon 
{
public:
  G4Quasmon(const G4QContent projQC, const G4int targPDG, const G4LorentzVector proj4M,
            const G4LorentzVector targ4Mom, G4int nOfParts);
  G4Quasmon(const G4int projPDG, const G4int targPDG, const G4LorentzVector proj4M,
            const G4LorentzVector targ4Momp, G4int nOfParts);
  G4Quasmon(const G4Quasmon &right);                                  // Quasmon duplication

  ~G4Quasmon();

  int operator==(const G4Quasmon &right) const;
  int operator!=(const G4Quasmon &right) const;

  //Selectors
  G4LorentzVector Get4Momentum() const;
  //General
  G4QHadronVector *  Fragment();
  // Static functions
  static void SetParameters(G4double temperature, G4double ssin2g, G4double etaetap);

private:  
  G4QHadronVector    HadronizeQuasmon();
  G4QParticleVector* InitQuasmonEnvironment(G4int nOfParts); // nOfParts<0 kills the CHIPS World
  G4QParticle*       GetPDGParticle(G4int PDGCode);
  G4QParticle*       GetQParticle(G4int QCode);
  G4double           GetRandomMass(G4int PDGCode, G4double maxM);
  G4double           CoulombBarrier(const G4double& tZ, const G4double& tA, const G4double& cZ,
                                    const G4double& cA);
  void               InitQuasmon(const G4QContent projQC,      const G4int targPDG,
                                 const G4LorentzVector proj4M, const G4LorentzVector targ4Mom,
                                 G4int nOfParts);
  void               ModifyInMatterCandidates();
  void               EvaporateResidual();
  void               FillNEnvInVector();
  void               AntyPDG(G4QHadron*);                     // @@ move to G4QHadron !!
  void               InitCandidateVector(G4int maxMes, G4int maxBar, G4int maxClust);
  void               CalculateNumberOfQPartons(G4double qMass);
  void               PrepareClusters();
  void               PrepareCandidates(G4int j);
  void               PrepareInteractionProbabilities(const G4QContent& projQC);
  void               CalculateHadronizationProbabilities(G4double kQ, G4double kLS, G4int j);
  void               FillHadronVector(G4QHadron* qHadron);
  G4int              RandomPoisson(G4double meanValue);
  G4double           GetQPartonMomentum(G4double mMinResidual2, G4double mCandidate2);
  G4bool             DecayOutHadron(G4QHadron* qHadron);
  G4ThreeVector      RndmDir();

// Body
private:
  G4QParticleVector* qWorld;          // Pointer to the Vector of Particles of the CHIPS World 
  // Static Parameters
  static G4double    Temperature;     // Quasmon Temperature
  static G4double    SSin2Gluons;     // Percent of ssbar sea in a constituen gluon
  static G4double    EtaEtaprime;     // Part of eta-prime in all etas
  // Hadronic input
  G4LorentzVector    q4Mom;           // 4-momentum of the Quasmon
  G4QContent         valQ;            // Quark Content of Quasmon
  G4double           addPhoton;       // Additional energy of soft photon
  G4double           momPhoton;       // Additional momentum of virtual/real photon
  // Output hadrons
  G4QHadronVector    theQHadrons;     // Vector of generated secondary hadrons
  // Internal working parameters
  G4int              nBarClust;       // Maximum barion number of clusters (Calc. @ Interaction)
  G4int              qZ;              // a#of "protons" in the Quasmon        (***delete***)
  G4int              qN;              // a#of "neutrons" in the Quasmon
  G4int              qS;              // a#of "lambdas" in the Quasmon
  G4int              nOfQ;            // number of quark-partons just to accelerate @@ ??
  G4QCandidateVector theQCandidates;  // Vector of possible secondary hadrons (***delete***)
  G4QNucleus         theEnvironment;  // Initial Nucleus & later Residual Nuclear Environment
  G4double           f2all;           // Ratio of free nucleons to free+dense nucleons
  G4double           rEP;             // E+p for the Residual Coloured Quasmon im LS
  G4double           rMo;             // p for the Residual Coloured Quasmon im LS
};

inline int G4Quasmon::operator==(const G4Quasmon &right) const {return this == &right;}
inline int G4Quasmon::operator!=(const G4Quasmon &right) const {return this != &right;}
inline G4LorentzVector G4Quasmon::Get4Momentum() const {return q4Mom;}
inline G4QParticle*    G4Quasmon::GetQParticle(G4int QCode) {return (*qWorld)[QCode];}
inline G4QParticle*    G4Quasmon::GetPDGParticle(G4int PDGCode)
                                  {return (*qWorld)[G4QPDGCode(PDGCode).GetQCode()];}
inline G4double        G4Quasmon::GetRandomMass(G4int PDGCode, G4double maxM)
                                  {return G4QHadron(GetPDGParticle(PDGCode),maxM).GetMass();}
inline void            G4Quasmon::AntyPDG(G4QHadron* hadr)
                                  {if(hadr->TestRealNeutral()) hadr->NegPDGCode();}

#endif





