// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QCandidate.cc,v 1.14 2001-10-04 20:00:22 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -----------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QCandidate ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Quasmon initiated Candidates used by CHIPS Model
// ------------------------------------------------------------------

//#define debug

#include "G4QCandidate.hh"
#include "g4std/algorithm"

G4QCandidate::G4QCandidate() : relativeProbability(0.),integralProbability(0.) {};

G4QCandidate::G4QCandidate(G4int PDGcode) :
  relativeProbability(0.),integralProbability(0.)
{
  G4LorentzVector cur4Mom(0.,0.,0.,0.);
  G4QPDGCode QPDG(PDGcode);
  SetQPDG(QPDG);
  G4double vacMass  = QPDG.GetMass();
  cur4Mom.setE(vacMass);
  Set4Momentum(cur4Mom);
  SetQC(QPDG.GetQuarkContent());
}

G4QCandidate::G4QCandidate(const G4QCandidate& right)
{
  Set4Momentum         (right.Get4Momentum());
  SetQPDG              (right.GetQPDG());
  SetQC                (right.GetQC());
  SetNFragments        (right.GetNFragments());
  possible            = right.possible;
  parPossible         = right.parPossible;
  kMin                = right.kMin;
  denseProbability    = right.denseProbability;
  preProbability      = right.preProbability;
  relativeProbability = right.relativeProbability;
  integralProbability = right.integralProbability;
  secondRelProbability= right.secondRelProbability;
  secondIntProbability= right.secondIntProbability;
  // thePClusters
  G4int nParCl        = right.thePClusters.size();
  if(nParCl) for(G4int ip=0; ip<nParCl; ip++)
  {
    G4QParentCluster* curPC = new G4QParentCluster(right.thePClusters[ip]);
    thePClusters.push_back(curPC);
  }
}

G4QCandidate::G4QCandidate(G4QCandidate* right)
{
  Set4Momentum         (right->Get4Momentum());
  SetQPDG              (right->GetQPDG());
  SetQC                (right->GetQC());
  SetNFragments        (right->GetNFragments());
  possible            = right->possible;
  parPossible         = right->parPossible;
  kMin                = right->kMin;
  denseProbability    = right->denseProbability;
  preProbability      = right->preProbability;
  relativeProbability = right->relativeProbability;
  integralProbability = right->integralProbability;
  secondRelProbability= right->secondRelProbability;
  secondIntProbability= right->secondIntProbability;
  // thePClusters
  G4int nParCl        = right->thePClusters.size();
  if(nParCl) for(G4int ip=0; ip<nParCl; ip++)
  {
    G4QParentCluster* curPC = new G4QParentCluster(right->thePClusters[ip]);
    thePClusters.push_back(curPC);
  }
}

G4QCandidate::~G4QCandidate()
{
#ifdef debug
  G4cout<<"~G4QCandidate: before thePClusters nC="<<thePClusters.entries()<<G4endl;
#endif
  G4std::for_each(thePClusters.begin(), thePClusters.end(), DeleteQParentCluster());
#ifdef debug
  G4cout<<"~G4QCandidate: === DONE ==="<<G4endl;
#endif
}

// Assignment operator
const G4QCandidate& G4QCandidate::operator=(const G4QCandidate &right)
{
  Set4Momentum         (right.Get4Momentum());
  SetQPDG              (right.GetQPDG());
  SetQC                (right.GetQC());
  SetNFragments        (right.GetNFragments());
  possible            = right.possible;
  parPossible         = right.parPossible;
  kMin                = right.kMin;
  denseProbability    = right.denseProbability;
  preProbability      = right.preProbability;
  relativeProbability = right.relativeProbability;
  integralProbability = right.integralProbability;
  secondRelProbability= right.secondRelProbability;
  secondIntProbability= right.secondIntProbability;
  // thePClusters (Vector)
  G4int nParCl        = right.thePClusters.size();
  if(nParCl) for(G4int ip=0; ip<nParCl; ip++)
  {
    G4QParentCluster* curPC = new G4QParentCluster(right.thePClusters[ip]);
    thePClusters.push_back(curPC);
  }

  return *this;
}




