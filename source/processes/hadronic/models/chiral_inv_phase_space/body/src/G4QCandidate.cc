// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QCandidate.cc,v 1.3 2000-08-16 07:33:00 hpw Exp $
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

G4QCandidate::G4QCandidate(const G4QCandidate &right)
{
  Set4Momentum         (right.Get4Momentum());
  SetQPDG              (right.GetQPDG());
  SetQC                (right.GetQC());
  SetNFragments        (right.GetNFragments());
  possible            = right.possible;
  kMin                = right.kMin;
  preProbability      = right.preProbability;
  relativeProbability = right.relativeProbability;
  integralProbability = right.integralProbability;
  secondRelProbability= right.secondRelProbability;
  secondIntProbability= right.secondIntProbability;
}

G4QCandidate::~G4QCandidate() {thePClusters.clearAndDestroy();}

const G4QCandidate& G4QCandidate::operator=(const G4QCandidate &right)
{
  Set4Momentum         (right.Get4Momentum());
  SetQPDG              (right.GetQPDG());
  SetQC                (right.GetQC());
  SetNFragments        (right.GetNFragments());
  possible            = right.possible;
  kMin                = right.kMin;
  preProbability      = right.preProbability;
  relativeProbability = right.relativeProbability;
  integralProbability = right.integralProbability;
  secondRelProbability= right.secondRelProbability;
  secondIntProbability= right.secondIntProbability;

  return *this;
}




