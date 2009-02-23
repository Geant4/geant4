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
// $Id: G4QCandidate.cc,v 1.35 2009-02-23 09:49:24 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QCandidate ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Quasmon initiated Candidates used by CHIPS Model
// ------------------------------------------------------------------
// Short description: A candidate for hadronization. The candidates
// (hadrons or nuclear fragments) are competative, each quark of a
// Quasmon select which candidate to use for hadronization
// ------------------------------------------------------------------

//#define debug

#include "G4QCandidate.hh"
#include <algorithm>

G4QCandidate::G4QCandidate() : 
  G4QHadron(),possible(false),parPossible(false),kMin(0),denseProbability(0.),
  preProbability(0.),relativeProbability(0.),integralProbability(0.),
  secondRelProbability(0.),secondIntProbability(0.),EBMass(0.),NBMass(0.)
                            
{
}

G4QCandidate::G4QCandidate(G4int PDGcode) :
  G4QHadron(PDGcode),possible(false),parPossible(false),kMin(0),denseProbability(0.),
  preProbability(0.),relativeProbability(0.),integralProbability(0.),
  secondRelProbability(0.),secondIntProbability(0.),EBMass(0.),NBMass(0.)
{
#ifdef debug
  G4cout<<"G4QCandidate::Constructor: PDG="<<PDGcode<<G4endl;
#endif
  G4LorentzVector cur4Mom(0.,0.,0.,0.);
  G4QPDGCode QPDG(PDGcode);
#ifdef debug
  G4cout<<"G4QCandidate::Constructor: QPDG="<<QPDG<<G4endl;
#endif
  SetQPDG(QPDG);
  G4double vacMass=QPDG.GetMass();
#ifdef debug
  G4cout<<"G4QCandidate::Constructor: M="<<vacMass<<G4endl;
#endif
  cur4Mom.setE(vacMass);
  Set4Momentum(cur4Mom);
  SetQC(QPDG.GetQuarkContent());
}

G4QCandidate::G4QCandidate(const G4QCandidate& right) : 
  G4QHadron(&right)
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
  EBMass = right.EBMass;
  NBMass = right.NBMass;
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
  EBMass = right->EBMass;
  NBMass = right->NBMass;
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
  G4cout<<"~G4QCandidate: before thePClusters nC="<<thePClusters.size()<<G4endl;
#endif
  std::for_each(thePClusters.begin(), thePClusters.end(), DeleteQParentCluster());
#ifdef debug
  G4cout<<"~G4QCandidate: === DONE ==="<<G4endl;
#endif
}

// Assignment operator
const G4QCandidate& G4QCandidate::operator=(const G4QCandidate &right)
{
  if(this != &right)                          // Beware of self assignment
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
    EBMass = right.EBMass;
    NBMass = right.NBMass;
    // thePClusters (Vector)
    G4int nParCl        = right.thePClusters.size();
    if(nParCl) for(G4int ip=0; ip<nParCl; ip++)
    {
      G4QParentCluster* curPC = new G4QParentCluster(right.thePClusters[ip]);
      thePClusters.push_back(curPC);
    }
  }
  return *this;
}
