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

#include "G4DNASamplingTable.hh"
#include "G4EmParameters.hh"
#include "Randomize.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

#include <vector>
#include <fstream>
#include <sstream>


G4DNASamplingTable::G4DNASamplingTable(std::size_t npoints)
{
  fPrimaryEnergy.reserve(npoints);
  fSecEnergy.reserve(npoints);
  for (G4int i=0; i<5; ++i) { (fPDF[i]).reserve(npoints); }
}

G4DNASamplingTable::~G4DNASamplingTable()
{
  for (auto & p : fSecEnergy) { delete p; }
  for (G4int i=0; i<5; ++i) {
    for (auto & p : fPDF[i]) { delete p; }
  }
}

void G4DNASamplingTable::LoadData(const G4String& fname, G4double factE,
				  G4double fact, G4bool verbose)
{
  std::ostringstream ost;
  ost << G4EmParameters::Instance()->GetDirLEDATA() << "/" << fname;
  std::ifstream fin(ost.str().c_str());
  if (!fin.is_open()) {
    G4ExceptionDescription ed;
    ed << "File <" << ost.str().c_str() << "> is not opened!";
    G4Exception("G4DNASamplingTable::LoadDifferential ", "em0003",
		FatalException, ed, "");
    return;
  }

  G4double t, e, sig;
  G4double e0{0.0};
  G4int ntmax{0};
  G4int nt{0};
  std::vector<G4double>* v = nullptr;
  std::vector<G4double>* vPDF[5];
  for (;;) {
    fin >> e;
    if (fin.eof()) { break; }
    if (e != e0 || nullptr == v) {
      fPrimaryEnergy.push_back(e*factE);
      e0 = e;
      ++fNpoints;
      v = new std::vector<G4double>;
      fSecEnergy.push_back(v);
      for (G4int i=0; i<5; ++i) {
        vPDF[i] = new std::vector<G4double>;
        (fPDF[i]).push_back(vPDF[i]);
      }
      ntmax = std::max(ntmax, nt);
      nt = 0;
    }
    fin >> t;
    v->push_back(t*factE);
    ++nt;
    for (G4int i=0; i<5; ++i) {
      fin >> sig;
      sig *= fact;
      (vPDF[i])->push_back(sig);
    }
    if (fin.eof()) { break; }
  }
  if (verbose) {
    G4cout << "G4DNASamplingTable::LoadData from file:" << G4endl;
    G4cout << fname << G4endl;
    G4cout << "    Nenergy= " << fNpoints << " NmaxT= " << ntmax << G4endl;
  }
  if (fNpoints > 0) { --fNpoints; }
}

G4double G4DNASamplingTable::GetValue(G4double ekinPrimary,
				      G4double ekinSec, G4int shell) const
{
  std::vector<G4double>* e1{nullptr};
  std::vector<G4double>* e2{nullptr};
  std::vector<G4double>* s1{nullptr};
  std::vector<G4double>* s2{nullptr};
  G4int idx = GetIndex(fPrimaryEnergy, ekinPrimary);
  if (idx == -1) {
    e1 = fSecEnergy[0];
    s1 = (fPDF[shell])[0];
  } else if (idx > fNpoints) {
    e1 = fSecEnergy[fNpoints];
    s1 = (fPDF[shell])[fNpoints];
  } else {
    e1 = fSecEnergy[idx];
    s1 = (fPDF[shell])[idx];
    e2 = fSecEnergy[idx + 1];
    s2 = (fPDF[shell])[idx + 1];
  }
  // edge cases
  G4double res1 = VecInterpolation(e1, s1, ekinSec);
  if (nullptr == e2) { return res1; }

  // ordinary case
  G4double res2 = VecInterpolation(e2, s2, ekinSec);
  G4double res = Interpolate(fPrimaryEnergy[idx], fPrimaryEnergy[idx + 1],
			     ekinPrimary, res1, res2);
  return res;
}

G4int G4DNASamplingTable::GetIndex(const std::vector<G4double>& v, G4double x) const
{
  G4int idx;
  if (x <= v[0]) { idx = -1; }
  else if (x >= v.back()) { idx = (G4int)v.size(); }
  else {
    std::size_t i = std::upper_bound(v.cbegin(), v.cend(), x) - v.cbegin() - 1;
    idx = (G4int)i;
  }
  return idx;
}

G4double G4DNASamplingTable::VecInterpolation(const std::vector<G4double>* ener,
					      const std::vector<G4double>* val,
					      G4double e) const
{
  G4int idx = GetIndex(*ener, e);
  G4double res;
  if (idx == -1) { res = (*val)[0]; }
  else if (e >= ener->back()) { res = val->back(); }
  else {
    res = Interpolate((*ener)[idx], (*ener)[idx + 1], e, (*val)[idx], (*val)[idx + 1]);
  }
  return res;
}

G4double G4DNASamplingTable::Interpolate(G4double e1, G4double e2, G4double e,
		                         G4double xs1, G4double xs2) const
{
  G4double res;
  // special case
  if (e1 == e2) {
    res = 0.5 * (xs1 + xs2);

    // Log-log interpolation by default
  } else if (e1 > 0.0 && e2 > 0.0 && xs1 > 0.0 && xs2 > 0.0) {
    G4double y = G4Log(xs1) + G4Log(e/e1) * G4Log(xs2/xs1)/G4Log(e2/e1);
    res = G4Exp(y);

    // Lin-Log interpolation
  } else if (xs1 > 0.0 && xs2 > 0.0) {
    G4double y = G4Log(xs1) + (e - e1) * G4Log(xs2/xs1)/(e2 - e1);
    res = G4Exp(y);

    // Lin-Lin interpolation
  } else { 
    res = xs1 + (e - e1) * (xs2 - xs1)/(e2 - e1);
  }
  return res;
}

G4double
G4DNASamplingTable::SampleCumulative(G4double ekinPrimary, G4int shell) const
{
  std::vector<G4double>* e1{nullptr};
  std::vector<G4double>* e2{nullptr};
  std::vector<G4double>* s1{nullptr};
  std::vector<G4double>* s2{nullptr};
  G4int idx = GetIndex(fPrimaryEnergy, ekinPrimary);
  if (idx == -1) {
    e1 = fSecEnergy[0];
    s1 = (fPDF[shell])[0];
  } else if (idx > fNpoints) {
    e1 = fSecEnergy[fNpoints];
    s1 = (fPDF[shell])[fNpoints];
  } else {
    e1 = fSecEnergy[idx];
    s1 = (fPDF[shell])[idx];
    e2 = fSecEnergy[idx + 1];
    s2 = (fPDF[shell])[idx + 1];
  }
  G4double q = G4UniformRand();

  // edge cases
  G4double res1 = VecInterpolation(s1, e1, q);
  if (nullptr == e2) { return res1; }

  // ordinary case
  G4double res2 = VecInterpolation(s2, e2, q);
  G4double res = Interpolate(fPrimaryEnergy[idx], fPrimaryEnergy[idx + 1],
			     ekinPrimary, res1, res2);
  return res;
}
