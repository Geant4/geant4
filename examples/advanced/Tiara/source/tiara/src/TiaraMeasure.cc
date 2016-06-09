//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: TiaraMeasure.cc,v 1.3 2003/06/25 09:13:09 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//

#include "TiaraMeasure.hh"

TiaraMeasure::TiaraMeasure() : 
  fEntries(0),
  fSum(0),
  fSumSquared(0)
{}

TiaraMeasure::~TiaraMeasure()
{}

void TiaraMeasure::Xin(G4double v) {
  fEntries++;
  fSum+=v;
  fSumSquared+=(v*v);  
};

G4int TiaraMeasure::GetEntries() const {
  return fEntries;
}

G4double TiaraMeasure::GetMean() const {
  G4double mean(0);
  if (fEntries>0) {
    mean = fSum / fEntries;
  } else {
    G4cout << "TiaraMeasure::GetMean(): fEntries<=0" << G4endl;
  }
  return mean;
}

G4double TiaraMeasure::GetVariance() const {
  G4double variance(-1);
  G4double mean(GetMean());
  if (fEntries>1) {
    variance = (static_cast<double>(fEntries)/(fEntries - 1)) * 
      (fSumSquared/fEntries - mean*mean);
  } else {
    G4cout << "TiaraMeasure::GetVariance(): fEntries<=0" << G4endl;
  }
  return variance;
  
}

G4double TiaraMeasure::GetSum() const {
  return fSum;
}

G4double TiaraMeasure::GetSumSquared() const {
  return fSumSquared;
}
