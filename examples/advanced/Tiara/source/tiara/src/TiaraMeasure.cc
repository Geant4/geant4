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
// $Id: TiaraMeasure.cc,v 1.4 2006/06/29 15:45:17 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
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
