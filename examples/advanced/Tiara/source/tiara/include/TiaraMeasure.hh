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
// $Id: TiaraMeasure.hh,v 1.3 2003/06/25 09:12:47 gunter Exp $
// GEANT4 tag $Name: geant4-05-02-patch-01 $
//
// ----------------------------------------------------------------------
//
// Class TiaraMeasure
//

#ifndef TiaraMeasure_hh
#define TiaraMeasure_hh TiaraMeasure_hh

#include "globals.hh"

class TiaraMeasure {
public:
  TiaraMeasure();
  ~TiaraMeasure();
  void Xin(G4double v);
  G4int GetEntries() const;
  G4double GetMean() const;
  G4double GetVariance() const;
  G4double GetSum() const;
  G4double GetSumSquared() const;
private:
  G4int fEntries;
  G4double fSum;
  G4double fSumSquared;
};

#endif

