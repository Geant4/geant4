// $Id: TiaraMeasure.hh,v 1.2 2003-06-16 17:06:46 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

