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
