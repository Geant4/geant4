#include "G4ChipsBaseXS.hh"

G4ChipsBaseXS::G4ChipsBaseXS() {}

G4ChipsBaseXS::~G4ChipsBaseXS() {}

G4double G4ChipsBaseXS::EquLinearFit(G4double X, G4int N, G4double X0, G4double DX, G4double* Y)
{
  if(DX<=0. || N<2)
    {
      G4cerr<<"***G4ChipsBaseXS::EquLinearFit: DX="<<DX<<", N="<<N<<G4endl;
      return Y[0];
    }
  
  G4int    N2=N-2;
  G4double d=(X-X0)/DX;
  G4int         j=static_cast<int>(d);
  if     (j<0)  j=0;
  else if(j>N2) j=N2;
  d-=j; // excess
  G4double yi=Y[j];
  G4double sigma=yi+(Y[j+1]-yi)*d;
  
  return sigma;
}

G4double G4ChipsBaseXS::GetChipsCrossSection(G4double /*mom*/, G4int /*Z*/, G4int /*n*/, G4int /*pdg*/) {return 0;}
