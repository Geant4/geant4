#include "G4PiData.hh"
// by J.P Wellisch, Sun Sep 15 2002.

G4PiData::    
G4PiData(const G4double * aT, const G4double * aIn, const G4double * anE, G4int nP)
{
  G4int i=0;
  for(i=0; i<nP; i++)
  {
    G4std::pair<G4double, G4double> x;
    x.first=aT[i]*millibarn;
    x.second=aIn[i]*millibarn;
    G4std::pair<G4double, G4std::pair<G4double, G4double > > aP;
    aP.first=anE[i]*GeV;
    aP.second=x;
    push_back(aP);
  }
}

G4bool G4PiData::
AppliesTo(G4double kineticEnergy)
{
  G4bool result = true;
  if(kineticEnergy>back().first) result = false;
  return result;
}

G4double G4PiData::
ReactionXSection(G4double kineticEnergy)
{
  G4double result = 0;
  G4PiData::iterator it=begin();
  while(it!=end()&&kineticEnergy>(*it).first) {it++;}
  if(it==end()) G4Exception("G4PiData used outside validity range");
  if(it==begin()) it++;
  G4double x1,x2,e1,e2;
  e1=(*(it-1)).first;
  x1=(*(it-1)).second.second;
  e2=(*(it)).first;
  x2=(*(it)).second.second;
  result = G4std::max(0., x1 + (kineticEnergy-e1)*(x2-x1)/(e2-e1));
  return result;
}

G4double G4PiData::
ElasticXSection(G4double kineticEnergy)
{
  G4double result = 0;
  G4PiData::iterator it=begin();
  while(it!=end()&&kineticEnergy>(*it).first) {it++;}
  if(it==end()) G4Exception("G4PiData used outside validity range");
  if(it==begin()) it++;
  G4double x1,x2,e1,e2;
  e1=(*(it-1)).first;
  x1=(*(it-1)).second.first - (*(it-1)).second.second;
  e2=(*(it)).first;
  x2=(*(it)).second.first - (*(it)).second.second;
  result = G4std::max(0., x1 + (kineticEnergy-e1)*(x2-x1)/(e2-e1));
  return result;
}
