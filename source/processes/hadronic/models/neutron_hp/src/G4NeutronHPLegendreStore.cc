// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPLegendreStore.hh"
#include "G4NeutronHPVector.hh"
#include "G4NeutronHPInterpolator.hh"
#include "G4NeutronHPFastLegendre.hh"
#include "Randomize.hh"
#include <iostream.h>

  G4NeutronHPLegendreStore::G4NeutronHPLegendreStore(G4int n)
  {
    theCoeff = new G4NeutronHPLegendreTable[n];
    nEnergy = n;
  }
  
  G4NeutronHPLegendreStore::~G4NeutronHPLegendreStore()
  {
    delete [] theCoeff;
  }

G4double G4NeutronHPLegendreStore::SampleMax (G4double anEnergy)
{
  G4double result;
  
  G4int i0, i1, i2, i3, i4;
  G4int low, high;
  G4NeutronHPFastLegendre theLeg;
  for (i0=0; i0<nEnergy; i0++)
  {
    high = i0;
    if(theCoeff[i0].GetEnergy()>anEnergy) break;
  }
  low = max(0, high-1);
  G4NeutronHPInterpolator theInt;
  G4double lim=0.005;
  G4double x, x1, x2, y1, y2, y;
  x = anEnergy;
  x1 = theCoeff[low].GetEnergy();
  x2 = theCoeff[high].GetEnergy();
  G4double theNorm = 0;
  G4double try01=0, try02=0;
  G4double max1, max2, costh;
  max1 = 0; max2 = 0;
  G4int l;
  for(i0=0; i0<601; i0++)
  {
    costh = G4double(i0-300)/300.;
    try01 = 0;
    for(l=0; l<theCoeff[low].GetNumberOfPoly() ; l++)
    {
      try01 += (2.*l+1)/2.*theCoeff[low].GetCoeff(l)*theLeg.Evaluate(l, costh);
    } 
    if(try01>max1) max1=try01;
    try02 = 0;
    for(l=0; l<theCoeff[high].GetNumberOfPoly() ; l++)
    {
      try02 += (2.*l+1)/2.*theCoeff[high].GetCoeff(l)*theLeg.Evaluate(l, costh);
    }
    if(try02>max2) max2=try02;
  } 
  theNorm = theInt.Interpolate(theManager.GetScheme(high), x, x1, x2, max1, max2);
  
  G4double value, random;
  G4double v1, v2;
  do
  {
    v1 = 0;
    v2 = 0;
    result = 2.*G4UniformRand()-1.;
    for(l=0; l<theCoeff[low].GetNumberOfPoly() ; l++)
    {
      G4double legend = theLeg.Evaluate(l, result); // @@@ done to avoid optimization error on SUN
      v1 += (2.*l+1)/2.*theCoeff[low].GetCoeff(l)*legend;
    } 
    for(l=0; l<theCoeff[high].GetNumberOfPoly() ; l++)
    {
      G4double legend = theLeg.Evaluate(l, result); // @@@ done to avoid optimization error on SUN
      v2 += (2.*l+1)/2.*theCoeff[high].GetCoeff(l)*legend;
    } 
    value = theInt.Interpolate(theManager.GetScheme(high), x, x1, x2, v1, v2);
    random = G4UniformRand();
  }
  while(random>value/theNorm);
  
  return result;
}


G4double G4NeutronHPLegendreStore::SampleElastic (G4double anEnergy)
{
  G4double result;
  
  G4int i0, i1, i2, i3, i4;
  G4int low, high;
  G4NeutronHPFastLegendre theLeg;
  for (i0=0; i0<nEnergy; i0++)
  {
    high = i0;
    if(theCoeff[i0].GetEnergy()>anEnergy) break;
  }
  low = max(0, high-1);
  G4NeutronHPInterpolator theInt;
  G4double lim=0.005;
  G4double x, x1, x2, y1, y2, y;
  x = anEnergy;
  x1 = theCoeff[low].GetEnergy();
  x2 = theCoeff[high].GetEnergy();
  G4double theNorm = 0;
  G4double try01=0, try02=0, try11=0, try12=0;
  G4double try1, try2;
  G4int l;
  for(l=0; l<theCoeff[low].GetNumberOfPoly(); l++)
  {
    try01 += (2.*l+1)/2.*theCoeff[low].GetCoeff(l)*theLeg.Evaluate(l, -1.);
    try11 += (2.*l+1)/2.*theCoeff[low].GetCoeff(l)*theLeg.Evaluate(l, +1.);
  } 
  for(l=0; l<theCoeff[high].GetNumberOfPoly(); l++)
  {
    try02 += (2.*l+1)/2.*theCoeff[high].GetCoeff(l)*theLeg.Evaluate(l, -1.);
    try12 += (2.*l+1)/2.*theCoeff[high].GetCoeff(l)*theLeg.Evaluate(l, +1.);
  } 
  try1 = theInt.Interpolate(theManager.GetScheme(high), x, x1, x2, try01, try02);
  try2 = theInt.Interpolate(theManager.GetScheme(high), x, x1, x2, try11, try12);
  theNorm = max(try1, try2);
  
  G4double value, random;
  G4double v1, v2;
  do
  {
    v1 = 0;
    v2 = 0;
    result = 2.*G4UniformRand()-1.;
    for(l=0; l<theCoeff[low].GetNumberOfPoly() ; l++)
    {
      G4double legend = theLeg.Evaluate(l, result); // @@@ done to avoid optimization error on SUN
      v1 += (2.*l+1)/2.*theCoeff[low].GetCoeff(l)*legend;
    } 
    for(l=0; l<theCoeff[high].GetNumberOfPoly() ; l++)
    {
      G4double legend = theLeg.Evaluate(l, result); // @@@ done to avoid optimization error on SUN
      v2 += (2.*l+1)/2.*theCoeff[high].GetCoeff(l)*legend;
    } 
    value = theInt.Interpolate(theManager.GetScheme(high), x, x1, x2, v1, v2);
    random = G4UniformRand();
  }
  while(random>value/theNorm);
  
  return result;
}

G4double G4NeutronHPLegendreStore::Sample (G4double energy) // still in interpolation; do not use
{
  G4int i0, i1, i2, i3, i4;
  G4int low, high;
//  G4cout << "G4NeutronHPLegendreStore::Sample "<<energy<<" "<<energy<<" "<<nEnergy<<endl;
  for (i0=0; i0<nEnergy; i0++)
  {
//     G4cout <<"theCoeff["<<i0<<"].GetEnergy() = "<<theCoeff[i0].GetEnergy()<<endl;
    high = i0;
    if(theCoeff[i0].GetEnergy()>energy) break;
  }
  low = max(0, high-1);
//  G4cout << "G4NeutronHPLegendreStore::Sample high, low: "<<high<<", "<<low<<endl;
  G4NeutronHPVector theBuffer;
  G4NeutronHPInterpolator theInt;
  G4double lim=0.005;
  G4double x1, x2, y1, y2, y;
  x1 = theCoeff[low].GetEnergy();
  x2 = theCoeff[high].GetEnergy();
//  G4cout << "the xes "<<x1<<" "<<x2<<endl;
  G4double costh=0;
  for(i0=0; i0<601; i0++)
  {
    costh = G4double(i0-300)/300.;
    y1 = Integrate(low, costh);
    y2 = Integrate(high, costh);
    y = theInt.Interpolate(theManager.GetScheme(high), energy, x1, x2, y1, y2);
    theBuffer.SetData(i0, costh, y);
//     G4cout << "Integration "<<low<<" "<<costh<<" "<<y1<<" "<<y2<<" "<<y<<endl;
  }
  G4double rand = G4UniformRand();
  G4int it;
  for (i0=1; i0<601; i0++)
  {
    it = i0;
    if(rand < theBuffer.GetY(i0)/theBuffer.GetY(600)) break;
//    G4cout <<"sampling now "<<i0<<" "
//         << theBuffer.GetY(i0)<<" "
//         << theBuffer.GetY(600)<<" "
//         << rand<<" "
//         << theBuffer.GetY(i0)/theBuffer.GetY(600)<<endl;;
  }
  if(it==601) it=600;
//  G4cout << "G4NeutronHPLegendreStore::Sample it "<<rand<<" "<<it<<endl;
  G4double norm = theBuffer.GetY(600);
  if(norm==0) return -DBL_MAX; 
  x1 = theBuffer.GetY(it)/norm;
  x2 = theBuffer.GetY(it-1)/norm;
  y1 = theBuffer.GetX(it);
  y2 = theBuffer.GetX(it-1);
//  G4cout << "G4NeutronHPLegendreStore::Sample x y "<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<endl;
  return theInt.Interpolate(theManager.GetScheme(high), rand, x1, x2, y1, y2);
}

G4double G4NeutronHPLegendreStore::Integrate(G4int k, G4double costh) // still in interpolation; not used anymore
{
  G4double result=0;
  G4NeutronHPFastLegendre theLeg;
//  G4cout <<"the COEFFS "<<k<<" ";
//  G4cout <<theCoeff[k].GetNumberOfPoly()<<" ";
  for(G4int l=0; l<theCoeff[k].GetNumberOfPoly() ; l++)
  {
    result += theCoeff[k].GetCoeff(l)*theLeg.Integrate(l, costh);
//    G4cout << theCoeff[k].GetCoeff(l)<<" ";
  } 
//  G4cout <<endl;
  return result;
}
