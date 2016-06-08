// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPKallbachMannSyst.hh" 
#include "Randomize.hh" 

G4double G4NeutronHPKallbachMannSyst::Sample(G4double anEnergy)
{
  G4double result;
  
  G4double zero = GetKallbachZero(anEnergy);
  if(zero>1) zero=1.;
  if(zero<-1)zero=-1.;
  G4double max = Kallbach(zero, anEnergy);
  double upper = Kallbach(1., anEnergy);
  double lower = Kallbach(-1., anEnergy);
  if(upper>max) max=upper;
  if(lower>max) max=lower;
  G4double value, random;
  do
  {
    result = 2.*G4UniformRand()-1;  
    value = Kallbach(result, anEnergy)/max;
    random = G4UniformRand();
  }
  while(random>value);
  
  return result;
}

G4double G4NeutronHPKallbachMannSyst::Kallbach(G4double cosTh, G4double anEnergy)
{
  // Kallbach-Mann systematics without normalization.
  G4double result;
  G4double theX = A(anEnergy)*cosTh;
  result = 0.5*(exp( theX)*(1+theCompoundFraction)
               +exp(-theX)*(1-theCompoundFraction));
  return result;
}

G4double G4NeutronHPKallbachMannSyst::GetKallbachZero(G4double anEnergy)
{
  G4double result;
  result = 0.5 * (1./A(anEnergy)) * log((1-theCompoundFraction)/(1+theCompoundFraction));
  return result;
}

G4double G4NeutronHPKallbachMannSyst::A(G4double anEnergy)
{ 
  G4double result;
  G4double C1 = 0.04/MeV;
  G4double C2 = 1.8E-6/(MeV*MeV*MeV);
  G4double C3 = 6.7E-7/(MeV*MeV*MeV*MeV);
  
    G4double epsa = anEnergy*theTargetMass/(theTargetMass+theIncidentMass);
    G4int Ac = theTargetA+1; 
    G4int Nc = Ac - theTargetZ;
    G4int AA = theTargetA;
    G4int ZA = theTargetZ;
    G4double ea = epsa+SeparationEnergy(Ac, Nc, AA, ZA);
    G4double Et1 = 130*MeV;
    G4double R1 = G4std::min(ea, Et1);
    // theProductEnergy is still in CMS!!!
    G4double epsb = theProductEnergy*(theProductMass+theResidualMass)/theResidualMass;
    G4int AB = theResidualA;
    G4int ZB = theResidualZ;
    G4double eb = epsb+SeparationEnergy(Ac, Nc, AB, ZB );
  G4double X1 = R1*eb/ea; 
    G4double Et3 = 41*MeV;
    G4double R3 = G4std::min(ea, Et3);
  G4double X3 = R3*eb/ea;
  G4double Ma = 1;
  G4double mb;
  G4int productA = theTargetA+1-theResidualA;
  G4int productZ = theTargetZ-theResidualZ;
  if(productZ==0)
  {
    mb = 0.5;
  }
  else if(productZ==1)
  {
    mb = 1;
  }
  else if(productZ==2)
  {
    mb = 2;
    if(productA==3) mb=1;
  }
  else
  {
    G4Exception("Severe error in the sampling of Kallbach-Mann Systematics");
  }
  
  result = C1*X1 + C2*pow(X1, 3.) + C3*Ma*mb*pow(X3, 4.);
  return result;
}

G4double G4NeutronHPKallbachMannSyst::SeparationEnergy(G4int Ac, G4int Nc, G4int AA, G4int ZA)
{
  G4double result;
  G4int NA = AA-ZA;
  G4int Zc = Ac-Nc;
  result = 15.68*(Ac-AA);
  result += -28.07*((Nc-Zc)*(Nc-Zc)/Ac - (NA-ZA)*(NA-ZA)/AA);
  result += -18.56*(pow(Ac, 2./3.) - pow(AA, 2./3.));
  result +=  33.22*((Nc-Zc)*(Nc-Zc)/pow(Ac, 4./3.) - (NA-ZA)*(NA-ZA)/pow(AA, 4./3.));
  result += -0.717*(Zc*Zc/pow(Ac,1./3.)-ZA*ZA/pow(AA,1./3.));
  result +=  1.211*(Zc*Zc/Ac-ZA*ZA/AA);
  G4double totalBinding;
  G4int productA = theTargetA+1-theResidualA;
  G4int productZ = theTargetZ-theResidualZ;
  if(productZ==0&&productA==1) totalBinding=0;
  if(productZ==1&&productA==1) totalBinding=0;
  if(productZ==1&&productA==2) totalBinding=2.22;
  if(productZ==1&&productA==3) totalBinding=8.48;
  if(productZ==2&&productA==3) totalBinding=7.72;
  if(productZ==2&&productA==4) totalBinding=28.3;
  result += -totalBinding;
  result *= MeV;
  return result;
}
