// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
#include "G4NeutronHPMadlandNixSpectrum.hh"
  
  G4double G4NeutronHPMadlandNixSpectrum::Madland(G4double aSecEnergy, G4double tm)
  {
    G4double result;
    G4double energy = aSecEnergy/eV;
    G4double EF;
    
    EF = theAvarageKineticPerNucleonForLightFragments/eV;
    G4double lightU1 = sqrt(energy)-sqrt(EF);
    lightU1 *= lightU1/tm;
    G4double lightU2 = sqrt(energy)+sqrt(EF);
    lightU2 *= lightU2/tm;
    G4double lightTerm=0;
    if(theAvarageKineticPerNucleonForLightFragments>1*eV)
    {
      lightTerm  = pow(lightU2, 1.5)*E1(lightU2);
      lightTerm -= pow(lightU1, 1.5)*E1(lightU1);
      lightTerm += Gamma15(lightU2)-Gamma15(lightU1);
      lightTerm /= 3.*sqrt(tm*EF);
    }
    
    EF = theAvarageKineticPerNucleonForHeavyFragments/eV;
    G4double heavyU1 = sqrt(energy)-sqrt(EF);
    heavyU1 *= heavyU1/tm;
    G4double heavyU2 = sqrt(energy)+sqrt(EF);
    heavyU2 *= heavyU2/tm;
    G4double heavyTerm=0	;
    if(theAvarageKineticPerNucleonForHeavyFragments> 1*eV)
    {
      heavyTerm  = pow(heavyU2, 1.5)*E1(heavyU2);
      heavyTerm -= pow(heavyU1, 1.5)*E1(heavyU1);
      heavyTerm += Gamma15(heavyU2)-Gamma15(heavyU1);
      heavyTerm /= 3.*sqrt(tm*EF);
    }
    
    result = 0.5*(lightTerm+heavyTerm);
    
    return result;
  }

  G4double G4NeutronHPMadlandNixSpectrum::Sample(G4double anEnergy) 
  {
    G4bool Done = false;
    G4double tm = theMaxTemp.GetY(anEnergy);
    G4double last=0, buff, current = 100*MeV;
    G4double precision = 0.001;
    G4double newValue = 0., oldValue=0., diff=0.;
    G4double random = G4UniformRand();
    
    do
    {
      oldValue = newValue;
      newValue = FissionIntegral(tm, current);
      if(newValue < random) 
      {
        buff = current;
	current+=abs(current-last)/2.;
	last = buff;
        if(current>190*MeV) G4Exception("Madland-Nix Spectrum has not converged in sampling");
      }
      else
      {
        buff = current;
        current-=abs(current-last)/2.;
	last = buff;
      }
    }
    while (abs(oldValue-newValue)>precision*newValue);
    return current;
  }

  G4double G4NeutronHPMadlandNixSpectrum::
  GIntegral(G4double tm, G4double anEnergy, G4double aMean)
  {
    if(aMean<1*eV) return 0;
    G4double b = anEnergy/eV;
    G4double sb = sqrt(b);
    G4double EF = aMean/eV;
    
    G4double alpha = sqrt(tm); 
    G4double beta = sqrt(EF);
    G4double A = EF/tm;
    G4double B =  (sb+beta)*(sb+beta)/tm;
    G4double Ap = A;
    G4double Bp = (sb-beta)*(sb-beta)/tm;
    
    G4double result;
    G4double alpha2 = alpha*alpha;
    G4double alphabeta = alpha*beta;
    if(b<EF)
    {
      result =
      (
        (0.4*alpha2*pow(B,2.5) - 0.5*alphabeta*B*B)*E1(B) -  
        (0.4*alpha2*pow(A,2.5) - 0.5*alphabeta*A*A)*E1(A) 
      )
       -
      (
        (0.4*alpha2*pow(Bp,2.5) + 0.5*alphabeta*Bp*Bp)*E1(Bp) -  
        (0.4*alpha2*pow(Ap,2.5) + 0.5*alphabeta*Ap*Ap)*E1(Ap) 
      )
      +
      (
        (alpha2*B-2*alphabeta*sqrt(B))*Gamma15(B)  - 
        (alpha2*A-2*alphabeta*sqrt(A))*Gamma15(A) 
      )
      -
      (
        (alpha2*Bp-2*alphabeta*sqrt(Bp))*Gamma15(Bp) -
        (alpha2*Ap-2*alphabeta*sqrt(Ap))*Gamma15(Ap)
      )
      - 0.6*alpha2*(Gamma25(B) - Gamma25(A) - Gamma25(Bp) + Gamma25(Ap))
      - 1.5*alphabeta*(exp(-B)*(1+B) - exp(-A)*(1+A) + exp(-Bp)*(1+Bp) + exp(-Ap)*(1+Ap)) ;
    }
    else
    {
      result =
      (
        (0.4*alpha2*pow(B,2.5) - 0.5*alphabeta*B*B)*E1(B) -  
        (0.4*alpha2*pow(A,2.5) - 0.5*alphabeta*A*A)*E1(A) 
      );
       result -=
      (
        (0.4*alpha2*pow(Bp,2.5) + 0.5*alphabeta*Bp*Bp)*E1(Bp) -  
        (0.4*alpha2*pow(Ap,2.5) + 0.5*alphabeta*Ap*Ap)*E1(Ap) 
      );
      result +=
      (
        (alpha2*B-2*alphabeta*sqrt(B))*Gamma15(B)  - 
        (alpha2*A-2*alphabeta*sqrt(A))*Gamma15(A) 
      );
      result -=
      (
        (alpha2*Bp+2*alphabeta*sqrt(Bp))*Gamma15(Bp) -
        (alpha2*Ap+2*alphabeta*sqrt(Ap))*Gamma15(Ap)
      );
      result -= 0.6*alpha2*(Gamma25(B) - Gamma25(A) - Gamma25(Bp) + Gamma25(Ap));
      result -= 1.5*alphabeta*(exp(-B)*(1+B) - exp(-A)*(1+A) + exp(-Bp)*(1+Bp) + exp(-Ap)*(1+Ap) - 2.) ;
    }
    result = result / (3.*sqrt(tm*EF));
    return result;
  }
