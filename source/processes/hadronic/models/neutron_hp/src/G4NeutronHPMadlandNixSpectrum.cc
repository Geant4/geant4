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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.

#include "G4NeutronHPMadlandNixSpectrum.hh"
#include "G4SystemOfUnits.hh"

  G4double G4NeutronHPMadlandNixSpectrum::Madland(G4double aSecEnergy, G4double tm)
  {
    G4double result;
    G4double energy = aSecEnergy/eV;
    G4double EF;
    
    EF = theAvarageKineticPerNucleonForLightFragments/eV;
    G4double lightU1 = std::sqrt(energy)-std::sqrt(EF);
    lightU1 *= lightU1/tm;
    G4double lightU2 = std::sqrt(energy)+std::sqrt(EF);
    lightU2 *= lightU2/tm;
    G4double lightTerm=0;
    if(theAvarageKineticPerNucleonForLightFragments>1*eV)
    {
      lightTerm  = std::pow(lightU2, 1.5)*E1(lightU2);
      lightTerm -= std::pow(lightU1, 1.5)*E1(lightU1);
      lightTerm += Gamma15(lightU2)-Gamma15(lightU1);
      lightTerm /= 3.*std::sqrt(tm*EF);
    }
    
    EF = theAvarageKineticPerNucleonForHeavyFragments/eV;
    G4double heavyU1 = std::sqrt(energy)-std::sqrt(EF);
    heavyU1 *= heavyU1/tm;
    G4double heavyU2 = std::sqrt(energy)+std::sqrt(EF);
    heavyU2 *= heavyU2/tm;
    G4double heavyTerm=0	;
    if(theAvarageKineticPerNucleonForHeavyFragments> 1*eV)
    {
      heavyTerm  = std::pow(heavyU2, 1.5)*E1(heavyU2);
      heavyTerm -= std::pow(heavyU1, 1.5)*E1(heavyU1);
      heavyTerm += Gamma15(heavyU2)-Gamma15(heavyU1);
      heavyTerm /= 3.*std::sqrt(tm*EF);
    }
    
    result = 0.5*(lightTerm+heavyTerm);
    
    return result;
  }

  G4double G4NeutronHPMadlandNixSpectrum::Sample(G4double anEnergy) 
  {
    G4double tm = theMaxTemp.GetY(anEnergy);
    G4double last=0, buff, current = 100*MeV;
    G4double precision = 0.001;
    G4double newValue = 0., oldValue=0.;
    G4double random = G4UniformRand();
    
    do
    {
      oldValue = newValue;
      newValue = FissionIntegral(tm, current);
      if(newValue < random) 
      {
        buff = current;
	current+=std::abs(current-last)/2.;
	last = buff;
        if(current>190*MeV) throw G4HadronicException(__FILE__, __LINE__, "Madland-Nix Spectrum has not converged in sampling");
      }
      else
      {
        buff = current;
        current-=std::abs(current-last)/2.;
	last = buff;
      }
    }
    while (std::abs(oldValue-newValue)>precision*newValue);
    return current;
  }

  G4double G4NeutronHPMadlandNixSpectrum::
  GIntegral(G4double tm, G4double anEnergy, G4double aMean)
  {
    if(aMean<1*eV) return 0;
    G4double b = anEnergy/eV;
    G4double sb = std::sqrt(b);
    G4double EF = aMean/eV;
    
    G4double alpha = std::sqrt(tm); 
    G4double beta = std::sqrt(EF);
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
        (0.4*alpha2*std::pow(B,2.5) - 0.5*alphabeta*B*B)*E1(B) -  
        (0.4*alpha2*std::pow(A,2.5) - 0.5*alphabeta*A*A)*E1(A) 
      )
       -
      (
        (0.4*alpha2*std::pow(Bp,2.5) + 0.5*alphabeta*Bp*Bp)*E1(Bp) -  
        (0.4*alpha2*std::pow(Ap,2.5) + 0.5*alphabeta*Ap*Ap)*E1(Ap) 
      )
      +
      (
        (alpha2*B-2*alphabeta*std::sqrt(B))*Gamma15(B)  - 
        (alpha2*A-2*alphabeta*std::sqrt(A))*Gamma15(A) 
      )
      -
      (
        (alpha2*Bp-2*alphabeta*std::sqrt(Bp))*Gamma15(Bp) -
        (alpha2*Ap-2*alphabeta*std::sqrt(Ap))*Gamma15(Ap)
      )
      - 0.6*alpha2*(Gamma25(B) - Gamma25(A) - Gamma25(Bp) + Gamma25(Ap))
      - 1.5*alphabeta*(std::exp(-B)*(1+B) - std::exp(-A)*(1+A) + std::exp(-Bp)*(1+Bp) + std::exp(-Ap)*(1+Ap)) ;
    }
    else
    {
      result =
      (
        (0.4*alpha2*std::pow(B,2.5) - 0.5*alphabeta*B*B)*E1(B) -  
        (0.4*alpha2*std::pow(A,2.5) - 0.5*alphabeta*A*A)*E1(A) 
      );
       result -=
      (
        (0.4*alpha2*std::pow(Bp,2.5) + 0.5*alphabeta*Bp*Bp)*E1(Bp) -  
        (0.4*alpha2*std::pow(Ap,2.5) + 0.5*alphabeta*Ap*Ap)*E1(Ap) 
      );
      result +=
      (
        (alpha2*B-2*alphabeta*std::sqrt(B))*Gamma15(B)  - 
        (alpha2*A-2*alphabeta*std::sqrt(A))*Gamma15(A) 
      );
      result -=
      (
        (alpha2*Bp+2*alphabeta*std::sqrt(Bp))*Gamma15(Bp) -
        (alpha2*Ap+2*alphabeta*std::sqrt(Ap))*Gamma15(Ap)
      );
      result -= 0.6*alpha2*(Gamma25(B) - Gamma25(A) - Gamma25(Bp) + Gamma25(Ap));
      result -= 1.5*alphabeta*(std::exp(-B)*(1+B) - std::exp(-A)*(1+A) + std::exp(-Bp)*(1+Bp) + std::exp(-Ap)*(1+Ap) - 2.) ;
    }
    result = result / (3.*std::sqrt(tm*EF));
    return result;
  }
