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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: TritonEvapProb.cc,v 1.1 2003-10-08 12:32:20 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "TritonEvapProb.hh"

G4double TritonEvapProb::CoulombBarrier(G4Fragment& fragm)
{
  G4double rc = 2.173*(1+0.006103*(fragm.GetZ()-1))/(1+0.009443*(fragm.GetZ()-1))*fermi;
  rc *= pow(fragm.GetA()-3,1./3.);
  return (fragm.GetZ()-1)*eplus*eplus*MeV*MeV/coulomb/coulomb/(1.5*fermi*pow(fragm.GetA()-3,1./3.)+2.1*fermi);
}

TritonEvapProb::TritonEvapProb() : MyEvapProb(3,1,6) // A,Z,Gamma
{
    std::vector<G4double>::size_type NumExcitedStatesEnergy = 31+1;
    std::vector<G4int>::size_type NumExcitedStatesSpin = 31+1;
    ExcitEnergies.reserve(NumExcitedStatesEnergy);
    ExcitSpins.reserve(NumExcitedStatesSpin);
    ExcitEnergies.insert(ExcitEnergies.begin(),NumExcitedStatesEnergy,0.0);
    ExcitSpins.insert(ExcitSpins.begin(),NumExcitedStatesSpin,0.0);

	
    ExcitEnergies[15] = 6.26*MeV;
    ExcitEnergies[17] = 3.59*MeV;
    ExcitEnergies[18] = 6.76*MeV;
    ExcitEnergies[20] = 6.34*MeV;
    ExcitEnergies[23] = 7.34*MeV;
    ExcitEnergies[25] = 5.11*MeV;
    ExcitEnergies[26] = 7.57*MeV;
    ExcitEnergies[28] = 7.28*MeV;
    ExcitEnergies[31] = 4.46*MeV;

    ExcitSpins[15] = 5;
    ExcitSpins[17] = 5;
    ExcitSpins[18] = 10;
    ExcitSpins[20] = 2;
    ExcitSpins[23] = 5;
    ExcitSpins[25] = 5;
    ExcitSpins[26] = 8;
    ExcitSpins[28] = 8;
    ExcitSpins[31] = 3;

    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);	
}

TritonEvapProb::TritonEvapProb(const TritonEvapProb &right)
{
    G4Exception("G4TritonEvaporationProbability::copy_constructor meant to not be accessable");
}




const TritonEvapProb & TritonEvapProb::operator=(const TritonEvapProb &right)
{
    G4Exception("G4TritonEvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool TritonEvapProb::operator==(const TritonEvapProb &right) const
{
    return false;
}

G4bool TritonEvapProb::operator!=(const TritonEvapProb &right) const
{
    return true;
}

#include "globals.hh"

G4double TritonEvapProb::CalcBetaParam(const G4Fragment& fragment) const
{
  G4double Common = 1.5*fermi*pow(fragment.GetA()-3,1./3)+2.1*fermi;
  G4double z = fragment.GetZ()-1;
  Common = -z*eplus*eplus*MeV*MeV/coulomb/coulomb /Common;
  if(z<=10)
    return Common*0.054*z;
  else if(z<=20)
    return Common*(0.54 + (0.70-0.54)*(z-10)*0.1);
  else if(z<=30)
    return Common*(0.70 + (0.80-0.70)*(z-20)*0.1);
  else if(z<50)
    return Common*(0.80 + (0.89-0.80)*(z-30)*0.05);
  else if(z<70) return Common*(0.89 +(0.92 - 0.89)*(z-50)*0.05);
  return Common*0.92;
}


G4double TritonEvapProb::CCoeficient(const G4double aZ) const
{
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // 
    // const G4int size = 5;
    // G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
    // G4double Cp[5] = { 0.50, 0.28, 0.20, 0.15, 0.10};
    // C for triton is equal to C for protons divided by 3
    G4double C = 0.0;
	
    if (aZ >= 70) {
	C = 0.10;
    } else {
      if(aZ<10) return 0.5*aZ/30.;
      else if(aZ<20) return (0.5 + (0.28-0.5)*(aZ-10)*0.1)/3.;
      else if(aZ<30) return (0.28 + (0.20-0.28)*(aZ-20)*0.1)/3.;
      else if(aZ<50) return (0.20 + (0.15-0.20)*(aZ-30)*0.05)/3.;
      else return (0.15 + (0.10-0.15)*(aZ-50)*0.05)/3.;
      //	C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
    }
	
    return C/3.0;
	
}
