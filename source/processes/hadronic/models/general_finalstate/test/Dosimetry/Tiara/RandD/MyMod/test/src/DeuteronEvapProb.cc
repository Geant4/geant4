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
// $Id: DeuteronEvapProb.cc,v 1.1 2003-10-08 12:32:19 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "DeuteronEvapProb.hh"

G4double DeuteronEvapProb::CoulombBarrier(G4Fragment& fragm)
{
  G4double rc = 2.173*(1+0.006103*(fragm.GetZ()-1))/(1+0.009443*(fragm.GetZ()-1))*fermi;
  rc *= pow(fragm.GetA()-2,1./3.);
  return (fragm.GetZ()-1)*eplus*eplus*MeV*MeV/coulomb/coulomb/(1.5*fermi*pow(fragm.GetA()-2,1./3.)+1.89*fermi);
}

DeuteronEvapProb::DeuteronEvapProb() : MyEvapProb(2,1,6) // A,Z,Gamma
{
    std::vector<G4double>::size_type NumExcitedStatesEnergy = 31+1;
    std::vector<G4int>::size_type NumExcitedStatesSpin = 31+1;
    ExcitEnergies.reserve(NumExcitedStatesEnergy);
    ExcitSpins.reserve(NumExcitedStatesSpin);
    ExcitEnergies.insert(ExcitEnergies.begin(),NumExcitedStatesEnergy,0.0);
    ExcitSpins.insert(ExcitSpins.begin(),NumExcitedStatesSpin,0.0);


    ExcitEnergies[15] = 6.18*MeV;
    ExcitEnergies[17] = 2.15*MeV;
    ExcitEnergies[18] = 5.02*MeV;
    ExcitEnergies[19] = 2.65*MeV;
    ExcitEnergies[20] = 4.80*MeV;
    ExcitEnergies[22] = 3.85*MeV;
    ExcitEnergies[23] = 6.96*MeV;
    ExcitEnergies[25] = 4.92*MeV;
    ExcitEnergies[26] = 7.22*MeV;
    ExcitEnergies[27] = 0.40*MeV;
    ExcitEnergies[28] = 6.83*MeV;
    ExcitEnergies[29] = 7.12*MeV;
    ExcitEnergies[30] = 3.84*MeV;
    ExcitEnergies[31] = 3.92*MeV;

    ExcitSpins[15] = 1;
    ExcitSpins[17] = 3;
    ExcitSpins[18] = 4;
    ExcitSpins[19] = 4;
    ExcitSpins[20] = 4;
    ExcitSpins[22] = 6;
    ExcitSpins[23] = 6;
    ExcitSpins[25] = 1;
    ExcitSpins[26] = 10;
    ExcitSpins[27] = 3;
    ExcitSpins[28] = 10;
    ExcitSpins[29] = 3;
    ExcitSpins[30] = 6;
    ExcitSpins[31] = 5;
	
    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);	
	
}
DeuteronEvapProb::DeuteronEvapProb(const DeuteronEvapProb &right)
{
    G4Exception("G4DeuteronEvaporationProbability::copy_constructor meant to not be accessable");
}




const DeuteronEvapProb & DeuteronEvapProb::operator=(const DeuteronEvapProb &right)
{
    G4Exception("G4DeuteronEvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool DeuteronEvapProb::operator==(const DeuteronEvapProb &right) const
{
    return false;
}

G4bool DeuteronEvapProb::operator!=(const DeuteronEvapProb &right) const
{
    return true;
}

#include "globals.hh"

G4double DeuteronEvapProb::CalcBetaParam(const G4Fragment& fragment) const
{
  G4double Common = 1.5*fermi*pow(fragment.GetA()-2,1./3) + 1.89*fermi;
  G4double z = fragment.GetZ()-1;
  Common = -z*eplus*eplus*MeV*MeV/coulomb/coulomb /Common;
  if(z<=10)
    return Common*0.048*z;
  else if(z<20)
    return Common*(0.48 + (0.64-0.48)*(z-10)*0.1);
  else if(z<30)
    return Common*(0.64 + (0.74-0.64)*(z-20)*0.1);
  else if(z<50)
    return Common*(0.74 + (0.83-0.74)*(z-30)*0.05);
  else if(z<70)
    return Common*(0.83 + (0.86-0.83)*(z-50)*0.05);
  return Common*0.86;
}

G4double DeuteronEvapProb::CCoeficient(const G4double aZ) const
{
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // 
    // const G4int size = 5;
    // G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
    // G4double Cp[5] = { 0.50, 0.28, 0.20, 0.15, 0.10};
    // C for deuteron is equal to C for protons divided by 2
    G4double C = 0.0;
	
    if (aZ >= 70) {
	C = 0.10;
    } else {
      if(aZ < 10) return 25*aZ*0.1;
      else if(aZ<20) return 25 + (14-25)*(aZ-10)*0.1;
      else if(aZ<30) return 14 + (10-14)*(aZ-20)*0.1;
      else if(aZ<50) return 10 + (7.5-10)*(aZ-30)*0.05;
      else return 7.5+(5-7.5)*(aZ-50)*0.05;
      //	C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
    }
	
    return C/2.0;
	
}
