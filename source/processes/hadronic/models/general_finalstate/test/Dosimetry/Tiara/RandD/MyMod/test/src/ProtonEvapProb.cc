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
// $Id: ProtonEvapProb.cc,v 1.1 2003-10-08 12:32:20 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "ProtonEvapProb.hh"

G4double ProtonEvapProb::CoulombBarrier(G4Fragment& fragm)
{
  G4double rc = 2.173*(1+0.006104*(fragm.GetZ()-1))/(1+0.009443*(fragm.GetZ()-1))*fermi;
  rc *= pow(fragm.GetA()-1,1./3.);
  return (fragm.GetZ()-1)*eplus*eplus*MeV*MeV/coulomb/coulomb/(1.5*fermi*pow(fragm.GetA()-1,1./3.)+0.88*fermi);
}

ProtonEvapProb::ProtonEvapProb() : MyEvapProb(1,1,2) // A,Z,Gamma
{
    std::vector<G4double>::size_type NumExcitedStatesEnergy = 31+1;
    std::vector<G4int>::size_type NumExcitedStatesSpin = 31+1;
    ExcitEnergies.reserve(NumExcitedStatesEnergy);
    ExcitSpins.reserve(NumExcitedStatesSpin);
    ExcitEnergies.insert(ExcitEnergies.begin(),NumExcitedStatesEnergy,0.0);
    ExcitSpins.insert(ExcitSpins.begin(),NumExcitedStatesSpin,0.0);



    ExcitEnergies[15] = 5.96*MeV;
    ExcitEnergies[17] = 1.74*MeV;
    ExcitEnergies[18] = 4.44*MeV;
    ExcitEnergies[19] = 1.67*MeV;
    ExcitEnergies[20] = 4.32*MeV;
    ExcitEnergies[22] = 3.68*MeV;
    ExcitEnergies[23] = 6.69*MeV;
    ExcitEnergies[25] = 3.95*MeV;
    ExcitEnergies[26] = 6.32*MeV;
    ExcitEnergies[27] = 0.30*MeV;
    ExcitEnergies[28] = 6.18*MeV;
    ExcitEnergies[29] = 6.92*MeV;
    ExcitEnergies[30] = 3.06*MeV;
    ExcitEnergies[31] = 3.57*MeV;


    ExcitSpins[15] = 8;
    ExcitSpins[17] = 1;
    ExcitSpins[18] = 6;
    ExcitSpins[19] = 5;
    ExcitSpins[20] = 6;
    ExcitSpins[22] = 4;
    ExcitSpins[23] = 8;
    ExcitSpins[25] = 3;
    ExcitSpins[26] = 4;
    ExcitSpins[27] = 7;
    ExcitSpins[28] = 4;
    ExcitSpins[29] = 5;
    ExcitSpins[30] = 2;
    ExcitSpins[31] = 10;
	
    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);	
	
}

ProtonEvapProb::ProtonEvapProb(const ProtonEvapProb &right)
{
    G4Exception("G4ProtonEvaporationProbability::copy_constructor meant to not be accessable");
}

const ProtonEvapProb & ProtonEvapProb::operator=(const ProtonEvapProb &right)
{
    G4Exception("G4ProtonEvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool ProtonEvapProb::operator==(const ProtonEvapProb &right) const
{
    return false;
}

G4bool ProtonEvapProb::operator!=(const ProtonEvapProb &right) const
{
    return true;
}
#include "globals.hh"

G4double ProtonEvapProb::CalcBetaParam(const G4Fragment& fragment) const
{
  G4double Common = 1.5*fermi*pow(fragment.GetA()-1,1./3)+0.88*fermi;
  G4double z = fragment.GetZ()-1;
  Common = -z*eplus*eplus*MeV*MeV/coulomb/coulomb /Common;
  if(z<=10)
    return Common*0.042*z;
  else if(z<=20)
    return Common*(0.42 + (0.58-0.42)*(z-10)*0.1);
  else if(z<=30)
    return Common*(0.58 + (0.68-0.58)*(z-20)*0.1);
  else if(z<=50)
    return Common*(0.68 + (0.77-0.68)*(z-30)*0.05);
  else if(z<70) return Common*(0.77 + (0.8-0.77)*(z-50)*0.05);
  return Common*0.8;
}

G4double ProtonEvapProb::CCoeficient(const G4double aZ) const
{
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // 
    // const G4int size = 5;
    // G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
    // G4double Cp[5] = { 0.50, 0.28, 0.20, 0.15, 0.10};
    G4double C = 0.0;
	
    if (aZ >= 70) {
	C = 0.10;
    } else {
      if(aZ<10) return 0.50*aZ/10;
      else if(aZ<20) return 0.50 + (0.28-0.50)*(aZ-10)*0.1;
      else if(aZ<30) return 0.28 + (0.20-0.28)*(aZ-20)*0.1;
      else if(aZ<50) return 0.20 + (0.15-0.20)*(aZ-30)*0.05;
      else return 0.15 + (0.10-0.15)*(aZ-50)*0.05;
      //	C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
    }
	
    return C;
	
}
