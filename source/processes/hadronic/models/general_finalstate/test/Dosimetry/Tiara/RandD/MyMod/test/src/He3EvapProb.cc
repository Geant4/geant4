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
// $Id: He3EvapProb.cc,v 1.1 2003-10-08 12:32:19 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "He3EvapProb.hh"

G4double He3EvapProb::CoulombBarrier(G4Fragment& fragm)
{
  G4double rc = 2.173*(1+0.006103*2*(fragm.GetZ()-2))/(1+0.009443*2*(fragm.GetZ()-2))*fermi;
  return 2*(fragm.GetZ()-2)*eplus*eplus*MeV*MeV/coulomb/coulomb/(1.5*fermi*pow(fragm.GetA()-3,1./3.)+2.1*fermi);
}

He3EvapProb::He3EvapProb() : MyEvapProb(3,2,6) // A,Z,Gamma
{
    std::vector<G4double>::size_type NumExcitedStatesEnergy = 31+1;
    std::vector<G4int>::size_type NumExcitedStatesSpin = 31+1;
    ExcitEnergies.reserve(NumExcitedStatesEnergy);
    ExcitSpins.reserve(NumExcitedStatesSpin);
    ExcitEnergies.insert(ExcitEnergies.begin(),NumExcitedStatesEnergy,0.0);
    ExcitSpins.insert(ExcitSpins.begin(),NumExcitedStatesSpin,0.0);


	
    ExcitEnergies[18] = 7.29*MeV;
    ExcitEnergies[20] = 6.48*MeV;
    ExcitEnergies[25] = 5.69*MeV;
    ExcitEnergies[26] = 8.31*MeV;
    ExcitEnergies[31] = 5.10*MeV;

    ExcitSpins[18] = 6;
    ExcitSpins[20] = 8;
    ExcitSpins[25] = 3;
    ExcitSpins[26] = 2;
    ExcitSpins[31] = 7;	
	
    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);	
}

He3EvapProb::He3EvapProb(const He3EvapProb &right)
{
    G4Exception("G4He3EvaporationProbability::copy_constructor meant to not be accessable");
}




const He3EvapProb & He3EvapProb::operator=(const He3EvapProb &right)
{
    G4Exception("G4He3EvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool He3EvapProb::operator==(const He3EvapProb &right) const
{
    return false;
}

G4bool He3EvapProb::operator!=(const He3EvapProb &right) const
{
    return true;
}

#include "globals.hh"

G4double He3EvapProb::CalcBetaParam(const G4Fragment& fragment) const
{
  G4double Common = 1.5*fermi*pow(fragment.GetA()-3,1./3.)+2.1*fermi;
  G4double z = fragment.GetZ()-2;
  Common = -2*z*eplus*eplus*MeV*MeV/coulomb/coulomb /Common;
  if(z<=10)
    return Common*0.62*z/10;
  else if(z<=20)
    return Common*(0.62 + (0.76-0.62)*(z-10)*0.1);
  else if(z<=30)
    return Common*(0.76 + (0.85-0.76)*(z-20)*0.1);
  else if(z<=50) return Common*(0.85 + (0.91-0.85)*(z-30)*0.05);
  else if(z<70) return Common*(0.91 + (0.92-0.91)*(z-50)*0.05);
  return Common*0.92;
}

G4double He3EvapProb::CCoeficient(const G4double aZ) const
{
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // 
    // const G4int size = 5;
    // G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
    //	G4double Calpha[5] = { 0.10, 0.10, 0.10, 0.08, 0.06};
    // C for He3 is equal to C for alpha times 4/3
    G4double C = 0.0;
	
	
    if (aZ <= 30) {
	C = 0.10;
    } else if (aZ <= 50) {
	C = 0.1 + -((aZ-50.)/20.)*0.02;
    } else if (aZ < 70) {
	C = 0.08 + -((aZ-70.)/20.)*0.02;
    } else {
	C = 0.06;
    }
    return C*(4.0/3.0);
}

