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
// $Id: AlphaEvapProb.cc,v 1.1 2003-10-08 12:32:19 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "AlphaEvapProb.hh"

G4double AlphaEvapProb::CoulombBarrier(G4Fragment& fragm)
{
  G4double rc = 2.173*(1+0.006103*2*(fragm.GetZ()-2))/(1+0.009443*2*(fragm.GetZ()-2));
  return 2*(fragm.GetZ()-2)*eplus*eplus*MeV*MeV/coulomb/coulomb/(1.5*fermi*pow(fragm.GetA()-4,1./3.)+2.4*fermi);
}

AlphaEvapProb::AlphaEvapProb() :
    MyEvapProb(4,2,4) // A,Z,Gamma
{
    //  const G4int NumExcitedStates = 31+1;
    std::vector<G4double>::size_type NumExcitedStatesEnergy = 31+1;
    std::vector<G4int>::size_type NumExcitedStatesSpin = 31+1;
    ExcitEnergies.reserve(NumExcitedStatesEnergy);
    ExcitSpins.reserve(NumExcitedStatesSpin);
    ExcitEnergies.insert(ExcitEnergies.begin(),NumExcitedStatesEnergy,0.0);
    ExcitSpins.insert(ExcitSpins.begin(),NumExcitedStatesSpin,0.0);
//    for (G4int i = 0; i < NumExcitedStates; i++) {
//      ExcitEnergies(i) = 0.0;
//      ExcitSpins(i) = 0;
//    }

    ExcitEnergies[18] = 7.98*MeV;
    ExcitEnergies[20] = 6.90*MeV;
    ExcitEnergies[25] = 5.83*MeV;
    ExcitEnergies[26] = 8.57*MeV;
    ExcitEnergies[31] = 5.33*MeV;

    ExcitSpins[18] = 4;
    ExcitSpins[20] = 6;
    ExcitSpins[25] = 7;
    ExcitSpins[26] = 4;
    ExcitSpins[31] = 13;
	
    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);		
}


AlphaEvapProb::AlphaEvapProb(const AlphaEvapProb &right)
{
    G4Exception("G4AlphaEvaporationProbability::copy_constructor meant to not be accessable");
}




const AlphaEvapProb & AlphaEvapProb::operator=(const AlphaEvapProb &right)
{
    G4Exception("G4AlphaEvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool AlphaEvapProb::operator==(const AlphaEvapProb &right) const
{
    return false;
}

G4bool AlphaEvapProb::operator!=(const AlphaEvapProb &right) const
{
    return true;
}
#include "globals.hh"
G4double AlphaEvapProb::CalcBetaParam(const G4Fragment& fragment) const
{
  G4double z = fragment.GetZ()-2;
  G4double Common = 1.5*fermi*pow(fragment.GetA()-4,1./3) + 2.4*fermi;
  Common = -2*z*eplus*eplus*MeV*MeV/coulomb/coulomb /Common;
  if(z<=10)
    return Common*0.068*z;
  else if(z<=20)
    return Common*(0.68 + (0.82-0.68)*(z-10)*0.1);
  else if(z<=30)
    return Common*(0.82 +(0.91-0.82)*(z-20)*0.1);
  else if(z<50)
    return Common*(0.91 + (0.97-0.91)*(z-30)*0.05);
  else if(z<70)
    return Common*(0.97 + (0.98-0.97)*(z-50)*0.05);
  return Common*0.98;
}

G4double AlphaEvapProb::CCoeficient(const G4double aZ) const
{
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // 
    // const G4int size = 5;
    // G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
    //	G4double Calpha[5] = { 0.10, 0.10, 0.10, 0.08, 0.06};
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
    return C;
}
