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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ElectronNucleusDataSet.cc,v 1.2 2001-07-11 10:03:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Physics class: ElectronNucleusDataSet for gamma+A cross sections
// M.V. Kossov, ITEP(Moscow), 10-SEP-00
// 

#include "G4ElectronNucleusDataSet.hh"

// The main member function giving the gamma-A cross section (E in GeV, CS in mb)
G4double G4ElectronNucleusDataSet::GetCrossSection(const G4DynamicParticle* aPart,
                                                       const G4Element* anEle)
{
}
