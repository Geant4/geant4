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
//
//---------------------------------------------------------------
//
// G4LowEnergyEmProcessSubType.hh
//
// Class Description:
//   This is an enumerator to define sub-type of electro-magnetic
//   low-energy processes
//
// Creation date: 23.09.2021
// Modifications:
//
//---------------------------------------------------------------

#ifndef G4LowEnergyEmProcessSubType_h
#define G4LowEnergyEmProcessSubType_h 1

enum G4LowEnergyEmProcessSubType 
{ 
  fLowEnergyElastic = 51,
  fLowEnergyExcitation = 52,
  fLowEnergyIonisation = 53,
  fLowEnergyVibrationalExcitation = 54,
  fLowEnergyAttachment = 55,
  fLowEnergyChargeDecrease = 56,
  fLowEnergyChargeIncrease = 57,
  fLowEnergyElectronSolvation = 58,
  fLowEnergyMolecularDecay = 59,
  fLowEnergyTransportation = 60,
  fLowEnergyBrownianTransportation = 61,
  fLowEnergyDoubleIonisation = 62,
  fLowEnergyDoubleCap = 63,
  fLowEnergyIoniTransfer = 64,
  fLowEnergyStaticMol = 65,
  fLowEnergyScavenger = 66
};

#endif
