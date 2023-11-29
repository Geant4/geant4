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
//
//
//---------------------------------------------------------------
//
// G4HadronicProcessType.hh
//
// Class Description:
//   Enumerator to define hadronic process sub-type
//
//---------------------------------------------------------------

#ifndef G4HadronicProcessType_h
#define G4HadronicProcessType_h 1

enum G4HadronicProcessType
{
  fHadronElastic =    111,
  fNeutronGeneral =   116,
  fHadronInelastic =  121,
  fCapture =          131,
  fMuAtomicCapture =  132,
  fFission =          141,
  fHadronAtRest =     151,
  fLeptonAtRest =     152,
  fChargeExchange =   161,
  fRadioactiveDecay = 210,
  fEMDissociation =   310
};
#endif
