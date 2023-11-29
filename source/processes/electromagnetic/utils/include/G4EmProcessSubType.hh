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
// G4EmProcessSubType.hh
//
// Class Description:
//   This is an enumerator to define sub-type of electro-magnetic
//   processes
//
// Creation date: 23.09.2008
// Modifications:
//
//---------------------------------------------------------------

#ifndef G4EmProcessSubType_h
#define G4EmProcessSubType_h 1

enum G4EmProcessSubType 
{ 
  fCoulombScattering = 1, 
  fIonisation = 2, 
  fBremsstrahlung = 3, 
  fPairProdByCharged = 4,
  fAnnihilation = 5, 
  fAnnihilationToMuMu = 6,
  fAnnihilationToHadrons = 7,
  fNuclearStopping = 8,
  fElectronGeneralProcess = 9,

  fMultipleScattering = 10, 
  
  fRayleigh = 11,
  fPhotoElectricEffect = 12,
  fComptonScattering = 13,
  fGammaConversion = 14,
  fGammaConversionToMuMu = 15,
  fGammaGeneralProcess = 16,
  fPositronGeneralProcess = 17,
  fAnnihilationToTauTau = 18,
 
  fCerenkov = 21,
  fScintillation = 22,
  fSynchrotronRadiation = 23,
  fTransitionRadiation = 24,

  fSurfaceReflection = 25,
  fDarkBremsstrahlung = 40,
  fMuonPairProdByCharged = 49
};

#endif
