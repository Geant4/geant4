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
//---------------------------------------------------------------
//
// G4EmSecondaryParticleType.hh
//
// Class Description:
//   This is an enumerator to define type of secondary according
//   to G4PhysicsModelCatalog
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 30.09.2021
//
//---------------------------------------------------------------

#ifndef G4EmSecondaryParticleType_h
#define G4EmSecondaryParticleType_h 1

enum G4EmSecondaryParticleType
  {
    _EM = 10000,

    _DeltaElectron = 10010,
    _DeltaEBelowCut = 10011,
    _PhotoElectron = 10012,
    _ComptonElectron = 10013,
    _TripletElectron = 10014,

    _Bremsstrahlung = 10020,
    // Legacy name for compatibility with Geant4 11.0 and patch01.
    _Bremsstruhlung = 10020,
    _SplitBremsstrahlung = 10021,
    _ComptonGamma = 10022,
    _Annihilation = 10023,
    _TripletGamma = 10024,
    _GammaGammaEntanglement = 10025,

    _PairProduction = 10030,

    _Fluorescence = 10040,
    _GammaPIXE = 10041,

    _AugerElectron = 10050,
    _ePIXE = 10051,

    _IonRecoil = 10060
  };

#endif


