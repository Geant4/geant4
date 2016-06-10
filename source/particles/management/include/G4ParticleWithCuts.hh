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
// $Id: G4ParticleWithCuts.hh 67971 2013-03-13 10:13:24Z gcosmo $
//
// 
// --------------------------------------------------------
//      GEANT 4 class header file
//
//      History: 
//       first implementation, based on object model of Hisaya Kurashige, 
//       calculation of Range Table is based on implementeation for Muon 
//                                           by L.Urban, 10 May 1996
//       added  RestoreCuts  H.Kurashige 09 Mar. 2001
//       introduced material dependent range cuts   08 Sep. 2001
//       restructuring for Cuts per Region  by Hisaya    07 Oct.2002 
//       change to dummy class              by Hisaya    11 MAr.2003 
// ----------------------------------------------------------------
// Class Description
//  Dummy to be removed in future 
//

#ifndef G4ParticleWithCuts_h
#define G4ParticleWithCuts_h 1
#include "G4ParticleDefinition.hh" 

typedef G4ParticleDefinition G4ParticleWithCuts;

#endif











