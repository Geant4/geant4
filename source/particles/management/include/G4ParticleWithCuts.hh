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
// $Id: G4ParticleWithCuts.hh,v 1.17 2003/03/17 00:51:03 kurasige Exp $
// GEANT4 tag $Name: geant4-05-01 $
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











