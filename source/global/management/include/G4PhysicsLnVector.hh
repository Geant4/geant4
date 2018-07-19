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
// $Id: G4PhysicsLnVector.hh 98864 2016-08-15 11:53:26Z gcosmo $
//
// 
//--------------------------------------------------------------------
//      GEANT 4 class header file
//
//  G4PhysicsLnVector.hh
//
//  Class description:
//
//    A physics vector which has values of energy-loss, cross-section, 
//    and other physics values of a particle in matter in a given 
//    range of the energy, momentum, etc. The scale of energy/momentum
//    bins is natural logarithmic.
//
//  History:
//    27 Apr. 1999, M.G. Pia: Created, copying from G4PhysicsLogVector 
//    11 Nov. 2000, H.Kurashige : Use STL vector for dataVector and binVector
//    16 Aug. 2011  H.Kurashige : Move dBin, baseBin to the base class
//    02 Oct. 2013  V.Ivanchenko : Remove FindBinLocation method
//
//--------------------------------------------------------------------

#ifndef G4PhysicsLnVector_h
#define G4PhysicsLnVector_h 1

#include "G4PhysicsLogVector.hh"

typedef G4PhysicsLogVector G4PhysicsLnVector;

#endif
