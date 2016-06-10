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
// $Id: G4PreCompoundParameters.cc 68028 2013-03-13 13:48:15Z gcosmo $
//
// by V. Lara
//
// 18.08.2010 V.Ivanchenko make this class as a standard singleton
//

#include "G4PreCompoundParameters.hh"
#include "G4SystemOfUnits.hh"

G4PreCompoundParameters::G4PreCompoundParameters() 
{
  fLevelDensity = 0.10/MeV;
  fR0 = 1.5*fermi;
  fTransitions_r0 = 0.6*fermi;
  fFermiEnergy = 35.0*MeV; 
}

G4PreCompoundParameters::~G4PreCompoundParameters() 
{}


