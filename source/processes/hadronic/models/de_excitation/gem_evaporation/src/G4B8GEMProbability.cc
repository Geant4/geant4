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
// $Id: G4B8GEMProbability.cc 87017 2014-11-21 16:26:26Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4B8GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4B8GEMProbability::G4B8GEMProbability() :
  G4GEMProbability(8,5,2.0) // A,Z,Spin
{
    ExcitEnergies.push_back(769.5*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(fPlanck/(40.0*keV));
    
    ExcitEnergies.push_back(2320.0*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(fPlanck/(350.0*keV));
    
    ExcitEnergies.push_back(10619.0*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(fPlanck/(60.0*keV));
}

G4B8GEMProbability::~G4B8GEMProbability()
{}
