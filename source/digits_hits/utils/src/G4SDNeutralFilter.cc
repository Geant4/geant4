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
// $Id: G4SDNeutralFilter.cc,v 1.1 2005/11/16 23:04:04 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// G4VSensitiveDetector
#include "G4SDNeutralFilter.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

////////////////////////////////////////////////////////////////////////////////
// class description:
//
//  This is the class of a filter to be associated with a
// sensitive detector. 
//
// Created: 2005-11-14  Tsukasa ASO.
// 
///////////////////////////////////////////////////////////////////////////////

G4SDNeutralFilter::G4SDNeutralFilter(G4String name)
  :G4VSDFilter(name)
{;}

G4SDNeutralFilter::~G4SDNeutralFilter()
{;}

G4bool G4SDNeutralFilter::Accept(const G4Step* aStep) const
{
  if (aStep->GetPreStepPoint()->GetCharge()== 0. ) return TRUE;
  return FALSE;
}

