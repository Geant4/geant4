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

#ifndef G4SDNeutralFilter_h
#define G4SDNeutralFilter_h 1

class G4Step;
class G4NeutralDefinition;
#include "globals.hh"
#include "G4VSDFilter.hh"

////////////////////////////////////////////////////////////////////////////////
// class description:
//
//  This is the class of a filter to be associated with a
// sensitive detector.
//
//  This filter accepts neutral tracks.
//
// Created: 2005-11-14  Tsukasa ASO.
//
///////////////////////////////////////////////////////////////////////////////

class G4SDNeutralFilter : public G4VSDFilter
{
 public:
  G4SDNeutralFilter(G4String name);
  ~G4SDNeutralFilter() override = default;

  G4bool Accept(const G4Step*) const override;
};

#endif
