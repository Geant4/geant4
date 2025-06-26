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
/*
 * G4PhysChemIO.cc
 *
 *  Created on: 3 févr. 2017
 *      Author: matkara
 */

#include "G4PhysChemIO.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

using namespace std;

//------------------------------------------------------------------------------

namespace G4PhysChemIO{
  
FormattedText::FormattedText(){
  fEventID = -1;
  fFileInitialized = false;
}
  
//------------------------------------------------------------------------------
  
FormattedText::~FormattedText(){
  CloseFile();
}
  
//------------------------------------------------------------------------------
  
void FormattedText::InitializeFile()
{
  if(fFileInitialized) return;
  
  fOfstream << std::setprecision(6) << std::scientific;
  fOfstream << setw(11) << left << "#Parent ID" << setw(10) << "Molecule"
  << setw(14) << "Elec Modif" << setw(13) << "Energy (eV)"
  << setw(22) << "X pos of parent [nm]" << setw(22)
  << "Y pos of parent [nm]" << setw(22) << "Z pos of parent [nm]"
  << setw(14) << "X pos [nm]" << setw(14) << "Y pos [nm]"
  << setw(14) << "Z pos [nm]" << G4endl<< setw(21) << "#"
  << setw(13) << "1)io/ex=0/1"
  << G4endl
  << setw(21) << "#"
  << setw(13) << "2)level=0...5"
  << G4endl;
  
  fFileInitialized = true;
}

//------------------------------------------------------------------------------

void FormattedText::WriteInto(const G4String& output,
                              ios_base::openmode mode)
{
  fOfstream.open(output.data(), mode);
  fFileInitialized = false;
}

//------------------------------------------------------------------------------

void FormattedText::AddEmptyLineInOutputFile()
{
  if(fFileInitialized) fOfstream << G4endl;
}

//------------------------------------------------------------------------------

void FormattedText::CloseFile()
{
  if (!fFileInitialized) return;
  
  if (fOfstream.is_open())
  {
    fOfstream.close();
  }
}

//------------------------------------------------------------------------------

void FormattedText::CreateWaterMolecule(G4int modification,
                                        G4int electronicLevel,
                                        G4double energy,
                                        const G4Track* theIncomingTrack)
{
  if(!fFileInitialized) InitializeFile();
  
  fOfstream << setw(11) << left << theIncomingTrack->GetTrackID()
  << setw(10) << "H2O" << left << modification << internal
  << ":" << right << electronicLevel << left << setw(11) << ""
  << std::setprecision(2) << std::fixed << setw(13)
  << energy / eV << std::setprecision(6) << std::scientific
  << setw(22)
  << (theIncomingTrack->GetPosition().x()) / nanometer
  << setw(22)
  << (theIncomingTrack->GetPosition().y()) / nanometer
  << setw(22)
  << (theIncomingTrack->GetPosition().z()) / nanometer
  << G4endl;
}

//------------------------------------------------------------------------------

void FormattedText::CreateSolvatedElectron(const G4Track* theIncomingTrack,
                                           G4ThreeVector* finalPosition)
{
  if(!fFileInitialized) InitializeFile();
  
  fOfstream << setw(11) << theIncomingTrack->GetTrackID() << setw(10)
  << "e_aq" << setw(14) << -1 << std::setprecision(2)
  << std::fixed << setw(13)
  << theIncomingTrack->GetKineticEnergy() / eV
  << std::setprecision(6) << std::scientific << setw(22)
  << (theIncomingTrack->GetPosition().x()) / nanometer
  << setw(22)
  << (theIncomingTrack->GetPosition().y()) / nanometer
  << setw(22)
  << (theIncomingTrack->GetPosition().z()) / nanometer;
  
  if (finalPosition != nullptr)
  {
    fOfstream << setw(14) << (finalPosition->x()) / nanometer << setw(14)
    << (finalPosition->y()) / nanometer << setw(14)
    << (finalPosition->z()) / nanometer;
  }
  
  fOfstream << G4endl;
}

}
