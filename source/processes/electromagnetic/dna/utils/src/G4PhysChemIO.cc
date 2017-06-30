/*
 * G4PhysChemIO.cc
 *
 *  Created on: 3 févr. 2017
 *      Author: matkara
 */

#include "G4PhysChemIO.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VAnalysisManager.hh"

using namespace std;

//------------------------------------------------------------------------------

namespace G4PhysChemIO{
  
FormattedText::FormattedText(){
  fRunID = -1;
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

void FormattedText::AddEmptyLineInOuputFile()
{
  if(fFileInitialized) fOfstream << G4endl;
}

//------------------------------------------------------------------------------

void FormattedText::CloseFile()
{
  if (fFileInitialized == false) return;
  
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
  
  if (finalPosition != 0)
  {
    fOfstream << setw(14) << (finalPosition->x()) / nanometer << setw(14)
    << (finalPosition->y()) / nanometer << setw(14)
    << (finalPosition->z()) / nanometer;
  }
  
  fOfstream << G4endl;
}

//------------------------------------------------------------------------------
//
// Using G4analysis
//

G4Analysis::G4Analysis(G4VAnalysisManager* analysisManager):
fpAnalysisManager(analysisManager)
{
  fFileInitialized = false;
  fNtupleID = -1;
}

//------------------------------------------------------------------------------

G4Analysis::~G4Analysis()
{
  fpAnalysisManager = 0;
}

//------------------------------------------------------------------------------

void G4Analysis::InitializeFile()
{
  if (fFileInitialized) return;
  
  fNtupleID = fpAnalysisManager->CreateNtuple("PhysChem","PhysChem");
  fpAnalysisManager->CreateNtupleIColumn(fNtupleID, "ParentID");
  fpAnalysisManager->CreateNtupleSColumn(fNtupleID, "Molecule");
  
  //----------------------------------------------------------------------------
  // valid for H2O only
  fpAnalysisManager->CreateNtupleIColumn(fNtupleID, "ElectronicModif");
  // ionization = 0 / excitation = 1 / diss att = 2
  fpAnalysisManager->CreateNtupleIColumn(fNtupleID, "level");
  // valid for ion and exc only
  fpAnalysisManager->CreateNtupleDColumn(fNtupleID, "Energy_eV");
  // valid for ion and exc only
  
  //----------------------------------------------------------------------------
  fpAnalysisManager->CreateNtupleDColumn(fNtupleID, "x_parent_nm");
  fpAnalysisManager->CreateNtupleDColumn(fNtupleID, "y_parent_nm");
  fpAnalysisManager->CreateNtupleDColumn(fNtupleID, "z_parent_nm");
  fpAnalysisManager->CreateNtupleDColumn(fNtupleID, "x_nm");
  fpAnalysisManager->CreateNtupleDColumn(fNtupleID, "y_nm");
  fpAnalysisManager->CreateNtupleDColumn(fNtupleID, "z_nm");
  fpAnalysisManager->FinishNtuple(fNtupleID);
  
  fFileInitialized = true;
}

//------------------------------------------------------------------------------

void G4Analysis::WriteInto(const G4String& output,
                           ios_base::openmode)
{
  fpAnalysisManager->OpenFile(output);
  fFileInitialized = false;
}

//------------------------------------------------------------------------------

void G4Analysis::CloseFile()
{
//  fpAnalysisManager->Write();
//  fpAnalysisManager->CloseFile();
}

//------------------------------------------------------------------------------

void G4Analysis::CreateWaterMolecule(G4int modification,
                                     G4int electronicLevel,
                                     G4double energy,
                                     const G4Track* theIncomingTrack)
{
  if(!fFileInitialized) InitializeFile();
  
  // parent ID
  fpAnalysisManager->FillNtupleIColumn(fNtupleID, 0,
                                       theIncomingTrack->GetTrackID());
  
  // molecule type
  fpAnalysisManager->FillNtupleSColumn(fNtupleID, 1, "H2O");
  
  //----------------------------------------------------------------------------
  // valid for H2O only
  
  // electronic modif
  fpAnalysisManager->FillNtupleIColumn(fNtupleID, 2, modification);
  // ionization = 0 / excitation = 1 / diss att = 2
  fpAnalysisManager->FillNtupleIColumn(fNtupleID, 3, electronicLevel);
  fpAnalysisManager->FillNtupleDColumn(fNtupleID, 4, energy / eV);
  
  //----------------------------------------------------------------------------
  const G4ThreeVector& parentPos = theIncomingTrack->GetPosition();
  
  fpAnalysisManager->FillNtupleDColumn(fNtupleID,5,(parentPos.x())/nanometer);
  fpAnalysisManager->FillNtupleDColumn(fNtupleID,6,(parentPos.y())/nanometer);
  fpAnalysisManager->FillNtupleDColumn(fNtupleID,7,(parentPos.z())/nanometer);
  
  fpAnalysisManager->FillNtupleDColumn(fNtupleID,8,(parentPos.x())/nanometer);
  fpAnalysisManager->FillNtupleDColumn(fNtupleID,9,(parentPos.y())/nanometer);
  fpAnalysisManager->FillNtupleDColumn(fNtupleID,10,(parentPos.z())/nanometer);
  fpAnalysisManager->AddNtupleRow(fNtupleID);
}

//------------------------------------------------------------------------------

void G4Analysis::CreateSolvatedElectron(const G4Track* electronTrack,
                                        G4ThreeVector* finalPosition)
{
  if(!fFileInitialized) InitializeFile();
  
  // parent ID
  fpAnalysisManager->FillNtupleIColumn(fNtupleID, 0,
                                       electronTrack->GetTrackID());
  
  // molecule type
  fpAnalysisManager->FillNtupleSColumn(fNtupleID, 1, "e_aq");
  
  //----------------------------------------------------------------------------
  // valid for H2O only
  
  // electronic modif
  fpAnalysisManager->FillNtupleIColumn(fNtupleID, 2, -1); // electronic modif
  fpAnalysisManager->FillNtupleIColumn(fNtupleID, 3, -1); // electronic level
  fpAnalysisManager->FillNtupleDColumn(fNtupleID, 4,
                                       electronTrack->GetKineticEnergy() / eV);
  
  //----------------------------------------------------------------------------
  const G4ThreeVector& parentPos = electronTrack->GetPosition();
  const double i_nm = 1./nanometer;
  
  fpAnalysisManager->FillNtupleDColumn(fNtupleID,5, parentPos.x() *i_nm);
  fpAnalysisManager->FillNtupleDColumn(fNtupleID,6, parentPos.y() *i_nm);
  fpAnalysisManager->FillNtupleDColumn(fNtupleID,7, parentPos.z() *i_nm);
  
  if (finalPosition != 0)
  {
    fpAnalysisManager->FillNtupleDColumn(fNtupleID,8, finalPosition->x()*i_nm);
    fpAnalysisManager->FillNtupleDColumn(fNtupleID,9, finalPosition->y()*i_nm);
    fpAnalysisManager->FillNtupleDColumn(fNtupleID,10, finalPosition->z()*i_nm);
  }
  else
  {
    fpAnalysisManager->FillNtupleDColumn(fNtupleID,8, parentPos.x() *i_nm);
    fpAnalysisManager->FillNtupleDColumn(fNtupleID,9, parentPos.y() *i_nm);
    fpAnalysisManager->FillNtupleDColumn(fNtupleID,10, parentPos.z() *i_nm);
  }
  
  fpAnalysisManager->AddNtupleRow(fNtupleID);
}
}
