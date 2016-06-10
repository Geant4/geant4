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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//      For information related to this code contact:
//
//      File name:     G4XAqmTotal
//
//      Author:        
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// Additive Quark Model total cross section 
// (H.J. Lipkin and F. Scheck, Phys.Rev. 16 (1966) 71
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"
#include "G4XAqmTotal.hh"
#include "G4KineticTrack.hh"
#include "G4ParticleDefinition.hh"


// Validity range of this cross-section
const G4double G4XAqmTotal::_lowLimit = 0.;
const G4double G4XAqmTotal::_highLimit = DBL_MAX;

G4XAqmTotal::G4XAqmTotal()
{ 
  // As a first approximation the model is assumed to be valid over 
  // the entire energy range
}


G4XAqmTotal::~G4XAqmTotal()
{ }


G4bool G4XAqmTotal::operator==(const G4XAqmTotal &right) const
{
  return (this == (G4XAqmTotal *) &right);
}


G4bool G4XAqmTotal::operator!=(const G4XAqmTotal &right) const
{
  return (this != (G4XAqmTotal *) &right);
}


G4double G4XAqmTotal::CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const
{
  G4double sigma = 0.;

  // Get strangeness content
  const G4ParticleDefinition* def1 = trk1.GetDefinition();
  G4int sTrk1 = def1->GetQuarkContent(3) + def1->GetAntiQuarkContent(3);
  const G4ParticleDefinition* def2 = trk2.GetDefinition();
  G4int sTrk2 = def2->GetQuarkContent(3) + def2->GetAntiQuarkContent(3);
  
  // Get non-strange quark content
  G4int qTrk1 = def1->GetQuarkContent(1) + def1->GetAntiQuarkContent(1) +
                def1->GetQuarkContent(2) + def1->GetAntiQuarkContent(2) +
                def1->GetQuarkContent(4) + def1->GetAntiQuarkContent(4) +
                def1->GetQuarkContent(5) + def1->GetAntiQuarkContent(5) +
                def1->GetQuarkContent(6) + def1->GetAntiQuarkContent(6);

  G4int qTrk2 = def2->GetQuarkContent(1) + def2->GetAntiQuarkContent(1) +
                def2->GetQuarkContent(2) + def2->GetAntiQuarkContent(2) +
                def2->GetQuarkContent(4) + def2->GetAntiQuarkContent(4) +
                def2->GetQuarkContent(5) + def2->GetAntiQuarkContent(5) +
                def2->GetQuarkContent(6) + def2->GetAntiQuarkContent(6);

  G4double sRatio1 = 0.;
  if (qTrk1 != 0) sRatio1 = sTrk1 / qTrk1;
  
  G4double sRatio2 = 0.;
  if (qTrk2 != 0) sRatio2 = sTrk2 / qTrk2;
  
  // Calculate the number of colliding mesons
  G4int nMesons = 0;
  G4int nQ1 = sTrk1 + qTrk1;

  if (nQ1 == 2) nMesons++;
  G4int nQ2 = sTrk2 + qTrk2;
  if (nQ2 == 2) nMesons++;

  // Cross-section (units to be checked!)
  sigma = 40. * G4Pow::GetInstance()->powN((2.0/3.0),nMesons) * (1. - 0.4 * sRatio1) * (1. - 0.4 * sRatio2) * millibarn;

  return sigma;
}


G4String G4XAqmTotal::Name() const
{
  G4String name("AqmTotalCrossSection");
  return name;
}



G4bool G4XAqmTotal::IsValid(G4double e) const
{
  G4bool answer = InLimits(e,_lowLimit,_highLimit);

  return answer;
}
