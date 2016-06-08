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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//      For information related to this code contact:
//
//      File name:     G4XAqmElastic
//
//      Author:        
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// Additive Quark Model Elastic cross section 
// (H.J. Lipkin and F. Scheck, Phys.Rev. 16 (1966) 71
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include "G4XAqmElastic.hh"
#include "G4XAqmTotal.hh"
#include "G4KineticTrack.hh"
#include "G4ParticleDefinition.hh"


// Validity range of this cross-section
const G4double G4XAqmElastic::_lowLimit = 0.;
const G4double G4XAqmElastic::_highLimit = DBL_MAX;

G4XAqmElastic::G4XAqmElastic()
{ 
  // As a first approximation the model is assumed to be valid over 
  // the entire energy range
}


G4XAqmElastic::~G4XAqmElastic()
{ }


G4bool G4XAqmElastic::operator==(const G4XAqmElastic &right) const
{
  return (this == (G4XAqmElastic *) &right);
}


G4bool G4XAqmElastic::operator!=(const G4XAqmElastic &right) const
{
  return (this != (G4XAqmElastic *) &right);
}


G4double G4XAqmElastic::CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const
{
  G4double sigma = 0.;

  // Reference to be checked!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  const G4double coeff = 0.39;
  const G4double param = 1.5;

  G4XAqmTotal aqmTotal;

  G4double sigmaTot = aqmTotal.CrossSection(trk1,trk2);
  sigma = coeff * pow(sigmaTot,param);  

  // Verify that elastic cross section < total cross section

  if (sigma > sigmaTot) 
    G4Exception("G4XAqmElastic::CrossSection - elastic cross section greater than total");

  return sigma;
}


G4String G4XAqmElastic::Name() const
{
  G4String name("AqmElasticCrossSection");
  return name;
}


G4bool G4XAqmElastic::IsValid(G4double e) const
{
  G4bool answer = InLimits(e,_lowLimit,_highLimit);

  return answer;
}
