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
//      File name:     G4XnpElasticLowE
//
//      Author:        
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// Neutron-Proton elastic cross section 
// Linear interpolation from data table
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4SystemOfUnits.hh"
#include "G4XnpElasticLowE.hh"
#include "G4KineticTrack.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLnVector.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

const G4double G4XnpElasticLowE::_lowLimit = 0.;
const G4double G4XnpElasticLowE::_highLimit = 3.*GeV;

// Low energy limit of the cross-section table (in GeV)
// Units are assigned when filling the PhysicsVector
const G4double G4XnpElasticLowE::_eMinTable = 1.8964808;
const G4double G4XnpElasticLowE::_eStepLog = 0.01;

// Cross-sections in mb
// Units are assigned when filling the PhysicsVector
const G4int G4XnpElasticLowE::_tableSize = 101;
const G4double G4XnpElasticLowE::_sigmaTable[101] = 
{ 
  1500.0,  //GF  was 0.0
  248.20, 93.38, 55.26, 44.50, 41.33, 38.48, 37.20, 35.98,
  35.02, 34.47, 32.48, 30.76, 29.46, 28.53, 27.84, 27.20,
  26.53, 25.95, 25.59, 25.46, 25.00, 24.49, 24.08, 23.86,
  23.17, 22.70, 21.88, 21.48, 20.22, 19.75, 18.97, 18.39,
  17.98, 17.63, 17.21, 16.72, 16.68, 16.58, 16.42, 16.22,
  15.98, 15.71, 15.42, 15.14, 14.87, 14.65, 14.44, 14.26,
  14.10, 13.95, 13.80, 13.64, 13.47, 13.29, 13.09, 12.89,
  12.68, 12.47, 12.27, 12.06, 11.84, 11.76, 11.69, 11.60,
  11.50, 11.41, 11.29, 11.17, 11.06, 10.93, 10.81, 10.68,
  10.56, 10.44, 10.33, 10.21, 10.12, 10.03,  9.96,  9.89,
  9.83,  9.80,  9.77,  9.75,  9.74,  9.74,  9.74,  9.76,
  9.73,  9.70,  9.68,  9.65,  9.63,  9.60,  9.57,  9.55,
  9.52,  9.49,  9.46,  9.43
};


G4XnpElasticLowE::G4XnpElasticLowE() 
{ 
  // Cross-sections are available in the range (_eMin,_eMax)

  _eMin = _eMinTable * GeV;
  _eMin = G4Exp(G4Log(_eMinTable)-_eStepLog)*GeV;
  _eMax = G4Exp(G4Log(_eMinTable) + _tableSize * _eStepLog) * GeV;

  // Protections: validity limits must be compatible with available data

  if (_eMin < _lowLimit)
    throw G4HadronicException(__FILE__, __LINE__, "G4XnpElasticLowE::G4XnpElasticLowE - Low energy limit not valid");
    
  if (_highLimit > _eMax)
    throw G4HadronicException(__FILE__, __LINE__, "G4XnpElasticLowE::G4XnpElasticLowE - High energy limit not valid");
    
  _sigma = new G4PhysicsLnVector(_eMin,_eMax,_tableSize);
  G4int i;
  for (i=0; i<_tableSize; i++)
    {
      G4double value = _sigmaTable[i] * millibarn;
      _sigma->PutValue(i,value);
    }
}


G4XnpElasticLowE::~G4XnpElasticLowE()
{
   delete _sigma;
   _sigma = 0;
}


G4bool G4XnpElasticLowE::operator==(const G4XnpElasticLowE &right) const
{
  return (this == (G4XnpElasticLowE *) &right);
}


G4bool G4XnpElasticLowE::operator!=(const G4XnpElasticLowE &right) const
{
  return (this != (G4XnpElasticLowE *) &right);
}


G4double G4XnpElasticLowE::CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const
{
  G4double sigma = 0.;
  G4double sqrtS = (trk1.Get4Momentum() + trk2.Get4Momentum()).mag();
  G4bool dummy = false;

  const G4ParticleDefinition* proton = G4Proton::ProtonDefinition();
  const G4ParticleDefinition* neutron = G4Neutron::NeutronDefinition();

  const G4ParticleDefinition* def1 = trk1.GetDefinition();
  const G4ParticleDefinition* def2 = trk2.GetDefinition();
  if ( (def1 == proton && def2 == neutron) ||
       (def1 == neutron && def2 == proton) )
  {
      if (sqrtS >= _eMin && sqrtS <= _eMax)
      {
	  sigma = _sigma->GetValue(sqrtS,dummy);
      } else if ( sqrtS < _eMin )
      {
          sigma = _sigma->GetValue(_eMin,dummy);
      }
  }

  return sigma;
}

void G4XnpElasticLowE::Print() const
{
  // Dump the cross-section table
  G4cout << Name() << "Cross-section table: " << G4endl;
  G4bool dummy = false;
  G4int i;

  for (i=0; i<_tableSize; i++)
    {
      G4double e = _sigma->GetLowEdgeEnergy(i) / GeV;
      G4double sigma = _sigma->GetValue(e,dummy) / millibarn;
      G4cout << i << ") e = " << e << " GeV ---- Cross section = " << sigma << " mb " << G4endl;
    }

  G4VCrossSectionSource::Print();
}


G4String G4XnpElasticLowE::Name() const
{
  G4String name("npElasticLowE");
  return name;
}


G4bool G4XnpElasticLowE::IsValid(G4double e) const
{
  G4bool answer = InLimits(e,_lowLimit,_highLimit);

  return answer;
}


