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
//      File name:     G4XnpTotalLowE
//
//      Author:        
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// Neutron-Proton total  cross section 
// Linear interpolation from data table
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4SystemOfUnits.hh"
#include "G4XnpTotalLowE.hh"
#include "G4KineticTrack.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLnVector.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

const G4double G4XnpTotalLowE::_lowLimit = 0.;
const G4double G4XnpTotalLowE::_highLimit = 3.*GeV;

// Low energy limit of the cross-section table (in GeV)
// Units are assigned when filling the PhysicsVector
const G4double G4XnpTotalLowE::_eMinTable = 1.8964808;
const G4double G4XnpTotalLowE::_eStepLog = 0.01;

// Cross-sections in mb
// Units are assigned when filling the PhysicsVector
const G4int G4XnpTotalLowE::_tableSize = 101;
const G4double G4XnpTotalLowE::_sigmaTable[101] = 
{ 
  1500.0,  
  248.20, 93.38, 55.26, 44.50, 41.33, 38.48, 37.20, 35.98,
  35.02, 34.47, 34.37, 34.67, 35.23, 35.97, 36.75, 37.37,
  37.77, 38.03, 38.40, 38.83, 39.26, 39.67, 40.06, 40.45,
  40.79, 41.06, 41.31, 41.52, 41.70, 41.81, 41.87, 41.98,
  42.12, 42.29, 42.55, 42.82, 43.01, 43.12, 43.16, 43.14,
  43.06, 42.95, 42.81, 42.67, 42.54, 42.45, 42.38, 42.33,
  42.30, 42.29, 42.28, 42.26, 42.24, 42.21, 42.17, 42.14,
  42.10, 42.07, 42.06, 42.05, 42.04, 42.03, 42.02, 42.00,
  41.97, 41.94, 41.89, 41.84, 41.79, 41.73, 41.67, 41.61,
  41.55, 41.49, 41.44, 41.38, 41.34, 41.31, 41.29, 41.28,
  41.27, 41.28, 41.30, 41.33, 41.36, 41.40, 41.44, 41.49,
  41.50, 41.51, 41.51, 41.51, 41.52, 41.51, 41.51, 41.50,
  41.50, 41.49, 41.47, 41.46
};


G4XnpTotalLowE::G4XnpTotalLowE() 
{ 
  // Cross-sections are available in the range (_eMin,_eMax)

  _eMin = _eMinTable * GeV;
  _eMin = G4Exp(G4Log(_eMinTable)-_eStepLog)*GeV;
  _eMax = G4Exp(G4Log(_eMinTable) + _tableSize * _eStepLog) * GeV;

  // Protections: validity limits must be compatible with available data
//  @@GF  this ought to be _lowLimit < _eMin
  if (_eMin < _lowLimit)
    throw G4HadronicException(__FILE__, __LINE__, "G4XnpTotalLowE::G4XnpTotalLowE - Low energy limit not valid");
    
  if (_highLimit > _eMax)
    throw G4HadronicException(__FILE__, __LINE__, "G4XnpTotalLowE::G4XnpTotalLowE - High energy limit not valid");
    
  _sigma = new G4PhysicsLnVector(_eMin,_eMax,_tableSize);
  G4int i;
  for (i=0; i<_tableSize; i++)
    {
      G4double value = _sigmaTable[i] * millibarn;
      _sigma->PutValue(i,value);
    }
}


G4XnpTotalLowE::~G4XnpTotalLowE()
{
  if (_sigma) delete _sigma;
  _sigma=0;
}


G4bool G4XnpTotalLowE::operator==(const G4XnpTotalLowE &right) const
{
  return (this == (G4XnpTotalLowE *) &right);
}


G4bool G4XnpTotalLowE::operator!=(const G4XnpTotalLowE &right) const
{
  return (this != (G4XnpTotalLowE *) &right);
}


G4double G4XnpTotalLowE::CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const
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

void G4XnpTotalLowE::Print() const
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


G4String G4XnpTotalLowE::Name() const
{
  G4String name("NNTotalLowE");
  return name;
}


G4bool G4XnpTotalLowE::IsValid(G4double e) const
{
  G4bool answer = InLimits(e,_lowLimit,_highLimit);

  return answer;
}


