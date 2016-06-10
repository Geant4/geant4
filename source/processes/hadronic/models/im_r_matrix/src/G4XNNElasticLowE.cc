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
#include "globals.hh"
#include "G4ios.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4SystemOfUnits.hh"
#include "G4XNNElasticLowE.hh"
#include "G4KineticTrack.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLnVector.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

const G4double G4XNNElasticLowE::_lowLimit = 0.;
const G4double G4XNNElasticLowE::_highLimit = 3.*GeV;

// Low energy limit of the cross-section table (in GeV)
// Units are assigned when filling the PhysicsVector
const G4double G4XNNElasticLowE::_eMinTable = 1.8964808;
const G4double G4XNNElasticLowE::_eStepLog = 0.01;

// Cross-sections in mb
// Units are assigned when filling the PhysicsVector

const G4int G4XNNElasticLowE::tableSize = 101;

const G4double G4XNNElasticLowE::ppTable[101] = 
{ 
  60.00, //was 0.
  33.48, 26.76, 25.26, 24.55, 23.94, 23.77, 23.72, 23.98,
  25.48, 27.52, 27.72, 27.21, 25.80, 26.00, 24.32, 23.81,
  24.37, 24.36, 23.13, 22.43, 21.71, 21.01, 20.83, 20.74,
  20.25, 20.10, 20.59, 20.04, 20.83, 20.84, 21.07, 20.83,
  20.79, 21.88, 21.15, 20.92, 19.00, 18.60, 17.30, 17.00,
  16.70, 16.50, 16.20, 15.80, 15.57, 15.20, 15.00, 14.60,
  14.20, 14.00, 13.80, 13.60, 13.40, 13.20, 13.00, 12.85,
  12.70, 12.60, 12.50, 12.40, 12.30, 12.20, 12.10, 12.00,
  11.90, 11.80, 11.75, 11.70, 11.64, 11.53, 11.41, 11.31,
  11.22, 11.13, 11.05, 10.97, 10.89, 10.82, 10.75, 10.68,
  10.61, 10.54, 10.48, 10.41, 10.35, 10.28, 10.22, 10.16,
  10.13, 10.10, 10.08, 10.05, 10.02,  9.99,  9.96,  9.93,
  9.90,  9.87,  9.84,  9.80       
};

const G4double G4XNNElasticLowE::npTable[101] = 
{ 
  1500.00, // was 0.
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


G4XNNElasticLowE::G4XNNElasticLowE() 
{ 
  // Cross-sections are available in the range (_eMin,_eMax)

  _eMin = _eMinTable * GeV;
  _eMax = G4Exp(G4Log(_eMinTable) + tableSize * _eStepLog) * GeV;
  if (_eMin < _lowLimit)
    throw G4HadronicException(__FILE__, __LINE__, "G4XNNElasticLowE::G4XNNElasticLowE - Low energy limit not valid");    
  if (_highLimit > _eMax)
    throw G4HadronicException(__FILE__, __LINE__, "G4XNNElasticLowE::G4XNNElasticLowE - High energy limit not valid");    
  G4PhysicsVector* pp = new G4PhysicsLnVector(_eMin,_eMax,tableSize);

  _eMin = G4Exp(G4Log(_eMinTable)-_eStepLog)*GeV;
  if (_eMin < _lowLimit)
    throw G4HadronicException(__FILE__, __LINE__, "G4XNNElasticLowE::G4XNNElasticLowE - Low energy limit not valid");
  G4PhysicsVector* np = new G4PhysicsLnVector(_eMin,_eMax,tableSize);

  G4int i;
  for (i=0; i<tableSize; i++)
    {
      G4double value = ppTable[i] * millibarn;
      pp->PutValue(i,value);
      value = npTable[i] * millibarn;
      np->PutValue(i,value);
    }
  xMap[G4Proton::ProtonDefinition()] = pp;
  xMap[G4Neutron::NeutronDefinition()] = np;
}


G4XNNElasticLowE::~G4XNNElasticLowE()
{
  delete xMap[G4Proton::ProtonDefinition()];
  delete xMap[G4Neutron::NeutronDefinition()];
}


G4bool G4XNNElasticLowE::operator==(const G4XNNElasticLowE &right) const
{
  return (this == (G4XNNElasticLowE *) &right);
}


G4bool G4XNNElasticLowE::operator!=(const G4XNNElasticLowE &right) const
{

  return (this != (G4XNNElasticLowE *) &right);
}


G4double G4XNNElasticLowE::CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const  
{
  G4double sigma = 0.;
  G4double sqrtS = (trk1.Get4Momentum() + trk2.Get4Momentum()).mag();
  G4bool dummy = false;

  const G4ParticleDefinition * key = FindKeyParticle(trk1,trk2);

  typedef std::map <const G4ParticleDefinition *, G4PhysicsVector*, std::less<const G4ParticleDefinition *> > StringPhysMap;

  if (xMap.find(key)!= xMap.end())
    {

      StringPhysMap::const_iterator iter;
      for (iter = xMap.begin(); iter != xMap.end(); ++iter)
	{
	  const G4ParticleDefinition * str = (*iter).first;
          if (str == key)
	    {
	      G4PhysicsVector* physVector = (*iter).second; 
	      //     G4PhysicsVector* physVector = xMap[key];
	      if (sqrtS >= _eMin && sqrtS <= _eMax)
		{
		  sigma = physVector->GetValue(sqrtS,dummy);
		} else if ( sqrtS < _eMin )
                {
                  sigma = physVector->GetValue(_eMin,dummy);
                }
		//G4cout << " sqrtS / sigma " << sqrtS/GeV << " / " <<
		//          sigma/millibarn << G4endl;
	    }
	}
    }
  return sigma;
}


void G4XNNElasticLowE::Print() const
{
  // Dump the pp cross-section table

  G4cout << Name() << ", pp cross-section: " << G4endl;

  G4bool dummy = false;
  G4int i;
  const G4ParticleDefinition * key = G4Proton::ProtonDefinition();
  G4PhysicsVector* pp = 0;

  typedef std::map <const G4ParticleDefinition *, G4PhysicsVector*, std::less<const G4ParticleDefinition *> > StringPhysMap;
  StringPhysMap::const_iterator iter;

  for (iter = xMap.begin(); iter != xMap.end(); ++iter)
    {
      const G4ParticleDefinition * str = (*iter).first;
      if (str == key)
	{
	  pp = (*iter).second; 
	}
    }
  
  if (pp != 0)
    {	
      for (i=0; i<tableSize; i++)
	{
	  G4double e = pp->GetLowEdgeEnergy(i);
	  G4double sigma = pp->GetValue(e,dummy) / millibarn;
	  G4cout << i << ") e = " << e / GeV << " GeV ---- Cross section = " << sigma << " mb " << G4endl;
	}
    }
  
  // Dump the np cross-section table

  G4cout << Name() << ", np cross-section: " << G4endl;

  key = G4Neutron::NeutronDefinition();
  G4PhysicsVector* np = 0;
  for (iter = xMap.begin(); iter != xMap.end(); ++iter)
    {
      const G4ParticleDefinition * str = (*iter).first;
      if (str == key)
	{
	  np = (*iter).second; 
	}
    }
  
  //  G4PhysicsVector* np = xMap[G4Neutron::NeutronDefinition()->GetParticleName()];
  
  if (np != 0)
    {	
      for (i=0; i<tableSize; i++)
	{
	  G4double e = np->GetLowEdgeEnergy(i);
	  G4double sigma = np->GetValue(e,dummy) / millibarn;
	  G4cout << i << ") e = " << e / GeV << " GeV ---- Cross section = " << sigma << " mb " << G4endl;
	}
    }
  G4VCrossSectionSource::Print();
}


G4String G4XNNElasticLowE::Name() const
{
  G4String name("NNElasticLowE");
  return name;
}



G4bool G4XNNElasticLowE::IsValid(G4double e) const
{
  G4bool answer = InLimits(e,_lowLimit,_highLimit);

  return answer;
}


