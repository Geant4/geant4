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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPData.hh"
#include "G4PhysicsFreeVector.hh"

  G4ParticleHPData::G4ParticleHPData(G4ParticleDefinition* projectile )
    : theProjectile(projectile)
  {
    //const char* theDataDirVariable;
    if( projectile == G4Neutron::Neutron() ) {
        theDataDirVariable = "G4NEUTRONHPDATA";
    }else if( projectile == G4Proton::Proton() ) {
      theDataDirVariable = "G4PROTONHPDATA";
    }else if( projectile == G4Deuteron::Deuteron() ) {
      theDataDirVariable = "G4DEUTERONHPDATA";
    }else if( projectile == G4Triton::Triton() ) {
      theDataDirVariable = "G4TRITONHPDATA";
    }else if( projectile == G4He3::He3() ) {
      theDataDirVariable = "G4HE3HPDATA";
    }else if( projectile == G4Alpha::Alpha() ) {
      theDataDirVariable = "G4ALPHAHPDATA";
    }

    numEle = (G4int)G4Element::GetNumberOfElements();
    for ( G4int i=0 ; i<numEle ; ++i )
    {
      theData.push_back ( new G4ParticleHPElementData );
    }
    for (G4int i=0; i<numEle; ++i)
    {
      (*theData[i]).Init((*(G4Element::GetElementTable()))[i], projectile, theDataDirVariable);
    }
  }
  
  G4ParticleHPData::~G4ParticleHPData()
  {
     for (auto it = theData.cbegin() ; it != theData.cend() ; ++it) delete *it;
     theData.clear();
  }
  
  G4ParticleHPData * G4ParticleHPData::Instance(G4ParticleDefinition* projectile)
  {
    static G4ThreadLocal G4ParticleHPData *theCrossSectionData_G4MT_TLS_ = nullptr ;
    if ( !theCrossSectionData_G4MT_TLS_ ) theCrossSectionData_G4MT_TLS_ = new G4ParticleHPData(projectile);  
    G4ParticleHPData &theCrossSectionData = *theCrossSectionData_G4MT_TLS_;
    return &theCrossSectionData;
  } 

  G4PhysicsVector * G4ParticleHPData::DoPhysicsVector(G4ParticleHPVector * theVector)
  {
    G4int len = theVector->GetVectorLength();
    if(len==0) return new G4PhysicsFreeVector(0, 0., 0.);
    G4double emin = theVector->GetX(0);
    G4double emax = theVector->GetX(len-1);

    G4PhysicsFreeVector * theResult = new G4PhysicsFreeVector(len, emin, emax);
    for (G4int i=0; i<len; ++i)
    {
      theResult->PutValues(i, theVector->GetX(i), theVector->GetY(i));
    }
    return theResult;
  }

void G4ParticleHPData::addPhysicsVector()
{
   for ( G4int i = numEle; i < (G4int)G4Element::GetNumberOfElements() ; ++i )
   {
      theData.push_back ( new G4ParticleHPElementData );
      (*theData[i]).Init((*(G4Element::GetElementTable()))[i], theProjectile, theDataDirVariable);
   }
   numEle = (G4int)G4Element::GetNumberOfElements();
}
