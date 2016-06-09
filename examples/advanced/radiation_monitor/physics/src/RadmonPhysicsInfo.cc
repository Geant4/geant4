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
// File name:     RadmonPhysicsInfo.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsInfo.cc,v 1.3 2006/06/29 16:18:51 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//

#include "RadmonPhysicsInfo.hh"
#include "G4UnitsTable.hh"

                                                RadmonPhysicsInfo :: RadmonPhysicsInfo(const RadmonPhysicsInfo & copy)
:
 minEnergy(copy.minEnergy),
 maxEnergy(copy.maxEnergy),
 name(copy.name),
 particle(copy.particle)
{
}
 



  
RadmonPhysicsInfo &                             RadmonPhysicsInfo :: operator=(const RadmonPhysicsInfo & copy)
{
 minEnergy=copy.minEnergy;
 maxEnergy=copy.maxEnergy;
 name=copy.name;
 particle=copy.particle;
 
 return (*this);
}



   

G4bool                                          RadmonPhysicsInfo :: CollidesWith(const RadmonPhysicsInfo & other) const
{
 if (name.empty() || particle==0 || minEnergy>maxEnergy ||
     other.name.empty() || other.particle==0 || other.minEnergy>other.maxEnergy)
  return false;
  
 return (name==other.name && particle->GetParticleName()==other.particle->GetParticleName() && maxEnergy>=other.minEnergy && minEnergy<=other.maxEnergy);
}



   
 
std::ostream &                                  operator<<(std::ostream & out, const RadmonPhysicsInfo & info)
{
 if (info.GetProcessName().empty())
  out << "UNDEFINED";
 else 
  out << info.GetProcessName();
 
 if (info.GetParticleDefinition()==0)
  out << "(UNDEFINED)";
 else
  out << info.GetParticleDefinition()->GetParticleName();
 
 if (info.GetMinEnergy()>info.GetMaxEnergy())
  out << " [UNDEFINED]";
 else if (info.GetMinEnergy()==info.GetMaxEnergy())
 {
  if (info.GetMinEnergy()==0.*eV)
   out << " [REST]";
  else
   out << " [" << G4BestUnit(info.GetMinEnergy(), "Energy") << ']';
 }
 else
  out << " [" << G4BestUnit(info.GetMinEnergy(), "Energy") << ", " << G4BestUnit(info.GetMaxEnergy(), "Energy") << ']';

 return out;
}
 
