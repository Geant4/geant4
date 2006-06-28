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
// File name:     RadmonPhysicsInfo.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsInfo.cc,v 1.2 2006-06-28 13:56:17 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
 
