//
// File name:     RadmonPhysicsInfo.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsInfo.cc,v 1.1 2005-11-10 08:14:10 capra Exp $
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
 
