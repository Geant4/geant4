#include "expdata.h"

ClassImp(expdata);


expdata::expdata(const expdata & v) 
{
  targetA = v.targetA;
  targetZ = v.targetZ;
  targetSymbol = v. targetSymbol;
  projectileA = v.projectileA;
  projectileZ = v.projectileZ;
  projectileEnergy = v.projectileEnergy;
  projectileSymbol = v.projectileSymbol;;
  entry = v.entry;
  reaction = v.reaction;
  particleA = v.particleA;
  particleZ = v.particleZ;
  particleSymbol = v.particleSymbol;
  CM = v.CM;
  cutoff = v.cutoff;
}

const expdata & expdata::operator=(const expdata & v)
{
  targetA = v.targetA;
  targetZ = v.targetZ;
  targetSymbol = v. targetSymbol;
  projectileA = v.projectileA;
  projectileZ = v.projectileZ;
  projectileEnergy = v.projectileEnergy;
  projectileSymbol = v.projectileSymbol;;
  entry = v.entry;
  reaction = v.reaction;
  particleA = v.particleA;
  particleZ = v.particleZ;
  particleSymbol = v.particleSymbol;
  CM = v.CM;
  cutoff = v.cutoff;
  return *this;
}

  void expdata::ShowYourSelf(std::ostringstream & os) const 
{
  // Target data
  os << "+-------------+\n"
     << "|   Target    |  " << targetSymbol << "(" << targetA 
     << "," << targetZ << ")\n"
     << "+-------------+\n\n";
  // Projectile data
  os << "+-------------+\n"
     << "|  Projectile |  " << projectileSymbol << "(" << projectileA 
     << "," << projectileZ << ") at " << projectileEnergy << " MeV" << '\n'
     << "+-------------+\n\n";
  // Ejectile data
  os << "+-------------+\n"
     << "|   Ejectile  |  " << particleSymbol << "(" << particleA 
     << "," << particleZ << ") cross section";
  if (this->IsCM()) os << " in CM";
  if (this->ThereIsCutoff()) os << " cutoff = " << this->GetCutoff() << " MeV";
  os << '\n';
  os << "+-------------+\n\n";
  // Experimental data  
  
  os << "Data comes from EXFOR database:\n\n" 
     << '\t' << "ENTRY:    " << entry << '\n'
     << '\t' << "REACTION: " << reaction << '\n';
  return;
}
