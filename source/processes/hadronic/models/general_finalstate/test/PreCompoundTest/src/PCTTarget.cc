#include "PCTTarget.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include "PCTOptions.hh"


PCTTarget::PCTTarget(const PCTOptions& opt)
  : theParticle(0), theMaterial(0), lastIsotope(-1)
{
  if ( opt.IsNaturalTarget() ) this->Initialize(opt.GetTargetMaterial(),opt.GetTargetZ());
  else this->Initialize(opt.GetTargetA(),opt.GetTargetZ());
}

void PCTTarget::Initialize(const G4int a, const G4int z)
{
  theParticle = const_cast<G4IonTable*>(G4ParticleTable::GetParticleTable()->GetIonTable())
    ->GetIon(z,a);
  if (!theParticle)
    {
      G4cerr << "PCTTarget::PCTTarget(): Can not find particle A = "
	     << a << " Z = " << z << '\n';
      theParticle = 0;
    }
  else 
    {
      ParticleMomentum.setPx(0.0);
      ParticleMomentum.setPy(0.0);
      ParticleMomentum.setPy(0.0);
      G4double mass = theParticle->GetPDGMass();
      ParticleMomentum.setE(mass);
    }
}

void PCTTarget::Initialize(const std::vector< std::pair<G4int,G4double> >& mat, const G4int z)
{
    theMaterial = new std::vector< std::pair<G4double,G4ParticleDefinition*> >;
    theMaterial->reserve(mat.size());
    MaterialMomentum.reserve(mat.size());
    std::vector< std::pair<G4int,G4double> >::const_iterator it;
    G4double norma = mat.back().second;
    for (it = mat.begin(); it != mat.end(); ++it)
    {
	G4int a = it->first;
	G4double pcent = it->second;
	G4ParticleDefinition * part =  const_cast<G4IonTable*>(G4ParticleTable::GetParticleTable()->GetIonTable())
	  ->GetIon(z,a);
	if (!part)
	{
	    G4cerr << "PCTTarget::PCTTarget(): Can not find ion A = "
		   << a << " Z = " << z << '\n';
	}
	else 
	{
	    MaterialMomentum.push_back(G4LorentzVector(0.0,0.0,0.0,part->GetPDGMass()));
	    theMaterial->push_back(std::make_pair(pcent/norma,part));
	}
    }
    if (theMaterial->size() != MaterialMomentum.size() || theMaterial->size() != mat.size())
    {
	delete theMaterial;
	theMaterial = 0;
    } 
}


PCTTarget::~PCTTarget()
{
    if (theMaterial) delete theMaterial;
}


void PCTTarget::SetParticleMomentum(const G4ThreeVector& p)
{
    ParticleMomentum.setVect(p);
    G4double mass = theParticle->GetPDGMass();
    ParticleMomentum.setE(sqrt(p.mag2()+mass*mass));
    return;
}

void PCTTarget::SetMaterialMomentum(const G4ThreeVector& p)
{
    G4double pmag2 = p.mag2();
    std::vector< G4LorentzVector >::iterator it;
    for (it = MaterialMomentum.begin(); it != MaterialMomentum.end(); ++it)
    {
#ifdef G4NO_ISO_VECDIST
	std::vector<G4LorentzVector>::difference_type n = 0;
	std::distance(MaterialMomentum.begin(),it,n);
#else         
	G4int n = std::distance(MaterialMomentum.begin(),it);
#endif
	it->setVect(p);
	G4double mass = (theMaterial->operator[](n)).second->GetPDGMass();
	it->setVect(sqrt(pmag2+mass*mass));
    }
    return;
}

G4int PCTTarget::GetZ() const
{
    if (theParticle)
    {
	return G4int(theParticle->GetPDGCharge()+0.5);
    }
    else if (theMaterial)
    {
	return G4int((theMaterial->begin()->second)->GetPDGCharge()+0.5);
    }
    return -1;
}


G4int PCTTarget::GetA() 
{
  if (theParticle)
    {
      return theParticle->GetBaryonNumber();
    }
  else if (theMaterial)
    {
      G4double rnd = G4UniformRand();
      std::vector< std::pair<G4double,G4ParticleDefinition*> >::iterator it;
      for (it = theMaterial->begin(); it != theMaterial->end(); ++it)
	{
	  if (rnd < it->first) 
	    {
#ifdef G4NO_ISO_VECDIST
	      std::vector<G4LorentzVector>::difference_type n = 0;
	      std::distance(theMaterial->begin(),it,n);
#else         
	      G4int n = std::distance(theMaterial->begin(),it);
#endif
	      lastIsotope = n;
	      return (it->second)->GetBaryonNumber();
	    }
	}
    }
  return -1;
}
