#include "PCTProjectile.hh"


#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"

PCTProjectile::PCTProjectile(const G4String partName)
    : dir(NoDirection), momentum(0.0,0.0,0.0) 
{
    theParticle = G4ParticleTable::GetParticleTable()->FindParticle(partName);
    if (!theParticle)
    {
	G4cerr << "PCTProjectile::PCTProjectile(): Can not find particle "
	       << partName << '\n';
    }
}

PCTProjectile::PCTProjectile(const G4int a, const G4int z)
    : dir(NoDirection), momentum(0.0,0.0,0.0) 
{
    theParticle = const_cast<G4IonTable*>(G4ParticleTable::GetParticleTable()->GetIonTable())->GetIon(z,a);
    if (!theParticle)
    {
	G4cerr << "PCTProjectile::PCTProjectile(): Can not find ion  A = "
	       << a << " Z = " << z << '\n';
    }
}

G4bool PCTProjectile::operator==(const PCTProjectile& right) const
{
    return (theParticle == right.theParticle && momentum == right.momentum);
}

G4bool PCTProjectile::operator!=(const PCTProjectile& right) const
{
    return  (!(this->operator==(right)));
}


void PCTProjectile::SetDirection(const PCTProjectileDirection aDirection)
{
    if (dir != aDirection) 
    {
	dir = aDirection;
	this->ResetMomentum();
    }
    return;
}


G4bool PCTProjectile::SetKineticEnergy(const G4double KineticEnergy)
{
    G4bool succes = true;
    if (dir == NoDirection || !theParticle) 
    {
	succes = false;
    }
    else 
    {
	G4double Mass = theParticle->GetPDGMass();
	G4double MomentumMagnitude = sqrt(KineticEnergy*(KineticEnergy+2.0*Mass));

	if (dir == XaxisDirection) 
	{
	    momentum.setPx(MomentumMagnitude);
	    momentum.setPy(0.0);
	    momentum.setPz(0.0);
	}
	else if (dir == YaxisDirection)
	{
	    momentum.setPx(0.0);
	    momentum.setPy(MomentumMagnitude);
	    momentum.setPz(0.0);
	}
	else if (dir == ZaxisDirection)
	{
	    momentum.setPx(0.0);
	    momentum.setPy(0.0);
	    momentum.setPz(MomentumMagnitude);
	}
	else 
	{
	    G4double CosTheta = 1.0 - 2.0*G4UniformRand();
	    G4double SinTheta = sqrt(1.0 - CosTheta*CosTheta);
	    G4double Phi = twopi*G4UniformRand();
	    momentum.setPx(MomentumMagnitude*cos(Phi)*SinTheta);
	    momentum.setPy(MomentumMagnitude*sin(Phi)*SinTheta);
	    momentum.setPz(MomentumMagnitude*CosTheta);
	}
	
	momentum.setE(Mass+KineticEnergy);
    }
    return succes;
}


    //   G4double protonBindingEnergy = G4NucleiProperties::GetMassExcess(1,1) +
    //   G4NucleiProperties::GetMassExcess(A0,Z0) -
    //   G4NucleiProperties::GetMassExcess(A0+1,Z0+1);
