//
// GridParticleGun
//
// Declaration of grid particle gun
//
// Grid gun: particles travel parallel to vector "direction".
//           Each particle has a double index (ij) where (0 <= i < n), (0 <= j < m).
//           The starting point of particle ij is:
//
//                  startingPoint = origin + (i/n)*g1 + (j/m)*g2
//
//           where g1 and g2 are directions and origin a position.

#ifndef GridParticleGun_hh
#define GridParticleGun_hh

#include "globals.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryVertex.hh"
#include "G4ParticleMomentum.hh"
 
class G4Event;
class GridParticleGunMessenger;
 
class GridParticleGun : public G4VPrimaryGenerator
{
	public:
	GridParticleGun();
	GridParticleGun( G4ThreeVector direction,
			 G4ThreeVector origin,
			 G4ThreeVector grid1,
			 G4ThreeVector grid2,
			 G4int n1, G4int n2 );
	~GridParticleGun();
	
	void GeneratePrimaryVertex( G4Event *evt );
	
	inline void SetDirection( G4ThreeVector newDirection) { direction = newDirection; }
	inline G4ThreeVector GetDirection() const { return direction; }
	
	inline void SetOrigin( G4ThreeVector newOrigin) { origin = newOrigin; }
	inline G4ThreeVector GetOrigin() const { return origin; }

        inline void SetGrid1( G4ThreeVector newGrid1 ) { grid1 = newGrid1; }
	inline G4ThreeVector GetGrid1() const { return grid1; }

        inline void SetGrid2( G4ThreeVector newGrid2 ) { grid2 = newGrid2; }
	inline G4ThreeVector GetGrid2() const { return grid2; }
	
	inline void SetN1( G4int newN1 ) { n1 = newN1; };
	inline G4int GetN1() const { return n1; }
	
	inline void SetN2( G4int newN2 ) { n2 = newN2; };
	inline G4int GetN2() const { return n2; }


	private:
	void SetDefaults();
	
	protected:
	G4int n1, n2;
	G4ThreeVector direction, origin, grid1, grid2;
	G4ParticleDefinition	*particleType;
	
	private:
	GridParticleGunMessenger *messenger;
};

#endif
