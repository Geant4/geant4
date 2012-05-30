
//  A parameterisation that describes a series of boxes along Z
//    The boxes have equal width, & their lengths are a linear equation.
//    They are spaced an equal distance apart, starting from given location.
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef BEAMTESTCELLPARAMETRISATION_H
#define BEAMTESTCELLPARAMETRISATION_H 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"

class G4VPhysicalVolume;
class G4Box;

// Dummy declarations to get rid of warnings ...
class G4Trd;
class G4Trap;
class G4Cons;
class G4Orb;
class G4Sphere;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Polycone;
class G4Polyhedra;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class BeamTestCellParameterisation : public G4VPVParameterisation
{ 
	public:

		BeamTestCellParameterisation(G4int    NoChambers, 
				G4double startZ, 
				G4double spacing,
				G4double widthChamber);

		virtual	~BeamTestCellParameterisation();

		void ComputeTransformation (const G4int copyNo,
				G4VPhysicalVolume* physVol) const;

	//	void ComputeDimensions (G4Box & trackerLayer, const G4int copyNo,
	//			const G4VPhysicalVolume* physVol) const;

	private:  // Dummy declarations to get rid of warnings ...

		/*void ComputeDimensions (G4Trd&,const G4int,const G4VPhysicalVolume*) const {}
		void ComputeDimensions (G4Trap&,const G4int,const G4VPhysicalVolume*) const {}
		void ComputeDimensions (G4Cons&,const G4int,const G4VPhysicalVolume*) const {}
		void ComputeDimensions (G4Sphere&,const G4int,const G4VPhysicalVolume*) const {}
		void ComputeDimensions (G4Orb&,const G4int,const G4VPhysicalVolume*) const {}
		void ComputeDimensions (G4Torus&,const G4int,const G4VPhysicalVolume*) const {}
		void ComputeDimensions (G4Para&,const G4int,const G4VPhysicalVolume*) const {}
		void ComputeDimensions (G4Hype&,const G4int,const G4VPhysicalVolume*) const {}
		void ComputeDimensions (G4Tubs&,const G4int,const G4VPhysicalVolume*) const {}
		void ComputeDimensions (G4Polycone&,const G4int,const G4VPhysicalVolume*) const {}
		void ComputeDimensions (G4Polyhedra&,const G4int,const G4VPhysicalVolume*) const {}
		*/
	private:

		G4int    fNoChambers;   
		G4double fStartZ;
		G4double fHalfWidth;    //  The half-width of each tracker chamber
		G4double fSpacing;      //  The distance between the chambers' center
    public:
        void SetNumberChambers(G4int val){ fNoChambers = val; }
        void SetStartZ(G4double val) { fStartZ = val; }
        void SetWidth(G4double val) { fHalfWidth = 0.5*val; }
        void SetSpacing(G4double val) { fSpacing = val; }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

