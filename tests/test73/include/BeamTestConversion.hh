// work out the position and direction of beam given intial p and pT

#ifndef BEAMTESTCONVERSION_HPP
#define BEAMTESTCONVERSION_HPP

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "BeamTestParameters.hh"


class BeamTestConversion
{

	public:

		// Constructor
		BeamTestConversion();
		BeamTestConversion(G4double p, G4double pt/*Parameters* parameter*/);
		// Destructor
		~BeamTestConversion();
		
		// Method functoins
		//void Angle();
		//double length() const { return len; }
		
		//geters and setters
		G4double getPT() const { return pT; }
		G4double getP() const { return p; }

		void setPT(const G4double pT_);
		void setP(const G4double p_);

		// statics
		static const double Pi;
		
	private:
		// data members
		G4double pT;
		G4double p;

};

// Free functions
G4double pZ(const BeamTestConversion& a);
G4double xzAngle(const BeamTestConversion& a);
G4ThreeVector PositionVec(const BeamTestConversion& a);
// arctangent
//double arctan(double y, double x);

#endif

