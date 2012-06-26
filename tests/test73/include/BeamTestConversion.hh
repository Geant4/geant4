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

