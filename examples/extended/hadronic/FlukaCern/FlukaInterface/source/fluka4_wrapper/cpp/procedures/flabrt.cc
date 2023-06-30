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
// Kills the FLUKA run.
//
// Author: G.Hugo, 01 August 2022
//
// ***************************************************************************
#ifdef G4_USE_FLUKA

#include "flabrt.h"


extern "C" {
	// OWN HELPER FUNCTION
	extern void flabrw_(const char* cOut, const int& lenOut, const char* cMes, const int& lenMes);
}


namespace fluka_cpp_wrapper {

	void flabrt(const std::string& CHROUT, const std::string& CHMESS) {
		const int lenOut = CHROUT.size();
		const int lenMes = CHMESS.size();

		flabrw_(CHROUT.c_str(), lenOut, CHMESS.c_str(), lenMes);
	}

}


#endif // G4_USE_FLUKA
