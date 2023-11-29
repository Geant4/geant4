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
#ifdef G4_USE_FLUKA


#include "string_conversion.h"


namespace cpp_to_fortran {

	std::vector<char> convertStringNoSizeNeeded(const std::string& cString) {
		std::vector<char> fString(cString.c_str(), cString.c_str() + cString.size());
		return fString;
	}

}


namespace fortran_to_cpp {

	std::string convertString(const char* fString, const int fStringLength) {

		std::string cString(fStringLength, ' ');
		for (int charIndex = 0; charIndex < fStringLength - 1; ++charIndex) {
			cString[charIndex] = *(fString + charIndex);
		}
		return cString;
	}

}


#endif // G4_USE_FLUKA
