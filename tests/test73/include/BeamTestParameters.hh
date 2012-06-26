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
// $Id: Parameters.h,v 1.14 2010/09/29 07:59:26 paula Exp $
#ifndef BEAMTESTPARAMETERS_H 
#define BEAMTESTPARAMETERS_H 1

// Include files
#include<iostream>
#include<ctype.h>
#include<string>
#include<list>
#include<map>

/** @class Parameters Parameters.h
 *
 */

class Parameters {
	public: 
		/// Standard constructor
		Parameters( ); 
		int readCommandLine(int,char**);
		//int readConditions();
		bool exists(std::string);
		void help();

		virtual ~Parameters( ); ///< Destructor
		
		double pTransverse;
		double momentum;
		double zThickness;
		std::string filename;
        int numberOfChambers;
        double chamberSpacing;
};
#endif // BEAMTESTPARAMETERS_H
