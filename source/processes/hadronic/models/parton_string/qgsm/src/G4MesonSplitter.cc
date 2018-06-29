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
#include "G4MesonSplitter.hh"
#include "Randomize.hh"

G4bool G4MesonSplitter::SplitMeson(G4int PDGcode, G4int* aEnd, G4int* bEnd)
{
	G4bool result = true;
	G4int absPDGcode = std::abs(PDGcode);
	if (absPDGcode >= 1000) return false;
	if(absPDGcode == 22)                   // For gamma -> 4 (u ubar) + 1 (d dbar)
	{
		G4int it=1;
		if(G4UniformRand()<0.8) it++;   // Uzhi Oct. 2016 0.5 -> 0.8
		*aEnd = it;
		*bEnd = -it;
	}
	else
	{
		G4int heavy =  absPDGcode/100;
		G4int light = (absPDGcode%100)/10;
		G4int anti  = 1 - 2*(std::max(heavy, light)%2);
		if (PDGcode < 0 ) anti = -anti;
		heavy *=  anti;
		light *= -anti;
		if ( anti < 0) G4SwapObj(&heavy, &light);
		*aEnd = heavy;
		*bEnd = light;
	}
	return result;
}
