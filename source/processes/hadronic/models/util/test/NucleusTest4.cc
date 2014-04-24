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
//
// $Id$
//
// ------------------------------------------------------------
//      GEANT 4 file
//
//
//             by Gunter Folger, June 1998.
//       class exercising G4Nucleus class.
//         testing nucleus creation - performance
// ------------------------------------------------------------

#include "G4Fancy3DNucleus.hh"



#include "G4StableIsotopes.hh"
#include "G4NucleiPropertiesTable.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

#include "G4ProtonField.hh"
#include "G4NeutronField.hh"


int main()
{

	G4Fancy3DNucleus nucleus;
	
// 	
// 	     
//
	
	G4StableIsotopes theIso;
		

	for (int Z=1; Z<93; Z++ )
	{
	   G4cout <<G4endl<< G4endl<< "new Element " ;
	   G4cout << theIso.GetName(Z) << G4endl;
	   G4cout << "         with " << theIso.GetNumberOfIsotopes(Z) ;
	   G4cout << " Isotopes"<< G4endl;
	   for (G4int iso=0; iso<theIso.GetNumberOfIsotopes(Z); iso++)
	   {
	       G4double massnumber=
	          theIso.GetIsotopeNucleonCount(theIso.GetFirstIsotope(Z)+iso);
              for ( G4int repeat=0; repeat< 500; ++repeat )
	      {
	           nucleus.Init(massnumber,(G4double) Z);
	      }
	   }   
	}
	return 0;
	
} 
	
	
	
