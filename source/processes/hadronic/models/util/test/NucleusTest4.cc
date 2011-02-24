//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: NucleusTest4.cc,v 1.1 2004-05-12 09:38:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
	
	
	
