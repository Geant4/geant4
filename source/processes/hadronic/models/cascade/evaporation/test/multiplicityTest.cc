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
// $Id: multiplicityTest.cc,v 1.2 2001-07-11 10:03:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Test program for G4 Bertini Evaporation.

#define G4VERBOSE 1

#include "G4LayeredNucleus.hh"
#include "G4BertiniEvaporation.hh"

int main(int argc, char *argv[])
{
  G4LayeredNucleus nucl;
  G4BertiniEvaporation bert;
  G4VParticleChange * pc;
  G4int A, Z;
  G4double E;

  G4int rounds;

  if ( argc == 4 )
    {
      sscanf( argv[1], "%d", &A);
      sscanf( argv[2], "%d", &Z);
      sscanf( argv[3], "%d", &rounds);
      G4cout << A << " " << Z << " " << rounds << endl;
    }
  
  else 
    {
      G4cout << "***************************" << endl
	     << "*  Parameters: A Z rounds  * " << endl
	     << "*  use 68 32 5000, takes  * " << endl
	     << "***************************" << endl << endl;
      return 0 ;
    } 

  bert.setVerboseLevel(0);
  
  for ( G4double energy = 1 ; energy < 150 ; energy += 0.5 )
    {
      G4int n=0, p=0, d=0, t=0, h3=0, h4=0, g=0;
      G4int m1=0, m2=0, m3=0;

      nucl.SetParameters( A, Z );
      nucl.AddExcitationEnergy( energy - nucl.GetEnergyDeposit() );

      for ( G4int j = 0 ; j < rounds ; j++)
	{
	  n=p=d=t=h3=h4=g=0;
	  pc = bert.BreakItUp( nucl );
	  
	  for ( G4int i = 0 ; i < pc->GetNumberOfSecondaries() ;  i++ )
	    {
	      char * name = pc->GetSecondary( i )->GetDefinition()->GetParticleName() ;
	      if ( strcmp ( name , "gamma" ) == 0 ) g++;
	      
	      delete pc->GetSecondary( i );
	    } // loop over particle change vector
	  
	  if ( pc->GetNumberOfSecondaries() - g == 1 ) m1++;
	  else if ( pc->GetNumberOfSecondaries() - g == 2 ) m2++;
	  else if ( pc->GetNumberOfSecondaries() - g == 3 ) m3++;

	  pc->Clear();
	  delete pc;
      
	}


      G4cout << energy << "\t"
	     << m1 << " "
	     << m2 << " "
	     << m3 << endl;

    } // energy loop
  
  return 0;
}
