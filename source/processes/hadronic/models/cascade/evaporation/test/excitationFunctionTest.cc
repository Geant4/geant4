// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: excitationFunctionTest.cc,v 1.1 2000-08-16 07:56:17 miheikki Exp $
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
      G4cout << " % excitation function for " << A << " " << Z << " " << rounds << endl;
    }
  
  else 
    {
      G4cout << "***************************" << endl
	     << "* Exc Fs for : n p pn 2n 2p p2n * " << endl
	     << "*  Parameters: A Z rounds  * " << endl
	     << "*  use 68 32 5000, takes  * " << endl
	     << "*  copper 63 (70%) -> zn 64,30 * " << endl
	     << "*  copper 65 (30%) -> zn 66,30 * " << endl
	     << "***************************" << endl << endl;
      return 0 ;
    } 

  G4cout << "% n p pn 2n 2p p2n alfan * " << endl;

  bert.setVerboseLevel(0);

  G4int ntemp=0, ptemp=0, dtemp=0, ttemp=0, h3temp=0, h4temp=0, gtemp=0;
  
  for ( G4double energy = 1 ; energy < 100 ; energy += 0.5 )
    {
      G4int n=0, p=0, d=0, t=0, h3=0, h4=0, g=0, pn=0, _2n=0, _2p=0, p2n=0, alfan=0;
      nucl.SetParameters( A, Z );
      nucl.AddExcitationEnergy( energy - nucl.GetEnergyDeposit() );
      n=p=d=t=h3=h4=g=0;

      for ( G4int j = 0 ; j < rounds ; j++)
	{
	  ntemp=ptemp=dtemp=ttemp=h3temp=h4temp=gtemp=0;
	  pc = bert.BreakItUp( nucl );
	  
	  for ( G4int i = 0 ; i < pc->GetNumberOfSecondaries() ;  i++ )
	    {
	      const char * name = pc->GetSecondary( i )->GetDefinition()->GetParticleName() ;
	      if ( strcmp ( name , "neutron" ) == 0 ) ntemp++;
	      else if ( strcmp ( name , "proton" ) == 0 ) ptemp++;
	      else if ( strcmp ( name , "deuteron" ) == 0 ) dtemp++;
	      else if ( strcmp ( name , "triton" ) == 0 ) ttemp++;
	      else if ( strcmp ( name , "He3" ) == 0 ) h3temp++;
	      else if ( strcmp ( name , "alpha" ) == 0 ) h4temp++;
	      else if ( strcmp ( name , "gamma" ) == 0 ) gtemp++;
	      
	      delete pc->GetSecondary( i );
	    } // loop over particle change vector

//  	  G4cout << "a" 
//  		 << ntemp << " " 
//  		 << ptemp << " " 
//  		 << dtemp << " " 
//  		 << ttemp << " " 
//  		 << h3temp << " " 
//  		 << h4temp << " " 
//  		 << gtemp << " "
//  		 << pc->GetNumberOfSecondaries() << endl;
	  
//  	  if ( pc->GetNumberOfSecondaries() - g == 1 ) m1++;
//  	  else if ( pc->GetNumberOfSecondaries() - g == 2 ) m2++;
//  	  else if ( pc->GetNumberOfSecondaries() - g == 3 ) m3++;
	  
	  if ( ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) n++;
	  if ( ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) p++;
	  if ( ( ntemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ||
	       ( dtemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) ) pn++;
	  if ( ntemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2n++;
	  if ( ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2p++;
	  if ( ( ntemp == 2 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ) p2n++;
	  if ( ntemp == 1 && h4temp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) alfan++;
	  //x  G4cout << "% n p pn 2n 2p p2n * " << en
	  pc->Clear();
	  delete pc;
      
	}

      G4cout << energy << "\t"
	     << n << " "
	     << p << " "
	     << pn << " "
	     << _2n << " "
	     << _2p << " "
	     << p2n << " "
	     << alfan << endl;
	  //	     << g << endl;

    } // energy loop
  
  return 0;
}
