// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ratioTest.cc,v 1.1 2000-08-16 07:56:18 miheikki Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Test program for G4 Bertini Evaporation.

#define G4VERBOSE 1

#include "G4LayeredNucleus.hh"
#include "G4BertiniEvaporation.hh"

int main(int argc, char *argv[])
{
  G4cout << "%*********************************" << endl
	 << "%*  Ratio test                   *" << endl 
	 << "%*  Output in Latex table format *" << endl 
	 << "%*********************************" << endl << endl;

  G4LayeredNucleus nucl;
  G4BertiniEvaporation bert;
  G4VParticleChange * pc;
  G4int A, Z;
  G4double E;

  G4int rounds=10000;
  bert.setVerboseLevel(0);

  G4int ntemp=0, ptemp=0, dtemp=0, ttemp=0, h3temp=0, h4temp=0, gtemp=0;
  G4int n=0, p=0, d=0, t=0, h3=0, h4=0, g=0, pn=0, _2n=0, _2p=0, alfan=0, p2n=0, _2pn=0;
  G4int j;

  G4cout << "% * Rounds : "<< rounds << endl;

  ////////////////////////// Se 74 30, 35 MeV

  A = 74;
  Z = 34;
  E = 35;
  
  nucl.SetParameters( A, Z );
  nucl.AddExcitationEnergy( E - nucl.GetEnergyDeposit() );
  n=p=d=t=h3=h4=g=0;

  for ( j = 0 ; j < rounds ; j++)
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
  
  G4cout << "Se$^{74}$     &  " << E << "   &   pn/2n   &   1.7  & " << (G4double) pn/((G4double)_2n)
	 << " \\\\ % " << pn << "  " << _2n << endl;

  ///////////////////////// Ge 68 32, E=20

  A = 68;
  Z = 32;
  E = 20;

  nucl.SetParameters( A, Z );
  nucl.AddExcitationEnergy( E - nucl.GetEnergyDeposit() );
  ntemp=0, ptemp=0, dtemp=0, ttemp=0, h3temp=0, h4temp=0, gtemp=0;
  n=p=d=t=h3=h4=g=pn=_2n=_2p=alfan=p2n=0;

  for ( j = 0 ; j < rounds ; j++)
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
  
  G4cout << "Ge$^{68}$   &  " << E << "   &   p/n     &   1.76 & " << (G4double) p/((G4double)n)
	 << " \\\\ % " << p << "  " << n << endl;

  ///////////////////////// Ge 68 32, E=35

  A = 68;
  Z = 32;
  E = 35;

  nucl.SetParameters( A, Z );
  nucl.AddExcitationEnergy( E - nucl.GetEnergyDeposit() );
  ntemp=0, ptemp=0, dtemp=0, ttemp=0, h3temp=0, h4temp=0, gtemp=0;
  n=p=d=t=h3=h4=g=pn=_2n=_2p=alfan=p2n=0;

  for ( j = 0 ; j < rounds ; j++)
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
  
  G4cout << "            &  " << E << "   &   pn/2n   &   8.4  & " << (G4double) pn/((G4double)_2n)
	 << " \\\\ % " << pn << "  " << _2n << endl;

  ///////////////////////// Ge 68 32, E=40

  A = 68;
  Z = 32;
  E = 40;

  nucl.SetParameters( A, Z );
  nucl.AddExcitationEnergy( E - nucl.GetEnergyDeposit() );
  ntemp=0, ptemp=0, dtemp=0, ttemp=0, h3temp=0, h4temp=0, gtemp=0;
  n=p=d=t=h3=h4=g=pn=_2n=_2p=alfan=p2n=_2pn=0;


  for ( j = 0 ; j < rounds ; j++)
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
      
      if ( ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) n++;
      if ( ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) p++;
      if ( ( ntemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ||
	   ( dtemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) ) pn++;
      if ( ntemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2n++;
      if ( ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2p++;
      if ( ( ntemp == 2 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ) p2n++;
      if ( ( ntemp == 1 && ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) )_2pn++;
      if ( ntemp == 1 && h4temp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) alfan++;
      //x  G4cout << "% n p pn 2n 2p p2n * " << en
      pc->Clear();
      delete pc;
      
    }
  
  G4cout << "            &  " << E << "   &   2pn/p2n &   5.8  & " << (G4double) _2pn/((G4double)p2n)
	 << " \\\\ % " << _2pn << "  " << p2n << endl;

  ///////////////////////// Ge 67 31, E=35

  A = 67;
  Z = 31;
  E = 35;

  nucl.SetParameters( A, Z );
  nucl.AddExcitationEnergy( E - nucl.GetEnergyDeposit() );
  ntemp=0, ptemp=0, dtemp=0, ttemp=0, h3temp=0, h4temp=0, gtemp=0;
  n=p=d=t=h3=h4=g=pn=_2n=_2p=alfan=p2n=_2pn=0;


  for ( j = 0 ; j < rounds ; j++)
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
      
      if ( ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) n++;
      if ( ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) p++;
      if ( ( ntemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ||
	   ( dtemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) ) pn++;
      if ( ntemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2n++;
      if ( ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2p++;
      if ( ( ntemp == 2 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ) p2n++;
      if ( ( ntemp == 1 && ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) )_2pn++;
      if ( ntemp == 1 && h4temp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) alfan++;
      //x  G4cout << "% n p pn 2n 2p p2n * " << en
      pc->Clear();
      delete pc;
      
    }
  
  G4cout << "Ga^{67}     &  " << E << "   &   pn/2n   &   3.3  & " << (G4double) pn/((G4double)_2n)
	 << " \\\\ % " << pn << "  " << _2n << endl;

 ///////////////////////// Ni 58 28

  A = 58;
  Z = 28;
  E = 25;

  nucl.SetParameters( A, Z );
  nucl.AddExcitationEnergy( E - nucl.GetEnergyDeposit() );
  ntemp=0, ptemp=0, dtemp=0, ttemp=0, h3temp=0, h4temp=0, gtemp=0;
  n=p=d=t=h3=h4=g=pn=_2n=_2p=alfan=p2n=_2pn=0;


  for ( j = 0 ; j < rounds ; j++)
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
      
      if ( ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) n++;
      if ( ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) p++;
      if ( ( ntemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ||
	   ( dtemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) ) pn++;
      if ( ntemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2n++;
      if ( ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2p++;
      if ( ( ntemp == 2 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ) p2n++;
      if ( ( ntemp == 1 && ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) )_2pn++;
      if ( ntemp == 1 && h4temp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) alfan++;
      //x  G4cout << "% n p pn 2n 2p p2n * " << en
      pc->Clear();
      delete pc;
      
    }
  
  G4cout << "Ni^{58}     &  " << E << "   &   p/n     &   3.8  & " << (G4double) p/((G4double)n)
	 << " \\\\ % " << p << "  " << n << endl;

 ///////////////////////// 

  ///////////////////////// Ni 58 28, 35

  A = 58;
  Z = 28;
  E = 35;

  nucl.SetParameters( A, Z );
  nucl.AddExcitationEnergy( E - nucl.GetEnergyDeposit() );
  ntemp=0, ptemp=0, dtemp=0, ttemp=0, h3temp=0, h4temp=0, gtemp=0;
  n=p=d=t=h3=h4=g=pn=_2n=_2p=alfan=p2n=_2pn=0;


  for ( j = 0 ; j < rounds ; j++)
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
      
      if ( ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) n++;
      if ( ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) p++;
      if ( ( ntemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ||
	   ( dtemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) ) pn++;
      if ( ntemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2n++;
      if ( ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2p++;
      if ( ( ntemp == 2 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ) p2n++;
      if ( ( ntemp == 1 && ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) )_2pn++;
      if ( ntemp == 1 && h4temp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) alfan++;
      //x  G4cout << "% n p pn 2n 2p p2n * " << en
      pc->Clear();
      delete pc;
      
    }
  
  G4cout << "Ni^{58}     &  " << E << "   &   pn/2n   &   67   & " << (G4double) pn/((G4double)_2n)
	 << " \\\\ % " << pn << "  " << _2n << endl;

  ///////////////////////// Ni 58 28 40

  A = 58;
  Z = 28;
  E = 40;

  nucl.SetParameters( A, Z );
  nucl.AddExcitationEnergy( E - nucl.GetEnergyDeposit() );
  ntemp=0, ptemp=0, dtemp=0, ttemp=0, h3temp=0, h4temp=0, gtemp=0;
  n=p=d=t=h3=h4=g=pn=_2n=_2p=alfan=p2n=_2pn=0;


  for ( j = 0 ; j < rounds ; j++)
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
      
      if ( ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) n++;
      if ( ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) p++;
      if ( ( ntemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ||
	   ( dtemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) ) pn++;
      if ( ntemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2n++;
      if ( ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2p++;
      if ( ( ntemp == 2 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ) p2n++;
      if ( ( ntemp == 1 && ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) )_2pn++;
      if ( ntemp == 1 && h4temp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) alfan++;
      //x  G4cout << "% n p pn 2n 2p p2n * " << en
      pc->Clear();
      delete pc;
      
    }
  
  G4cout << "Ni^{58}     &  " << E << "   &   2pn/p2n &   7.3  & " << (G4double) _2pn/((G4double)p2n)
	 << " \\\\ % " << _2pn << "  " << p2n << endl;

  ///////////////////////// Fe 54 26 35

  A = 54;
  Z = 26;
  E = 35;

  nucl.SetParameters( A, Z );
  nucl.AddExcitationEnergy( E - nucl.GetEnergyDeposit() );
  ntemp=0, ptemp=0, dtemp=0, ttemp=0, h3temp=0, h4temp=0, gtemp=0;
  n=p=d=t=h3=h4=g=pn=_2n=_2p=alfan=p2n=_2pn=0;


  for ( j = 0 ; j < rounds ; j++)
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
      
      if ( ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) n++;
      if ( ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) p++;
      if ( ( ntemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ||
	   ( dtemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) ) pn++;
      if ( ntemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2n++;
      if ( ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2p++;
      if ( ( ntemp == 2 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ) p2n++;
      if ( ( ntemp == 1 && ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) )_2pn++;
      if ( ntemp == 1 && h4temp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) alfan++;
      //x  G4cout << "% n p pn 2n 2p p2n * " << en
      pc->Clear();
      delete pc;
      
    }
  
  G4cout << "Fe^{54}     &  " << E << "   &   pn/2n   &   31   & " << (G4double) pn/((G4double)_2n)
	 << " \\\\ % " << pn << "  " << _2n << endl;

  /// 
  ///////////////////////// Mn 54 25 30

  A = 54;
  Z = 25;
  E = 30;

  nucl.SetParameters( A, Z );
  nucl.AddExcitationEnergy( E - nucl.GetEnergyDeposit() );
  ntemp=0, ptemp=0, dtemp=0, ttemp=0, h3temp=0, h4temp=0, gtemp=0;
  n=p=d=t=h3=h4=g=pn=_2n=_2p=alfan=p2n=_2pn=0;


  for ( j = 0 ; j < rounds ; j++)
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
      
      if ( ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) n++;
      if ( ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) p++;
      if ( ( ntemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ||
	   ( dtemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) ) pn++;
      if ( ntemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2n++;
      if ( ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2p++;
      if ( ( ntemp == 2 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ) p2n++;
      if ( ( ntemp == 1 && ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) )_2pn++;
      if ( ntemp == 1 && h4temp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) alfan++;
      //x  G4cout << "% n p pn 2n 2p p2n * " << en
      pc->Clear();
      delete pc;
      
    }
  
  G4cout << "Mn^{54}     &  " << E << "   &   pn/2n   &   <1   & " << (G4double) pn/((G4double)_2n)
	 << " \\\\ % " << pn << "  " << _2n << endl;

 ///////////////////////// Cr 50 24

  A = 50;
  Z = 24;
  E = 25;

  nucl.SetParameters( A, Z );
  nucl.AddExcitationEnergy( E - nucl.GetEnergyDeposit() );
  ntemp=0, ptemp=0, dtemp=0, ttemp=0, h3temp=0, h4temp=0, gtemp=0;
  n=p=d=t=h3=h4=g=pn=_2n=_2p=alfan=p2n=_2pn=0;


  for ( j = 0 ; j < rounds ; j++)
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
      
      if ( ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) n++;
      if ( ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) p++;
      if ( ( ntemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ||
	   ( dtemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) ) pn++;
      if ( ntemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2n++;
      if ( ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2p++;
      if ( ( ntemp == 2 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ) p2n++;
      if ( ( ntemp == 1 && ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) )_2pn++;
      if ( ntemp == 1 && h4temp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) alfan++;
      //x  G4cout << "% n p pn 2n 2p p2n * " << en
      pc->Clear();
      delete pc;
      
    }
  
  G4cout << "Cr^{50}     &  " << E << "   &   p/n     &   0.51 & " << (G4double) p/((G4double)n)
	 << " \\\\ % " << p << "  " << n << endl;

  ///////////////////////// Cr 50 24 35

  A = 50;
  Z = 24;
  E = 35;

  nucl.SetParameters( A, Z );
  nucl.AddExcitationEnergy( E - nucl.GetEnergyDeposit() );
  ntemp=0, ptemp=0, dtemp=0, ttemp=0, h3temp=0, h4temp=0, gtemp=0;
  n=p=d=t=h3=h4=g=pn=_2n=_2p=alfan=p2n=_2pn=0;


  for ( j = 0 ; j < rounds ; j++)
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
      
      if ( ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) n++;
      if ( ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) p++;
      if ( ( ntemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ||
	   ( dtemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 1 ) ) pn++;
      if ( ntemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2n++;
      if ( ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 2 ) _2p++;
      if ( ( ntemp == 2 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ntemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) ) p2n++;
      if ( ( ntemp == 1 && ptemp == 2 && pc->GetNumberOfSecondaries() - gtemp == 3 ) ||
	       ( dtemp == 1 && ptemp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) )_2pn++;
      if ( ntemp == 1 && h4temp == 1 && pc->GetNumberOfSecondaries() - gtemp == 2 ) alfan++;
      //x  G4cout << "% n p pn 2n 2p p2n * " << en
      pc->Clear();
      delete pc;
      
    }
  
  G4cout << "            &  " << E << "   &   pn/2n   &   46   & " << (G4double) pn/((G4double)_2n)
	 << " \\\\ % " << pn << "  " << _2n << endl;

  
  return 0;
}







