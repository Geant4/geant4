#include "g4std/fstream"
#include "g4std/iostream"
#include <stdlib.h>
#include "globals.hh"

int main()
{
  G4String theName = "lead.inelasticxsec.kumac";

  G4std::ifstream aDataSet(theName, G4std::ios::in);
  int count = 0;
  aDataSet >> count;
  double * ee = new double[count];
  double * xsec = new double[count];
  for(int i=0; i<count; i++)
  {
    aDataSet >> ee[i]>>xsec[i];
    ee[i]/=1000000;
  }
  cout << "ve/cr ee("<<count<<") r ";
  for( i=0; i<count; i++) 
  {
    if(i==30*(i/30)) cout <<" _"<<G4endl;
    cout << ee[i]<<" ";
  }
  cout << G4endl<<"ve/cr xsec("<<count<<") r ";
  for( i=0; i<count; i++) 
  { 
    if(i==30*(i/30)) cout <<" _"<<G4endl;
    cout << xsec[i]<<" ";
  }
  cout << G4endl<<"null -10 30 0 3"<<G4endl;
  cout << "ve/cr err("<<count<<") r "<<count<<"*0.01"<<G4endl;
  cout << "hpl/err ee xsec err err "<<count<<" 20 0.15"<<G4endl;
}
