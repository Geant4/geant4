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
