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
#include "Histo.hh"


int main(int argc,char* argv[])
{
  
  G4std::ifstream in(argv[1]),comp(argv[2]);

  odHisto hin,hcomp;
  char c,bla[80];
  int nhisto[2];
  
  in >> c >> bla >> nhisto[0];
  comp >> c >> bla >> nhisto[1];

  if(nhisto[0]!=nhisto[1]) exit(-1);

  for(int i=0;i<nhisto[0];i++)
    {
      hin.input(in);
      hcomp.input(comp);
      if(hin.compare(hcomp)>1.) exit(-1);
    }

  exit(0);
}
  
