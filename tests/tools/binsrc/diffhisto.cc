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
#include <fstream>
#include "Histo.hh"


int main(int argc,char* argv[])
{
  
  std::ifstream in(argv[1]),comp(argv[2]);

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
  
