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
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     Test30
//
//      Author:        V.Ivanchenko 
// 
//      Creation date: 12 March 2002
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
#include <string>
#include <iostream.h>
#include <stdlib.h>
#include <strstream.h>


int main(int argc, char** argv)
{
  int verbose = 1;

  // -------------------------------------------------------------------
  // Control on input

  if(argc < 2) {
    G4cout << "Input file is not specified! Exit" << G4endl;
    exit(1);
  }

  ifstream* fin = new ifstream();
  string fname = argv[1];
  fin->open(fname.c_str());
  if( !fin->is_open()) {
    G4cout << "Input file <" << fname << "> does not exist! Exit" << G4endl;
    exit(1);
  }

  /*
  ofstream* fout_a = new ofstream();
  string fname1 = "dsde_1.dat";
  fout_a->open(fname1.c_str(), std::ios::out|std::ios::trunc);
  ofstream* fout_b = new ofstream();
  string fname2 = "dsdtet_1.dat";
  fout_b->open(fname2.c_str(), std::ios::out|std::ios::trunc);
  */
  ofstream* fout_c = new ofstream();
  string fname3 = "dsdedtet_1.dat";
  fout_c->open(fname3.c_str(), std::ios::out|std::ios::trunc);

  //there can't be lines longer than nmax characters
  const int nmax = 200;
  char line[nmax]; 
  std::string line1, line2, word1, word2, word3;
  G4bool end = true;


  G4DataVector* angle = new G4DataVector();
  std::vector<G4DataVector*> cs;
  G4int counter = 0;
  G4double an;

  // main loop 

  do {

    // read next line

    counter++;
    for(int ii = 0; ii < nmax; ii++) {line[ii] = ' ';}
    fin->getline( line, nmax);
    line1 = std::string("");
    line1 = std::string(line, nmax);
    if(1 < verbose) {
      cout << "Next line # " << counter << ": " << line1 << endl;
      cout << "First symbols= " << line[0] << line[1] << line[2] << endl;
    }  

    // analize line contence

    // end of file
    if(fin->eof()) {
      end = false;
    } else if(line[0] == 'E' && line[1] == 'N' && line[2] == 'D') {
      end = false;
    } else if(line[0] == 'A' && line[1] == 'D' && line[2] == 'E') {

      G4DataVector* energy = new G4DataVector();
      G4DataVector* cross = new G4DataVector();
      (*fin) >> an;
      counter++;
      an *= degree;
      angle->push_back(an);
      (*fout_c) << "#####..New data..##### Theta(deg)= " << an/degree << G4endl;
      G4int nbin, i;
      std::string mystr;
      G4double e0, e1, e2, es1, es2, x, xs;
      (*fin) >> mystr >> i >> nbin;
      (*fin) >> mystr >> mystr >> mystr >> mystr >> mystr >> mystr;
      (*fin) >> mystr >> mystr >> mystr >> mystr >> mystr >> mystr;

      for (i=0; i<nbin; i++) {
        (*fin) >> e1 >> e2 >> es1 >> es2 >> x >> xs;
        x *= 1000.0;      
        xs *= 1000.0;      
        if(1 < verbose) {
          cout << "an= " << an/degree << " e1= " << e1 
               << " e2= " << e2 << " cross= " << x
               << " +- " << xs  
               << endl;
        }  

        e0 = 0.5*(e1 + e2);
        energy->push_back(e0);
        cross->push_back(x);
        (*fout_c) << "e(MeV)= " << e0/MeV 
                  << " cross(mb/sr)= " << x 
                  << " +- " << xs 
                  << endl;

      } 
    }
  } while (end);
  //  (*fout_a) << "#####..End..#####" << G4endl;
  //  (*fout_b) << "#####..End..#####" << G4endl;
  (*fout_c) << "#####..End..#####" << G4endl;
}








