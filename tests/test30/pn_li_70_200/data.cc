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
#include <fstream>
#include <iomanip>
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
  } else if(argc > 2) {
    verbose = 2;
  }

  ifstream* fin = new ifstream();
  string fname = argv[1];
  fin->open(fname.c_str());
  if( !fin->is_open()) {
    G4cout << "Input file <" << fname << "> does not exist! Exit" << G4endl;
    exit(1);
  }

  ofstream* fout_b = new ofstream();
  string fname2 = "dsdedtet.out";
  fout_b->open(fname2.c_str(), std::ios::out|std::ios::trunc);
  ofstream* fout_c = new ofstream();
  string fname3 = "dsdedtet.dat";
  fout_c->open(fname3.c_str(), std::ios::out|std::ios::trunc);

  (*fout_c) << "#####..Result.of.parcing..#######.. "
            << " Angle(degree)= 0" << G4endl;


  //there can't be lines longer than nmax characters
  const int nmax = 100;
  char line[nmax];
  std::string line1, line2, word1, word2, word3;
  G4bool end = true;
  int counter = 0;
  G4DataVector* energy = new G4DataVector();
  G4DataVector* cross = new G4DataVector();

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
    if(fin->eof() || line[0] == 'E' && line[1] == 'N' && line[2] == 'D') {
      end = false;
    } else if(line[0] == 'M' && line[1] == 'E' && line[2] == 'V') {
      end = false;

      double e0, e, s, s1;
      double ebeam = 0.0;

      for(int i=0; i<1080; i++) {
        (*fin) >> e0 >> e >> s >> s1;
        if(1 < verbose) {
          cout << "e0= " << e0 << " e= " << e << " cross= " << s << endl;
        }
	if(ebeam != e0) {
          if(1 < verbose) {
            cout << "New beam energy " << e0 << endl;
          }
	  if(ebeam > 0.0) {
	    int nbin = cross->size();
            (*fout_b) << "#####..Ebeam(MeV) = " << ebeam/MeV << G4endl;
            (*fout_b) << "ve/create X(" << nbin << ") R ";
            (*fout_c) << "#####..Ebeam(MeV) = " << ebeam/MeV << G4endl;

	    int j;
	    for(j=0; j<nbin; j++) {
              (*fout_b) << (*energy)[j] << " ";
              (*fout_c) << "e(MeV)= " << (*energy)[j]
                        << " cross(mb/MeV/sr)= " << (*cross)[j] << endl;
	    }
            (*fout_b) << endl;
            (*fout_b) << "ve/create Y(" << nbin << ") R ";
	    for(j=0; j<nbin; j++) {
              (*fout_b) << (*cross)[j] << " ";
	    }
            (*fout_b) << endl;
          }
	  if(e0 == 0.0) break;
	  ebeam = e0;
          energy->clear();
          cross->clear();
	}
	energy->push_back(e);
	cross->push_back(s);
      }
    }
  } while (end);
  (*fout_c) << "#####..End..#####" << G4endl;
}






