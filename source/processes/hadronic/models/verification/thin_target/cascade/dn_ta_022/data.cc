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
//      Author:        V.Ivanchenko
//
//      Creation date: 12 March 2002
//
//      Modifications: 14.11.03 new directory
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


  //there can't be lines longer than nmax characters
  const int nmax = 100;
  char line[nmax];
  std::string line1, line2, word1, word2, word3;
  G4bool end = true;
  double ebeam = 0.0;
  std::vector<int>    inn;
  std::vector<double> angle;
  std::vector<double> energy;
  std::vector<G4DataVector*> cross;
  std::vector<G4DataVector*> err;

  // main loop
  int counter = 0;

  do {
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

      double an, a, e, s, s1, e0;
      int i, j, jj, nbin;
      (*fin) >> ebeam >> e;
      (*fin) >> word1 >> word2 >> nbin;
      (*fin) >> word1 >> word2 >> word3 >> line1;
      (*fin) >> word1 >> word2 >> word3 >> line1;
      inn.clear();
      angle.clear();
      energy.clear();
      cross.clear();
      err.clear();
      angle.push_back(30.);
      angle.push_back(60.);
      angle.push_back(90.);
      angle.push_back(120.);
      angle.push_back(150.);
      for (i=0; i<5; i++) {
        cross.push_back(new G4DataVector());
        err.push_back(new G4DataVector());
      }
      an = 0.;
      e0 = 0.;
      j = 0;
      for(i=0; i<nbin; i++) {
        (*fin) >> e >> a >> s >> s1;
        if(1 < verbose) {
          cout << "an= " << a << "e= " << e << " cross= " << s << " +- " << s1 << endl;
        }
	if(e != e0) {
	  e0 = e;
	  energy.push_back(e);
	}
        for (j=0; j<5; j++) {if(abs(a - angle[j])<1.) break;}
	cross[j]->push_back(s);
	err[j]->push_back(s1);
      }

      for(i=0; i<5; i++) {
        int k = cross[i]->size();
	(*fout_b) << "#####..Ebeam(MeV)= " << ebeam/MeV << " Theta(deg)= " << angle[i] << G4endl;
        (*fout_b) << "ve/create X(" << k << ") R ";
        for(jj=0; jj<k; jj++) {
          (*fout_b) << energy[jj] << " ";
        }
        (*fout_b) << endl;
        (*fout_b) << "ve/create Y(" << k << ") R ";

        for(jj=0; jj<k; jj++) {
          (*fout_b) << (*cross[i])[jj] << " ";
	}
        (*fout_b) << endl;
        (*fout_b) << "ve/create Z(" << k << ") R ";
        for(jj=0; jj<k; jj++) {
          (*fout_b) << (*err[i])[jj] << " ";
	}
        (*fout_b) << endl;
      }
      (*fout_b) << endl;
    }
  } while (end);
  (*fout_b) << "#####..End..#####" << G4endl;
}






