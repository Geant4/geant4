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
            << G4endl;


  //there can't be lines longer than nmax characters
  const int nmax = 100;
  char line[nmax];
  std::string line1, line2, word1, word2, word3;
  G4bool end = true;
  double ebeam = 0.0;
  std::vector<int>    inn;
  std::vector<double> angle;
  std::vector<double> energy;
  std::vector<double> cross;
  std::vector<double> err;

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
      (*fin) >> ebeam;

      double an, a, e, s, s1;
      int i, j, jj, nbin;
      (*fin) >> word1 >> word2 >> nbin;
      (*fin) >> word1 >> word2 >> word3 >> line1;
      (*fin) >> word1 >> word2 >> word3 >> line1;
      inn.clear();
      angle.clear();
      energy.clear();
      cross.clear();
      err.clear();
      an = 0.;
      j = 0;
      for(i=0; i<nbin; i++) {
        (*fin) >> a >> e >> s >> s1;
	s1 *= 0.01*s;
        if(1 < verbose) {
          cout << "an= " << a << "e= " << e << " cross= " << s << " +- " << s1 << endl;
        }
	if(an != a) {
          if(an>0.) inn.push_back(j);
	  j = 0;
          (*fout_c) << "#####..Ebeam(MeV)= " << ebeam/MeV << " Theta(deg)= " << a << G4endl;
	  an = a;
 	  angle.push_back(a);
	}
	j++;
	energy.push_back(e);
	cross.push_back(s);
	err.push_back(s1);

        (*fout_c) << j << ".   e(MeV)= " << e
	          << "  cross(mb/sr/MeV)= "
	          << s << " +- " << s1 << endl;
      }
      inn.push_back(j);
      int jbin = inn.size();
      i = 0;
      for(j=0; j<jbin; j++) {
        int k = inn[j];
        (*fout_b) << "#####..Ebeam(MeV)= " << ebeam/MeV << " Theta(deg)= " << angle[j] << G4endl;
        (*fout_b) << "ve/create X(" << k << ") R ";
        for(jj=0; jj<k; jj++) {
          (*fout_b) << energy[i+jj] << " ";
        }
        (*fout_b) << endl;
        (*fout_b) << "ve/create Y(" << k << ") R ";

        for(jj=0; jj<k; jj++) {
          (*fout_b) << cross[i+jj] << " ";
	}
        (*fout_b) << endl;
        (*fout_b) << "ve/create Z(" << k << ") R ";
        for(jj=0; jj<k; jj++) {
          (*fout_b) << err[i+jj] << " ";
	}
        (*fout_b) << endl;
	i += k;
      }
      (*fout_b) << endl;
    }
  } while (end);
  (*fout_b) << "#####..End..#####" << G4endl;
  (*fout_c) << "#####..End..#####" << G4endl;
}






