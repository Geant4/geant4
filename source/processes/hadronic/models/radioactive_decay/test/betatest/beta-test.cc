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
#include <stdlib.h>
#include <string>
#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>

using namespace std;

#include "G4FermiFunction.cc"

int main ()
{
  int    Ap = 210;
  int    Zp = 83;
  double Xp = -14.801;

  int    Ad = 210;
  int    Zd = 84;
  double Xd = -15.963;

  double  Q = 1.1615;

  double Me = 0.511;

  double rd1, rd2, rd;
  double pd[3];
  double pmax = 0.0;
  double psum = 0.0;
  double e[3] = {0.0,0.0,0.0};
  double Ntotal = 100;
  string filename;
  char mode;

  cout << " Parent isotope (A,Z):" << endl;
  cin >>Ap>>Zp ;
  cout << " beta + or - decay:" << endl;
  cin>>mode;
  cout <<" The end point energy (MeV):" << endl;
  cin >>Q;
  cout <<" The daughter mass deficit (MeV):" << endl;
  cin >>Xd;
  cout << " Number of decays to be simulated:"<< endl;
  cin >>Ntotal;
  cout << " Name of the output file" << endl;
  cin>> filename;

  ofstream OutputFile2("test-beta.data", ios::out);

  Ad = Ap;
  Zd = Zp+1;
  if (mode == '+') Zd = -(Zp-1);
  double Md = fabs(Zd) * 938.2796 + (Ad-fabs(Zd)) * 939.5731 - Xd;
  G4FermiFunction* aFermiFunction = new G4FermiFunction (Ad, Zd);
  
  double FermiFN  = aFermiFunction->GetFFN (Q/Me);
  cout << " Q = " << Q << endl;
  cout<< " FermiFN =  " <<FermiFN<<endl;
  
  for (int i=0; i<Ntotal; i++) {
    if (!(i%100)) cout <<i <<endl;
    do {
      rd1 = rand();
      rd1 = rd1 / RAND_MAX;
      rd2 = rand();
      rd2 = rd2 / RAND_MAX;
      pmax = 0.0;
      psum = 0.0;

//    pd[0] = pow(rd2,1.0/3.0) * sqrt((Q+2.0*Me)*Q);
      pd[0] = sqrt(rd2) * sqrt((Q+2.0*Me)*Q);
      e[0]  = sqrt(pd[0]*pd[0] + Me*Me) - Me;
      if (pd[0] > pmax) pmax = pd[0];
      psum += pd[0];

//    e[1] = pow(rd1,1.0/3.0) * Q;
      e[1]  = sqrt(rd1) * Q;
      pd[1] = e[1];
      if (pd[1] > pmax) pmax = pd[1];
      psum += pd[1];
      
      e[2]  = Q - e[0] - e[1];
      if (e[2]>0.0) {
	pd[2] = sqrt(e[2]*e[2] + 2.0*e[2]*Md);
	if (pd[2] > pmax) pmax = pd[2];
	psum += pd[2];
      } else {
	pmax = Q;
	psum = Q;
      }
      e[0] = e[0]/0.511;
      double fermif = aFermiFunction->GetFF(e[0])/FermiFN;
      //    cout << fermif << endl;
      if ((rand() / RAND_MAX) > fermif) pmax = psum =  Q;
      
    } while (pmax > psum-pmax);
    
    OutputFile2 << setw(12) <<e[0]*511.
		<< setw(12) <<e[1]
		<< setw(12) <<e[2]
		<< endl;
  }
  OutputFile2.close();}






