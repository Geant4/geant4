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
// -------------------------------------------------------------------
// $Id: RunAction.cc,v 1.3 2008/04/10 12:02:57 sincerti Exp $
// -------------------------------------------------------------------

#include "G4SteppingManager.hh"
#include "G4Run.hh"
#include "G4Material.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>
#include <assert.h>

#include "RunAction.hh"

// MATRIX
#define MATRIX_BOUND_CHECK
#include "globals.hh"
#include <iostream>
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"
#include <fstream>
#include "G4ios.hh"
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::RunAction(DetectorConstruction* det)
:detector(det)
{   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginOfRunAction(const G4Run* /*aRun*/)
{
  // Cleaning result files
  system ("rm -rf  ./results/x.txt");
  system ("rm -rf  ./results/y.txt");
  system ("rm -rf  ./results/theta.txt");
  system ("rm -rf  ./results/phi.txt");
  system ("rm -rf  ./results/image.txt");
  system ("rm -rf  ./results/matrix.txt");
  system ("rm -rf  ./results/profile.txt");
  system ("rm -rf  ./results/grid.txt");
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndOfRunAction(const G4Run* /*aRun*/)
{

if (detector->GetCoef()==1)
{

	CLHEP::HepMatrix m,s;
	m = CLHEP::HepMatrix(32,32);
	s = CLHEP::HepMatrix(32,32);
	CLHEP::HepVector v;
	v = CLHEP::HepVector(32);

	float aa,ab,ac,ad,ae,af,ag,ah,ai,aj;
	float ba,bb,bc,bd,be,bf,bg,bh,bi,bj;
	float ca,cb,cc,cd,ce,cf,cg,ch,ci,cj;
	float da,db;
	float vv;

	// MATRIX READING

	int ncols;
	int nlines=1;

	FILE * fp1 = fopen("results/matrix.txt","r");
	while (1) 
	{
	ncols = fscanf(fp1,
      	"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
      	&aa,&ab,&ac,&ad,&ae,&af,&ag,&ah,&ai,&aj,
      	&ba,&bb,&bc,&bd,&be,&bf,&bg,&bh,&bi,&bj,
      	&ca,&cb,&cc,&cd,&ce,&cf,&cg,&ch,&ci,&cj,
	&da,&db     
	);

      	if (ncols<0) break;
	if (nlines>32) G4Exception("Try to read more than 32 lines in matrix file !");

	m(nlines,1)=aa;
      	m(nlines,2)=ab;
      	m(nlines,3)=ac;
      	m(nlines,4)=ad;
      	m(nlines,5)=ae;
      	m(nlines,6)=af;
      	m(nlines,7)=ag;
      	m(nlines,8)=ah;
      	m(nlines,9)=ai;
      	m(nlines,10)=aj;
      	m(nlines,11)=ba;
      	m(nlines,12)=bb;
      	m(nlines,13)=bc;
      	m(nlines,14)=bd;
      	m(nlines,15)=be;
      	m(nlines,16)=bf;
      	m(nlines,17)=bg;
      	m(nlines,18)=bh;
      	m(nlines,19)=bi;
      	m(nlines,20)=bj;
      	m(nlines,21)=ca;
      	m(nlines,22)=cb;
      	m(nlines,23)=cc;
      	m(nlines,24)=cd;
      	m(nlines,25)=ce;
      	m(nlines,26)=cf;
      	m(nlines,27)=cg;
      	m(nlines,28)=ch;
      	m(nlines,29)=ci;
      	m(nlines,30)=cj;
      	m(nlines,31)=da;
      	m(nlines,32)=db;

      	nlines++;
	
    	}
	fclose(fp1);

	// VECTOR READING


	G4cout << G4endl;
	G4cout << "===> NANOBEAM LINE INTRINSIC ABERRATION COEFFICIENTS (units of micrometer and mrad) :" << G4endl;
	G4cout << G4endl;

	int inv;
	nlines = 1;
	FILE * fp2 = fopen("results/x.txt","r");
	while (1) 
	{
	  ncols = fscanf(fp2,"%f", &vv);
          if (ncols<0) break;
          v(nlines)=vv;
          nlines++;
        }
	fclose(fp2);

	m.invert(inv);
	CLHEP::HepVector seb(32,0);
	seb=m*v;
	CLHEP::HepVector b;
	b=seb.sub(2,2);   G4cout << "<x|theta>=" << b << G4endl;
	b=seb.sub(8,8);   G4cout << "<x|theta*delta>=" << b << G4endl;
	b=seb.sub(10,10); G4cout << "<x|theta^3>=" << b << G4endl;
	b=seb.sub(12,12); G4cout << "<x|theta*phi^2>=" << b << G4endl;
	m.invert(inv);

	nlines = 1;
	FILE * fp3 = fopen("results/theta.txt","r");
	while (1) 
   	{
      	  ncols = fscanf(fp3,"%f", &vv);
          if (ncols<0) break;
          v(nlines)=vv;
          nlines++;
    	}
	fclose(fp3);

	m.invert(inv);
	seb = m*v;
	m.invert(inv);
	b=seb.sub(2,2); G4cout << "<x|x>=" << b << G4endl;

	nlines = 1;
	FILE * fp4 = fopen("results/y.txt","r");
	while (1) 
   	{
          ncols = fscanf(fp4,"%f", &vv);
          if (ncols<0) break;
          v(nlines)=vv;
          nlines++;
    	}
	fclose(fp4);
	m.invert(inv);
	seb=m*v;
	b=seb.sub(3,3);   G4cout << "<y|phi>=" << b << G4endl;
	b=seb.sub(9,9);   G4cout << "<y|phi*delta>=" << b << G4endl;
	b=seb.sub(11,11); G4cout << "<y|theta^2*phi>=" << b << G4endl;
	b=seb.sub(13,13); G4cout << "<y|phi^3>=" << b << G4endl;
	m.invert(inv);

	nlines = 1;
	FILE * fp5 = fopen("results/phi.txt","r");
	while (1) 
   	{
          ncols = fscanf(fp5,"%f", &vv);
          if (ncols<0) break;
          v(nlines)=vv;
          nlines++;
    	}
	fclose(fp5);
	m.invert(inv);
	seb = m*v;
	m.invert(inv);
	b=seb.sub(3,3); G4cout << "<y|y>=" << b << G4endl;

}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
