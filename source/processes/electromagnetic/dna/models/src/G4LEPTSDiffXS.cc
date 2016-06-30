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
// AMR Simplification /4
// read Diff XSection & Interpolate
#include <string.h>
#include <stdio.h>
#include <string>

#include <cmath>
#include "globals.hh"
#include <iostream> 
using namespace std;
#include "CLHEP/Units/PhysicalConstants.h"

#include "G4LEPTSDiffXS.hh"
#include "G4Exp.hh"

G4LEPTSDiffXS::G4LEPTSDiffXS(string file) {
  fileName = file;

  readDXS();
  BuildCDXS();
  //BuildCDXS(1.0, 0.5);
  NormalizeCDXS();
  InterpolateCDXS();
}



//DXS y KT
void G4LEPTSDiffXS::readDXS( ) {

  FILE   *fp;
  float data, data2;

  if ((fp=fopen(fileName.c_str(), "r"))==NULL){
    //G4cout << "Error reading " << fileName << G4endl;
    NumEn = 0;
    bFileFound = false;
    return;
  }

  bFileFound = true;

  //G4cout << "Reading2 " << fileName << G4endl;

  //NumAng = 181;
  fscanf(fp, "%d %d %s", &NumAng, &NumEn, DXSTypeName);
  if( !strcmp(DXSTypeName, "KTC") )     DXSType = 2;  // read DXS & calculate KT
  else if( !strcmp(DXSTypeName, "KT") ) DXSType = 1;  // read DXS & KT
  else                                  DXSType = 0;

  //  if( verboseLevel >= 1 ) G4cout << "Read DXS   (" << fileName  <<  ")\t NEg " << NumEn << " NAng " << NumAng
  //       << "DXSType " << DXSTypeName << " " << DXSType << G4endl;

  for (G4int eBin=1; eBin<=NumEn; eBin++){
    fscanf(fp,"%f ",&data);
    Eb[eBin] = (G4double)data;
  }


  //for (aBin=1;aBin<NumAng;aBin++){

  if(DXSType==1) {
    G4cout << "DXSTYpe 1" << G4endl;
    for (G4int aBin=0;aBin<NumAng;aBin++){
      fscanf(fp,"%f ",&data);
      DXS[0][aBin]=(G4double)data;
      for (G4int eBin=1;eBin<=NumEn;eBin++){
	fscanf(fp,"%f %f ",&data2, &data);
	DXS[eBin][aBin]=(G4double)data;
	KT[eBin][aBin]=(G4double)data2;
      }
    }
  }
  else {
    for(G4int aBin=0; aBin<NumAng; aBin++){
      for(G4int eBin=0; eBin<=NumEn; eBin++){
	fscanf(fp,"%f ",&data);
	DXS[eBin][aBin] = (G4double)data;
      }
    }
    for(G4int aBin=0; aBin<NumAng; aBin++){
      for(G4int eBin=1; eBin<=NumEn; eBin++){
	G4double A = DXS[0][aBin];                         // Angle
	G4double E = Eb[eBin];                             // Energy
	G4double p = sqrt(pow( (E/27.2/137),2) +2*E/27.2); // Momentum
	KT[eBin][aBin] = p *sqrt(2.-2.*cos(A*CLHEP::twopi/360.)); // Mom. Transfer
	//G4cout << "aEpKt " << aBin << " " << A << " E " << E << " p " << p << " KT "
	//   << KT[eBin][aBin] << " DXS " << DXS[eBin][aBin] << G4endl;
      }
    }
  }

  fclose(fp);
}



// CDXS from DXS
void G4LEPTSDiffXS::BuildCDXS(G4double E, G4double El) {

  for(G4int aBin=0;aBin<NumAng;aBin++) {
    for(G4int eBin=0;eBin<=NumEn;eBin++){
      CDXS[eBin][aBin]=0.0;
    }
  }

  for(G4int aBin=0;aBin<NumAng;aBin++)
    CDXS[0][aBin] = DXS[0][aBin];

  for (G4int eBin=1;eBin<=NumEn;eBin++){
    G4double sum=0.0;
    for (G4int aBin=0;aBin<NumAng;aBin++){
      sum += pow(DXS[eBin][aBin], (1.0-El/E) );
      CDXS[eBin][aBin]=sum;
    }
  }
}



// CDXS from DXS
void G4LEPTSDiffXS::BuildCDXS() {

  BuildCDXS(1.0, 0.0); // El = 0
}



// CDXS & DXS
void G4LEPTSDiffXS::NormalizeCDXS() {

  // Normalize:  1/area
  for (G4int eBin=1; eBin<=NumEn; eBin++){
    G4double area = CDXS[eBin][NumAng-1];
    //G4cout << eBin << " area = " << area << G4endl;

    for (G4int aBin=0; aBin<NumAng; aBin++) {
      CDXS[eBin][aBin] /= area;
      //DXS[eBin][aBin]  /= area;
    }
  }
}



//ICDXS from CDXS   & IKT from KT
void G4LEPTSDiffXS::InterpolateCDXS() {  // *10 angles, linear

  G4double eps = 1e-5;
  G4int ia = 0;

  for( G4int aBin=0; aBin<NumAng-1; aBin++) {
    G4double x1 = CDXS[0][aBin] + eps;
    G4double x2 = CDXS[0][aBin+1] + eps;
    G4double dx = (x2-x1)/100;

    //if( x1<10 || x1) dx = (x2-x1)/100;

    for( G4double x=x1; x < (x2-dx/10); x += dx) {
      for( G4int eBin=0; eBin<=NumEn; eBin++) {
	G4double y1 = CDXS[eBin][aBin];
	G4double y2 = CDXS[eBin][aBin+1];
	G4double z1 = KT[eBin][aBin];
	G4double z2 = KT[eBin][aBin+1];

	if( aBin==0 ) {
	  y1 /=100;
	  z1 /=100;
	}

	if( eBin==0 ) {	  //linear abscisa
	  ICDXS[eBin][ia] = (y1*(x2-x) + y2*(x-x1))/(x2-x1);
	}
	else {           //log-log ordenada
	  ICDXS[eBin][ia] = G4Exp( (log(y1)*log(x2/x)+log(y2)*log(x/x1))/log(x2/x1) );
	}

	IKT[eBin][ia] = (z1*(x2-x) + z2*(x-x1))/(x2-x1);
	//IKT[eBin][ia] = exp( (log(z1)*log(x2/x)+log(z2)*log(x/x1))/log(x2/x1) );
      }

      ia++;
    }

  }

  INumAng = ia;
}



// from ICDXS
#include "Randomize.hh"
G4double G4LEPTSDiffXS::SampleAngle(G4double Energy) {
  G4int  ii,jj,kk=0, Ebin;

  Ebin=1;
  for(ii=2; ii<=NumEn; ii++)
    if(Energy >= Eb[ii])
      Ebin=ii;
  if(Energy > Eb[NumEn]) Ebin=NumEn;
  else if(Energy > (Eb[Ebin]+Eb[Ebin+1])*0.5 ) Ebin++;

  //G4cout << "SampleAngle E " << Energy << " Ebin " << Ebin << " E[] " << Eb[Ebin] << G4endl;

  ii=0;
  jj=INumAng-1;
  G4double rnd=G4UniformRand();

  while ((jj-ii)>1){
    kk=(ii+jj)/2;
    G4double dxs = ICDXS[Ebin][kk];
    if (dxs < rnd) ii=kk;
    else           jj=kk;
  }


  //G4double x = ICDXS[0][jj];
  G4double x = ICDXS[0][kk] *CLHEP::twopi/360.;

  return(x);
}



G4double G4LEPTSDiffXS::SampleAngleEthylene(G4double E, G4double El) {

  BuildCDXS(E, El);
  NormalizeCDXS();
  InterpolateCDXS();

  return( SampleAngle(E) );
}



//Momentum Transfer formula
G4double G4LEPTSDiffXS::SampleAngleMT(G4double Energy, G4double Elost) {
  G4int  ii, jj, kk=0, Ebin, iMin, iMax;

  G4double Ei = Energy;
  G4double Ed = Energy - Elost;
  G4double Pi = sqrt( pow( (Ei/27.2/137),2) +2*Ei/27.2); //incidente
  G4double Pd = sqrt( pow( (Ed/27.2/137),2) +2*Ed/27.2); //dispersado
  G4double Kmin = Pi - Pd;
  G4double Kmax = Pi + Pd;

  if(Pd <= 1e-9 ) return (0.0);
 

  // locate Energy bin
  Ebin=1;
  for(ii=2; ii<=NumEn; ii++)
    if(Energy > Eb[ii]) Ebin=ii;
  if(Energy > Eb[NumEn]) Ebin=NumEn;
  else if(Energy > (Eb[Ebin]+Eb[Ebin+1])*0.5 ) Ebin++;
 
  //G4cout << "SampleAngle2 E " << Energy << " Ebin " << Ebin << " E[] " << Eb[Ebin] << G4endl;

  ii=0; jj=INumAng-1;
  while ((jj-ii)>1) {
    kk=(ii+jj)/2;
    if( IKT[Ebin][kk] < Kmin ) ii=kk;
    else                      jj=kk;
  }
  iMin = ii;

  ii=0; jj=INumAng-1;
  while ((jj-ii)>1) {
    kk=(ii+jj)/2;
    if( IKT[Ebin][kk] < Kmax ) ii=kk;
    else                      jj=kk;
  }
  iMax = ii;


  // r -> a + (b-a)*r = a*(1-r) + b*r
  G4double rnd = G4UniformRand();
  rnd = (1-rnd)*ICDXS[Ebin][iMin] + rnd*ICDXS[Ebin][iMax];
  //G4double rnd = (ICDXS[Ebin][iMax] - ICDXS[Ebin][iMin]) * G4UniformRand()
  //  +  ICDXS[Ebin][iMin];

  ii=0; jj=INumAng-1;
  while ((jj-ii)>1){
    kk=(ii+jj)/2;
    if( ICDXS[Ebin][kk] < rnd) ii=kk;
    else                      jj=kk;
  }

  //Sampled
  G4double KR = IKT[Ebin][kk];

  G4double co = (Pi*Pi + Pd*Pd - KR*KR) / (2*Pi*Pd); //cos ang. disp.
  if(co > 1) co =1;
  G4double x = acos(co); //*360/twopi;            //ang. dispers.

  // Elastic aprox.
  //x = 2*asin(KR/Pi/2)*360/twopi;

  return(x);
}



void G4LEPTSDiffXS::PrintDXS(G4int NE) {
// Debug
//#include <string>
//using namespace std;


  G4double dxs;

  G4cout << G4endl<< "DXS & CDXS: " << fileName << G4endl<< G4endl;

  for (G4int aBin=0; aBin<NumAng; aBin++) {
    if( aBin>0)
      dxs = (CDXS[NE][aBin] - CDXS[NE][aBin-1])/(CDXS[0][aBin] - CDXS[0][aBin-1]);
    else
      dxs = CDXS[NE][aBin];

    G4cout << CDXS[0][aBin] << " " << dxs << " " << CDXS[NE][aBin] << G4endl;
  }

  G4cout << G4endl<< "IDXS & ICDXS: " << fileName << G4endl<< G4endl;

  for (G4int aBin=0; aBin<INumAng; aBin++) {
    if( aBin>0)
      dxs = (ICDXS[NE][aBin] - ICDXS[NE][aBin-1])/(ICDXS[0][aBin] - ICDXS[0][aBin-1]);
    else
      dxs = ICDXS[NE][aBin];

    G4cout << ICDXS[0][aBin] << " " << dxs << " " << ICDXS[NE][aBin] << G4endl;
  }


  // if(jmGlobals->VerboseHeaders) {
  //   G4cout << G4endl << "dxskt1" << G4endl;
  //   for (G4int aBin=0;aBin<NumAng;aBin++){
  //     G4cout << DXS[0][aBin] << "\t" << DXS[1][aBin] << "\t" << DXS[2][aBin] << "\t"
  // 	     << CDXS[1][aBin] << "\t" << KT[12][aBin] << G4endl;
  //   }
  // }

}
