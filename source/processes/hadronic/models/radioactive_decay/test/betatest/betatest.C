{
/*
#include <stdlib.h>
#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
//#include "

int main ()
*/

int    Ap = 24;
int    Zp = 11;
double Xp = -8.4175;

int    Ad = 24;
int    Zd = 12;
double Xd = -13.9306;

double  Q = 1.0;

double Me = 0.511;

double Md = Zd * 938.2796 + (Ad-Zd) * 939.5731 - Xd;
double Mp = Md + Me + Q;

double rd1, rd2, rd;
double pd[3];
double pmax = 0.0;
double psum = 0.0;
double e[3] = {0.0,0.0,0.0};
int i;
ofstream OutputFile1("sampled_org.data", ios::out);
for ( i= 0; i < 10; i++) {
  cout << i << endl;
  do {
    rd1 = rand();
    rd1 = rd1 / RAND_MAX;
    rd2 = rand();
    rd2 = rd2 / RAND_MAX;
    if (rd2 > rd1)
      {
	rd  = rd1;
	rd1 = rd2;
	rd2 = rd;
      }
    pmax = 0.0;
    psum = 0.0;
    
    e[0] = rd2 * Q;
    pd[0] = sqrt(e[0]*e[0] + 2.0*e[0]*Me);
    if (pd[0] > pmax) pmax = pd[0];
    psum += pd[0];
    
    e[1] = (rd1-rd2) * Q;
    pd[1] = e[1];
    if (pd[1] > pmax) pmax = pd[1];
    psum += pd[1];
    
    e[2] = (1.-rd1) * Q;
    pd[2] = sqrt(e[2]*e[2] + 2.0*e[2]*Md);
    if (pd[2] > pmax) pmax = pd[2];
    psum += pd[2];
  } while (pmax > psum-pmax);
  OutputFile1 << setw(12) <<e[0]
	      << setw(12) <<e[1]
	      << setw(12) <<e[2]
	      << endl;
  }
OutputFile1.close();

ofstream OutputFile2("sampled_new.data", ios::out);

for ( i=0; i<100000; i++) {
  do {
    rd1 = rand();
    rd1 = rd1 / RAND_MAX;
      rd2 = rand();
      rd2 = rd2 / RAND_MAX;
//    if (rd2 > rd1)
//    {
//      rd  = rd1;
//      rd1 = rd2;
//      rd2 = rd;
//    }
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
    } while (pmax > psum-pmax);
    OutputFile2 << setw(12) <<e[0]
                << setw(12) <<e[1]
                << setw(12) <<e[2]
                << endl;
  }
OutputFile2.close();
}











