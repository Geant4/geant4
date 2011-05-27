#include <iostream>
#include <string>
#include <vector>
#include "globals.hh"
#include "G4Types.hh"
#include "G4Incl.hh"
#include "G4Abla.hh"
#include "G4AblaFission.hh"
#include "G4InclRandomNumbers.hh"
#include "G4Ranecu.hh"

// FORTRAN bindings
extern "C" {
  //common/rseed/iseed1,iseed2
  struct {
    long iseed1;
    long iseed2;
  } rseed_;
  //common/hazard/ial,IY
  struct {
    long ial;
    long iy[19];
  } hazard_;
  
  // COMMON /EC2SUB/ ECNZ                                              
//       REAL*8 ECNZ                                                       
//       DIMENSION ECNZ(0:153,0:98)                                       
  struct {
    double ecnz[154][99];
  } ec2sub_;

  // subroutines
  void ablainit_();
  //      SUBROUTINE DIRECT(ZPRF,A,EE,JPRF,PROBP,PROBN,PROBA,PROBF,PTOTL,   
  //     +		SN,SBP,SBA,ECN,ECP,ECA,BP,BA,INTTYPE,INUM,itest)        
  void direct_(double *ZPRF, double *A, double *EE, double *JPRF, double *PROBP, double *PROBN,
	       double *PROBA, double *PROBF, double *PTOTL, double *SN, double *SBP, double *SBA, double *ECN, double *ECP, double *ECA, double *BP, double * BA, int *INTTYPE, int *INUM, int *itest);
  void densniv_(double *A, double *Z,double *EE, double *ESOUS, double *DENS, double *BSHELL, double *BS, double *BK, double *TEMP, int *OPTSHP, int *OPTCOL, double *DEFBET);
  void fission_distri__(float *A, float *Z, float *E, float *A1, float *Z1, float *E1, float *A2, float *Z2, float *E2);
}

G4Ec2sub *ec2sub;

void printTable()
{
  std::cout <<"Frldm (shell effect) table comparison." << std::endl;
  std::cout <<"i \t j \t FORTRAN \t C++ \t diff "<< std::endl;
  for(int i = 0; i < 154; i++) {
    for(int j = 0; j < 99; j++) {
      std::cout << i << setw(5) << j << setw(5) << ec2sub_.ecnz[i][j] << setw(20) << ec2sub->ecnz[i][j] << setw(20) << ec2sub_.ecnz[i][j] - ec2sub->ecnz[i][j] << std::endl;
    }
  }
}

int main(int argc, char *argv[])
{
  ablainit_();
  G4Abla *abla;
  G4AblaFission *fission;
  G4Hazard *hazard;
  G4VarNtp *varntp;
  G4Volant *volant;
  G4Ranecu *rndm;

  hazard = new G4Hazard();
  varntp = new G4VarNtp();
  volant = new G4Volant();
    
    hazard->ial = 979678188;
    hazard_.ial = hazard->ial;
    hazard->igraine[0] = 3997;
    hazard->igraine[1] = 15573;
    hazard->igraine[2] = 9971;
    hazard->igraine[3] = 9821; 
    hazard->igraine[4] = 99233; 
    hazard->igraine[5] = 11167; 
    hazard->igraine[6] = 12399;
    hazard->igraine[7] = 11321; 
    hazard->igraine[8] = 9825;
    hazard->igraine[9] = 2587; 
    hazard->igraine[10] = 1775;
    hazard->igraine[11] = 56799; 
    hazard->igraine[12] = 1156;
    //  hazard->igraine[13] = 11207;
    hazard->igraine[13] = 38957; 
    hazard->igraine[14] = 35779; 
    hazard->igraine[15] = 10055; 
    hazard->igraine[16] = 76533; 
    hazard->igraine[17] = 33759;
    hazard->igraine[18] = 13227;
    for(int i = 0; i < 19; i++) {
      hazard_.iy[i] = hazard->igraine[i];
    }
    abla = new G4Abla(hazard, volant, varntp);
    abla->setVerboseLevel(4);
    abla->initEvapora();
    ec2sub = abla->getFrldmTable();
    printTable();

  return 0;
}
