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

  // subroutines
  void ablainit_();
  //      SUBROUTINE DIRECT(ZPRF,A,EE,JPRF,PROBP,PROBN,PROBA,PROBF,PTOTL,   
  //     +		SN,SBP,SBA,ECN,ECP,ECA,BP,BA,INTTYPE,INUM,itest)        
  void direct_(double *ZPRF, double *A, double *EE, double *JPRF, double *PROBP, double *PROBN,
	       double *PROBA, double *PROBF, double *PTOTL, double *SN, double *SBP, double *SBA, double *ECN, double *ECP, double *ECA, double *BP, double * BA, int *INTTYPE, int *INUM, int *itest);
  void densniv_(double *A, double *Z,double *EE, double *ESOUS, double *DENS, double *BSHELL, double *BS, double *BK, double *TEMP, int *OPTSHP, int *OPTCOL, double *DEFBET);
  void fission_distri__(float *A, float *Z, float *E, float *A1, float *Z1, float *E1, float *A2, float *Z2, float *E2);
}

class UnitTester {
public:
  UnitTester() {
  }

  ~UnitTester() {
  }

  void init() {
    ablainit_();
    rseed_.iseed1=666;
    rseed_.iseed2=777;
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
    //  G4Incl *incl = new G4Incl(hazard, calincl, ws, mat, varntp);
    //abla = new G4Abla(hazard, volant, varntp);
    //abla->setVerboseLevel(0);
    //abla->initEvapora();
    rndm = new G4Ranecu();
    fission = new G4AblaFission(hazard, rndm);
    fission->about();
  }
  
  void testDensniv() {
    double a = 208, z = 82, ee = 300.0, esous = 1.0, bshell = 1.0, dens = 0.0, bs = 1.0, bk = 1.0, temp = 0.0, defbet = 1.0;
    int optshp = 1, optcol = 1;

    double A, Z, EE, ESOUS, DENS, BSHELL, BS, BK, TEMP, DEFBET;
    int OPTSHP, OPTCOL;
    double cA, cZ, cEE, cESOUS, cDENS, cBSHELL, cBS, cBK, cTEMP, cDEFBET;
    int cOPTSHP, cOPTCOL;

    A = a; cA = a;
    Z = z; cZ = z;
    EE = ee; cEE = ee;
    ESOUS = esous; cESOUS = esous; 
    DENS = dens; cDENS = dens;
    BSHELL = bshell; cBSHELL = bshell;
    BS = bs; cBS = bs;
    BK = bk; cBK = bk;
    TEMP = temp; cTEMP = temp;
    DEFBET = defbet; cDEFBET = defbet;
    OPTSHP = optshp; cOPTSHP = optshp; 
    OPTCOL = optcol; cOPTCOL = optcol; 

    G4cout <<"a/C++ b/FORTRAN" << G4endl;
    for(G4int i = 0; i < 30; i++) {
      densniv_(&A, &Z, &EE, &ESOUS, &DENS, &BSHELL, &BS, &BK, &TEMP, &OPTSHP, &OPTCOL, &DEFBET);
      abla->densniv(cA, cZ, cEE, cESOUS, &cDENS, cBSHELL, cBS, cBK, &cTEMP, cOPTSHP, cOPTCOL, cDEFBET);
      EE += 100.0;
      cEE = EE;
      G4cout <<"temp: C++: " << cTEMP << " FORTRAN:  " << TEMP << G4endl;
      G4cout <<"dens: C++: " << cDENS << " FORTRAN:  " << DENS << G4endl;
    }
  }

  void testDirect() {
// void G4Abla::direct(G4double zprf, G4double a, G4double ee, G4double jprf, 
// 		    G4double *probp_par, G4double *probn_par, G4double *proba_par, 
// 		    G4double *probf_par, G4double *ptotl_par, G4double *sn_par,
// 		    G4double *sbp_par, G4double *sba_par, G4double *ecn_par, 
// 		    G4double *ecp_par,G4double *eca_par, G4double *bp_par,
// 		    G4double *ba_par, G4int inttype, G4int inum, G4int itest)
    G4double zprf = 82, a = 208, ee = 300.0, jprf = 25;
    G4double probp_par, probn_par, proba_par; 
    G4double probf_par, ptotl_par, sn_par;
    G4double sbp_par, sba_par, ecn_par;
    G4double ecp_par, eca_par, bp_par, ba_par;
    G4int inttype = 0, inum = 1, itest = 0;

    G4double cZPRF, cA, cEE, cJPRF;
    G4double cPROBP_PAR, cPROBN_PAR, cPROBA_PAR; 
    G4double cPROBF_PAR, cPTOTL_PAR, cSN_PAR;
    G4double cSBP_PAR, cSBA_PAR, cECN_PAR;
    G4double cECP_PAR, cECA_PAR, cBP_PAR, cBA_PAR;
    G4int cINTTYPE, cINUM, cITEST;

    G4double ZPRF, A, EE, JPRF;
    G4double PROBP_PAR, PROBN_PAR, PROBA_PAR; 
    G4double PROBF_PAR, PTOTL_PAR, SN_PAR;
    G4double SBP_PAR, SBA_PAR, ECN_PAR;
    G4double ECP_PAR, ECA_PAR, BP_PAR, BA_PAR;
    G4int INTTYPE, INUM, ITEST;
    ZPRF = zprf; cZPRF = zprf;
    A = a; cA = a;
    EE = ee; cEE = ee;
    JPRF = jprf; cJPRF = jprf;
    INTTYPE = inttype; cINTTYPE = inttype;
    INUM = inum; cINUM = inum;
    ITEST = itest; cITEST = itest;
    G4cout <<"Calling C++" << G4endl;
    abla->direct(cZPRF, cA, cEE, cJPRF, &cPROBP_PAR, &cPROBN_PAR, &cPROBA_PAR, 
		 &cPROBF_PAR, &cPTOTL_PAR, &cSN_PAR,
		 &cSBP_PAR, &cSBA_PAR, &cECN_PAR, 
		 &cECP_PAR, &cECA_PAR, &cBP_PAR,
		 &cBA_PAR, cINTTYPE, cINUM, cITEST);
    G4cout <<"Calling FORTRAN" << G4endl;
    direct_(&ZPRF, &A, &EE, &JPRF, &PROBP_PAR, &PROBN_PAR, &PROBA_PAR, 
		 &PROBF_PAR, &PTOTL_PAR, &SN_PAR,
		 &SBP_PAR, &SBA_PAR, &ECN_PAR, 
		 &ECP_PAR, &ECA_PAR, &BP_PAR,
		 &BA_PAR, &INTTYPE, &INUM, &ITEST);
    G4cout <<"a/C++ b/FORTRAN" << G4endl;
    G4cout <<"A = " << A << G4endl;
    G4cout <<"cPROBP = " << cPROBP_PAR << G4endl;
    G4cout <<"cPROBN = " << cPROBN_PAR << G4endl;
    G4cout <<"cPROBA = " << cPROBA_PAR << G4endl;
    G4cout <<"cPROBF = " << cPROBF_PAR << G4endl;
    G4cout <<"cPTOTL = " << cPTOTL_PAR << G4endl;
    G4cout <<"cSN_PAR = " << cSN_PAR << G4endl;
    G4cout <<"cSBP_PAR = " << cSBP_PAR << G4endl;
    G4cout <<"cSBA_PAR = " << cSBA_PAR << G4endl;
    G4cout <<"cECN_PAR = " << cECN_PAR << G4endl;
    G4cout <<"cECP_PAR = " << cECP_PAR << G4endl;
    G4cout <<"cECA_PAR = " << cECA_PAR << G4endl;
    G4cout <<"cBP_PAR = " << cBP_PAR << G4endl;
    G4cout <<"cBA_PAR = " << cBA_PAR << G4endl;
    G4cout <<"===========================================" << G4endl;
    G4cout <<"PROBP = " << PROBP_PAR << G4endl;
    G4cout <<"PROBN = " << PROBN_PAR << G4endl;
    G4cout <<"PROBA = " << PROBA_PAR << G4endl;
    G4cout <<"PROBF = " << PROBF_PAR << G4endl;
    G4cout <<"PTOTL = " << PTOTL_PAR << G4endl;
    G4cout <<"SN_PAR = " << SN_PAR << G4endl;
    G4cout <<"SBP_PAR = " << SBP_PAR << G4endl;
    G4cout <<"SBA_PAR = " << SBA_PAR << G4endl;
    G4cout <<"ECN_PAR = " << ECN_PAR << G4endl;
    G4cout <<"ECP_PAR = " << ECP_PAR << G4endl;
    G4cout <<"ECA_PAR = " << ECA_PAR << G4endl;
    G4cout <<"BP_PAR = " << BP_PAR << G4endl;
    G4cout <<"BA_PAR = " << BA_PAR << G4endl;
  }

  void testFissionDistri() {
    //    float a = 208, z = 82, e = 200.0, a1 = 0, z1 = 0, e1 = 0, a2 = 0, z2 = 0, e2 = 0;
    float a = 238, z = 92, e = 100.0, a1 = 0, z1 = 0, e1 = 0, a2 = 0, z2 = 0, e2 = 0;
    float A, Z, E, A1, Z1, E1, A2, Z2, E2;
    G4double cA, cZ, cE, cA1, cZ1, cE1, cV1, cA2, cZ2, cE2, cV2;
    A = a; cA = a;
    Z = z; cZ = z;
    E = e; cE = e;
    G4cout <<"A = " << a << " Z = " << z << " E = " << e << G4endl;
    //    G4cout.setw(8);

    G4cout << "i" << setw(4) 
	   <<  "A1" << setw(5) 
	   << "cA1" << setw(8)
	   << "Z1" << setw(5)
	   << "cZ1" << setw(12) 
	   << "E1" << setw(12)
	   << "cE1" << setw(10)
	   << "A2" << setw(5)
	   << "cA2" << setw(8)
	   << "Z2" << setw(5)
	   << "cZ2" << setw(8)
	   << "E2" << setw(8)
	   << "cE2" << setw(4)
	   << "dA1" << setw(4)
	   << "dZ1" << setw(8)
	   << "dE1" << setw(4)
	   << "dA2" << setw(4)
	   << "dZ2" << setw(8)
	   << "dE2" << G4endl;

    for(int i = 0; i < 3; i++) {
      fission_distri__(&A, &Z, &E, &A1, &Z1, &E1, &A2, &Z2, &E2);
      fission->fissionDistri(cA, cZ, cE, cA1, cZ1, cV1, cE1, cA2, cZ2, cE2, cV2);
      G4cout << i << setw(4) 
	     <<  A1 << setw(5) 
	     << cA1 << setw(8)
	     << Z1 << setw(5)
	     << cZ1 << setw(12) 
	     << E1 << setw(12)
	     << cE1 << setw(10)
	     << A2 << setw(5)
	     << cA2 << setw(8)
	     << Z2 << setw(5)
	     << cZ2 << setw(8)
	     << E2 << setw(8)
	     << cE2 << setw(4)
	     << (A1 - cA1) << setw(4)
	     << (Z1 - cZ1) << setw(8)
	     << (E1 - cE1) << setw(4)
	     << (A2 - cA2) << setw(4)
	     << (Z2 - cZ2) << setw(8)
	     << (E2 - cE2)
	     << G4endl;
      
//       G4cout <<"A = " << A << " Z = " << Z << " E = " << E << G4endl;
//       G4cout <<"A1 = " << A1 << " Z1 = " << Z1 << " E1 = " << E1 << G4endl;
//       G4cout <<"A2 = " << A2 << " Z2 = " << Z2 << " E2 = " << E2 << G4endl;
//       G4cout <<"========================================" << G4endl;
//       G4cout <<"cA1 = " << cA1 << " cZ1 = " << cZ1 << " cE1 = " << cE1 << " cV1 = " << cV1 << G4endl;
//       G4cout <<"cA2 = " << cA2 << " cZ2 = " << cZ2 << " cE2 = " << cE2 << " cV2 = " << cV2 << G4endl;
//       G4cout << G4endl;
    }
  }

private:
  G4Abla *abla;
  G4AblaFission *fission;
  G4Hazard *hazard;
  G4VarNtp *varntp;
  G4Volant *volant;
  G4Ranecu *rndm;
};

int main(int *argc, char **argv[])
{
  UnitTester *tester = new UnitTester();
  tester->init();
  //  tester-> testDensniv();
  tester-> testFissionDistri();
  //tester-> testDirect();
  return 0;
}
