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
//
// -----------------------------------------------------------------------------
//
// GEANT4 Class implementation file
//
// File name:     G4GoudsmitSaundersonTable
//
// Author:        Mihaly Novak / (Omrane Kadri)
//
// Creation date: 20.02.2009
//
// Class description:
//   Class to handle multiple scattering angular distributions precomputed by
//   using Kawrakow-Bielajew Goudsmit-Saunderson MSC model based on the screened
//   Rutherford DCS for elastic scattering of electrons/positrons [1,2]. This
//   class is used by G4GoudsmitSaundersonMscModel to sample the angular
//   deflection of electrons/positrons after travelling a given path.
//
// Modifications:
// 04.03.2009 V.Ivanchenko cleanup and format according to Geant4 EM style
// 26.08.2009 O.Kadri: avoiding unuseful calculations and optimizing the root
//                     finding parameter error's within SampleTheta method
// 08.02.2010 O.Kadri: reduce delared variables; reduce error of finding root
//                     in secant method
// 26.03.2010 O.Kadri: minimum of used arrays in computation within the dichotomie
//                     finding method the error was the lowest value of uvalues
// 12.05.2010 O.Kadri: changing of sqrt((b-a)*(b-a)) with fabs(b-a)
// 18.05.2015 M. Novak This class has been completely replaced (only the original
//            class name was kept; class description was also inserted):
//            A new version of Kawrakow-Bielajew Goudsmit-Saunderson MSC model
//            based on the screened Rutherford DCS for elastic scattering of
//            electrons/positrons has been introduced[1,2]. The corresponding MSC
//            angular distributions over a 2D parameter grid have been recomputed
//            and the CDFs are now stored in a variable transformed (smooth) form
//            together with the corresponding rational interpolation parameters.
//            The new version is several times faster, more robust and accurate
//            compared to the earlier version (G4GoudsmitSaundersonMscModel class
//            that use these data has been also completely replaced)
// 28.04.2017 M. Novak: New representation of the angular distribution data with
//            significantly reduced data size.
// 23.08.2017 M. Novak: Added funtionality to handle Mott-correction to the
//            base GS angular distributions and some other factors (screening
//            parameter, first and second moments) when Mott-correction is
//            activated in the GS-MSC model.
//
// References:
//   [1] A.F.Bielajew, NIMB, 111 (1996) 195-208
//   [2] I.Kawrakow, A.F.Bielajew, NIMB 134(1998) 325-336
//
// -----------------------------------------------------------------------------

#include "G4GoudsmitSaundersonTable.hh"


#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

#include "G4GSMottCorrection.hh"
#include "G4MaterialTable.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"

#include "G4String.hh"

#include <fstream>
#include <cstdlib>
#include <cmath>

#include <iostream>
#include <iomanip>

// perecomputed GS angular distributions, based on the Screened-Rutherford DCS
// are the same for e- and e+ so make sure we load them only onece
G4bool G4GoudsmitSaundersonTable::gIsInitialised = false;
//
std::vector<G4GoudsmitSaundersonTable::GSMSCAngularDtr*> G4GoudsmitSaundersonTable::gGSMSCAngularDistributions1;
std::vector<G4GoudsmitSaundersonTable::GSMSCAngularDtr*> G4GoudsmitSaundersonTable::gGSMSCAngularDistributions2;
//
std::vector<double> G4GoudsmitSaundersonTable::gMoliereBc;
std::vector<double> G4GoudsmitSaundersonTable::gMoliereXc2;


G4GoudsmitSaundersonTable::G4GoudsmitSaundersonTable(G4bool iselectron) {
  fIsElectron         = iselectron;
  // set initial values: final values will be set in the Initialize method
  fLogLambda0         = 0.;              // will be set properly at init.
  fLogDeltaLambda     = 0.;              // will be set properly at init.
  fInvLogDeltaLambda  = 0.;              // will be set properly at init.
  fInvDeltaQ1         = 0.;              // will be set properly at init.
  fDeltaQ2            = 0.;              // will be set properly at init.
  fInvDeltaQ2         = 0.;              // will be set properly at init.
  //
  fLowEnergyLimit     =   0.1*CLHEP::keV; // will be set properly at init.
  fHighEnergyLimit    = 100.0*CLHEP::MeV; // will be set properly at init.
  //
  fIsMottCorrection   = false;            // will be set properly at init.
  fIsPWACorrection    = false;            // will be set properly at init.
  fMottCorrection     = nullptr;
  //
  fNumSPCEbinPerDec   = 3;
}

G4GoudsmitSaundersonTable::~G4GoudsmitSaundersonTable() {
  for (std::size_t i=0; i<gGSMSCAngularDistributions1.size(); ++i) {
    if (gGSMSCAngularDistributions1[i]) {
      delete [] gGSMSCAngularDistributions1[i]->fUValues;
      delete [] gGSMSCAngularDistributions1[i]->fParamA;
      delete [] gGSMSCAngularDistributions1[i]->fParamB;
      delete gGSMSCAngularDistributions1[i];
    }
  }
  gGSMSCAngularDistributions1.clear();
  for (std::size_t i=0; i<gGSMSCAngularDistributions2.size(); ++i) {
    if (gGSMSCAngularDistributions2[i]) {
      delete [] gGSMSCAngularDistributions2[i]->fUValues;
      delete [] gGSMSCAngularDistributions2[i]->fParamA;
      delete [] gGSMSCAngularDistributions2[i]->fParamB;
      delete gGSMSCAngularDistributions2[i];
    }
  }
  gGSMSCAngularDistributions2.clear();
  if (fMottCorrection) {
    delete fMottCorrection;
    fMottCorrection = nullptr;
  }
  // clear scp correction data
  for (std::size_t imc=0; imc<fSCPCPerMatCuts.size(); ++imc) {
    if (fSCPCPerMatCuts[imc]) {
      fSCPCPerMatCuts[imc]->fVSCPC.clear();
      delete fSCPCPerMatCuts[imc];
    }
  }
  fSCPCPerMatCuts.clear();
  //
  gIsInitialised = false;
}

void G4GoudsmitSaundersonTable::Initialise(G4double lownergylimit, G4double highenergylimit) {
  fLowEnergyLimit     = lownergylimit;
  fHighEnergyLimit    = highenergylimit;
  G4double lLambdaMin = G4Log(gLAMBMIN);
  G4double lLambdaMax = G4Log(gLAMBMAX);
  fLogLambda0         = lLambdaMin;
  fLogDeltaLambda     = (lLambdaMax-lLambdaMin)/(gLAMBNUM-1.);
  fInvLogDeltaLambda  = 1./fLogDeltaLambda;
  fInvDeltaQ1         = 1./((gQMAX1-gQMIN1)/(gQNUM1-1.));
  fDeltaQ2            = (gQMAX2-gQMIN2)/(gQNUM2-1.);
  fInvDeltaQ2         = 1./fDeltaQ2;
  // load precomputed angular distributions and set up several values used during the sampling
  // these are particle independet => they go to static container: load them only onece
  if (!gIsInitialised) {
    // load pre-computed GS angular distributions (computed based on Screened-Rutherford DCS)
    LoadMSCData();
    gIsInitialised = true;
  }
  InitMoliereMSCParams();
  // Mott-correction: particle(e- or e+) dependet so init them
  if (fIsMottCorrection) {
    if (!fMottCorrection) {
      fMottCorrection = new G4GSMottCorrection(fIsElectron);
    }
    fMottCorrection->Initialise();
  }
  // init scattering power correction data; used only together with Mott-correction
  // (Moliere's parameters must be initialised before)
  if (fMottCorrection) {
    InitSCPCorrection();
  }
}


// samplig multiple scattering angles cos(theta) and sin(thata)
//  - including no-scattering, single, "few" scattering cases as well
//  - Mott-correction will be included if it was requested by the user (i.e. if fIsMottCorrection=true)
// lambdaval : s/lambda_el
// qval      : s/lambda_el G1
// scra      : screening parameter
// cost      : will be the smapled cos(theta)
// sint      : will be the smapled sin(theta)
// lekin     : logarithm of the current kinetic energy
// beta2     : the corresponding beta square
// matindx   : index of the current material
// returns true if it was msc
G4bool G4GoudsmitSaundersonTable::Sampling(G4double lambdaval, G4double qval, G4double scra, G4double &cost,
                                           G4double &sint, G4double lekin, G4double beta2, G4int matindx,
                                           GSMSCAngularDtr **gsDtr, G4int &mcekini, G4int &mcdelti,
                                           G4double &transfPar, G4bool isfirst) {
  G4double rand0 = G4UniformRand();
  G4double expn  = G4Exp(-lambdaval);
  //
  // no scattering case
  if (rand0<expn) {
    cost = 1.0;
    sint = 0.0;
    return false;
  }
  //
  // single scattering case : sample from the single scattering PDF
  // - Mott-correction will be included if it was requested by the user (i.e. if fIsMottCorrection=true)
  if (rand0<(1.+lambdaval)*expn) {
    // cost is sampled in SingleScattering()
    cost = SingleScattering(lambdaval, scra, lekin, beta2, matindx);
    // add protections
    if (cost<-1.0) cost = -1.0;
    if (cost>1.0)  cost =  1.0;
    // compute sin(theta) from the sampled cos(theta)
    G4double dum0 = 1.-cost;
    sint = std::sqrt(dum0*(2.0-dum0));
    return false;
  }
  //
  // handle this case:
  //      -lambdaval < 1 i.e. mean #elastic events along the step is < 1 but
  //       the currently sampled case is not 0 or 1 scattering. [Our minimal
  //       lambdaval (that we have precomputed, transformed angular distributions
  //       stored in a form of equally probabe intervalls together with rational
  //       interp. parameters) is 1.]
  //      -probability of having n elastic events follows Poisson stat. with
  //       lambdaval parameter.
  //      -the max. probability (when lambdaval=1) of having more than one
  //       elastic events is 0.2642411 and the prob of having 2,3,..,n elastic
  //       events decays rapidly with n. So set a max n to 10.
  //      -sampling of this cases is done in a one-by-one single elastic event way
  //       where the current #elastic event is sampled from the Poisson distr.
  if (lambdaval<1.0) {
    G4double prob, cumprob;
    prob = cumprob = expn;
    G4double curcost,cursint;
    // init cos(theta) and sin(theta) to the zero scattering values
    cost = 1.0;
    sint = 0.0;
    for (G4int iel=1; iel<10; ++iel) {
      // prob of having iel scattering from Poisson
      prob    *= lambdaval/(G4double)iel;
      cumprob += prob;
      //
      //sample cos(theta) from the singe scattering pdf:
      // - Mott-correction will be included if it was requested by the user (i.e. if fIsMottCorrection=true)
      curcost       = SingleScattering(lambdaval, scra, lekin, beta2, matindx);
      G4double dum0 = 1.-curcost;
      cursint       = dum0*(2.0-dum0); // sin^2(theta)
      //
      // if we got current deflection that is not too small
      // then update cos(theta) sin(theta)
      if (cursint>1.0e-20) {
        cursint         = std::sqrt(cursint);
        G4double curphi = CLHEP::twopi*G4UniformRand();
        cost            = cost*curcost-sint*cursint*std::cos(curphi);
        sint            = std::sqrt(std::max(0.0, (1.0-cost)*(1.0+cost)));
      }
      //
      // check if we have done enough scattering i.e. sampling from the Poisson
      if (rand0<cumprob) {
        return false;
      }
    }
    // if reached the max iter i.e. 10
    return false;
  }
  //
  // multiple scattering case with lambdavalue >= 1:
  //   - use the precomputed and transformed Goudsmit-Saunderson angular
  //     distributions to sample cos(theta)
  //   - Mott-correction will be included if it was requested by the user (i.e. if fIsMottCorrection=true)
  cost = SampleCosTheta(lambdaval, qval, scra, lekin, beta2, matindx, gsDtr, mcekini, mcdelti, transfPar, isfirst);
  // add protections
  if (cost<-1.0)  cost = -1.0;
  if (cost> 1.0)  cost =  1.0;
  // compute cos(theta) and sin(theta) from the sampled 1-cos(theta)
  G4double dum0 = 1.0-cost;
  sint = std::sqrt(dum0*(2.0-dum0));
  // return true if it was msc
  return true;
}


G4double G4GoudsmitSaundersonTable::SampleCosTheta(G4double lambdaval, G4double qval, G4double scra,
                                                   G4double lekin, G4double beta2, G4int matindx,
                                                   GSMSCAngularDtr **gsDtr, G4int &mcekini,G4int &mcdelti,
                                                   G4double &transfPar, G4bool isfirst) {
  G4double cost = 1.;
  // determine the base GS angular distribution if it is the first call (when sub-step sampling is used)
  if (isfirst) {
    *gsDtr = GetGSAngularDtr(scra, lambdaval, qval, transfPar);
  }
  // sample cost from the GS angular distribution (computed based on Screened-Rutherford DCS)
  cost = SampleGSSRCosTheta(*gsDtr, transfPar);
  // Mott-correction if it was requested by the user
  if (fIsMottCorrection && *gsDtr) {     // no Mott-correction in case of izotropic theta
    static const G4int nlooplim = 1000;
    G4int    nloop    =  0 ; // rejection loop counter
//    G4int    ekindx   = -1; // evaluate only in the first call
//    G4int    deltindx = -1 ; // evaluate only in the first call
    G4double val      = fMottCorrection->GetMottRejectionValue(lekin, beta2, qval, cost, matindx, mcekini, mcdelti);
    while (G4UniformRand()>val && ++nloop<nlooplim) {
      // sampling cos(theta)
      cost = SampleGSSRCosTheta(*gsDtr, transfPar);
      val  = fMottCorrection->GetMottRejectionValue(lekin, beta2, qval, cost, matindx, mcekini, mcdelti);
    };
  }
  return cost;
}


// returns with cost sampled from the GS angular distribution computed based on Screened-Rutherford DCS
G4double G4GoudsmitSaundersonTable::SampleGSSRCosTheta(const GSMSCAngularDtr *gsDtr, G4double transfpar) {
  // check if isotropic theta (i.e. cost is uniform on [-1:1])
  if (!gsDtr) {
    return 1.-2.0*G4UniformRand();
  }
  //
  // sampling form the selected distribution
  G4double ndatm1 = gsDtr->fNumData-1.;
  G4double delta  = 1.0/ndatm1;
  // determine lower cumulative bin inidex
  G4double rndm   = G4UniformRand();
  G4int indxl     = rndm*ndatm1;
  G4double  aval  = rndm-indxl*delta;
  G4double  dum0  = delta*aval;

  G4double  dum1  = (1.0+gsDtr->fParamA[indxl]+gsDtr->fParamB[indxl])*dum0;
  G4double  dum2  = delta*delta + gsDtr->fParamA[indxl]*dum0 + gsDtr->fParamB[indxl]*aval*aval;
  G4double sample = gsDtr->fUValues[indxl] +  dum1/dum2 *(gsDtr->fUValues[indxl+1]-gsDtr->fUValues[indxl]);
  // transform back u to cos(theta) :
  // this is the sampled cos(theta) = (2.0*para*sample)/(1.0-sample+para)
  return 1.-(2.0*transfpar*sample)/(1.0-sample+transfpar);
}


// determine the GS angular distribution we need to sample from: will set other things as well ...
G4GoudsmitSaundersonTable::GSMSCAngularDtr* G4GoudsmitSaundersonTable::GetGSAngularDtr(G4double scra,
                                                  G4double &lambdaval, G4double &qval, G4double &transfpar) {
  GSMSCAngularDtr *dtr = nullptr;
  G4bool first         = false;
  // isotropic cost above gQMAX2 (i.e. dtr stays nullptr)
  if (qval<gQMAX2) {
    G4int    lamIndx  = -1; // lambda value index
    G4int    qIndx    = -1; // lambda value index
    // init to second grid Q values
    G4int    numQVal  = gQNUM2;
    G4double minQVal  = gQMIN2;
    G4double invDelQ  = fInvDeltaQ2;
    G4double pIndxH   = 0.; // probability of taking higher index
    // check if first or second grid needs to be used
    if (qval<gQMIN2) {  // first grid
      first = true;
      // protect against qval<gQMIN1
      if (qval<gQMIN1) {
        qval   = gQMIN1;
        qIndx  = 0;
        //pIndxH = 0.;
      }
      // set to first grid Q values
      numQVal  = gQNUM1;
      minQVal  = gQMIN1;
      invDelQ  = fInvDeltaQ1;
    }
    // make sure that lambda = s/lambda_el is in [gLAMBMIN,gLAMBMAX)
    // lambda<gLAMBMIN=1 is already handeled before so lambda>= gLAMBMIN for sure
    if (lambdaval>=gLAMBMAX) {
      lambdaval = gLAMBMAX-1.e-8;
      lamIndx   = gLAMBNUM-1;
    }
    G4double lLambda  = G4Log(lambdaval);
    //
    // determine lower lambda (=s/lambda_el) index: linear interp. on log(lambda) scale
    if (lamIndx<0) {
      pIndxH  = (lLambda-fLogLambda0)*fInvLogDeltaLambda;
      lamIndx = (G4int)(pIndxH);        // lower index of the lambda bin
      pIndxH  = pIndxH-lamIndx;       // probability of taking the higher index distribution
      if (G4UniformRand()<pIndxH) {
        ++lamIndx;
      }
    }
    //
    // determine lower Q (=s/lambda_el G1) index: linear interp. on Q
    if (qIndx<0) {
      pIndxH  = (qval-minQVal)*invDelQ;
      qIndx   = (G4int)(pIndxH);        // lower index of the Q bin
      pIndxH  = pIndxH-qIndx;
      if (G4UniformRand()<pIndxH) {
        ++qIndx;
      }
    }
    // set indx
    G4int indx = lamIndx*numQVal+qIndx;
    if (first) {
      dtr = gGSMSCAngularDistributions1[indx];
    } else {
      dtr = gGSMSCAngularDistributions2[indx];
    }
    // dtr might be nullptr that indicates isotropic cot distribution because:
    // - if the selected lamIndx, qIndx correspond to L(=s/lambda_el) and Q(=s/lambda_el G1) such that G1(=Q/L) > 1
    //   G1 should always be < 1 and if G1 is ~1 -> the dtr is isotropic (this can only happen in case of the 2. grid)
    //
    // compute the transformation parameter
    if (lambdaval>10.0) {
      transfpar = 0.5*(-2.77164+lLambda*( 2.94874-lLambda*(0.1535754-lLambda*0.00552888) ));
    } else {
      transfpar = 0.5*(1.347+lLambda*(0.209364-lLambda*(0.45525-lLambda*(0.50142-lLambda*0.081234))));
    }
    transfpar *= (lambdaval+4.0)*scra;
  }
  // return with the selected GS angular distribution that we need to sample cost from (if nullptr => isotropic cost)
  return dtr;
}


void G4GoudsmitSaundersonTable::LoadMSCData() {
  const char* path = G4FindDataDir("G4LEDATA");
  if (!path) {
    G4Exception("G4GoudsmitSaundersonTable::LoadMSCData()","em0006",
		FatalException,
		"Environment variable G4LEDATA not defined");
    return;
  }
  //
  gGSMSCAngularDistributions1.resize(gLAMBNUM*gQNUM1,nullptr);
  const G4String str1 = G4String(path) + "/msc_GS/GSGrid_1/gsDistr_";
  for (G4int il=0; il<gLAMBNUM; ++il) {
    G4String fname = str1 + std::to_string(il);
    std::ifstream infile(fname,std::ios::in);
    if (!infile.is_open()) {
      G4String msgc = "Cannot open file: " + fname;
      G4Exception("G4GoudsmitSaundersonTable::LoadMSCData()","em0006",
	 	  FatalException, msgc.c_str());
      return;
    }
    for (G4int iq=0; iq<gQNUM1; ++iq) {
      auto gsd = new GSMSCAngularDtr();
      infile >> gsd->fNumData;
      gsd->fUValues = new G4double[gsd->fNumData]();
      gsd->fParamA  = new G4double[gsd->fNumData]();
      gsd->fParamB  = new G4double[gsd->fNumData]();
      G4double ddummy;
      infile >> ddummy; infile >> ddummy;
      for (G4int i=0; i<gsd->fNumData; ++i) {
        infile >> gsd->fUValues[i];
        infile >> gsd->fParamA[i];
        infile >> gsd->fParamB[i];
      }
      gGSMSCAngularDistributions1[il*gQNUM1+iq] = gsd;
    }
    infile.close();
  }
  //
  // second grid
  gGSMSCAngularDistributions2.resize(gLAMBNUM*gQNUM2,nullptr);
  const G4String str2 = G4String(path) + "/msc_GS/GSGrid_2/gsDistr_";
  for (G4int il=0; il<gLAMBNUM; ++il) {
    G4String fname = str2 + std::to_string(il);
    std::ifstream infile(fname,std::ios::in);
    if (!infile.is_open()) {
      G4String msgc = "Cannot open file: " + fname;
      G4Exception("G4GoudsmitSaundersonTable::LoadMSCData()","em0006",
	 	  FatalException, msgc.c_str());
      return;
    }
    for (G4int iq=0; iq<gQNUM2; ++iq) {
      G4int numData;
      infile >> numData;
      if (numData>1) {
        auto gsd = new GSMSCAngularDtr();
        gsd->fNumData = numData;
        gsd->fUValues = new G4double[gsd->fNumData]();
        gsd->fParamA  = new G4double[gsd->fNumData]();
        gsd->fParamB  = new G4double[gsd->fNumData]();
        double ddummy;
        infile >> ddummy; infile >> ddummy;
        for (G4int i=0; i<gsd->fNumData; ++i) {
          infile >> gsd->fUValues[i];
          infile >> gsd->fParamA[i];
          infile >> gsd->fParamB[i];
        }
        gGSMSCAngularDistributions2[il*gQNUM2+iq] = gsd;
      } else {
        gGSMSCAngularDistributions2[il*gQNUM2+iq] = nullptr;
      }
    }
    infile.close();
  }
}

// samples cost in single scattering based on Screened-Rutherford DCS
// (with Mott-correction if it was requested)
G4double G4GoudsmitSaundersonTable::SingleScattering(G4double /*lambdaval*/, G4double scra,
                                                     G4double lekin, G4double beta2,
                                                     G4int matindx) {
  G4double rand1 = G4UniformRand();
  // sample cost from the Screened-Rutherford DCS
  G4double cost  = 1.-2.0*scra*rand1/(1.0-rand1+scra);
  // Mott-correction if it was requested by the user
  if (fIsMottCorrection) {
    static const G4int nlooplim = 1000;  // rejection loop limit
    G4int    nloop    =  0 ; // loop counter
    G4int    ekindx   = -1 ; // evaluate only in the first call
    G4int    deltindx =  0 ; // single scattering case
    G4double q1       =  0.; // not used when deltindx = 0;
    // computing Mott rejection function value
    G4double val      = fMottCorrection->GetMottRejectionValue(lekin, beta2, q1, cost,
                                                               matindx, ekindx, deltindx);
    while (G4UniformRand()>val && ++nloop<nlooplim) {
      // sampling cos(theta) from the Screened-Rutherford DCS
      rand1 = G4UniformRand();
      cost  = 1.-2.0*scra*rand1/(1.0-rand1+scra);
      // computing Mott rejection function value
      val   = fMottCorrection->GetMottRejectionValue(lekin, beta2, q1, cost, matindx,
                                                     ekindx, deltindx);
    };
  }
  return cost;
}


void  G4GoudsmitSaundersonTable::GetMottCorrectionFactors(G4double logekin, G4double beta2,
                                                          G4int matindx, G4double &mcToScr,
                                                          G4double &mcToQ1, G4double &mcToG2PerG1) {
  if (fIsMottCorrection) {
    fMottCorrection->GetMottCorrectionFactors(logekin, beta2, matindx, mcToScr, mcToQ1, mcToG2PerG1);
  }
}


// compute material dependent Moliere MSC parameters at initialisation
void G4GoudsmitSaundersonTable::InitMoliereMSCParams() {
   const G4double const1   = 7821.6;      // [cm2/g]
   const G4double const2   = 0.1569;      // [cm2 MeV2 / g]
   const G4double finstrc2 = 5.325135453E-5; // fine-structure const. square

   G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
   // get number of materials in the table
   std::size_t numMaterials = theMaterialTable->size();
   // make sure that we have long enough vectors
   if(gMoliereBc.size()<numMaterials) {
     gMoliereBc.resize(numMaterials);
     gMoliereXc2.resize(numMaterials);
   }
   G4double xi   = 1.0;
   G4int    maxZ = 200;
   if (fIsMottCorrection || fIsPWACorrection) {
     // xi   = 1.0;  <= always set to 1 from now on
     maxZ = G4GSMottCorrection::GetMaxZet();
   }
   //
   for (std::size_t imat=0; imat<numMaterials; ++imat) {
     const G4Material*      theMaterial     = (*theMaterialTable)[imat];
     const G4ElementVector* theElemVect     = theMaterial->GetElementVector();
     const G4int            numelems        = (G4int)theMaterial->GetNumberOfElements();
     //
     const G4double*        theNbAtomsPerVolVect  = theMaterial->GetVecNbOfAtomsPerVolume();
     G4double               theTotNbAtomsPerVol   = theMaterial->GetTotNbOfAtomsPerVolume();
     //
     G4double zs = 0.0;
     G4double zx = 0.0;
     G4double ze = 0.0;
     G4double sa = 0.0;
     //
     for(G4int ielem = 0; ielem < numelems; ielem++) {
       G4double zet = (*theElemVect)[ielem]->GetZ();
       if (zet>maxZ) {
         zet = (G4double)maxZ;
       }
       G4double iwa  = (*theElemVect)[ielem]->GetN();
       G4double ipz  = theNbAtomsPerVolVect[ielem]/theTotNbAtomsPerVol;
       G4double dum  = ipz*zet*(zet+xi);
       zs           += dum;
       ze           += dum*(-2.0/3.0)*G4Log(zet);
       zx           += dum*G4Log(1.0+3.34*finstrc2*zet*zet);
       sa           += ipz*iwa;
     }
     G4double density = theMaterial->GetDensity()*CLHEP::cm3/CLHEP::g; // [g/cm3]
     //
     gMoliereBc[theMaterial->GetIndex()]  = const1*density*zs/sa*G4Exp(ze/zs)/G4Exp(zx/zs);  //[1/cm]
     gMoliereXc2[theMaterial->GetIndex()] = const2*density*zs/sa;  // [MeV2/cm]
     // change to Geant4 internal units of 1/length and energ2/length
     gMoliereBc[theMaterial->GetIndex()]  *= 1.0/CLHEP::cm;
     gMoliereXc2[theMaterial->GetIndex()] *= CLHEP::MeV*CLHEP::MeV/CLHEP::cm;
   }
}


// this method is temporary, will be removed/replaced with a more effictien solution after 10.3.ref09
G4double G4GoudsmitSaundersonTable::ComputeScatteringPowerCorrection(const G4MaterialCutsCouple *matcut, G4double ekin) {
  G4int    imc       = matcut->GetIndex();
  G4double corFactor = 1.0;
  if (!(fSCPCPerMatCuts[imc]->fIsUse) || ekin<=fSCPCPerMatCuts[imc]->fPrCut) {
    return corFactor;
  }
  // get the scattering power correction factor
  G4double lekin      = G4Log(ekin);
  G4double remaining  = (lekin-fSCPCPerMatCuts[imc]->fLEmin)*fSCPCPerMatCuts[imc]->fILDel;
  std::size_t lindx   = (std::size_t)remaining;
  remaining          -= lindx;
  std::size_t imax    = fSCPCPerMatCuts[imc]->fVSCPC.size()-1;
  if (lindx>=imax) {
    corFactor = fSCPCPerMatCuts[imc]->fVSCPC[imax];
  } else {
    corFactor = fSCPCPerMatCuts[imc]->fVSCPC[lindx] + remaining*(fSCPCPerMatCuts[imc]->fVSCPC[lindx+1]-fSCPCPerMatCuts[imc]->fVSCPC[lindx]);
  }
  return corFactor;
}


void G4GoudsmitSaundersonTable::InitSCPCorrection() {
  // get the material-cuts table
  G4ProductionCutsTable *thePCTable = G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t numMatCuts            = thePCTable->GetTableSize();
  // clear container if any
  for (std::size_t imc=0; imc<fSCPCPerMatCuts.size(); ++imc) {
    if (fSCPCPerMatCuts[imc]) {
      fSCPCPerMatCuts[imc]->fVSCPC.clear();
      delete fSCPCPerMatCuts[imc];
      fSCPCPerMatCuts[imc] = nullptr;
    }
  }
  //
  // set size of the container and create the corresponding data structures
  fSCPCPerMatCuts.resize(numMatCuts,nullptr);
  // loop over the material-cuts and create scattering power correction data structure for each
  for (G4int imc=0; imc<(G4int)numMatCuts; ++imc) {
    const G4MaterialCutsCouple *matCut =  thePCTable->GetMaterialCutsCouple(imc);
    // get e- production cut in the current material-cuts in energy
    G4double limit;
    G4double ecut;
    if (fIsElectron) {
      ecut  = (*(thePCTable->GetEnergyCutsVector(idxG4ElectronCut)))[matCut->GetIndex()];
      limit = 2.*ecut;
    } else {
      ecut  = (*(thePCTable->GetEnergyCutsVector(idxG4PositronCut)))[matCut->GetIndex()];
      limit = ecut;
    }
    G4double min = std::max(limit,fLowEnergyLimit);
    G4double max = fHighEnergyLimit;
    if (min>=max) {
      fSCPCPerMatCuts[imc] = new SCPCorrection();
      fSCPCPerMatCuts[imc]->fIsUse = false;
      fSCPCPerMatCuts[imc]->fPrCut = min;
      continue;
    }
    G4int numEbins       = fNumSPCEbinPerDec*G4lrint(std::log10(max/min));
    numEbins             = std::max(numEbins,3);
    G4double lmin        = G4Log(min);
    G4double ldel        = G4Log(max/min)/(numEbins-1.0);
    fSCPCPerMatCuts[imc] = new SCPCorrection();
    fSCPCPerMatCuts[imc]->fVSCPC.resize(numEbins,1.0);
    fSCPCPerMatCuts[imc]->fIsUse = true;
    fSCPCPerMatCuts[imc]->fPrCut = min;
    fSCPCPerMatCuts[imc]->fLEmin = lmin;
    fSCPCPerMatCuts[imc]->fILDel = 1./ldel;
    for (G4int ie=0; ie<numEbins; ++ie) {
      G4double ekin    = G4Exp(lmin+ie*ldel);
      G4double scpCorr = 1.0;
      // compute correction factor: I.Kawrakow NIMB 114(1996)307-326 (Eqs(32-37))
      if (ie>0) {
         G4double tau     = ekin/CLHEP::electron_mass_c2;
         G4double tauCut  = ecut/CLHEP::electron_mass_c2;
         // Moliere's screening parameter
         G4int    matindx = (G4int)matCut->GetMaterial()->GetIndex();
         G4double A       = GetMoliereXc2(matindx)/(4.0*tau*(tau+2.)*GetMoliereBc(matindx));
         G4double gr      = (1.+2.*A)*G4Log(1.+1./A)-2.;
         G4double dum0    = (tau+2.)/(tau+1.);
         G4double dum1    = tau+1.;
         G4double gm      = G4Log(0.5*tau/tauCut) + (1.+dum0*dum0)*G4Log(2.*(tau-tauCut+2.)/(tau+4.))
                            - 0.25*(tau+2.)*( tau+2.+2.*(2.*tau+1.)/(dum1*dum1))*
                              G4Log((tau+4.)*(tau-tauCut)/tau/(tau-tauCut+2.))
                            + 0.5*(tau-2*tauCut)*(tau+2.)*(1./(tau-tauCut)-1./(dum1*dum1));
         if (gm<gr) {
           gm = gm/gr;
         } else {
           gm = 1.;
         }
         G4double z0 = matCut->GetMaterial()->GetIonisation()->GetZeffective();
         scpCorr     = 1.-gm*z0/(z0*(z0+1.));
      }
      fSCPCPerMatCuts[imc]->fVSCPC[ie] = scpCorr;
    }
  }
}
