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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eDPWAElasticDCS
//
// Author:        Mihaly Novak
//
// Creation date: 02.07.2020
//
// Modifications:
//
//
// -------------------------------------------------------------------

#include "G4eDPWAElasticDCS.hh"

#include "G4Physics2DVector.hh"

#include "zlib.h"

//
// Global variables:
//
G4bool                G4eDPWAElasticDCS::gIsGridLoaded  = false;
G4String              G4eDPWAElasticDCS::gDataDirectory = "";
// final values of these variables will be set in LoadGrid() called by master
std::size_t           G4eDPWAElasticDCS::gNumEnergies   = 106;
std::size_t           G4eDPWAElasticDCS::gIndxEnergyLim =  35;
std::size_t           G4eDPWAElasticDCS::gNumThetas1    = 247;
std::size_t           G4eDPWAElasticDCS::gNumThetas2    = 128;
G4double              G4eDPWAElasticDCS::gLogMinEkin    = 1.0;
G4double              G4eDPWAElasticDCS::gInvDelLogEkin = 1.0;
// containers for grids: Ekin, mu(t)=0.5[1-cos(t)] and u(mu,A)=(A+1)mu/(mu+A)
std::vector<G4double> G4eDPWAElasticDCS::gTheEnergies(G4eDPWAElasticDCS::gNumEnergies);
std::vector<G4double> G4eDPWAElasticDCS::gTheMus1(G4eDPWAElasticDCS::gNumThetas1);
std::vector<G4double> G4eDPWAElasticDCS::gTheMus2(G4eDPWAElasticDCS::gNumThetas2);
std::vector<G4double> G4eDPWAElasticDCS::gTheU1(G4eDPWAElasticDCS::gNumThetas1);
std::vector<G4double> G4eDPWAElasticDCS::gTheU2(G4eDPWAElasticDCS::gNumThetas2);
// abscissas and weights of an 8 point Gauss-Legendre quadrature
// for numerical integration on [0,1]
const G4double        G4eDPWAElasticDCS::gXGL[] = {
  1.98550718E-02, 1.01666761E-01, 2.37233795E-01, 4.08282679E-01,
  5.91717321E-01, 7.62766205E-01, 8.98333239E-01, 9.80144928E-01
};
const G4double        G4eDPWAElasticDCS::gWGL[] = {
  5.06142681E-02, 1.11190517E-01, 1.56853323E-01, 1.81341892E-01,
  1.81341892E-01, 1.56853323E-01, 1.11190517E-01, 5.06142681E-02
};


// - iselectron   : data for e- (for e+ otherwise)
// - isrestricted : sampling of angular deflection on restricted interavl is
//                  required (i.e. in case of mixed-simulation models)
G4eDPWAElasticDCS::G4eDPWAElasticDCS(G4bool iselectron, G4bool isrestricted)
: fIsRestrictedSamplingRequired(isrestricted), fIsElectron(iselectron) {
  fDCS.resize(gMaxZ+1, nullptr);
  fDCSLow.resize(gMaxZ+1, nullptr);
  fSamplingTables.resize(gMaxZ+1, nullptr);
}


 // DTR
G4eDPWAElasticDCS::~G4eDPWAElasticDCS() {
  for (std::size_t i=0; i<fDCS.size(); ++i) {
    if (fDCS[i]) delete fDCS[i];
  }
  for (std::size_t i=0; i<fDCSLow.size(); ++i) {
    if (fDCSLow[i]) delete fDCSLow[i];
  }
  for (std::size_t i=0; i<fSamplingTables.size(); ++i) {
    if (fSamplingTables[i]) delete fSamplingTables[i];
  }
  // clear scp correction data
  for (std::size_t imc=0; imc<fSCPCPerMatCuts.size(); ++imc) {
    if (fSCPCPerMatCuts[imc]) {
      fSCPCPerMatCuts[imc]->fVSCPC.clear();
      delete fSCPCPerMatCuts[imc];
    }
  }
  fSCPCPerMatCuts.clear();
}


// initialise for a given 'iz' atomic number:
//  - nothing happens if it has already been initialised for that Z.
void G4eDPWAElasticDCS::InitialiseForZ(std::size_t iz) {
  if (!gIsGridLoaded) {
    LoadGrid();
  }
  LoadDCSForZ((G4int)iz);
  BuildSmplingTableForZ((G4int)iz);
}


// loads the kinetic energy and theta grids for the DCS data (first init step)
// should be called only by the master
void G4eDPWAElasticDCS::LoadGrid() {
  G4String fname = FindDirectoryPath() + "grid.dat";
  std::ifstream infile(fname.c_str());
  if (!infile.is_open()) {
    G4String msg  =
         "    Problem while trying to read " + fname + " file.\n"+
         "    G4LEDATA version should be G4EMLOW7.12 or later.\n";
    G4Exception("G4eDPWAElasticDCS::ReadCompressedFile","em0006",
                FatalException,msg.c_str());
    return;
  }
  // read size
  infile >> gNumEnergies;
  infile >> gNumThetas1;
  infile >> gNumThetas2;
  // read the grids
  // - energy in [MeV]
  G4double dum = 0.0;
  gTheEnergies.resize(gNumEnergies);
  for (std::size_t ie=0; ie<gNumEnergies; ++ie) {
    infile >> dum;
    gTheEnergies[ie] = G4Log(dum*CLHEP::MeV);
    if (gTheEnergies[ie]<G4Log(2.0*CLHEP::keV)) gIndxEnergyLim = ie;  // only for e-
  }
  ++gIndxEnergyLim;
  // store/set usefull logarithms of the kinetic energy grid
  gLogMinEkin    = gTheEnergies[0];
  gInvDelLogEkin = (gNumEnergies-1)/(gTheEnergies[gNumEnergies-1]-gTheEnergies[0]);
  // - theta1 in [deg.] (247): we store mu(theta) = 0.5[1-cos(theta)]
  gTheMus1.resize(gNumThetas1);
  gTheU1.resize(gNumThetas1);
  const double theA = 0.01;
  for (std::size_t it=0; it<gNumThetas1; ++it) {
    infile >> dum;
    gTheMus1[it] = 0.5*(1.0-std::cos(dum*CLHEP::degree));
    gTheU1[it]   = (theA+1.0)*gTheMus1[it]/(theA+gTheMus1[it]);
  }
  // - theta2 in [deg.] (128): we store mu(theta) = 0.5[1-cos(theta)]
  gTheMus2.resize(gNumThetas2);
  gTheU2.resize(gNumThetas2);
  for (std::size_t it=0; it<gNumThetas2; ++it) {
    infile >> dum;
    gTheMus2[it] = 0.5*(1.0-std::cos(dum*CLHEP::degree));
    gTheU2[it]   = (theA+1.0)*gTheMus2[it]/(theA+gTheMus2[it]);

  }
  infile.close();
  gIsGridLoaded = true;
}


// load DCS data for a given Z
void G4eDPWAElasticDCS::LoadDCSForZ(G4int iz) {
  // Check if it has already been done:
  if (fDCS[iz]) return;
  // Do it otherwise
  if (fIsElectron) {
    // e-
    // load the high energy part firt:
    // - with gNumThetas2 theta and gNumEnergies-gIndxEnergyLim energy values
    const std::size_t hNumEnergries = gNumEnergies-gIndxEnergyLim;
    auto v2DHigh = new G4Physics2DVector(gNumThetas2, hNumEnergries);
    v2DHigh->SetBicubicInterpolation(true);
    for (std::size_t it=0; it<gNumThetas2; ++it) {
      v2DHigh->PutX(it, gTheMus2[it]);
    }
    for (std::size_t ie=0; ie<hNumEnergries; ++ie) {
      v2DHigh->PutY(ie, gTheEnergies[gIndxEnergyLim+ie]);
    }
    std::ostringstream ossh;
    ossh << FindDirectoryPath() << "dcss/el/dcs_"<< iz<<"_h";
    std::istringstream finh(std::ios::in);
    ReadCompressedFile(ossh.str(), finh);
    G4double dum = 0.0;
    for (std::size_t it=0; it<gNumThetas2; ++it) {
      finh >> dum;
      for (std::size_t ie=0; ie<hNumEnergries; ++ie) {
        finh >> dum;
        v2DHigh->PutValue(it, ie, G4Log(dum*CLHEP::cm2/CLHEP::sr));
      }
    }
    // load the low energy part:
    // - with gNumThetas1 theta and gIndxEnergyLim+1 energy values (the +1 is
    //   for including the firts DCS from the higher part above for being
    //   able to perform interpolation between the high and low energy DCS set)
    auto v2DLow = new G4Physics2DVector(gNumThetas1, gIndxEnergyLim+1);
    v2DLow->SetBicubicInterpolation(true);
    for (std::size_t it=0; it<gNumThetas1; ++it) {
      v2DLow->PutX(it, gTheMus1[it]);
    }
    for (std::size_t ie=0; ie<gIndxEnergyLim+1; ++ie) {
      v2DLow->PutY(ie, gTheEnergies[ie]);
    }
    std::ostringstream ossl;
    ossl << FindDirectoryPath() << "dcss/el/dcs_"<< iz<<"_l";
    std::istringstream finl(std::ios::in);
    ReadCompressedFile(ossl.str(), finl);
    for (std::size_t it=0; it<gNumThetas1; ++it) {
      finl >> dum;
      for (std::size_t ie=0; ie<gIndxEnergyLim; ++ie) {
        finl >> dum;
        v2DLow->PutValue(it, ie, G4Log(dum*CLHEP::cm2/CLHEP::sr));
      }
    }
    // add the +1 part: interpolate the firts DCS from the high energy
    std::size_t ix = 0;
    std::size_t iy = 0;
    for (std::size_t it=0; it<gNumThetas1; ++it) {
      const G4double val = v2DHigh->Value(gTheMus1[it], gTheEnergies[gIndxEnergyLim], ix, iy);
      v2DLow->PutValue(it, gIndxEnergyLim, val);
    }
    // store
    fDCSLow[iz] = v2DLow;
    fDCS[iz]    = v2DHigh;
  } else {
    // e+
    auto v2D= new G4Physics2DVector(gNumThetas2, gNumEnergies);
    v2D->SetBicubicInterpolation(true);
    for (std::size_t it=0; it<gNumThetas2; ++it) {
      v2D->PutX(it, gTheMus2[it]);
    }
    for (std::size_t ie=0; ie<gNumEnergies; ++ie) {
      v2D->PutY(ie, gTheEnergies[ie]);
    }
    std::ostringstream oss;
    oss << FindDirectoryPath() << "dcss/pos/dcs_"<< iz;
    std::istringstream fin(std::ios::in);
    ReadCompressedFile(oss.str(), fin);
    G4double dum = 0.0;
    for (std::size_t it=0; it<gNumThetas2; ++it) {
      fin >> dum;
      for (std::size_t ie=0; ie<gNumEnergies; ++ie) {
        fin >> dum;
        v2D->PutValue(it, ie, G4Log(dum*CLHEP::cm2/CLHEP::sr));
      }
    }
    fDCS[iz]= v2D;
  }
}



// Computes the elastic, first and second cross sections for the given kinetic
// energy and target atom.
// Cross sections are zero ff ekin is below/above the kinetic energy grid
void G4eDPWAElasticDCS::ComputeCSPerAtom(G4int iz, G4double ekin, G4double& elcs,
                                         G4double& tr1cs, G4double& tr2cs,
                                         G4double mumin, G4double mumax) {
  // init all cross section values to zero;
  elcs  = 0.0;
  tr1cs = 0.0;
  tr2cs = 0.0;
  // make sure that mu(theta) = 0.5[1-cos(theta)] limits have appropriate vals
  mumin = std::max(0.0, std::min(1.0, mumin));
  mumax = std::max(0.0, std::min(1.0, mumax));
  if (mumin>=mumax) return;
  // make sure that kin. energy is within the available range (10 eV-100MeV)  
  const G4double lekin = std::max(gTheEnergies[0], std::min(gTheEnergies[gNumEnergies-1], G4Log(ekin)));
  // if the lower, denser in theta, DCS set should be used
  const G4bool isLowerGrid = (fIsElectron && lekin<gTheEnergies[gIndxEnergyLim]);
  const std::vector<G4double>& theMuVector = (isLowerGrid) ? gTheMus1    : gTheMus2;
  const G4Physics2DVector*     the2DDCS    = (isLowerGrid) ? fDCSLow[iz] : fDCS[iz];
  // find lower/upper mu bin of integration:
  // 0.0 <= mumin < 1.0 for sure here
  const std::size_t iMuStart = (mumin == 0.0) ? 0 : std::distance( theMuVector.begin(), std::upper_bound(theMuVector.begin(), theMuVector.end(), mumin) )-1 ;
  // 0.0 < mumax <= 1.0 for sure here
  const std::size_t iMuEnd   = (mumax == 1.0) ? theMuVector.size()-2 : std::distance( theMuVector.begin(), std::upper_bound(theMuVector.begin(), theMuVector.end(), mumax) )-1 ;
  // perform numerical integration of the DCS over the given [mumin, mumax]
  // interval (where mu(theta) = 0.5[1-cos(theta)]) to get the elastic, first
  std::size_t ix = 0;
  std::size_t iy = 0;
  for (std::size_t imu=iMuStart; imu<=iMuEnd; ++imu) {
    G4double elcsPar   = 0.0;
    G4double tr1csPar  = 0.0;
    G4double tr2csPar  = 0.0;
    const G4double low = (imu==iMuStart) ? mumin     : theMuVector[imu];
    const G4double del = (imu==iMuEnd)   ? mumax-low : theMuVector[imu+1]-low;
    ix = imu;
    for (std::size_t igl=0; igl<8; ++igl) {
      const double mu  = low + del*gXGL[igl];
      const double dcs = G4Exp(the2DDCS->Value(mu, lekin, ix, iy));
      elcsPar  += gWGL[igl]*dcs;             // elastic
      tr1csPar += gWGL[igl]*dcs*mu;          // first transport
      tr2csPar += gWGL[igl]*dcs*mu*(1.0-mu); // second transport
    }
    elcs  += del*elcsPar;
    tr1cs += del*tr1csPar;
    tr2cs += del*tr2csPar;
  }
  elcs  *=  2.0*CLHEP::twopi;
  tr1cs *=  4.0*CLHEP::twopi;
  tr2cs *= 12.0*CLHEP::twopi;
}


// data structure to store one sampling table: combined Alias + RatIn
// NOTE: when Alias is used, sampling on a resctricted interval is not possible
//       However, Alias makes possible faster sampling. Alias is used in case
//       of single scattering model while it's not used in case of mixed-model
//       when restricted interval sampling is needed. This is controlled by
//       the fIsRestrictedSamplingRequired flag (false by default).
struct OneSamplingTable {
  OneSamplingTable () = default;
  void SetSize(std::size_t nx, G4bool useAlias)  {
    fN = nx;
    // Alias
    if (useAlias) {
      fW.resize(nx);
      fI.resize(nx);
    }
    // Ratin
    fCum.resize(nx);
    fA.resize(nx);
    fB.resize(nx);
  }

  // members
  std::size_t           fN;            // # data points
  G4double              fScreenParA;   // the screening parameter
  std::vector<G4double> fW;
  std::vector<G4double> fCum;
  std::vector<G4double> fA;
  std::vector<G4double> fB;
  std::vector<G4int>    fI;
};


// loads sampling table for the given Z over the enrgy grid
void G4eDPWAElasticDCS::BuildSmplingTableForZ(G4int iz) {
  // Check if it has already been done:
  if (fSamplingTables[iz]) return;
  // Do it otherwise:
  // allocate space
  auto sTables = new std::vector<OneSamplingTable>(gNumEnergies);
  // read compressed sampling table data
  std::ostringstream oss;
  const G4String fname = fIsElectron ? "stables/el/" : "stables/pos/";
  oss << FindDirectoryPath() << fname << "stable_" << iz;
  std::istringstream fin(std::ios::in);
  ReadCompressedFile(oss.str(), fin);
  std::size_t ndata = 0;
  for (std::size_t ie=0; ie<gNumEnergies; ++ie) {
    OneSamplingTable& aTable = (*sTables)[ie];
    // #data in this table
    fin >> ndata;
    aTable.SetSize(ndata, !fIsRestrictedSamplingRequired);
    // the A screening parameter value used for transformation of mu to u
    fin >> aTable.fScreenParA;
    // load data: Alias(W,I) + RatIn(Cum, A, B)
    if (!fIsRestrictedSamplingRequired)  {
      for (std::size_t id=0; id<ndata; ++id) {
        fin >> aTable.fW[id];
      }
      for (std::size_t id=0; id<ndata; ++id) {
        fin >> aTable.fI[id];
      }
    }
    for (std::size_t id=0; id<ndata; ++id) {
      fin >> aTable.fCum[id];
    }
    for (std::size_t id=0; id<ndata; ++id) {
      fin >> aTable.fA[id];
    }
    for (std::size_t id=0; id<ndata; ++id) {
      fin >> aTable.fB[id];
    }
  }
  fSamplingTables[iz] = sTables;
}





// samples cos(theta) i.e. cosine of the polar angle of scattering in elastic
// interaction (Coulomb scattering) of the projectile (e- or e+ depending on
// fIsElectron) with kinetic energy of exp('lekin'), target atom with atomic
// muber of 'iz'. See the 'SampleCosineThetaRestricted' for obtain samples on
// a restricted inteval.
G4double
G4eDPWAElasticDCS::SampleCosineTheta(std::size_t iz, G4double lekin, G4double r1,
                                     G4double r2, G4double r3) {
  lekin = std::max(gTheEnergies[0], std::min(gTheEnergies[gNumEnergies-1], lekin));
  // determine the discrete ekin sampling table to be used:
  //  - statistical interpolation (i.e. linear) on log energy scale
  const G4double      rem = (lekin-gLogMinEkin)*gInvDelLogEkin;
  const auto            k = (std::size_t)rem;
  const std::size_t iekin = (r1 < rem-k) ? k+1 : k;
  // sample the mu(t)=0.5(1-cos(t))
  const double mu    = SampleMu(iz, iekin, r2, r3);
  return std::max(-1.0, std::min(1.0, 1.0-2.0*mu));
}


// samples cos(theta) i.e. cosine of the polar angle of scattering in elastic
// interaction (Coulomb scattering) of the projectile (e- or e+ depending on
// fIsElectron) with kinetic energy of exp('lekin'), target atom with atomic
// muber of 'iz'.
// The cosine theta will be in the [costMin, costMax] interval where costMin
// corresponds to a maximum allowed polar scattering angle thetaMax while
// costMin corresponds to minimum allowed polar scatterin angle thetaMin.
// See the 'SampleCosineTheta' for obtain samples on the entire [-1,1] range.
G4double
G4eDPWAElasticDCS::SampleCosineThetaRestricted(std::size_t iz, G4double lekin,
                                               G4double r1, G4double r2,
                                               G4double costMax, G4double costMin) {
  // costMin corresponds to mu-max while costMax to mu-min: mu(t)=0.5[1-cos(t)]
  lekin = std::max(gTheEnergies[0], std::min(gTheEnergies[gNumEnergies-1], lekin));
  // determine the discrete ekin sampling table to be used:
  //  - statistical interpolation (i.e. linear) on log energy scale
  const G4double      rem = (lekin-gLogMinEkin)*gInvDelLogEkin;
  const auto            k = (size_t)rem;
  const std::size_t iekin = (r1 < rem-k) ? k : k+1;
  // sample the mu(t)=0.5(1-cos(t))
  const G4double mu  = SampleMu(iz, iekin, r2, 0.5*(1.0-costMax), 0.5*(1.0-costMin));
  return std::max(-1.0, std::min(1.0, 1.0-2.0*mu));
}



G4double
G4eDPWAElasticDCS::SampleMu(std::size_t izet, std::size_t ie, G4double r1, G4double r2) {
  OneSamplingTable& rtn = (*fSamplingTables[izet])[ie];
  // get the lower index of the bin by using the alias part
  const G4double rest = r1 * (rtn.fN - 1);
  auto indxl          = (std::size_t)rest;
  const G4double dum0 = rest - indxl;
  if (rtn.fW[indxl] < dum0) indxl = rtn.fI[indxl];
  // sample value within the selected bin by using ratin based numerical inversion
  const G4double delta = rtn.fCum[indxl + 1] - rtn.fCum[indxl];
  const G4double aval  = r2 * delta;

  const G4double dum1  = (1.0 + rtn.fA[indxl] + rtn.fB[indxl]) * delta * aval;
  const G4double dum2  = delta * delta + rtn.fA[indxl] * delta * aval + rtn.fB[indxl] * aval * aval;
  const std::vector<G4double>& theUVect = (fIsElectron && ie<gIndxEnergyLim) ? gTheU1 : gTheU2;
  const G4double    u  = theUVect[indxl] + dum1 / dum2 * (theUVect[indxl + 1] - theUVect[indxl]);
  // transform back u to mu
  return rtn.fScreenParA*u/(rtn.fScreenParA+1.0-u);
}



G4double
G4eDPWAElasticDCS::FindCumValue(G4double u, const OneSamplingTable& stable,
                                const std::vector<G4double>& uvect) {
  const std::size_t iLow = std::distance( uvect.begin(), std::upper_bound(uvect.begin(), uvect.end(), u) )-1;
  const G4double    tau  = (u-uvect[iLow])/(uvect[iLow+1]-uvect[iLow]); // Note: I could store 1/(fX[iLow+1]-fX[iLow])
  const G4double    dum0 = (1.0+stable.fA[iLow]*(1.0-tau)+stable.fB[iLow]);
  const G4double    dum1 = 2.0*stable.fB[iLow]*tau;
  const G4double    dum2 = 1.0 - std::sqrt(std::max(0.0, 1.0-2.0*dum1*tau/(dum0*dum0)));
  return std::min(stable.fCum[iLow+1], std::max(stable.fCum[iLow], stable.fCum[iLow]+dum0*dum2*(stable.fCum[iLow+1]-stable.fCum[iLow])/dum1 ));
}

// muMin and muMax : no checks on these
G4double G4eDPWAElasticDCS::SampleMu(std::size_t izet, std::size_t ie, G4double r1,
                                     G4double muMin, G4double muMax) {
  const OneSamplingTable& rtn = (*fSamplingTables[izet])[ie];
  const G4double theA    = rtn.fScreenParA;
  //
  const std::vector<G4double>& theUVect = (fIsElectron && ie<gIndxEnergyLim) ? gTheU1 : gTheU2;
  const G4double xiMin   = (muMin > 0.0) ? FindCumValue((theA+1.0)*muMin/(theA+muMin), rtn, theUVect) : 0.0;
  const G4double xiMax   = (muMax < 1.0) ? FindCumValue((theA+1.0)*muMax/(theA+muMax), rtn, theUVect) : 1.0;
  //
  const G4double    xi   = xiMin+r1*(xiMax-xiMin); // a smaple within the range
  const std::size_t iLow = std::distance( rtn.fCum.begin(), std::upper_bound(rtn.fCum.begin(), rtn.fCum.end(), xi) )-1;
  const G4double   delta = rtn.fCum[iLow + 1] - rtn.fCum[iLow];
  const G4double   aval  = xi - rtn.fCum[iLow];

  const G4double dum1    = (1.0 + rtn.fA[iLow] + rtn.fB[iLow]) * delta * aval;
  const G4double dum2    = delta * delta + rtn.fA[iLow] * delta * aval + rtn.fB[iLow] * aval * aval;
  const G4double    u    = theUVect[iLow] + dum1 / dum2 * (theUVect[iLow + 1] - theUVect[iLow]);
  return theA*u/(theA+1.0-u);
}



// set the DCS data directory path
const G4String& G4eDPWAElasticDCS::FindDirectoryPath() {
  // check environment variable
  if (gDataDirectory.empty()) {
    const char* path = G4FindDataDir("G4LEDATA");
    if (path) {
      std::ostringstream ost;
      ost << path << "/dpwa/";
      gDataDirectory = ost.str();
    } else {
      G4Exception("G4eDPWAElasticDCS::FindDirectoryPath()","em0006",
                  FatalException,
                  "Environment variable G4LEDATA not defined");
    }
  }
  return gDataDirectory;
}



// uncompress one data file into the input string stream
void
G4eDPWAElasticDCS::ReadCompressedFile(G4String fname, std::istringstream &iss) {
  G4String *dataString = nullptr;
  G4String compfilename(fname+".z");
  // create input stream with binary mode operation and positioning at the end of the file
  std::ifstream in(compfilename, std::ios::binary | std::ios::ate);
  if (in.good()) {
     // get current position in the stream (was set to the end)
     std::size_t fileSize = in.tellg();
     // set current position being the beginning of the stream
     in.seekg(0,std::ios::beg);
     // create (zlib) byte buffer for the data
     Bytef *compdata = new Bytef[fileSize];
     while(in) {
        in.read((char*)compdata, fileSize);
     }
     // create (zlib) byte buffer for the uncompressed data
     uLongf complen    = (uLongf)(fileSize*4);
     Bytef *uncompdata = new Bytef[complen];
     while (Z_OK!=uncompress(uncompdata, &complen, compdata, fileSize)) {
        // increase uncompressed byte buffer
        delete[] uncompdata;
        complen   *= 2;
        uncompdata = new Bytef[complen];
     }
     // delete the compressed data buffer
     delete [] compdata;
     // create a string from the uncompressed data (will be deallocated by the caller)
     dataString = new G4String((char*)uncompdata, (long)complen);
     // delete the uncompressed data buffer
     delete [] uncompdata;
  } else {
    G4String msg  =
         "    Problem while trying to read " + fname + " data file.\n"+
         "    G4LEDATA version should be G4EMLOW7.12 or later.\n";
    G4Exception("G4eDPWAElasticDCS::ReadCompressedFile","em0006",
                FatalException,msg.c_str());
    return;
  }
  // create the input string stream from the data string
  if (dataString) {
    iss.str(*dataString);
    in.close();
    delete dataString;
  }
}




G4double
G4eDPWAElasticDCS::ComputeScatteringPowerCorrection(const G4MaterialCutsCouple *matcut,
                                                    G4double ekin) {
  const G4int imc  = matcut->GetIndex();
  G4double corFactor = 1.0;
  if (!(fSCPCPerMatCuts[imc]->fIsUse) || ekin<=fSCPCPerMatCuts[imc]->fPrCut) {
    return corFactor;
  }
  // get the scattering power correction factor
  const G4double lekin = G4Log(ekin);
  G4double remaining   = (lekin-fSCPCPerMatCuts[imc]->fLEmin)*fSCPCPerMatCuts[imc]->fILDel;
  std::size_t lindx    = (G4int)remaining;
  remaining           -= lindx;
  std::size_t imax     = fSCPCPerMatCuts[imc]->fVSCPC.size()-1;
  if (lindx>=imax) {
    corFactor = fSCPCPerMatCuts[imc]->fVSCPC[imax];
  } else {
    corFactor = fSCPCPerMatCuts[imc]->fVSCPC[lindx] + remaining*(fSCPCPerMatCuts[imc]->fVSCPC[lindx+1]-fSCPCPerMatCuts[imc]->fVSCPC[lindx]);
  }
  return corFactor;
}


void G4eDPWAElasticDCS::InitSCPCorrection(G4double lowEnergyLimit,
                                          G4double highEnergyLimit) {
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
    const G4Material* mat  = matCut->GetMaterial();
    // get e- production cut in the current material-cuts in energy
    const G4double ecut    = (*(thePCTable->GetEnergyCutsVector(idxG4ElectronCut)))[matCut->GetIndex()];
    const G4double limit   = fIsElectron ? 2.0*ecut : ecut;
    const G4double min     = std::max(limit,lowEnergyLimit);
    const G4double max     = highEnergyLimit;
    if (min>=max) {
      fSCPCPerMatCuts[imc] = new SCPCorrection();
      fSCPCPerMatCuts[imc]->fIsUse = false;
      fSCPCPerMatCuts[imc]->fPrCut = min;
      continue;
    }
    G4int numEbins         = fNumSPCEbinPerDec*G4lrint(std::log10(max/min));
    numEbins               = std::max(numEbins,3);
    const G4double lmin    = G4Log(min);
    const G4double ldel    = G4Log(max/min)/(numEbins-1.0);
    fSCPCPerMatCuts[imc]   = new SCPCorrection();
    fSCPCPerMatCuts[imc]->fVSCPC.resize(numEbins,1.0);
    fSCPCPerMatCuts[imc]->fIsUse = true;
    fSCPCPerMatCuts[imc]->fPrCut = min;
    fSCPCPerMatCuts[imc]->fLEmin = lmin;
    fSCPCPerMatCuts[imc]->fILDel = 1./ldel;
    // compute Moliere material dependet parameetrs
    G4double moliereBc     = 0.0;
    G4double moliereXc2    = 0.0;
    ComputeMParams(mat, moliereBc, moliereXc2);
    // compute scattering power correction over the enrgy grid
    for (G4int ie=0; ie<numEbins; ++ie) {
      const G4double ekin = G4Exp(lmin+ie*ldel);
      G4double scpCorr    = 1.0;
      // compute correction factor: I.Kawrakow, Med.Phys.24,505-517(1997)(Eqs(32-37)
      if (ie>0) {
         const G4double tau     = ekin/CLHEP::electron_mass_c2;
         const G4double tauCut  = ecut/CLHEP::electron_mass_c2;
         // Moliere's screening parameter
         const G4double A       = moliereXc2/(4.0*tau*(tau+2.)*moliereBc);
         const G4double gr      = (1.+2.*A)*G4Log(1.+1./A)-2.;
         const G4double dum0    = (tau+2.)/(tau+1.);
         const G4double dum1    = tau+1.;
         G4double gm  = G4Log(0.5*tau/tauCut) + (1.+dum0*dum0)*G4Log(2.*(tau-tauCut+2.)/(tau+4.))
                        - 0.25*(tau+2.)*( tau+2.+2.*(2.*tau+1.)/(dum1*dum1))*
                          G4Log((tau+4.)*(tau-tauCut)/tau/(tau-tauCut+2.))
                        + 0.5*(tau-2*tauCut)*(tau+2.)*(1./(tau-tauCut)-1./(dum1*dum1));
         if (gm<gr) {
           gm = gm/gr;
         } else {
           gm = 1.;
         }
         const G4double z0 = matCut->GetMaterial()->GetIonisation()->GetZeffective();
         scpCorr = 1.-gm*z0/(z0*(z0+1.));
      }
      fSCPCPerMatCuts[imc]->fVSCPC[ie] = scpCorr;
    }
  }
}


// compute material dependent Moliere MSC parameters at initialisation
void G4eDPWAElasticDCS::ComputeMParams(const G4Material* mat, G4double& theBc,
                                       G4double& theXc2) {
   const G4double const1   = 7821.6;      // [cm2/g]
   const G4double const2   = 0.1569;      // [cm2 MeV2 / g]
   const G4double finstrc2 = 5.325135453E-5; // fine-structure const. square
   //   G4double xi   = 1.0;
   const G4ElementVector* theElemVect     = mat->GetElementVector();
   const std::size_t      numelems        = mat->GetNumberOfElements();
   //
   const G4double*  theNbAtomsPerVolVect  = mat->GetVecNbOfAtomsPerVolume();
   G4double         theTotNbAtomsPerVol   = mat->GetTotNbOfAtomsPerVolume();
   //
   G4double zs = 0.0;
   G4double zx = 0.0;
   G4double ze = 0.0;
   G4double sa = 0.0;
   //
   for(std::size_t ielem = 0; ielem < numelems; ++ielem) {
     const G4double zet  = (*theElemVect)[ielem]->GetZ();
     const G4double iwa  = (*theElemVect)[ielem]->GetN();
     const G4double ipz  = theNbAtomsPerVolVect[ielem]/theTotNbAtomsPerVol;
     const G4double dum  = ipz*zet*(zet+1.0);
     zs           += dum;
     ze           += dum*(-2.0/3.0)*G4Log(zet);
     zx           += dum*G4Log(1.0+3.34*finstrc2*zet*zet);
     sa           += ipz*iwa;
   }
   const G4double density = mat->GetDensity()*CLHEP::cm3/CLHEP::g; // [g/cm3]
   //
   theBc  = const1*density*zs/sa*G4Exp(ze/zs)/G4Exp(zx/zs);  //[1/cm]
   theXc2 = const2*density*zs/sa;  // [MeV2/cm]
   // change to Geant4 internal units of 1/length and energ2/length
   theBc  *= 1.0/CLHEP::cm;
   theXc2 *= CLHEP::MeV*CLHEP::MeV/CLHEP::cm;
}

