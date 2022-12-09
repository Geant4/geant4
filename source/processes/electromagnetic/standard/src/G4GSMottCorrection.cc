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
// ----------------------------------------------------------------------------
//
//
// File name:     G4GSMottCorrection
//
// Author:        Mihaly Novak
//
// Creation date: 23.08.2017
//
// Modifications:
// 02.02.2018 M.Novak: fixed initialization of first moment correction.
//
// Class description: see the header file.
//
// -----------------------------------------------------------------------------

#include "G4GSMottCorrection.hh"

#include "G4PhysicalConstants.hh"
#include "zlib.h"
#include "Randomize.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Material.hh"
#include "G4ElementVector.hh"
#include "G4Element.hh"

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>


const std::string G4GSMottCorrection::gElemSymbols[] = {"H","He","Li","Be","B" ,
 "C" ,"N" ,"O" ,"F" ,"Ne","Na","Mg","Al","Si","P" , "S","Cl","Ar","K" ,"Ca","Sc",
 "Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb",
 "Sr","Y" ,"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I" ,
 "Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",
 "Yb","Lu","Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At",
 "Rn","Fr","Ra","Ac","Th","Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf"};

G4GSMottCorrection::G4GSMottCorrection(G4bool iselectron) : fIsElectron(iselectron) {
  // init grids related data member values
  fMaxEkin        = CLHEP::electron_mass_c2*(1./std::sqrt(1.-gMaxBeta2)-1.);
  fLogMinEkin     = G4Log(gMinEkin);
  fInvLogDelEkin  = (gNumEkin-gNumBeta2)/G4Log(gMidEkin/gMinEkin);
  G4double pt2    = gMidEkin*(gMidEkin+2.0*CLHEP::electron_mass_c2);
  fMinBeta2       = pt2/(pt2+CLHEP::electron_mass_c2*CLHEP::electron_mass_c2);
  fInvDelBeta2    = (gNumBeta2-1.)/(gMaxBeta2-fMinBeta2);
  fInvDelDelta    = (gNumDelta-1.)/gMaxDelta;
  fInvDelAngle    = gNumAngle-1.;
}


G4GSMottCorrection::~G4GSMottCorrection() {
  ClearMCDataPerElement();
  ClearMCDataPerMaterial();
}


void G4GSMottCorrection::GetMottCorrectionFactors(G4double logekin, G4double beta2, G4int matindx, G4double &mcToScr,
                                                  G4double &mcToQ1, G4double &mcToG2PerG1) {
  G4int    ekinIndxLow = 0;
  G4double remRfaction = 0.;
  if (beta2>=gMaxBeta2) {
    ekinIndxLow = gNumEkin - 1;
    // remRfaction = -1.
  } else if (beta2>=fMinBeta2) {  // linear interpolation on \beta^2
    remRfaction   = (beta2 - fMinBeta2) * fInvDelBeta2;
    ekinIndxLow   = (G4int)remRfaction;
    remRfaction  -= ekinIndxLow;
    ekinIndxLow  += (gNumEkin - gNumBeta2);
  } else if (logekin>=fLogMinEkin) {
    remRfaction   = (logekin - fLogMinEkin) * fInvLogDelEkin;
    ekinIndxLow   = (G4int)remRfaction;
    remRfaction  -= ekinIndxLow;
  } // the defaults otherwise i.e. use the lowest energy values when ekin is smaller than the minum ekin
  //
  DataPerEkin *perEkinLow  = fMCDataPerMaterial[matindx]->fDataPerEkin[ekinIndxLow];
  mcToScr      = perEkinLow->fMCScreening;
  mcToQ1       = perEkinLow->fMCFirstMoment;
  mcToG2PerG1  = perEkinLow->fMCSecondMoment;
  if (remRfaction>0.) {
    DataPerEkin *perEkinHigh = fMCDataPerMaterial[matindx]->fDataPerEkin[ekinIndxLow+1];
    mcToScr      += remRfaction*(perEkinHigh->fMCScreening    - perEkinLow->fMCScreening);
    mcToQ1       += remRfaction*(perEkinHigh->fMCFirstMoment  - perEkinLow->fMCFirstMoment);
    mcToG2PerG1  += remRfaction*(perEkinHigh->fMCSecondMoment - perEkinLow->fMCSecondMoment);
  }
}


// accept cost if rndm [0,1] < return value
double G4GSMottCorrection::GetMottRejectionValue(G4double logekin, G4double beta2, G4double q1, G4double cost,
                                                 G4int matindx, G4int &ekindx, G4int &deltindx) {
  G4double val   = 1.0;
  G4double delta = q1/(0.5+q1);
  // check if converged to 1 for all angles => accept cost
  if (delta>=gMaxDelta) {
    return val;
  }
  //
  // check if kinetic energy index needs to be determined
  if (ekindx<0) {
    G4int    ekinIndxLow  = 0;
    G4double probIndxHigh = 0.;  // will be the prob. of taking the ekinIndxLow+1 bin
    if (beta2>gMaxBeta2) {
      ekinIndxLow = gNumEkin - 1;
      // probIndxHigh = -1.
    } else if (beta2>=fMinBeta2) {    // linear interpolation on \beta^2
      probIndxHigh  = (beta2 - fMinBeta2) * fInvDelBeta2;
      ekinIndxLow   = (G4int)probIndxHigh;
      probIndxHigh -= ekinIndxLow;
      ekinIndxLow  += (gNumEkin - gNumBeta2);
    } else if (logekin>fLogMinEkin) { // linear interpolation on \ln(E_{kin})
      probIndxHigh  = (logekin - fLogMinEkin) * fInvLogDelEkin;
      ekinIndxLow   = (G4int)probIndxHigh;
      probIndxHigh -= ekinIndxLow;
    } // the defaults otherwise i.e. use the lowest energy values when ekin is smaller than the minum ekin
    //
    // check if need to take the higher ekin index
    if (G4UniformRand()<probIndxHigh) {
      ++ekinIndxLow;
    }
    // set kinetic energy grid index
    ekindx = ekinIndxLow;
  }
  // check if delta value index needs to be determined (note: in case of single scattering deltindx will be set to 0 by
  // by the caller but the ekindx will be -1: kinetic energy index is not known but the delta index is known)
  if (deltindx<0) {
    // note: delta is for sure < gMaxDelta at this point ( and minimum delta value is 0)
    G4double probIndxHigh = delta*fInvDelDelta;  // will be the prob. of taking the deltIndxLow+1 bin
    G4int    deltIndxLow  = (G4int)probIndxHigh;
    probIndxHigh         -= deltIndxLow;
    // check if need to take the higher delta index
    if (G4UniformRand()<probIndxHigh) {
      ++deltIndxLow;
    }
    // set the delta value grid index
    deltindx = deltIndxLow;
  }
  //
  // get the corresponding distribution
  DataPerDelta *perDelta  = fMCDataPerMaterial[matindx]->fDataPerEkin[ekindx]->fDataPerDelta[deltindx];
  //
  // determine lower index of the angular bin
  G4double ang         = std::sqrt(0.5*(1.-cost)); // sin(0.5\theta) in [0,1]
  G4double remRfaction = ang*fInvDelAngle;
  G4int    angIndx     = (G4int)remRfaction;
  remRfaction         -= angIndx;
  if (angIndx<gNumAngle-2) { // normal case: linear interpolation
    val          = remRfaction*(perDelta->fRejFuntion[angIndx+1]-perDelta->fRejFuntion[angIndx]) + perDelta->fRejFuntion[angIndx];
  } else {   // last bin
    G4double dum = ang-1.+1./fInvDelAngle;
    val          = perDelta->fSA + dum*(perDelta->fSB + dum*(perDelta->fSC + dum*perDelta->fSD));
  }
  return val;
}


void G4GSMottCorrection::Initialise() {
  // load Mott-correction data for each elements that belongs to materials that are used in the detector
  InitMCDataPerElement();
  // clrea Mott-correction data per material
  ClearMCDataPerMaterial();
  // initialise Mott-correction data for the materials that are used in the detector
  InitMCDataPerMaterials();
}


void G4GSMottCorrection::InitMCDataPerElement() {
  // do it only once
  if (fMCDataPerElement.size()<gMaxZet+1) {
    fMCDataPerElement.resize(gMaxZet+1,nullptr);
  }
  // loop over all materials, for those that are used check the list of elements and load data from file if the
  // corresponding data has not been loaded yet
  G4ProductionCutsTable *thePCTable = G4ProductionCutsTable::GetProductionCutsTable();
  G4int numMatCuts = (G4int)thePCTable->GetTableSize();
  for (G4int imc=0; imc<numMatCuts; ++imc) {
    const G4MaterialCutsCouple *matCut = thePCTable->GetMaterialCutsCouple(imc);
    if (!matCut->IsUsed()) {
      continue;
    }
    const G4Material      *mat      = matCut->GetMaterial();
    const G4ElementVector *elemVect = mat->GetElementVector();
    //
    std::size_t numElems = elemVect->size();
    for (std::size_t ielem=0; ielem<numElems; ++ielem) {
      const G4Element *elem = (*elemVect)[ielem];
      G4int izet = G4lrint(elem->GetZ());
      if (izet>gMaxZet) {
        izet = gMaxZet;
      }
      if (!fMCDataPerElement[izet]) {
        LoadMCDataElement(elem);
      }
    }
  }
}


void G4GSMottCorrection::InitMCDataPerMaterials() {
  // prepare size of the container
  std::size_t numMaterials = G4Material::GetNumberOfMaterials();
  if (fMCDataPerMaterial.size()!=numMaterials) {
    fMCDataPerMaterial.resize(numMaterials);
  }
  // init. Mott-correction data for the Materials that are used in the geometry
  G4ProductionCutsTable *thePCTable = G4ProductionCutsTable::GetProductionCutsTable();
  G4int numMatCuts = (G4int)thePCTable->GetTableSize();
  for (G4int imc=0; imc<numMatCuts; ++imc) {
    const G4MaterialCutsCouple *matCut = thePCTable->GetMaterialCutsCouple(imc);
    if (!matCut->IsUsed()) {
      continue;
    }
    const G4Material *mat = matCut->GetMaterial();
    if (!fMCDataPerMaterial[mat->GetIndex()]) {
      InitMCDataMaterial(mat);
    }
  }
}


// it's called only if data has not been loaded for this element yet
void G4GSMottCorrection::LoadMCDataElement(const G4Element *elem) {
  // allocate memory
  G4int izet = elem->GetZasInt();
  if (izet>gMaxZet) {
    izet = gMaxZet;
  }
  auto perElem = new DataPerMaterial();
  AllocateDataPerMaterial(perElem);
  fMCDataPerElement[izet]  = perElem;
  //
  // load data from file
  const char* tmppath = G4FindDataDir("G4LEDATA");
  if (!tmppath) {
    G4Exception("G4GSMottCorrection::LoadMCDataElement()","em0006",
		FatalException,
		"Environment variable G4LEDATA not defined");
    return;
  }
  std::string path(tmppath);
  if (fIsElectron) {
    path += "/msc_GS/MottCor/el/";
  } else {
    path += "/msc_GS/MottCor/pos/";
  }
  std::string fname = path+"rej_"+gElemSymbols[izet-1];
  std::istringstream infile(std::ios::in);
  ReadCompressedFile(fname, infile);
  // check if file is open !!!
  for (G4int iek=0; iek<gNumEkin; ++iek) {
    DataPerEkin *perEkin = perElem->fDataPerEkin[iek];
    // 1. get the 3 Mott-correction factors for the current kinetic energy
    infile >> perEkin->fMCScreening;
    infile >> perEkin->fMCFirstMoment;
    infile >> perEkin->fMCSecondMoment;
    // 2. load each data per delta:
    for (G4int idel=0; idel<gNumDelta; ++idel) {
      DataPerDelta *perDelta = perEkin->fDataPerDelta[idel];
      // 2./a.  : first the rejection function values
      for (G4int iang=0; iang<gNumAngle; ++iang) {
        infile >> perDelta->fRejFuntion[iang];
      }
      // 2./b. : then the 4 spline parameter for the last bin
      infile >> perDelta->fSA;
      infile >> perDelta->fSB;
      infile >> perDelta->fSC;
      infile >> perDelta->fSD;
    }
  }
}

// uncompress one data file into the input string stream
void G4GSMottCorrection::ReadCompressedFile(std::string fname, std::istringstream &iss) {
  std::string *dataString = nullptr;
  std::string compfilename(fname+".z");
  // create input stream with binary mode operation and positioning at the end of the file
  std::ifstream in(compfilename, std::ios::binary | std::ios::ate);
  if (in.good()) {
     // get current position in the stream (was set to the end)
     std::streamoff fileSize = in.tellg();
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
     dataString = new std::string((char*)uncompdata, (long)complen);
     // delete the uncompressed data buffer
     delete [] uncompdata;
  } else {
    std::string msg = "  Problem while trying to read " + compfilename + " data file.\n";
    G4Exception("G4GSMottCorrection::ReadCompressedFile","em0006", FatalException,msg.c_str());
    return;
  }
  // create the input string stream from the data string
  if (dataString) {
    iss.str(*dataString);
    in.close();
    delete dataString;
  }
}


void G4GSMottCorrection::InitMCDataMaterial(const G4Material *mat) {
  constexpr G4double const1   = 7821.6;      // [cm2/g]
  constexpr G4double const2   = 0.1569;      // [cm2 MeV2 / g]
  constexpr G4double finstrc2 = 5.325135453E-5; // fine-structure const. square

  G4double constFactor        = CLHEP::electron_mass_c2*CLHEP::fine_structure_const/0.88534;
  constFactor                *= constFactor;  // (mc^2)^2\alpha^2/( C_{TF}^2)
  // allocate memory
  auto perMat = new DataPerMaterial();
  AllocateDataPerMaterial(perMat);
  fMCDataPerMaterial[mat->GetIndex()] = perMat;
  //
  const G4ElementVector* elemVect           = mat->GetElementVector();
  const G4int            numElems           = (G4int)mat->GetNumberOfElements();
  const G4double*        nbAtomsPerVolVect  = mat->GetVecNbOfAtomsPerVolume();
  G4double               totNbAtomsPerVol   = mat->GetTotNbOfAtomsPerVolume();
  //
  // 1. Compute material dependent part of Moliere's b_c \chi_c^2
  //    (with \xi=1 (i.e. total sub-threshold scattering power correction)
  G4double moliereBc  = 0.0;
  G4double moliereXc2 = 0.0;
  G4double zs         = 0.0;
  G4double ze         = 0.0;
  G4double zx         = 0.0;
  G4double sa         = 0.0;
  G4double xi         = 1.0;
  for (G4int ielem=0; ielem<numElems; ++ielem) {
    G4double zet = (*elemVect)[ielem]->GetZ();
    if (zet>gMaxZet) {
      zet = (G4double)gMaxZet;
    }
    G4double iwa  = (*elemVect)[ielem]->GetN();
    G4double ipz  = nbAtomsPerVolVect[ielem]/totNbAtomsPerVol;
    G4double dum  = ipz*zet*(zet+xi);
    zs           += dum;
    ze           += dum*(-2.0/3.0)*G4Log(zet);
    zx           += dum*G4Log(1.0+3.34*finstrc2*zet*zet);
    sa           += ipz*iwa;
  }
  G4double density = mat->GetDensity()*CLHEP::cm3/CLHEP::g; // [g/cm3]
  //
  moliereBc  = const1*density*zs/sa*G4Exp(ze/zs)/G4Exp(zx/zs);  //[1/cm]
  moliereXc2 = const2*density*zs/sa;  // [MeV2/cm]
  // change to Geant4 internal units of 1/length and energ2/length
  moliereBc  *= 1.0/CLHEP::cm;
  moliereXc2 *= CLHEP::MeV*CLHEP::MeV/CLHEP::cm;
  //
  // 2. loop over the kinetic energy grid
  for (G4int iek=0; iek<gNumEkin; ++iek) {
    // 2./a. set current kinetic energy and pt2 value
      G4double ekin          = G4Exp(fLogMinEkin+iek/fInvLogDelEkin);
      G4double pt2  = ekin*(ekin+2.0*CLHEP::electron_mass_c2);
      if (ekin>gMidEkin) {
        G4double b2   = fMinBeta2+(iek-(gNumEkin-gNumBeta2))/fInvDelBeta2;
        ekin = CLHEP::electron_mass_c2*(1./std::sqrt(1.-b2)-1.);
        pt2  = ekin*(ekin+2.0*CLHEP::electron_mass_c2);
      }
    // 2./b. loop over the elements at the current kinetic energy point
    for (G4int ielem=0; ielem<numElems; ++ielem) {
      const G4Element *elem = (*elemVect)[ielem];
      G4double zet  = elem->GetZ();
      if (zet>gMaxZet) {
        zet = (G4double)gMaxZet;
      }
      G4int izet         = G4lrint(zet);
      // xi should be one i.e. z(z+1) since total sub-threshold scattering power correction
      G4double nZZPlus1  = nbAtomsPerVolVect[ielem]*zet*(zet+1.0)/totNbAtomsPerVol;
      G4double Z23       = std::pow(zet,2./3.);
      //
      DataPerEkin *perElemPerEkin  = fMCDataPerElement[izet]->fDataPerEkin[iek];
      DataPerEkin *perMatPerEkin   = perMat->fDataPerEkin[iek];
      //
      // 2./b./(i) Add the 3 Mott-correction factors
      G4double mcScrCF = perElemPerEkin->fMCScreening;     // \kappa_i[1.13+3.76(\alpha Z_i)^2] with \kappa_i=scr_mc/scr_sr
      // compute the screening parameter correction factor (Z_i contribution to the material)
      // src_{mc} = C \exp\left[ \frac{ \sum_i n_i Z_i(Z_i+1)\ln[Z_{i}^{2/3}\kappa_i(1.13+3.76(\alpha Z_i)^2)] } {\sum_i n_i Z_i(Z_i+1)}
      // with C = \frac{(mc^2)^\alpha^2} {4(pc)^2 C_{TF}^2} = constFactor/(4*(pc)^2)
      // here we compute the \sum_i n_i Z_i(Z_i+1)\ln[Z_{i}^{2/3}\kappa_i(1.13+3.76(\alpha Z_i)^2)] part
      perMatPerEkin->fMCScreening += nZZPlus1*G4Log(Z23*mcScrCF);
      // compute the corrected screening parameter for the current Z_i and E_{kin}
      // src(Z_i)_{mc} = \frac{(mc^2)^\alpha^2 Z_i^{2/3}} {4(pc)^2 C_{TF}^2} \kappa_i[1.13+3.76(\alpha Z_i)^2]
      mcScrCF *= constFactor*Z23/(4.*pt2);
      // compute first moment correction factor
      // q1_{mc} = \frac{ \sum_i n_i Z_i(Z_i+1) A_i  B_i } {\sum_i n_i Z_i(Z_i+1)} \frac{1}{C}
      // where:
      // A_i(src(Z_i)_{mc}) = [\ln(1+1/src(Z_i)_{mc}) - 1/(1+src(Z_i)_{mc})]; where \sigma(Z_i)_{tr1}^(sr) = A_i(src(Z_i)_{mc}) [2\pi r_0 Z_i mc^2/(pc)\beta]^2
      // B_i = \beta_i \gamma_i with beta_i(Z_i) = \sigma(Z_i)_{tr1}^(PWA)/\sigma(Z_i,src(Z_i)_{mc})_{tr1}^(sr)
      // and \gamma_i = \sigma(Z_i)_{el}^(MC-DCS)/\sigma(Z_i,src(Z_i)_{mc})_{el}^(sr)
      // C(src_{mc}) = [\ln(1+1/src_{mc}) - 1/(1+src_{mc})]; where \sigma_{tr1}^(sr) = C(src_{mc}) [2\pi r_0 Z_i mc^2/(pc)\beta]^2
      // A_i x B_i is stored in file per e-/e+, E_{kin} and Z_i
      // here we compute the \sum_i n_i Z_i(Z_i+1) A_i  B_i part
      perMatPerEkin->fMCFirstMoment += nZZPlus1*(G4Log(1.+1./mcScrCF)-1./(1.+mcScrCF))*perElemPerEkin->fMCFirstMoment;
      // compute the second moment correction factor
      // [G2/G1]_{mc} = \frac{ \sum_i n_i Z_i(Z_i+1) A_i } {\sum_i n_i Z_i(Z_i+1)} \frac{1}{C}
      // with A_i(Z_i) = G2(Z_i)^{PWA}/G1(Z_i)^{PWA} and C=G2(Z_i,scr_{mc})^{sr}/G1(Z_i,scr_{mc})^{sr}}
      // here we compute the \sum_i n_i Z_i(Z_i+1) A_i part
      perMatPerEkin->fMCSecondMoment += nZZPlus1*perElemPerEkin->fMCSecondMoment;
      //
      // 2./b./(ii) Go for the rejection funtion part
      // I. loop over delta values
      for (G4int idel=0; idel<gNumDelta; ++idel) {
        DataPerDelta *perMatPerDelta  = perMatPerEkin->fDataPerDelta[idel];
        DataPerDelta *perElemPerDelta = perElemPerEkin->fDataPerDelta[idel];
        // I./a. loop over angles (i.e. the \sin(0.5\theta) values) and add the rejection function
        for (G4int iang=0; iang<gNumAngle; ++iang) {
          perMatPerDelta->fRejFuntion[iang] += nZZPlus1*perElemPerDelta->fRejFuntion[iang];
        }
        // I./b. get the last bin spline parameters and add them (a+bx+cx^2+dx^3)
        perMatPerDelta->fSA += nZZPlus1*perElemPerDelta->fSA;
        perMatPerDelta->fSB += nZZPlus1*perElemPerDelta->fSB;
        perMatPerDelta->fSC += nZZPlus1*perElemPerDelta->fSC;
        perMatPerDelta->fSD += nZZPlus1*perElemPerDelta->fSD;
      }
      //
      // 2./b./(iii) When the last element has been added:
      if (ielem==numElems-1) {
        //
        // 1. the remaining part of the sreening correction and divide the corrected screening par. with Moliere's one:
        //    (Moliere screening parameter = moliereXc2/(4(pc)^2 moliereBc) )
        G4double dumScr   = G4Exp(perMatPerEkin->fMCScreening/zs);
        perMatPerEkin->fMCScreening = constFactor*dumScr*moliereBc/moliereXc2;
        //
        // 2. the remaining part of the first moment correction and divide by the one computed by using the corrected
        //    screening parameter (= (mc^2)^\alpha^2/(4(pc)^2C_{TF}^2) dumScr
        G4double scrCorTed = constFactor*dumScr/(4.*pt2);
        G4double dum0      = G4Log(1.+1./scrCorTed);
        perMatPerEkin->fMCFirstMoment = perMatPerEkin->fMCFirstMoment/(zs*(dum0-1./(1.+scrCorTed)));
        //
        // 3. the remaining part of the second moment correction and divide by the one computed by using the corrected
        //    screening parameter
        G4double G2PerG1   =  3.*(1.+scrCorTed)*((1.+2.*scrCorTed)*dum0-2.)/((1.+scrCorTed)*dum0-1.);
        perMatPerEkin->fMCSecondMoment = perMatPerEkin->fMCSecondMoment/(zs*G2PerG1);
        //
        // 4. scale the maximum of the rejection function to unity and correct the last bin spline parameters as well
        // I. loop over delta values
        for (G4int idel=0; idel<gNumDelta; ++idel) {
          DataPerDelta *perMatPerDelta  = perMatPerEkin->fDataPerDelta[idel];
          G4double maxVal = -1.;
          // II. llop over angles
          for (G4int iang=0; iang<gNumAngle; ++iang) {
            if (perMatPerDelta->fRejFuntion[iang]>maxVal)
              maxVal = perMatPerDelta->fRejFuntion[iang];
          }
          for (G4int iang=0; iang<gNumAngle; ++iang) {
            perMatPerDelta->fRejFuntion[iang] /=maxVal;
          }
          perMatPerDelta->fSA /= maxVal;
          perMatPerDelta->fSB /= maxVal;
          perMatPerDelta->fSC /= maxVal;
          perMatPerDelta->fSD /= maxVal;
        }
      }
    }
  }
}


void G4GSMottCorrection::AllocateDataPerMaterial(DataPerMaterial *data) {
  data->fDataPerEkin = new DataPerEkin*[gNumEkin]();
  for (G4int iek=0; iek<gNumEkin; ++iek) {
    auto perEkin = new DataPerEkin();
    perEkin->fDataPerDelta = new DataPerDelta*[gNumDelta]();
    for (G4int idel=0; idel<gNumDelta; ++idel) {
      auto perDelta                = new DataPerDelta();
      perDelta->fRejFuntion        = new double[gNumAngle]();
      perEkin->fDataPerDelta[idel] = perDelta;
    }
    data->fDataPerEkin[iek] = perEkin;
  }
}

void G4GSMottCorrection::DeAllocateDataPerMaterial(DataPerMaterial *data) {
  for (G4int iek=0; iek<gNumEkin; ++iek) {
    DataPerEkin *perEkin = data->fDataPerEkin[iek]; //new DataPerEkin();
    for (G4int idel=0; idel<gNumDelta; ++idel) {
      DataPerDelta *perDelta = perEkin->fDataPerDelta[idel];
      delete [] perDelta->fRejFuntion;
      delete perDelta;
    }
    delete [] perEkin->fDataPerDelta;
    delete perEkin;
  }
  delete [] data->fDataPerEkin;
}


void G4GSMottCorrection::ClearMCDataPerElement() {
  for (std::size_t i=0; i<fMCDataPerElement.size(); ++i) {
    if (fMCDataPerElement[i]) {
      DeAllocateDataPerMaterial(fMCDataPerElement[i]);
      delete fMCDataPerElement[i];
    }
  }
  fMCDataPerElement.clear();
}

void G4GSMottCorrection::ClearMCDataPerMaterial() {
  for (std::size_t i=0; i<fMCDataPerMaterial.size(); ++i) {
    if (fMCDataPerMaterial[i]) {
      DeAllocateDataPerMaterial(fMCDataPerMaterial[i]);
      delete fMCDataPerMaterial[i];
    }
  }
  fMCDataPerMaterial.clear();
}
