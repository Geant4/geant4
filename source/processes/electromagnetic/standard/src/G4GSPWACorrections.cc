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
// File name:     G4GSPWACorrections
//
// Author:        Mihaly Novak
//
// Creation date: 17.10.2017
//
// Modifications:
// 02.02.2018 M.Novak: fixed initialization of first moment correction.
//
// Class description: see the header file.
//
// -----------------------------------------------------------------------------

#include "G4GSPWACorrections.hh"

#include "G4PhysicalConstants.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Material.hh"
#include "G4ElementVector.hh"
#include "G4Element.hh"


const std::string G4GSPWACorrections::gElemSymbols[] = {"H","He","Li","Be","B" ,
 "C" ,"N" ,"O" ,"F" ,"Ne","Na","Mg","Al","Si","P" , "S","Cl","Ar","K" ,"Ca","Sc",
 "Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb",
 "Sr","Y" ,"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I" ,
 "Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",
 "Yb","Lu","Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At",
 "Rn","Fr","Ra","Ac","Th","Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf"};

G4GSPWACorrections::G4GSPWACorrections(G4bool iselectron) : fIsElectron(iselectron) {
  // init grids related data member values
  fMaxEkin        = CLHEP::electron_mass_c2*(1./std::sqrt(1.-gMaxBeta2)-1.);
  fLogMinEkin     = G4Log(gMinEkin);
  fInvLogDelEkin  = (gNumEkin-gNumBeta2)/G4Log(gMidEkin/gMinEkin);
  G4double pt2    = gMidEkin*(gMidEkin+2.0*CLHEP::electron_mass_c2);
  fMinBeta2       = pt2/(pt2+CLHEP::electron_mass_c2*CLHEP::electron_mass_c2);
  fInvDelBeta2    = (gNumBeta2-1.)/(gMaxBeta2-fMinBeta2);
}


G4GSPWACorrections::~G4GSPWACorrections() {
  ClearDataPerElement();
  ClearDataPerMaterial();
}


void  G4GSPWACorrections::GetPWACorrectionFactors(G4double logekin, G4double beta2, G4int matindx,
                                                  G4double &corToScr, G4double &corToQ1, G4double &corToG2PerG1) {
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
  DataPerMaterial *data = fDataPerMaterial[matindx];
  corToScr      = data->fCorScreening[ekinIndxLow];
  corToQ1       = data->fCorFirstMoment[ekinIndxLow];
  corToG2PerG1  = data->fCorSecondMoment[ekinIndxLow];
  if (remRfaction>0.) {
    corToScr      += remRfaction*(data->fCorScreening[ekinIndxLow+1]    - data->fCorScreening[ekinIndxLow]);
    corToQ1       += remRfaction*(data->fCorFirstMoment[ekinIndxLow+1]  - data->fCorFirstMoment[ekinIndxLow]);
    corToG2PerG1  += remRfaction*(data->fCorSecondMoment[ekinIndxLow+1] - data->fCorSecondMoment[ekinIndxLow]);
  }
}


void  G4GSPWACorrections::Initialise() {
  // load PWA correction data for each elements that belongs to materials that are used in the detector
  InitDataPerElement();
  // clear  PWA correction data per material
  ClearDataPerMaterial();
  // initialise PWA correction data for the materials that are used in the detector
  InitDataPerMaterials();
}


void G4GSPWACorrections::InitDataPerElement() {
  // do it only once
  if (fDataPerElement.size()<gMaxZet+1) {
    fDataPerElement.resize(gMaxZet+1,nullptr);
  }
  // loop over all materials, for those that are used check the list of elements and load data from file if the
  // corresponding data has not been loaded yet
  G4ProductionCutsTable *thePCTable = G4ProductionCutsTable::GetProductionCutsTable();
  G4int numMatCuts = (G4int)thePCTable->GetTableSize();
  for (G4int imc=0; imc<numMatCuts; ++imc) {
    const G4MaterialCutsCouple *matCut =  thePCTable->GetMaterialCutsCouple(imc);
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
      if (!fDataPerElement[izet]) {
        LoadDataElement(elem);
      }
    }
  }
}


void G4GSPWACorrections::InitDataPerMaterials() {
  // prepare size of the container
  std::size_t numMaterials = G4Material::GetNumberOfMaterials();
  if (fDataPerMaterial.size()!=numMaterials) {
    fDataPerMaterial.resize(numMaterials);
  }
  // init. PWA correction data for the Materials that are used in the geometry
  G4ProductionCutsTable *thePCTable = G4ProductionCutsTable::GetProductionCutsTable();
  G4int numMatCuts = (G4int)thePCTable->GetTableSize();
  for (G4int imc=0; imc<numMatCuts; ++imc) {
    const G4MaterialCutsCouple *matCut =  thePCTable->GetMaterialCutsCouple(imc);
    if (!matCut->IsUsed()) {
      continue;
    }
    const G4Material *mat = matCut->GetMaterial();
    if (!fDataPerMaterial[mat->GetIndex()]) {
      InitDataMaterial(mat);
    }
  }
}


// it's called only if data has not been loaded for this element yet
void G4GSPWACorrections::LoadDataElement(const G4Element *elem) {
  // allocate memory
  G4int izet = elem->GetZasInt();
  if (izet>gMaxZet) {
    izet = gMaxZet;
  }
  // load data from file
  const char* tmppath = G4FindDataDir("G4LEDATA");
  if (!tmppath) {
    G4Exception("G4GSPWACorrection::LoadDataElement()","em0006",
		FatalException,
		"Environment variable G4LEDATA not defined");
    return;
  }
  std::string path(tmppath);
  if (fIsElectron) {
    path += "/msc_GS/PWACor/el/";
  } else {
    path += "/msc_GS/PWACor/pos/";
  }
  std::string  fname = path+"cf_"+gElemSymbols[izet-1];
  std::ifstream infile(fname,std::ios::in);
  if (!infile.is_open()) {
    std::string msg = "  Problem while trying to read " + fname + " data file.\n";
    G4Exception("G4GSPWACorrection::LoadDataElement","em0006", FatalException,msg.c_str());
    return;
  }
  // allocate data structure
  auto perElem = new DataPerMaterial();
  perElem->fCorScreening.resize(gNumEkin,0.0);
  perElem->fCorFirstMoment.resize(gNumEkin,0.0);
  perElem->fCorSecondMoment.resize(gNumEkin,0.0);
  fDataPerElement[izet]  = perElem;
  G4double dum0;
  for (G4int iek=0; iek<gNumEkin; ++iek) {
    infile >> dum0;
    infile >> perElem->fCorScreening[iek];
    infile >> perElem->fCorFirstMoment[iek];
    infile >> perElem->fCorSecondMoment[iek];
  }
  infile.close();
}


void G4GSPWACorrections::InitDataMaterial(const G4Material *mat) {
  constexpr G4double const1   = 7821.6;      // [cm2/g]
  constexpr G4double const2   = 0.1569;      // [cm2 MeV2 / g]
  constexpr G4double finstrc2 = 5.325135453E-5; // fine-structure const. square

  G4double constFactor        = CLHEP::electron_mass_c2*CLHEP::fine_structure_const/0.88534;
  constFactor                *= constFactor;  // (mc^2)^2\alpha^2/( C_{TF}^2)
  // allocate memory
  auto perMat = new DataPerMaterial();
  perMat->fCorScreening.resize(gNumEkin,0.0);
  perMat->fCorFirstMoment.resize(gNumEkin,0.0);
  perMat->fCorSecondMoment.resize(gNumEkin,0.0);
  fDataPerMaterial[mat->GetIndex()] = perMat;
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
      G4int izet = G4lrint(zet);
      // loaded PWA corrections for the current element
      DataPerMaterial *perElem  = fDataPerElement[izet];
      //
      // xi should be one i.e. z(z+1) since total sub-threshold scattering power correction
      G4double nZZPlus1  = nbAtomsPerVolVect[ielem]*zet*(zet+1.0)/totNbAtomsPerVol;
      G4double Z23       = std::pow(zet,2./3.);
      //
      // 2./b./(i) Add the 3 PWA correction factors
      G4double mcScrCF = perElem->fCorScreening[iek];     // \kappa_i[1.13+3.76(\alpha Z_i)^2] with \kappa_i=scr_mc/scr_sr
      // compute the screening parameter correction factor (Z_i contribution to the material)
      // src_{mc} = C \exp\left[ \frac{ \sum_i n_i Z_i(Z_i+1)\ln[Z_{i}^{2/3}\kappa_i(1.13+3.76(\alpha Z_i)^2)] } {\sum_i n_i Z_i(Z_i+1)}
      // with C = \frac{(mc^2)^\alpha^2} {4(pc)^2 C_{TF}^2} = constFactor/(4*(pc)^2)
      // here we compute the \sum_i n_i Z_i(Z_i+1)\ln[Z_{i}^{2/3}\kappa_i(1.13+3.76(\alpha Z_i)^2)] part
      perMat->fCorScreening[iek] += nZZPlus1*G4Log(Z23*mcScrCF);
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
      perMat->fCorFirstMoment[iek] += nZZPlus1*(G4Log(1.+1./mcScrCF)-1./(1.+mcScrCF))*perElem->fCorFirstMoment[iek];
      // compute the second moment correction factor
      // [G2/G1]_{mc} = \frac{ \sum_i n_i Z_i(Z_i+1) A_i } {\sum_i n_i Z_i(Z_i+1)} \frac{1}{C}
      // with A_i(Z_i) = G2(Z_i)^{PWA}/G1(Z_i)^{PWA} and C=G2(Z_i,scr_{mc})^{sr}/G1(Z_i,scr_{mc})^{sr}}
      // here we compute the \sum_i n_i Z_i(Z_i+1) A_i part
      perMat->fCorSecondMoment[iek] += nZZPlus1*perElem->fCorSecondMoment[iek];
      //
      // 2./b./(ii) When the last element has been added:
      if (ielem==numElems-1) {
        //
        // 1. the remaining part of the sreening correction and divide the corrected screening par. with Moliere's one:
        //    (Moliere screening parameter = moliereXc2/(4(pc)^2 moliereBc) )
        G4double dumScr   = G4Exp(perMat->fCorScreening[iek]/zs);
        perMat->fCorScreening[iek] = constFactor*dumScr*moliereBc/moliereXc2;
        //
        // 2. the remaining part of the first moment correction and divide by the one computed by using the corrected
        //    screening parameter (= (mc^2)^\alpha^2/(4(pc)^2C_{TF}^2) dumScr
        G4double scrCorTed = constFactor*dumScr/(4.*pt2);
        G4double dum0      = G4Log(1.+1./scrCorTed);
        perMat->fCorFirstMoment[iek] = perMat->fCorFirstMoment[iek]/(zs*(dum0-1./(1.+scrCorTed)));
        //
        // 3. the remaining part of the second moment correction and divide by the one computed by using the corrected
        //    screening parameter
        G4double G2PerG1   =  3.*(1.+scrCorTed)*((1.+2.*scrCorTed)*dum0-2.)/((1.+scrCorTed)*dum0-1.);
        perMat->fCorSecondMoment[iek] = perMat->fCorSecondMoment[iek]/(zs*G2PerG1);
      }
    }
  }
}



void G4GSPWACorrections::ClearDataPerElement() {
  for (std::size_t i=0; i<fDataPerElement.size(); ++i) {
    if (fDataPerElement[i]) {
      fDataPerElement[i]->fCorScreening.clear();
      fDataPerElement[i]->fCorFirstMoment.clear();
      fDataPerElement[i]->fCorSecondMoment.clear();
      delete fDataPerElement[i];
    }
  }
  fDataPerElement.clear();
}


void G4GSPWACorrections::ClearDataPerMaterial() {
  for (std::size_t i=0; i<fDataPerMaterial.size(); ++i) {
    if (fDataPerMaterial[i]) {
      fDataPerMaterial[i]->fCorScreening.clear();
      fDataPerMaterial[i]->fCorFirstMoment.clear();
      fDataPerMaterial[i]->fCorSecondMoment.clear();
      delete fDataPerMaterial[i];
    }
  }
  fDataPerMaterial.clear();
}
