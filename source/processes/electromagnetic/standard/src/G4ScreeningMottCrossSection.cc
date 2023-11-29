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
//      G4ScreeningMottCrossSection.cc
//-------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:    G4ScreeningMottCrossSection
//
// Author:      Cristina Consolandi
//
// Creation date: 20.10.2011
//
// Modifications:
// 27-05-2012 Added Analytic Fitting to the Mott Cross Section by means of G4MottCoefficients class.
//
//
// Class Description:
//	Computation of electron Coulomb Scattering Cross Section.
//	Suitable for high energy electrons and light target materials.
//
//      Reference:
//      M.J. Boschini et al.
//     "Non Ionizing Energy Loss induced by Electrons in the Space Environment"
//      Proc. of the 13th International Conference on Particle Physics and Advanced Technology
//      (13th ICPPAT, Como 3-7/10/2011), World Scientific (Singapore).
//	Available at: http://arxiv.org/abs/1111.4042v4
//
//      1) Mott Differential Cross Section Approximation:
//	   For Target material up to Z=92 (U):
//         As described in http://arxiv.org/abs/1111.4042v4
//         par. 2.1 , eq. (16)-(17)
//         Else (Z>92):
//	   W. A. McKinley and H. Fashbach, Phys. Rev. 74, (1948) 1759.
//      2) Screening coefficient:
//      vomn G. Moliere, Z. Naturforsh A2 (1947), 133-145; A3 (1948), 78.
//      3) Nuclear Form Factor:
//      A.V. Butkevich et al. Nucl. Instr. Meth. A488 (2002), 282-294.
//
// ---------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ScreeningMottCrossSection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Proton.hh"
#include "G4LossTableManager.hh"
#include "G4NucleiProperties.hh"
#include "G4Element.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4ThreeVector.hh"
#include "G4Pow.hh"
#include "G4MottData.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static G4double invlog10 = 1./std::log(10.);

static const G4double angle[DIMMOTT]={
1e-07,1.02329e-07,1.04713e-07,1.07152e-07,1.09648e-07,1.12202e-07,1.14815e-07,1.1749e-07,1.20226e-07,1.23027e-07,1.25893e-07,1.28825e-07,1.31826e-07,1.34896e-07,1.38038e-07,1.41254e-07,1.44544e-07,1.47911e-07,1.51356e-07,1.54882e-07,1.58489e-07,1.62181e-07,1.65959e-07,1.69824e-07,1.7378e-07,1.77828e-07,1.8197e-07,1.86209e-07,1.90546e-07,1.94984e-07,1.99526e-07,2.04174e-07,2.0893e-07,2.13796e-07,2.18776e-07,2.23872e-07,2.29087e-07,2.34423e-07,2.39883e-07,2.45471e-07,2.51189e-07,2.5704e-07,2.63027e-07,2.69153e-07,2.75423e-07,2.81838e-07,2.88403e-07,2.95121e-07,3.01995e-07,3.0903e-07,3.16228e-07,3.23594e-07,3.31131e-07,3.38844e-07,3.46737e-07,3.54813e-07,3.63078e-07,3.71535e-07,3.80189e-07,3.89045e-07,3.98107e-07,4.0738e-07,4.16869e-07,4.2658e-07,4.36516e-07,4.46684e-07,4.57088e-07,4.67735e-07,4.7863e-07,4.89779e-07,5.01187e-07,5.12861e-07,5.24807e-07,5.37032e-07,5.49541e-07,5.62341e-07,5.7544e-07,5.88844e-07,6.0256e-07,6.16595e-07,6.30957e-07,6.45654e-07,6.60693e-07,6.76083e-07,6.91831e-07,7.07946e-07,7.24436e-07,7.4131e-07,7.58578e-07,7.76247e-07,7.94328e-07,8.12831e-07,8.31764e-07,8.51138e-07,8.70964e-07,8.91251e-07,9.12011e-07,9.33254e-07,9.54993e-07,9.77237e-07,1e-06,1.02329e-06,1.04713e-06,1.07152e-06,1.09648e-06,1.12202e-06,1.14815e-06,1.1749e-06,1.20226e-06,1.23027e-06,1.25893e-06,1.28825e-06,1.31826e-06,1.34896e-06,1.38038e-06,1.41254e-06,1.44544e-06,1.47911e-06,1.51356e-06,1.54882e-06,1.58489e-06,1.62181e-06,1.65959e-06,1.69824e-06,1.7378e-06,1.77828e-06,1.8197e-06,1.86209e-06,1.90546e-06,1.94984e-06,1.99526e-06,2.04174e-06,2.0893e-06,2.13796e-06,2.18776e-06,2.23872e-06,2.29087e-06,2.34423e-06,2.39883e-06,2.45471e-06,2.51189e-06,2.5704e-06,2.63027e-06,2.69153e-06,2.75423e-06,2.81838e-06,2.88403e-06,2.95121e-06,3.01995e-06,3.0903e-06,3.16228e-06,3.23594e-06,3.31131e-06,3.38844e-06,3.46737e-06,3.54813e-06,3.63078e-06,3.71535e-06,3.80189e-06,3.89045e-06,3.98107e-06,4.0738e-06,4.16869e-06,4.2658e-06,4.36516e-06,4.46684e-06,4.57088e-06,4.67735e-06,4.7863e-06,4.89779e-06,5.01187e-06,5.12861e-06,5.24807e-06,5.37032e-06,5.49541e-06,5.62341e-06,5.7544e-06,5.88844e-06,6.0256e-06,6.16595e-06,6.30957e-06,6.45654e-06,6.60693e-06,6.76083e-06,6.91831e-06,7.07946e-06,7.24436e-06,7.4131e-06,7.58578e-06,7.76247e-06,7.94328e-06,8.12831e-06,8.31764e-06,8.51138e-06,8.70964e-06,8.91251e-06,9.12011e-06,9.33254e-06,9.54993e-06,9.77237e-06,1e-05,1.02329e-05,1.04713e-05,1.07152e-05,1.09648e-05,1.12202e-05,1.14815e-05,1.1749e-05,1.20226e-05,1.23027e-05,1.25893e-05,1.28825e-05,1.31826e-05,1.34896e-05,1.38038e-05,1.41254e-05,1.44544e-05,1.47911e-05,1.51356e-05,1.54882e-05,1.58489e-05,1.62181e-05,1.65959e-05,1.69824e-05,1.7378e-05,1.77828e-05,1.8197e-05,1.86209e-05,1.90546e-05,1.94984e-05,1.99526e-05,2.04174e-05,2.0893e-05,2.13796e-05,2.18776e-05,2.23872e-05,2.29087e-05,2.34423e-05,2.39883e-05,2.45471e-05,2.51189e-05,2.5704e-05,2.63027e-05,2.69153e-05,2.75423e-05,2.81838e-05,2.88403e-05,2.95121e-05,3.01995e-05,3.0903e-05,3.16228e-05,3.23594e-05,3.31131e-05,3.38844e-05,3.46737e-05,3.54813e-05,3.63078e-05,3.71535e-05,3.80189e-05,3.89045e-05,3.98107e-05,4.0738e-05,4.16869e-05,4.2658e-05,4.36516e-05,4.46684e-05,4.57088e-05,4.67735e-05,4.7863e-05,4.89779e-05,5.01187e-05,5.12861e-05,5.24807e-05,5.37032e-05,5.49541e-05,5.62341e-05,5.7544e-05,5.88844e-05,6.0256e-05,6.16595e-05,6.30957e-05,6.45654e-05,6.60693e-05,6.76083e-05,6.91831e-05,7.07946e-05,7.24436e-05,7.4131e-05,7.58578e-05,7.76247e-05,7.94328e-05,8.12831e-05,8.31764e-05,8.51138e-05,8.70964e-05,8.91251e-05,9.12011e-05,9.33254e-05,9.54993e-05,9.77237e-05,0.0001,0.000102329,0.000104713,0.000107152,0.000109648,0.000112202,0.000114815,0.00011749,0.000120226,0.000123027,0.000125893,0.000128825,0.000131826,0.000134896,0.000138038,0.000141254,0.000144544,0.000147911,0.000151356,0.000154882,0.000158489,0.000162181,0.000165959,0.000169824,0.00017378,0.000177828,0.00018197,0.000186209,0.000190546,0.000194984,0.000199526,0.000204174,0.00020893,0.000213796,0.000218776,0.000223872,0.000229087,0.000234423,0.000239883,0.000245471,0.000251189,0.00025704,0.000263027,0.000269153,0.000275423,0.000281838,0.000288403,0.000295121,0.000301995,0.00030903,0.000316228,0.000323594,0.000331131,0.000338844,0.000346737,0.000354813,0.000363078,0.000371535,0.000380189,0.000389045,0.000398107,0.00040738,0.000416869,0.00042658,0.000436516,0.000446684,0.000457088,0.000467735,0.00047863,0.000489779,0.000501187,0.000512861,0.000524807,0.000537032,0.000549541,0.000562341,0.00057544,0.000588844,0.00060256,0.000616595,0.000630957,0.000645654,0.000660693,0.000676083,0.000691831,0.000707946,0.000724436,0.00074131,0.000758578,0.000776247,0.000794328,0.000812831,0.000831764,0.000851138,0.000870964,0.000891251,0.000912011,0.000933254,0.000954993,0.000977237,0.001,0.00102329,0.00104713,0.00107152,0.00109648,0.00112202,0.00114815,0.0011749,0.00120226,0.00123027,0.00125893,0.00128825,0.00131826,0.00134896,0.00138038,0.00141254,0.00144544,0.00147911,0.00151356,0.00154882,0.00158489,0.00162181,0.00165959,0.00169824,0.0017378,0.00177828,0.0018197,0.00186209,0.00190546,0.00194984,0.00199526,0.00204174,0.0020893,0.00213796,0.00218776,0.00223872,0.00229087,0.00234423,0.00239883,0.00245471,0.00251189,0.0025704,0.00263027,0.00269153,0.00275423,0.00281838,0.00288403,0.00295121,0.00301995,0.0030903,0.00316228,0.00323594,0.00331131,0.00338844,0.00346737,0.00354813,0.00363078,0.00371535,0.00380189,0.00389045,0.00398107,0.0040738,0.00416869,0.0042658,0.00436516,0.00446684,0.00457088,0.00467735,0.0047863,0.00489779,0.00501187,0.00512861,0.00524807,0.00537032,0.00549541,0.00562341,0.0057544,0.00588844,0.0060256,0.00616595,0.00630957,0.00645654,0.00660693,0.00676083,0.00691831,0.00707946,0.00724436,0.0074131,0.00758578,0.00776247,0.00794328,0.00812831,0.00831764,0.00851138,0.00870964,0.00891251,0.00912011,0.00933254,0.00954993,0.00977237,0.01,0.0102329,0.0104713,0.0107152,0.0109648,0.0112202,0.0114815,0.011749,0.0120226,0.0123027,0.0125893,0.0128825,0.0131826,0.0134896,0.0138038,0.0141254,0.0144544,0.0147911,0.0151356,0.0154882,0.0158489,0.0162181,0.0165959,0.0169824,0.017378,0.0177828,0.018197,0.0186209,0.0190546,0.0194984,0.0199526,0.0204174,0.020893,0.0213796,0.0218776,0.0223872,0.0229087,0.0234423,0.0239883,0.0245471,0.0251189,0.025704,0.0263027,0.0269153,0.0275423,0.0281838,0.0288403,0.0295121,0.0301995,0.030903,0.0316228,0.0323594,0.0331131,0.0338844,0.0346737,0.0354813,0.0363078,0.0371535,0.0380189,0.0389045,0.0398107,0.040738,0.0416869,0.042658,0.0436516,0.0446684,0.0457088,0.0467735,0.047863,0.0489779,0.0501187,0.0512861,0.0524807,0.0537032,0.0549541,0.0562341,0.057544,0.0588844,0.060256,0.0616595,0.0630957,0.0645654,0.0660693,0.0676083,0.0691831,0.0707946,0.0724436,0.074131,0.0758578,0.0776247,0.0794328,0.0812831,0.0831764,0.0851138,0.0870964,0.0891251,0.0912011,0.0933254,0.0954993,0.0977237,0.1,0.102329,0.104713,0.107152,0.109648,0.112202,0.114815,0.11749,0.120226,0.123027,0.125893,0.128825,0.131826,0.134896,0.138038,0.141254,0.144544,0.147911,0.151356,0.154882,0.158489,0.162181,0.165959,0.169824,0.17378,0.177828,0.18197,0.186209,0.190546,0.194984,0.199526,0.204174,0.20893,0.213796,0.218776,0.223872,0.229087,0.234423,0.239883,0.245471,0.251189,0.25704,0.263027,0.269153,0.275423,0.281838,0.288403,0.295121,0.301995,0.30903,0.316228,0.323594,0.331131,0.338844,0.346737,0.354813,0.363078,0.371535,0.380189,0.389045,0.398107,0.40738,0.416869,0.42658,0.436516,0.446684,0.457088,0.467735,0.47863,0.489779,0.501187,0.512861,0.524807,0.537032,0.549541,0.562341,0.57544,0.588844,0.60256,0.616595,0.630957,0.645654,0.660693,0.676083,0.691831,0.707946,0.724436,0.74131,0.758578,0.776247,0.794328,0.812831,0.831764,0.851138,0.870964,0.891251,0.912011,0.933254,0.954993,0.977237,1,1.02329,1.04713,1.07152,1.09648,1.12202,1.14815,1.1749,1.20226,1.23027,1.25893,1.28825,1.31826,1.34896,1.38038,1.41254,1.44544,1.47911,1.51356,1.54882,1.58489,1.62181,1.65959,1.69824,1.7378,1.77828,1.8197,1.86209,1.90546,1.94984,1.99526,2.04174,2.0893,2.13796,2.18776,2.23872,2.29087,2.34423,2.39883,2.45471,2.51189,2.5704,2.63027,2.69153,2.75423,2.81838,2.88403,2.95121,3.01995,3.0903};

//using namespace std;

G4ScreeningMottCrossSection::G4ScreeningMottCrossSection():
   cosThetaMin(1.0),
   cosThetaMax(-1.0),
   alpha(fine_structure_const),
   htc2(hbarc_squared),
   e2(electron_mass_c2*classic_electr_radius)
{
  fTotalCross=0;

  fNistManager = G4NistManager::Instance();
  fG4pow = G4Pow::GetInstance();
  particle=nullptr;

  spin = mass = mu_rel = 0.0;
  tkinLab = momLab2 = invbetaLab2 = 0.0;
  tkin = mom2 = invbeta2 = beta = gamma = 0.0;

  targetMass = As = etag = ecut = 0.0;

  targetZ = targetA = 0;

  cosTetMinNuc = cosTetMaxNuc = 0.0;
}

//....Ooooo0ooooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ScreeningMottCrossSection::~G4ScreeningMottCrossSection() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ScreeningMottCrossSection::Initialise(const G4ParticleDefinition* p,
                                             G4double cosThetaLim)
{
  SetupParticle(p);
  tkin = mom2 = 0.0;
  ecut = etag = DBL_MAX;
  particle = p;
  cosThetaMin = cosThetaLim;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ScreeningMottCrossSection::SetupKinematic(G4double ekin, G4int Z)
{
  //...Target
  const G4int iz = std::min(92, Z);
  const G4int ia = G4lrint(fNistManager->GetAtomicMassAmu(iz));

  targetZ = iz;
  targetA = ia;
  targetMass = G4NucleiProperties::GetNuclearMass(ia, iz);

  //G4cout<<"......... targetA "<< targetA <<G4endl;
  //G4cout<<"......... targetMass "<< targetMass/MeV <<G4endl;

  // incident particle lab
  tkinLab = ekin;
  momLab2 = tkinLab*(tkinLab + 2.0*mass);
  invbetaLab2 = 1.0 +  mass*mass/momLab2;

  const G4double etot = tkinLab + mass;
  const G4double ptot = std::sqrt(momLab2);
  const G4double m12  = mass*mass;

  // relativistic reduced mass from publucation
  // A.P. Martynenko, R.N. Faustov, Teoret. mat. Fiz. 64 (1985) 179

  //incident particle & target nucleus
  const G4double Ecm = std::sqrt(m12 + targetMass*targetMass + 2.0*etot*targetMass);
  mu_rel = mass*targetMass/Ecm;
  const G4double momCM= ptot*targetMass/Ecm;
  // relative system
  mom2 = momCM*momCM;
  const G4double x = mu_rel*mu_rel/mom2;
  invbeta2 = 1.0 + x;
  tkin = momCM*std::sqrt(invbeta2) - mu_rel;//Ekin of mu_rel
  const G4double  beta2 = 1./invbeta2;
  beta = std::sqrt(beta2) ;
  const G4double gamma2= invbeta2/x;
  gamma = std::sqrt(gamma2);

  //Thomas-Fermi screening length
  const G4double alpha2 = alpha*alpha;
  const G4double aU = 0.88534*CLHEP::Bohr_radius/fG4pow->Z13(targetZ);
  const G4double twoR2 = aU*aU;
  As = 0.25*htc2*(1.13 + 3.76*targetZ*targetZ*invbeta2*alpha2)/(twoR2*mom2);

  //Integration Angles definition
  cosTetMinNuc = cosThetaMin;
  cosTetMaxNuc = cosThetaMax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ScreeningMottCrossSection::FormFactor2ExpHof(G4double sin2t2)
{
  G4double M=targetMass;
  G4double E=tkinLab;
  G4double Etot=E+mass;
  G4double Tmax=2.*M*E*(E+2.*mass)/(mass*mass+M*M+2.*M*Etot);
  G4double T=Tmax*sin2t2;
  G4double q2=T*(T+2.*M);
  q2/=htc2;//1/cm2
  G4double RN=1.27e-13*G4Exp(fG4pow->logZ(targetA)*0.27)*CLHEP::cm;
  G4double xN= (RN*RN*q2);
  G4double den=(1.+xN/12.);
  G4double FN=1./(den*den);
  G4double form2=(FN*FN);

  return form2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ScreeningMottCrossSection::FormFactor2Gauss(G4double sin2t2)
{
  G4double M=targetMass;
  G4double E=tkinLab;
  G4double Etot=E+mass;
  G4double Tmax=2.*M*E*(E+2.*mass)/(mass*mass+M*M+2.*M*Etot);
  G4double T=Tmax*sin2t2;
  G4double q2=T*(T+2.*M);
  q2/=htc2;//1/cm2
  G4double RN=1.27e-13*G4Exp(fG4pow->logZ(targetA)*0.27)*CLHEP::cm;
  G4double xN= (RN*RN*q2);

  G4double expo=(-xN/6.);
  G4double FN=G4Exp(expo);

  G4double form2=(FN*FN);

  return form2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ScreeningMottCrossSection::FormFactor2UniformHelm(G4double sin2t2)
{
  G4double M=targetMass;
  G4double E=tkinLab;
  G4double Etot=E+mass;
  G4double Tmax=2.*M*E*(E+2.*mass)/(mass*mass+M*M+2.*M*Etot);
  G4double T=Tmax*sin2t2;
  G4double q2=T*(T+2.*M);
  q2=q2/(htc2*0.01);//1/cm2

  G4double q=std::sqrt(q2);
  G4double R0=1.2E-13*fG4pow->Z13(targetA);
  G4double R1=2.0E-13;

  G4double x0=q*R0;
  G4double F0=(3./fG4pow->powN(x0,3))*(std::sin(x0)-x0*std::cos(x0));

  G4double x1=q*R1;
  G4double F1=(3./fG4pow->powN(x1,3))*(std::sin(x1)-x1*std::cos(x1));

  G4double F=F0*F1;

  G4double form2=(F*F);

  return form2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ScreeningMottCrossSection::McFcorrection(G4double sin2t2)
{
  const G4double sint = std::sqrt(sin2t2);
  return 1.-beta*beta*sin2t2 + targetZ*alpha*beta*pi*sint*(1.-sint);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ScreeningMottCrossSection::RatioMottRutherford(G4double angles)
{
  return RatioMottRutherfordCosT(std::sqrt(1. - std::cos(angles)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ScreeningMottCrossSection::RatioMottRutherfordCosT(G4double fcost)
{
  G4double R(0.);
  const G4double shift = 0.7181228;
  G4double beta0 = beta - shift;

  G4double b0 = 1.0;
  G4double b[6];
  for(G4int i=0; i<6; ++i) {
    b[i] = b0;
    b0 *= beta0;
  }

  b0 = 1.0;
  for(G4int j=0; j<=4; ++j) {
    G4double a = 0.0;
    for(G4int k=0; k<=5; ++k){
      a += fMottCoef[targetZ][j][k]*b[k];
    }
    R += a*b0;
    b0 *= fcost;
  }
  return R;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ScreeningMottCrossSection::GetTransitionRandom()
{
  G4double x = G4Log(tkinLab)*invlog10;
  G4double res(0.), x0(1.0);
  for(G4int i=0; i<11; ++i) {
    res += fPRM[targetZ][i]*x0;
    x0 *= x;
  }
  //G4cout << "Z= " << targetZ << " x0= " << x0 << " res= " << res << G4endl;
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4ScreeningMottCrossSection::DifferentialXSection(G4int idx, G4int form)
{
  const G4double ang = angle[idx];
  const G4double z1 = 1.0 - std::cos(ang);
  G4double step;
  if(0 == idx) { 
    step = 0.5*(angle[1] + angle[0]); 
  } else if(DIMMOTT == idx + 1) { 
    step = CLHEP::pi - 0.5*(angle[DIMMOTT-2] + angle[DIMMOTT-1]);
  } else {
    step = 0.5*(angle[idx+1] - angle[idx-1]);  
  }

  G4double F2 = 1.;
  switch (form) {
  case 1:
    F2 = FormFactor2ExpHof(z1*0.5);
    break;
  case 2:
    F2 = FormFactor2Gauss(z1*0.5);
    break;
  case 3:
    F2 = FormFactor2UniformHelm(z1*0.5);
    break;
  }

  const G4double R = RatioMottRutherfordCosT(std::sqrt(z1));

  G4double den  = 2.*As + z1;
  G4double func = 1./(den*den);

  G4double fatt = (G4double)targetZ/(mu_rel*gamma*beta*beta);
  G4double sigma= e2*e2*fatt*fatt*func;
  G4double pi2sintet = CLHEP::twopi*std::sqrt(z1*(2. - z1));

  G4double Xsec = std::max(pi2sintet*F2*R*sigma*step, 0.);
  return Xsec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4ScreeningMottCrossSection::NuclearCrossSection(G4int form, G4int fast)
{
  fTotalCross=0.;
  if(cosTetMaxNuc >= cosTetMinNuc) return fTotalCross;
  if(0 == cross.size()) { cross.resize(DIMMOTT, 0.0); } 

  //G4cout<<"MODEL: "<<fast<<G4endl;
  //************ PRECISE COMPUTATION
  if(fast == 0){

    for(G4int i=0; i<DIMMOTT; ++i){
      G4double y = DifferentialXSection(i, form);
      fTotalCross += y;
      cross[i] = fTotalCross;
      if(fTotalCross*1.e-9 > y) {
	for(G4int j=i+1; j<DIMMOTT; ++j) {
	  cross[j] = fTotalCross;
	}
        break;
      }
    }
    //**************** FAST COMPUTATION
  } else if(fast == 1) {

    G4double p0 = electron_mass_c2*classic_electr_radius;
    G4double coeff  = twopi*p0*p0;

    G4double fac = coeff*targetZ*targetZ*invbeta2/mom2;

    G4double x  = 1.0 - cosTetMinNuc;
    G4double x1 = x + 2*As;

    // scattering with nucleus
    fTotalCross = fac*(cosTetMinNuc - cosTetMaxNuc)/
      (x1*(1.0 - cosTetMaxNuc + 2*As));
  }
  /*
  G4cout << "Energy(MeV): " << tkinLab/CLHEP::MeV
         << " Total Cross(mb): " << fTotalCross
         << " form= " << form << " fast= " << fast << G4endl;
  */
  return fTotalCross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4ScreeningMottCrossSection::GetScatteringAngle(G4int form, G4int fast)
{
  G4double scattangle = 0.0;
  G4double r = G4UniformRand();
  //************ PRECISE COMPUTATION
  if(fast == 0){
    //G4cout << "r= " << r << " tot= " << fTotalCross << G4endl;
    r *= fTotalCross;
    for(G4int i=0; i<DIMMOTT; ++i){
      //G4cout << i << ". r= " << r << " xs= " << cross[i] << G4endl;
      if(r <= cross[i]) {
        scattangle = ComputeAngle(i, r);
	break;  
      }
    }

    //**************** FAST COMPUTATION
  } else if(fast == 1) {

    G4double limit = GetTransitionRandom();
    if(limit > 0.0) {
      G4double Sz = 2*As;
      G4double x = Sz-((Sz*(2+Sz))/(Sz+2*limit))+1.0;
      //G4cout << "limit= " << limit << " Sz= " << Sz << " x= " << x << G4endl; 
      G4double angle_limit = (std::abs(x) < 1.0) ? std::acos(x) : 0.0;
      //G4cout<<"ANGLE LIMIT: "<<angle_limit<<"  x= " << x << G4endl;

      if(r > limit && angle_limit != 0.0) {
	x = Sz-((Sz*(2+Sz))/(Sz+2*r))+1.0;
	scattangle = (x >= 1.0) ? 0.0 : ((x <= -1.0) ? pi : std::acos(x));
	//G4cout<<"FAST: scattangle= "<< scattangle <<G4endl;
      }
    } else {
      //G4cout<<"SLOW: total= "<<fTotalCross<< G4endl;
      r *= fTotalCross;
      G4double tot = 0.0;
      for(G4int i=0; i<DIMMOTT; ++i) {
	G4double xs = DifferentialXSection(i, form);
        tot += xs;
        cross[i] = tot;
	if(r <= tot) {
	  scattangle = ComputeAngle(i, r);
	  break;  
	}
      }
    }
  }
  //************************************************
  //G4cout<<"Energy(MeV): "<<tkinLab/MeV<<" SCATTANGLE: "<<scattangle<<G4endl;
  return scattangle;
}

G4double G4ScreeningMottCrossSection::ComputeAngle(G4int i, G4double& r)
{
  G4double x1, x2, y;
  if(0 == i) {
    x1 = 0.0;
    x2 = 0.5*(angle[0] + angle[1]);
    y = cross[0];
  } else if(i == DIMMOTT-1) {
    x1 = 0.5*(angle[i] + angle[i-1]);
    x2 = CLHEP::pi;
    y  = cross[i] - cross[i-1];
    r -= cross[i-1];
  } else {
    x1 = 0.5*(angle[i] + angle[i-1]);
    x2 = 0.5*(angle[i] + angle[i+1]);
    y  = cross[i] - cross[i-1];
    r -= cross[i-1];
  }
  //G4cout << i << ". r= " << r << " y= " << y 
  //	 << " x1= " << " x2= " << x2 << G4endl;
  return x1 + (x2 - x1)*r/y;
}
