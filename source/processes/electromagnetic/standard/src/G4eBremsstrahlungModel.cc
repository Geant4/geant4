//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eBremsstrahlungModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
// 
// Creation date: 03.01.2002
//
// Modifications: 
//
// Class Description: 
//
// 
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eBremsstrahlungModel.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eBremsstrahlungModel::G4eBremsstrahlungModel(const G4ParticleDefinition* p) 
  : G4VEmModel(),
  particle(0),
  highKinEnergy(100.*TeV),
  lowKinEnergy(1.0*keV),
  minThreshold(1.0*keV),
  probsup(1.0),
  MigdalConstant(classic_electr_radius*electron_Compton_length*electron_Compton_length/pi),
  LPMconstant(fine_structure_const*electron_mass_c2*electron_mass_c2/(8.*pi*hbarc)),
  isElectron(true),
  theLPMflag(true),
  oldMaterial(0)
{
  if(p) SetParticle(p);
  partialSumSigma.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eBremsstrahlungModel::~G4eBremsstrahlungModel() 
{
  size_t n = partialSumSigma.size();
  if(n > 0) {
    for(size_t i=0; i<n; i++) {
      delete partialSumSigma[i];
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremsstrahlungModel::SetParticle(const G4ParticleDefinition* p) 
{
  particle = p;
  if(p == G4Electron::Electron()) isElectron = true;
  else                            isElectron = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungModel::HighEnergyLimit(const G4ParticleDefinition* p,
                                            const G4Material*) 
{
  return highKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4double G4eBremsstrahlungModel::LowEnergyLimit(const G4ParticleDefinition* p,
                                           const G4Material*) 
{
  return lowKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungModel::MinEnergyCut(const G4ParticleDefinition* p,
                                              const G4Material*) 
{
  return minThreshold;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4bool G4eBremsstrahlungModel::IsInCharge(const G4ParticleDefinition* p,
	      		                  const G4Material*) 
{
  return (p == G4Electron::Electron() || p == G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungModel::ComputeDEDX(const G4Material* material,
                                             const G4ParticleDefinition* p,
                                                   G4double kineticEnergy,
                                                   G4double cutEnergy) 
{
  if(!particle) SetParticle(p);
  if(kineticEnergy < lowKinEnergy) return 0.0;

  const G4double thigh = 100.*GeV;
  const G4double xhigh = log(thigh/electron_mass_c2);

  G4double cut = G4std::min(cutEnergy, kineticEnergy);

  G4double x, rate, loss;
  const G4double factorHigh = 36./(1450.*GeV);
  const G4double coef1 = -0.5;
  const G4double coef2 = 2./9.;
                            
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector();

  G4double totalEnergy = kineticEnergy + electron_mass_c2;
  G4double dedx = 0.0;

  //  loop for elements in the material
  for (size_t i=0; i<material->GetNumberOfElements(); i++) {

    G4double Z     = (*theElementVector)[i]->GetZ();
    G4double natom = theAtomicNumDensityVector[i];

    // loss for MinKinEnergy<KineticEnergy<=100 GeV
    if (kineticEnergy <= thigh) {

      x = log(totalEnergy/electron_mass_c2);
      loss = ComputeBremLoss(Z, kineticEnergy, cut, x) ;
      if (!isElectron) loss *= PositronCorrFactorLoss(Z, kineticEnergy, cut);   
                  
    // extrapolation for KineticEnergy>100 GeV
    } else if(cut < thigh) {
                                 
      loss = ComputeBremLoss(Z, thigh, cut, xhigh) ;
      if (!isElectron) loss *= PositronCorrFactorLoss(Z, thigh, cut) ;   
      rate = cut/kineticEnergy;
      loss *= (1. + coef1*rate + coef2*rate*rate);
      rate = cut/thigh;
      loss /= (1.+coef1*rate+coef2*rate*rate);

    } else {
                                 
      loss = ComputeBremLoss(Z, thigh, 0.5*thigh, xhigh) ;
      if (!isElectron) loss *= PositronCorrFactorLoss(Z, thigh, 0.5*thigh) ;   
      rate = cut/kineticEnergy;
      loss *= (1. + coef1*rate + coef2*rate*rate);
      loss *= cut*factorHigh;
    }
    loss *= natom;

    G4double kp2 = MigdalConstant*totalEnergy*totalEnergy
                 * (material->GetElectronDensity()) ;

    // now compute the correction due to the supression(s)
    G4double kmin = 1.*eV;
    G4double kmax = cut;

    if (kmax > kmin) {

      G4double floss = 0.;
      G4int nmax = 100;
      G4int nn;
      G4double vmin=log(kmin);
      G4double vmax=log(kmax) ;
      nn = int(nmax*(vmax-vmin)/(log(highKinEnergy)-vmin)) ;
      G4double u,fac,c,v,dv ;
      dv = (vmax-vmin)/nn ;
      v  = vmin-dv ;
      if(nn > 0) {

        for(G4int n=0; n<=nn; n++) {

          v += dv;  
          u = exp(v);               
          fac = u*SupressionFunction(material,kineticEnergy,u);
	  fac *= probsup*(u*u/(u*u+kp2))+1.-probsup;
          if ((n==0)||(n==nn)) c=0.5;
          else                 c=1. ;
          fac   *= c ;
          floss += fac ;
        }
        floss *=dv/(kmax-kmin); 

      } else {
        floss = 1.;
      }
      if(floss > 1.) floss = 1.;
      // correct the loss
      loss *= floss;
    }
    dedx += loss;
  }
  if(dedx < 0.) dedx = 0.;
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungModel::ComputeBremLoss(G4double Z, G4double T,
                                                 G4double Cut, G4double x)

// compute loss due to soft brems 
{
  static const G4double beta=1.0, ksi=2.0;
  static const G4double clossh = 0.254 , closslow = 1./3. , alosslow = 1. ;
  static const G4double Tlim= 10.*MeV ;

  static const G4double xlim = 1.2 ;
  static const G4int NZ = 8 ;
  static const G4int Nloss = 11 ;
  static const G4double ZZ[NZ] =
        {2.,4.,6.,14.,26.,50.,82.,92.};
  static const G4double coefloss[NZ][Nloss] = {
  // Z=2
 { 0.98916,        0.47564,        -0.2505,       -0.45186,        0.14462,
   0.21307,      -0.013738,      -0.045689,     -0.0042914,      0.0034429,
   0.00064189},

  // Z=4
 { 1.0626,        0.37662,       -0.23646,       -0.45188,        0.14295,
   0.22906,      -0.011041,      -0.051398,     -0.0055123,      0.0039919,
   0.00078003},
  // Z=6
 { 1.0954,          0.315,       -0.24011,       -0.43849,        0.15017,
   0.23001,      -0.012846,      -0.052555,     -0.0055114,      0.0041283,
   0.00080318},

  // Z=14
 { 1.1649,        0.18976,       -0.24972,       -0.30124,         0.1555,
   0.13565,      -0.024765,      -0.027047,    -0.00059821,      0.0019373,
   0.00027647},

  // Z=26
 { 1.2261,        0.14272,       -0.25672,       -0.28407,        0.13874,
   0.13586,      -0.020562,      -0.026722,    -0.00089557,      0.0018665,
   0.00026981},

  // Z=50
 { 1.3147,       0.020049,       -0.35543,       -0.13927,        0.17666,
   0.073746,      -0.036076,      -0.013407,      0.0025727,     0.00084005,
  -1.4082e-05},

  // Z=82
 { 1.3986,       -0.10586,       -0.49187,     -0.0048846,        0.23621,
   0.031652,      -0.052938,     -0.0076639,      0.0048181,     0.00056486,
  -0.00011995},

  // Z=92
 { 1.4217,         -0.116,       -0.55497,      -0.044075,        0.27506,
   0.081364,      -0.058143,      -0.023402,      0.0031322,      0.0020201,
   0.00017519}

    } ;
  static G4double aaa = 0.414;
  static G4double bbb = 0.345;
  static G4double ccc = 0.460;

  G4int iz = 0;
  G4double delz = 1.e6;
  for (G4int ii=0; ii<NZ; ii++)
    {
      G4double dz = abs(Z-ZZ[ii]); 
      if(dz < delz)  { 
        iz = ii; 
        delz = dz;
      }
    }

  G4double xx = log10(T);
  G4double fl = 1.;
  
  if (xx <= xlim)
    {
      fl = coefloss[iz][Nloss-1];
      for (G4int j=Nloss-2; j>=0; j--) fl = fl*xx+coefloss[iz][j];
      if (fl < 0.) fl = 0.;
    }

  G4double loss;
  G4double E = T+electron_mass_c2 ;

  loss = Z*(Z+ksi)*E*E/(T+E)*exp(beta*log(Cut/T))*(2.-clossh*exp(log(Z)/4.));
  if (T <= Tlim) loss /= exp(closslow*log(Tlim/T)); 
  if( T <= Cut)  loss *= exp(alosslow*log(T/Cut));

  //  correction 
  loss *= (aaa+bbb*T/Tlim)/(1.+ccc*T/Tlim);
  loss *= fl;
  loss /= Avogadro; 

  return loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungModel::PositronCorrFactorLoss(G4double Z,
                                 G4double kineticEnergy, G4double cut)

//calculates the correction factor for the energy loss due to bremsstrahlung for positrons
//the same correction is in the (discrete) bremsstrahlung 

{
  static const G4double K = 132.9416*eV ;
  static const G4double a1=4.15e-1, a3=2.10e-3, a5=54.0e-5 ;

  G4double x   = log(kineticEnergy/(K*Z*Z)), x2 = x*x, x3 = x2*x;
  G4double eta = 0.5+atan(a1*x+a3*x3+a5*x3*x2)/pi;
  G4double e0  = cut/kineticEnergy;
  
  G4double factor = 0.0;
  if (e0 < 1.0) { 
    factor=log(1.-e0)/eta; 
    factor=exp(factor);
  }  
  factor = eta*(1.-factor)/e0;

  return factor;
}
      
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungModel::CrossSection(const G4Material* material,
                                              const G4ParticleDefinition* p,
                                                    G4double kineticEnergy,
                                                    G4double cutEnergy,
                                                    G4double maxEnergy) 
{
  if(!particle) SetParticle(p);
  G4double cross = 0.0;
  G4double tmax = G4std::min(maxEnergy, kineticEnergy);
  G4double cut  = G4std::max(cutEnergy, minThreshold);
  if(cut >= tmax) return cross;

  const G4ElementVector* theElementVector = material->GetElementVector() ;
  const G4double* theAtomNumDensityVector = material->GetAtomicNumDensityVector();

  if(material != oldMaterial) {
    oldMaterial = material;
    ComputePartialSumSigma(material, 0.5*highKinEnergy,  
                           G4std::min(cutEnergy, 0.25*highKinEnergy));
  }

  for (size_t i=0; i<material->GetNumberOfElements(); i++) {
             
    cross += theAtomNumDensityVector[i] * CrossSectionPerAtom(kineticEnergy,
             (*theElementVector)[i]->GetZ(), cut);
    if(tmax < kineticEnergy) {
      cross -= theAtomNumDensityVector[i] * CrossSectionPerAtom(kineticEnergy,
             (*theElementVector)[i]->GetZ(), tmax);
    }
  }       
           
  // now compute the correction due to the supression(s)

  G4double kmax = tmax;
  G4double kmin = cut;

  G4double totalEnergy = kineticEnergy+electron_mass_c2 ;
  G4double kp2 = MigdalConstant*totalEnergy*totalEnergy*(material->GetElectronDensity());

  G4double fsig = 0.;
  G4int nmax = 100;
  G4double vmin=log(kmin);
  G4double vmax=log(kmax) ;
  G4int nn = int(nmax*(vmax-vmin)/(log(highKinEnergy)-vmin));
  G4double u,fac,c,v,dv,y ;
  dv = (vmax-vmin)/nn ;
  v  = vmin-dv ;
  if(nn > 0) {

      for(G4int n=0; n<=nn; n++) {

        v += dv;  
        u = exp(v);              
        fac = SupressionFunction(material, kineticEnergy, u);
        y = u/kmax;
        fac *= (4.-4.*y+3.*y*y)/3.;
        fac *= probsup*(u*u/(u*u+kp2))+1.-probsup;
        if ((n==0)||(n==nn)) c=0.5;
        else                 c=1. ;
        fac  *= c;
        fsig += fac;
      }
      y = kmin/kmax ;
      fsig *=dv/(-4.*log(y)/3.-4.*(1.-y)/3.+0.5*(1.-y*y));

  } else {

      fsig = 1.;
  }
  if (fsig > 1.) fsig = 1.;

  // correct the cross section
  cross *= fsig;
  
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungModel::CrossSectionPerAtom(G4double kineticEnergy, 
                                                     G4double Z, G4double cut)
 
// Calculates the cross section per atom in GEANT4 internal units.
//
 
{
  G4double cross = 0.0 ;
  if ( kineticEnergy < 1*keV || kineticEnergy < cut) return cross;

  static const G4double ksi=2.0, alfa=1.00;
  static const G4double csigh = 0.127, csiglow = 0.25, asiglow = 0.020*MeV ;
  static const G4double Tlim = 10.*MeV ;

  static const G4double xlim = 1.2 ;
  static const G4int NZ = 8 ;
  static const G4int Nsig = 11 ;
  static const G4double ZZ[NZ] =
        {2.,4.,6.,14.,26.,50.,82.,92.} ;
  static const G4double coefsig[NZ][Nsig] = {
  // Z=2
  { 0.4638,        0.37748,        0.32249,      -0.060362,      -0.065004,
   -0.033457,      -0.004583,       0.011954,      0.0030404,     -0.0010077,
   -0.00028131},

  // Z=4
  { 0.50008,        0.33483,        0.34364,      -0.086262,      -0.055361,
   -0.028168,     -0.0056172,       0.011129,      0.0027528,    -0.00092265,
   -0.00024348},

  // Z=6
  { 0.51587,        0.31095,        0.34996,       -0.11623,      -0.056167,
   -0.0087154,     0.00053943,      0.0054092,     0.00077685,    -0.00039635,
   -6.7818e-05},

  // Z=14
  { 0.55058,        0.25629,        0.35854,      -0.080656,      -0.054308,
   -0.049933,    -0.00064246,       0.016597,      0.0021789,      -0.001327,
   -0.00025983},

  // Z=26
  { 0.5791,        0.26152,        0.38953,       -0.17104,      -0.099172,
    0.024596,       0.023718,     -0.0039205,     -0.0036658,     0.00041749,
    0.00023408},

  // Z=50
  { 0.62085,        0.27045,        0.39073,       -0.37916,       -0.18878,
    0.23905,       0.095028,      -0.068744,      -0.023809,      0.0062408,
    0.0020407},

  // Z=82
  { 0.66053,        0.24513,        0.35404,       -0.47275,       -0.22837,
    0.35647,        0.13203,        -0.1049,      -0.034851,      0.0095046,
    0.0030535},

  // Z=92
  { 0.67143,        0.23079,        0.32256,       -0.46248,       -0.20013,
    0.3506,        0.11779,        -0.1024,      -0.032013,      0.0092279,
    0.0028592}

    } ;

  G4int iz = 0 ;
  G4double delz = 1.e6 ;
  for (G4int ii=0; ii<NZ; ii++)
  {
    if(abs(Z-ZZ[ii]) < delz)
    {
      iz = ii ;
      delz = abs(Z-ZZ[ii]);
    }
  }

  G4double xx = log10(kineticEnergy) ;
  G4double fs = 1. ;
  
  if (xx <= xlim) {

    fs = coefsig[iz][Nsig-1] ;
    for (G4int j=Nsig-2; j>=0; j--) {

      fs = fs*xx+coefsig[iz][j] ;
    }
    if(fs < 0.) fs = 0.;
  }

  cross = Z*(Z+ksi)*(1.-csigh*exp(log(Z)/4.))*pow(log(kineticEnergy/cut),alfa);

  if (kineticEnergy <= Tlim)
     cross *= exp(csiglow*log(Tlim/kineticEnergy))*(1.+asiglow/(sqrt(Z)*kineticEnergy));

  if (!isElectron)
     cross *= PositronCorrFactorSigma(Z, kineticEnergy, cut);

  cross *= fs/Avogadro ;

  if (cross < 0.) cross = 0.;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4double G4eBremsstrahlungModel::PositronCorrFactorSigma( G4double Z,
                                 G4double kineticEnergy, G4double cut)
 
//Calculates the correction factor for the total cross section of the positron bremsstrahl.
// Eta is the ratio of positron to electron energy loss by bremstrahlung. 
// A parametrized formula from L. Urban is used to estimate eta. It is a fit to the results
// of L. Kim & al: Phys Rev. A33,3002 (1986)
 
{
  static const G4double K = 132.9416*eV;
  static const G4double a1 = 4.15e-1, a3 = 2.10e-3, a5 = 54.0e-5;

  G4double x    = log(kineticEnergy/(K*Z*Z));
  G4double x2 = x*x;
  G4double x3 = x2*x;
  G4double eta  = 0.5 + atan(a1*x + a3*x3 + a5*x3*x2)/pi ;
  G4double alfa = (1. - eta)/eta;
  return eta*pow((1. - cut/kineticEnergy), alfa);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremsstrahlungModel::ComputePartialSumSigma(const G4Material* material,
                                                          G4double kineticEnergy,
                                                          G4double cut)

// Build the table of cross section per element. The table is built for MATERIALS.
// This table is used by DoIt to select randomly an element in the material. 
{
  size_t index = material->GetIndex();
  G4int nElements = material->GetNumberOfElements();
  const G4ElementVector* theElementVector = material->GetElementVector(); 
  const G4double* theAtomNumDensityVector = material->GetAtomicNumDensityVector();

  G4DataVector* dv;

  if (index >= partialSumSigma.size()) {

    dv = new G4DataVector();
    partialSumSigma.push_back(dv);

  } else {

    dv = partialSumSigma[index];
    dv->clear();
  }

  G4double cross = 0.0;

  for (G4int i=0; i<nElements; i++ ) {
             
    cross += theAtomNumDensityVector[i] * CrossSectionPerAtom(kineticEnergy,
              (*theElementVector)[i]->GetZ(), cut);
    dv->push_back(cross);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4std::vector<G4DynamicParticle*>* G4eBremsstrahlungModel::SampleSecondary(
                             const G4Material* material,
                             const G4DynamicParticle* dp,
                                   G4double tmin,
                                   G4double maxEnergy) 
// The emitted gamma energy is sampled using a parametrized formula from L. Urban.
// This parametrization is derived from :
//    cross-section values of Seltzer and Berger for electron energies 1 keV - 10 GeV,
//    screened Bethe Heilter differential cross section above 10 GeV,
//    Migdal corrections in both case. 
//  Seltzer & Berger: Nim B 12:95 (1985)
//  Nelson, Hirayama & Rogers: Technical report 265 SLAC (1985)
//  Migdal: Phys Rev 103:1811 (1956); Messel & Crawford: Pergamon Press (1970)
//     
// A modified version of the random number techniques of Butcher & Messel is used 
//    (Nuc Phys 20(1960),15).
{
  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double tmax = G4std::min(maxEnergy, kineticEnergy);
  if(tmin >= tmax) return 0;
  G4ThreeVector momentum = dp->GetMomentum();

//
// GEANT4 internal units.
// 
  static const G4double
     ah10 = 4.67733E+00, ah11 =-6.19012E-01, ah12 = 2.02225E-02,
     ah20 =-7.34101E+00, ah21 = 1.00462E+00, ah22 =-3.20985E-02,
     ah30 = 2.93119E+00, ah31 =-4.03761E-01, ah32 = 1.25153E-02;

  static const G4double
     bh10 = 4.23071E+00, bh11 =-6.10995E-01, bh12 = 1.95531E-02,
     bh20 =-7.12527E+00, bh21 = 9.69160E-01, bh22 =-2.74255E-02,
     bh30 = 2.69925E+00, bh31 =-3.63283E-01, bh32 = 9.55316E-03;

  static const G4double
     al00 =-2.05398E+00, al01 = 2.38815E-02, al02 = 5.25483E-04,
     al10 =-7.69748E-02, al11 =-6.91499E-02, al12 = 2.22453E-03,
     al20 = 4.06463E-02, al21 =-1.01281E-02, al22 = 3.40919E-04;

  static const G4double
     bl00 = 1.04133E+00, bl01 =-9.43291E-03, bl02 =-4.54758E-04,
     bl10 = 1.19253E-01, bl11 = 4.07467E-02, bl12 =-1.30718E-03,
     bl20 =-1.59391E-02, bl21 = 7.27752E-03, bl22 =-1.94405E-04;

  G4double gammaEnergy;
  G4bool LPMOK = false;

  // select randomly one element constituing the material  
  const G4Element* anElement = SelectRandomAtom(material);

  // Extract Z factors for this Element
  G4double lnZ = 3.*(anElement->GetIonisation()->GetlogZ3());
  G4double FZ = lnZ* (4.- 0.55*lnZ);
  G4double ZZ = anElement->GetIonisation()->GetZZ3();

  // limits of the energy sampling
  G4double totalEnergy = kineticEnergy + electron_mass_c2;
  G4double xmin     = tmin/kineticEnergy;
  G4double xmax     = tmax/kineticEnergy;
  G4double kappa    = log(xmax)/log(xmin);
  G4double epsilmin = tmin/totalEnergy;
  G4double epsilmax = tmax/totalEnergy;

  // Migdal factor
  G4double MigdalFactor = (material->GetElectronDensity())*MigdalConstant
                        / (epsilmax*epsilmax);

  G4double x, epsil, greject, migdal, grejmax, q;
  G4double U  = log(kineticEnergy/electron_mass_c2);
  G4double U2 = U*U;

  //
  //  sample the energy rate of the emitted gamma for electron kinetic energy > 1 MeV
  //

  do {
   if (kineticEnergy > 1.*MeV) 
     {
       // parameters
       G4double ah1 = ah10 + ZZ* (ah11 + ZZ* ah12),
                ah2 = ah20 + ZZ* (ah21 + ZZ* ah22),
                ah3 = ah30 + ZZ* (ah31 + ZZ* ah32);

       G4double bh1 = bh10 + ZZ* (bh11 + ZZ* bh12),
                bh2 = bh20 + ZZ* (bh21 + ZZ* bh22),
                bh3 = bh30 + ZZ* (bh31 + ZZ* bh32);

       G4double ah = 1.   + (ah1*U2 + ah2*U + ah3) / (U2*U);
       G4double bh = 0.75 + (bh1*U2 + bh2*U + bh3) / (U2*U);

       // limit of the screening variable
       G4double screenfac =
       136.*electron_mass_c2/((anElement->GetIonisation()->GetZ3())*totalEnergy);
       G4double screenmin = screenfac*epsilmin/(1.-epsilmin);

       // Compute the maximum of the rejection function
       G4double F1 = G4std::max(ScreenFunction1(screenmin) - FZ ,0.);
       G4double F2 = G4std::max(ScreenFunction2(screenmin) - FZ ,0.);
       grejmax = (F1 - epsilmin* (F1*ah - bh*epsilmin*F2))/(42.392 - FZ);

       // sample the energy rate of the emitted Gamma 
       G4double screenvar;

       do {
             q = G4UniformRand();
             x = pow(xmin, q + kappa*(1.0 - q));  
             epsil = x*kineticEnergy/totalEnergy;
             screenvar = screenfac*epsil/(1-epsil);
             F1 = G4std::max(ScreenFunction1(screenvar) - FZ ,0.);
             F2 = G4std::max(ScreenFunction2(screenvar) - FZ ,0.);
             migdal = (1. + MigdalFactor)/(1. + MigdalFactor/(x*x));
             greject = migdal*(F1 - epsil* (ah*F1 - bh*epsil*F2))/(42.392 - FZ);      

        }  while( greject < G4UniformRand()*grejmax );

    }

   else
     {  
       // sample the energy rate of the emitted gamma for electron kinetic energy < 1 MeV
       //
       // parameters
       G4double al0 = al00 + ZZ* (al01 + ZZ* al02),
                al1 = al10 + ZZ* (al11 + ZZ* al12),
                al2 = al20 + ZZ* (al21 + ZZ* al22);
 
       G4double bl0 = bl00 + ZZ* (bl01 + ZZ* bl02),
                bl1 = bl10 + ZZ* (bl11 + ZZ* bl12),
                bl2 = bl20 + ZZ* (bl21 + ZZ* bl22);
 
       G4double al = al0 + al1*U + al2*U2;
       G4double bl = bl0 + bl1*U + bl2*U2;

       // Compute the maximum of the rejection function
       grejmax = G4std::max(1. + xmin* (al + bl*xmin), 1.+al+bl);
       G4double xm = -al/(2.*bl);
       if ((xmin < xm)&&(xm < 1.)) grejmax = G4std::max(grejmax, 1.+ xm* (al + bl*xm));

       // sample the energy rate of the emitted Gamma 

       do {  
             q = G4UniformRand();
             x = pow(xmin, q + kappa*(1.0 - q));  
             migdal = (1. + MigdalFactor)/(1. + MigdalFactor/(x*x));  
             greject = migdal*(1. + x* (al + bl*x));
        }  while( greject < G4UniformRand()*grejmax );
     }
  
    gammaEnergy = x*kineticEnergy; 

    if(theLPMflag)
      {
     // take into account the supression due to the LPM effect
      if (G4UniformRand() <= SupressionFunction(material,kineticEnergy,gammaEnergy))
        LPMOK = true ;
      }
    else LPMOK = true ;

  } while (!LPMOK) ;


  //protection: DO NOT PRODUCE a gamma with energy 0. !
  if (gammaEnergy <= 0.) return 0; 

  //
  //  angles of the emitted gamma. ( Z - axis along the parent particle)
  //
  //  universal distribution suggested by L. Urban (Geant3 manual (1993) Phys211),
  //  derived from Tsai distribution (Rev Mod Phys 49,421(1977))

  G4double u;
  const G4double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

  if (9./(9.+d) > G4UniformRand()) u = - log(G4UniformRand()*G4UniformRand())/a1 ;
     else                          u = - log(G4UniformRand()*G4UniformRand())/a2 ;

  G4double theta = u*electron_mass_c2/totalEnergy;

  G4double sint = sin(theta);
 
  G4double phi = twopi * G4UniformRand() ; 

  G4ThreeVector gammaDirection(sint*cos(phi),sint*sin(phi), cos(theta));
  gammaDirection.rotateUz(momentum);

  // create G4DynamicParticle object for the Gamma 
  G4DynamicParticle* g = new G4DynamicParticle(G4Gamma::Gamma(),
                                               gammaEnergy,
                                               gammaDirection);

  G4std::vector<G4DynamicParticle*>* vdp = new G4std::vector<G4DynamicParticle*>;
  vdp->push_back(g);

  return vdp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4Element* G4eBremsstrahlungModel::SelectRandomAtom(
           const G4Material* material) const
{
  // select randomly 1 element within the material

  G4int nElements = material->GetNumberOfElements();
  const G4ElementVector* theElementVector = material->GetElementVector();
  if(1 == nElements) return (*theElementVector)[0];
  else if(1 > nElements) return 0;  

  G4DataVector* dv = partialSumSigma[material->GetIndex()];
  G4double rval = G4UniformRand()*((*dv)[nElements-1]);
  for (G4int i=0; i<nElements; i++) {
    if (rval <= (*dv)[i]) return (*theElementVector)[i];
  }
  return (*theElementVector)[nElements-1];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungModel::SupressionFunction(const G4Material* material,
                                 G4double kineticEnergy, G4double gammaEnergy)
{
  // supression due to the LPM effect+polarisation of the medium/
  // supression due to the polarisation alone

  G4double TotalEnergy,TotalEnergySquare,LPMEnergy,LPMGammaEnergyLimit,
           LPMGammaEnergyLimit2,GammaEnergySquare,sp,s2lpm,supr,w,splim,Cnorm ;

  TotalEnergy = kineticEnergy+electron_mass_c2 ;
  TotalEnergySquare = TotalEnergy*TotalEnergy ;

  LPMEnergy = LPMconstant*(material->GetRadlen()) ;
  LPMGammaEnergyLimit = TotalEnergySquare/LPMEnergy ;
  GammaEnergySquare = gammaEnergy*gammaEnergy ;

  LPMGammaEnergyLimit2 = LPMGammaEnergyLimit*LPMGammaEnergyLimit;
  G4double electronDensity = material->GetElectronDensity();
  splim = LPMGammaEnergyLimit2/(LPMGammaEnergyLimit2+MigdalConstant*TotalEnergySquare*
                                   electronDensity) ;
  w = 1.+1./splim ;
  Cnorm = 2./(sqrt(w*w+4.)-w) ;

  sp = GammaEnergySquare/(GammaEnergySquare+MigdalConstant*TotalEnergySquare*
                                     electronDensity) ;
  if (theLPMflag)
    {
     s2lpm = LPMEnergy*gammaEnergy/TotalEnergySquare;
     if (s2lpm < 1.)
       {
        if ((1.-sp) < 1.e-6) w = s2lpm*(3.-sp);
        else                 w = s2lpm*(1.+1./sp);
        supr = Cnorm*(sqrt(w*w+4.*s2lpm)-w)/2. ;
       }
     else supr = sp;
    }
  else supr = sp;

  supr /= sp;
  return supr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


