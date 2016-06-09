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
// $Id: G4MuBremsstrahlungModel.cc,v 1.18 2005/08/18 14:38:55 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4MuBremsstrahlungModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 24.06.2002
//
// Modifications:
//
// 04-12-02 Change G4DynamicParticle constructor in PostStepDoIt (V.Ivanchenko)
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 24-01-03 Fix for compounds (V.Ivanchenko)
// 27-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 10-02-04 Add lowestKinEnergy (V.Ivanchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 03-08-05 Angular correlations according to PRM (V.Ivantchenko)
//

//
// Class Description:
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4MuBremsstrahlungModel.hh"
#include "G4Gamma.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleChangeForLoss.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// static members
//
G4double G4MuBremsstrahlungModel::zdat[]={1.,4.,13.,29.,92.};
G4double G4MuBremsstrahlungModel::adat[]={1.01,9.01,26.98,63.55,238.03};
G4double G4MuBremsstrahlungModel::tdat[]={1.e3,1.e4,1.e5,1.e6,1.e7,1.e8,1.e9,1.e10};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4MuBremsstrahlungModel::G4MuBremsstrahlungModel(const G4ParticleDefinition* p,
                                                 const G4String& nam)
  : G4VEmModel(nam),
  particle(0),
  lowestKinEnergy(1.0*GeV),
  minThreshold(1.0*keV),
  nzdat(5),
  ntdat(8),
  NBIN(1000),
  cutFixed(0.98*keV),
  samplingTablesAreFilled(false)
{
  theGamma = G4Gamma::Gamma();

  if(p) SetParticle(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuBremsstrahlungModel::~G4MuBremsstrahlungModel()
{
  size_t n = partialSumSigma.size();
  if(n > 0) {
    for(size_t i=0; i<n; i++) {
      delete partialSumSigma[i];
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuBremsstrahlungModel::MinEnergyCut(const G4ParticleDefinition*,
                                               const G4MaterialCutsCouple*)
{
  return minThreshold;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuBremsstrahlungModel::SetParticle(const G4ParticleDefinition* p)
{
  if(!particle) {
    particle = p;
    mass = particle->GetPDGMass();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuBremsstrahlungModel::Initialise(const G4ParticleDefinition* p,
                                         const G4DataVector& cuts)
{
  if(p) SetParticle(p);

  highKinEnergy = HighEnergyLimit();

  G4double fixedEnergy = 0.5*highKinEnergy;
//  G4double fixedEnergy = 500000.*TeV;

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  if(theCoupleTable) {
    G4int numOfCouples = theCoupleTable->GetTableSize();

    for (size_t ii=0; ii<partialSumSigma.size(); ii++){
      G4DataVector* a=partialSumSigma[ii];
      if ( a )  delete a;    
    } 
    partialSumSigma.clear();
    if(numOfCouples>0) {
      for (G4int i=0; i<numOfCouples; i++) {
	const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
	const G4Material* material = couple->GetMaterial();
	G4DataVector* dv = ComputePartialSumSigma(material, fixedEnergy,cuts[i]);
	partialSumSigma.push_back(dv);
      }
    }
  }
  if(!samplingTablesAreFilled) MakeSamplingTables();
  if(pParticleChange)
    fParticleChange = reinterpret_cast<G4ParticleChangeForLoss*>(pParticleChange);
  else
    fParticleChange = new G4ParticleChangeForLoss();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuBremsstrahlungModel::ComputeDEDXPerVolume(
					      const G4Material* material,
                                              const G4ParticleDefinition*,
                                                    G4double kineticEnergy,
                                                    G4double cutEnergy)
{
  G4double dedx = 0.0;
  if (kineticEnergy <= lowestKinEnergy) return dedx;

  G4double tmax = kineticEnergy;
  G4double cut  = min(cutEnergy,tmax);

  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector();

  //  loop for elements in the material
  for (size_t i=0; i<material->GetNumberOfElements(); i++) {

    G4double Z = (*theElementVector)[i]->GetZ();
    G4double A = (*theElementVector)[i]->GetA()/(g/mole) ;

    G4double loss = ComputMuBremLoss(Z, A, kineticEnergy, cut);

    dedx += loss*theAtomicNumDensityVector[i];
  }
  if(dedx < 0.) dedx = 0.;
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuBremsstrahlungModel::ComputMuBremLoss(G4double Z, G4double A,
                                                   G4double tkin, G4double cut)
{
  G4double totalEnergy = mass + tkin;
  G4double ak1 = 0.05;
  G4int    k2=5;
  G4double xgi[]={0.03377,0.16940,0.38069,0.61931,0.83060,0.96623};
  G4double wgi[]={0.08566,0.18038,0.23396,0.23396,0.18038,0.08566};
  G4double loss = 0.;

  G4double vcut = cut/totalEnergy;
  G4double vmax = tkin/totalEnergy;

  G4double aaa = 0.;
  G4double bbb = vcut;
  if(vcut>vmax) bbb=vmax ;
  G4int kkk = (G4int)((bbb-aaa)/ak1)+k2 ;
  G4double hhh=(bbb-aaa)/float(kkk) ;

  G4double aa = aaa;
  for(G4int l=0; l<kkk; l++)
  {
    for(G4int i=0; i<6; i++)
    {
      G4double ep = (aa + xgi[i]*hhh)*totalEnergy;
      loss += ep*wgi[i]*ComputeDMicroscopicCrossSection(tkin, Z, A, ep);
    }
    aa += hhh;
  }

  loss *=hhh*totalEnergy ;

  return loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuBremsstrahlungModel::ComputeMicroscopicCrossSection(
                                           G4double tkin,
                                           G4double Z,
                                           G4double A,
                                           G4double cut)
{
  G4double totalEnergy = tkin + mass;
  G4double ak1 = 2.3;
  G4int    k2  = 4;
  G4double xgi[]={0.03377,0.16940,0.38069,0.61931,0.83060,0.96623};
  G4double wgi[]={0.08566,0.18038,0.23396,0.23396,0.18038,0.08566};
  G4double cross = 0.;

  if(cut >= tkin) return cross;

  G4double vcut = cut/totalEnergy;
  G4double vmax = tkin/totalEnergy;

  G4double aaa = log(vcut);
  G4double bbb = log(vmax);
  G4int    kkk = (G4int)((bbb-aaa)/ak1)+k2 ;
  G4double hhh = (bbb-aaa)/float(kkk);

  G4double aa = aaa;

  for(G4int l=0; l<kkk; l++)
  {
    for(G4int i=0; i<6; i++)
    {
      G4double ep = exp(aa + xgi[i]*hhh)*totalEnergy;
      cross += ep*wgi[i]*ComputeDMicroscopicCrossSection(tkin, Z, A, ep);
    }
    aa += hhh;
  }

  cross *=hhh;

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuBremsstrahlungModel::ComputeDMicroscopicCrossSection(
                                           G4double tkin,
                                           G4double Z,
                                           G4double A,
                                           G4double gammaEnergy)
//  differential cross section
{
  static const G4double sqrte=sqrt(exp(1.)) ;
  static const G4double bh=202.4,bh1=446.,btf=183.,btf1=1429. ;
  static const G4double rmass=mass/electron_mass_c2 ;
  static const G4double cc=classic_electr_radius/rmass ;
  static const G4double coeff= 16.*fine_structure_const*cc*cc/3. ;

  G4double dxsection = 0.;

  if( gammaEnergy > tkin) return dxsection ;

  G4double E = tkin + mass ;
  G4double v = gammaEnergy/E ;
  G4double delta = 0.5*mass*mass*v/(E-gammaEnergy) ;
  G4double rab0=delta*sqrte ;

  G4double z13 = exp(-log(Z)/3.) ;
  G4double dn  = 1.54*exp(0.27*log(A)) ;

  G4double b,b1,dnstar ;

  if(Z<1.5)
  {
    b=bh;
    b1=bh1;
    dnstar=dn ;
  }
  else
  {
    b=btf;
    b1=btf1;
    dnstar = exp((1.-1./Z)*log(dn)) ;
  }

  // nucleus contribution logarithm
  G4double rab1=b*z13;
  G4double fn=log(rab1/(dnstar*(electron_mass_c2+rab0*rab1))*
              (mass+delta*(dnstar*sqrte-2.))) ;
  if(fn <0.) fn = 0. ;
  // electron contribution logarithm
  G4double epmax1=E/(1.+0.5*mass*rmass/E) ;
  G4double fe=0.;
  if(gammaEnergy<epmax1)
  {
    G4double rab2=b1*z13*z13 ;
    fe=log(rab2*mass/((1.+delta*rmass/(electron_mass_c2*sqrte))*
                              (electron_mass_c2+rab0*rab2))) ;
    if(fe<0.) fe=0. ;
  }

  dxsection = coeff*(1.-v*(1. - 0.75*v))*Z*(fn*Z + fe)/gammaEnergy;

  return dxsection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuBremsstrahlungModel::CrossSectionPerVolume(
					       const G4Material* material,
                                               const G4ParticleDefinition*,
                                                     G4double kineticEnergy,
                                                     G4double cutEnergy,
                                                     G4double maxEnergy)
{
  G4double cross = 0.0;
  if (cutEnergy >= maxEnergy || kineticEnergy <= lowestKinEnergy) return cross;
  
  G4double tmax = min(maxEnergy, kineticEnergy);
  G4double cut  = min(cutEnergy, tmax);

  const G4ElementVector* theElementVector = material->GetElementVector() ;
  const G4double* theAtomNumDensityVector = material->GetAtomicNumDensityVector();

  for (size_t i=0; i<material->GetNumberOfElements(); i++) {

    G4double Z = (*theElementVector)[i]->GetZ();
    G4double A = (*theElementVector)[i]->GetA()/(g/mole) ;

    G4double cr = ComputeMicroscopicCrossSection(kineticEnergy, Z, A, cut);

    if(tmax < kineticEnergy) {
      cr -= ComputeMicroscopicCrossSection(kineticEnergy, Z, A, tmax);
    }
    cross += theAtomNumDensityVector[i] * cr;
  }

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DataVector* G4MuBremsstrahlungModel::ComputePartialSumSigma(
                                       const G4Material* material,
                                             G4double kineticEnergy,
                                             G4double cut)

// Build the table of cross section per element. The table is built for MATERIALS.
// This table is used by DoIt to select randomly an element in the material.
{
  G4int nElements = material->GetNumberOfElements();
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = material->GetAtomicNumDensityVector();

  G4DataVector* dv = new G4DataVector();

  G4double cross = 0.0;

  for (G4int i=0; i<nElements; i++ ) {

    G4double Z = (*theElementVector)[i]->GetZ();
    G4double A = (*theElementVector)[i]->GetA()/(g/mole) ;
    cross += theAtomNumDensityVector[i] * ComputeMicroscopicCrossSection(kineticEnergy,
             Z, A, cut);
    dv->push_back(cross);
  }
  return dv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuBremsstrahlungModel::MakeSamplingTables()
{

  G4double AtomicNumber,AtomicWeight,KineticEnergy,
           TotalEnergy,Maxep;

  for (G4int iz=0; iz<nzdat; iz++)
   {
     AtomicNumber = zdat[iz];
     AtomicWeight = adat[iz]*g/mole ;

     for (G4int it=0; it<ntdat; it++)
     {
       KineticEnergy = tdat[it];
       TotalEnergy = KineticEnergy + mass;
       Maxep = KineticEnergy ;

       G4double CrossSection = 0.0 ;

       // calculate the differential cross section
       // numerical integration in
       //  log ...............
       G4double c = log(Maxep/cutFixed) ;
       G4double ymin = -5. ;
       G4double ymax = 0. ;
       G4double dy = (ymax-ymin)/NBIN ;

       G4double y = ymin - 0.5*dy ;
       G4double yy = ymin - dy ;
       G4double x = exp(y);
       G4double fac = exp(dy);
       G4double dx = exp(yy)*(fac - 1.0);

       for (G4int i=0 ; i<NBIN; i++)
       {
         y += dy ;
         x *= fac;
         dx*= fac;
         G4double ep = cutFixed*exp(c*x) ;

         CrossSection += ep*dx*ComputeDMicroscopicCrossSection(
                                                 KineticEnergy,AtomicNumber,
                                                 AtomicWeight,ep) ;
         ya[i]=y ;
         proba[iz][it][i] = CrossSection ;

       }

       proba[iz][it][NBIN] = CrossSection ;
       ya[NBIN] = 0. ;   //   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if(CrossSection > 0.)
       {
         for(G4int ib=0; ib<=NBIN; ib++)
         {
           proba[iz][it][ib] /= CrossSection ;
         }
       }
     }
   }
  samplingTablesAreFilled = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

vector<G4DynamicParticle*>* G4MuBremsstrahlungModel::SampleSecondaries(
                             const G4MaterialCutsCouple* couple,
                             const G4DynamicParticle* dp,
                                   G4double tmin,
                                   G4double maxEnergy)
{

  G4double kineticEnergy = dp->GetKineticEnergy();
  // check against insufficient energy
  G4double tmax = min(kineticEnergy, maxEnergy);
  if(tmin >= tmax) return 0;

  static const G4double ysmall = -100. ;
  static const G4double ytablelow = -5. ;

  G4ParticleMomentum partDirection = dp->GetMomentumDirection();

  // select randomly one element constituing the material
  const G4Element* anElement = SelectRandomAtom(couple);

  G4double totalEnergy   = kineticEnergy + mass;
  G4double totalMomentum = sqrt(kineticEnergy*(kineticEnergy + 2.0*mass));

  G4double dy = 5./G4float(NBIN);

  // This sampling should be checked!!! VI
  G4double ymin=log(log(tmin/cutFixed)/log(tmax/cutFixed));

  if(ymin < ysmall) return 0;

  //  sampling using tables

  G4double v,x,y ;
  G4int iy;
  // select sampling table ;
  G4double lnZ = log(anElement->GetZ()) ;
  G4double delmin = 1.e10 ;
  G4double del ;
  G4int izz = 0;
  G4int itt = 0;
  G4int NBINminus1;
  NBINminus1 = NBIN-1 ;
  for (G4int iz=0; iz<nzdat; iz++)
  {
    del = std::abs(lnZ-log(zdat[iz])) ;
    if(del<delmin)
    {
       delmin=del ;
       izz=iz ;
    }
  }

  delmin = 1.e10 ;
  for (G4int it=0; it<ntdat; it++)
  {
    del = std::abs(log(tmax)-log(tdat[it])) ;
    if(del<delmin)
    {
      delmin=del;
      itt=it ;
    }
  }
  G4int iymin = G4int((ymin+5.)/dy+0.5) ;

  do {
    if(ymin < ytablelow)
    {
      y = ymin + G4UniformRand()*(ytablelow-ymin) ;
    }
    else
    {
      G4double r = G4UniformRand() ;

      iy = iymin-1 ;
      delmin = proba[izz][itt][NBINminus1]-proba[izz][itt][iymin] ;
      do {
         iy += 1 ;
      } while ((r > (proba[izz][itt][iy]-proba[izz][itt][iymin])/delmin)
                 &&(iy < NBINminus1)) ;

      //sampling is Done uniformly in y in the bin
      y = ya[iy] + G4UniformRand() * ( ya[iy+1] - ya[iy] ) ;
    }

    x = exp(y) ;

    v = cutFixed*exp(x*log(tmax/cutFixed)) ;

  } while ( v <= 0.);

  // create G4DynamicParticle object for the Gamma
  G4double gEnergy = v;

  // sample angle
  G4double gam  = totalEnergy/mass;
  G4double rmax = gam*min(1.0, totalEnergy/gEnergy - 1.0);
  rmax *= rmax;
  x = G4UniformRand()*rmax/(1.0 + rmax);

  G4double theta = sqrt(x/(1.0 - x))/gam;
  G4double sint  = sin(theta);
  G4double phi   = twopi * G4UniformRand() ;
  G4double dirx  = sint*cos(phi), diry = sint*sin(phi), dirz = cos(theta) ;

  G4ThreeVector gDirection(dirx, diry, dirz);
  gDirection.rotateUz(partDirection);

  partDirection *= totalMomentum;
  partDirection -= gEnergy*gDirection;
  partDirection = partDirection.unit();

  // primary change
  kineticEnergy -= gEnergy;
  fParticleChange->SetProposedKineticEnergy(kineticEnergy);
  fParticleChange->SetProposedMomentumDirection(partDirection);

  // save secondary
  G4DynamicParticle* aGamma = new G4DynamicParticle(theGamma,gDirection,gEnergy);
  vector<G4DynamicParticle*>* vdp = new vector<G4DynamicParticle*>;
  vdp->push_back(aGamma);

  return vdp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4Element* G4MuBremsstrahlungModel::SelectRandomAtom(
           const G4MaterialCutsCouple* couple) const
{
  // select randomly 1 element within the material

  const G4Material* material = couple->GetMaterial();
  G4int nElements = material->GetNumberOfElements();
  const G4ElementVector* theElementVector = material->GetElementVector();
  if(1 == nElements) return (*theElementVector)[0];
  else if(1 > nElements) return 0;

  G4DataVector* dv = partialSumSigma[couple->GetIndex()];
  G4double rval = G4UniformRand()*((*dv)[nElements-1]);
  for (G4int i=0; i<nElements; i++) {
    if (rval <= (*dv)[i]) return (*theElementVector)[i];
  }
  return (*theElementVector)[nElements-1];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
