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
// File name:     G4MuPairProductionModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 24.06.2002
//
// Modifications:
//
// 04-12-02 Change G4DynamicParticle constructor in PostStep (V.Ivanchenko)
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 24-01-03 Fix for compounds (V.Ivanchenko)
// 27-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 Add model (V.Ivanchenko)

//
// Class Description:
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4MuPairProductionModel.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ProductionCutsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// static members
//
G4double G4MuPairProductionModel::zdat[]={1.,4.,13.,29.,92.};
G4double G4MuPairProductionModel::adat[]={1.01,9.01,26.98,63.55,238.03};
G4double G4MuPairProductionModel::tdat[]={1.e3,1.e4,1.e5,1.e6,1.e7,1.e8,1.e9,1.e10};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuPairProductionModel::G4MuPairProductionModel(const G4ParticleDefinition* p,
                                                 const G4String& nam)
  : G4VEmModel(nam),
  minPairEnergy(4.*electron_mass_c2),
  highKinEnergy(1000000.*TeV),
  lowKinEnergy(minPairEnergy),
  nzdat(5),
  ntdat(8),
  NBIN(1000),
  samplingTablesAreFilled(false)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuPairProductionModel::~G4MuPairProductionModel()
{
  size_t n = partialSumSigma.size();
  if(n > 0) {
    for(size_t i=0; i<n; i++) {
      delete partialSumSigma[i];
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuPairProductionModel::HighEnergyLimit(const G4ParticleDefinition*)
{
  return highKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuPairProductionModel::LowEnergyLimit(const G4ParticleDefinition*)
{
  return lowKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuPairProductionModel::MinEnergyCut(const G4ParticleDefinition*,
                                               const G4MaterialCutsCouple* couple)
{

  size_t index = couple->GetIndex();
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();

  G4double eCut = (*(theCoupleTable->GetEnergyCutsVector(1)))[index];
  G4double pCut = (*(theCoupleTable->GetEnergyCutsVector(2)))[index];
  G4double x = 2.0*electron_mass_c2 + eCut + pCut;
  if(x < minPairEnergy) x = minPairEnergy;
/*
  if(eCut < highKinEnergy && pCut < highKinEnergy) {
    x +=  eCut + pCut;
  } else {
    x = 0.5*highKinEnergy;
  }
  */
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4MuPairProductionModel::IsInCharge(const G4ParticleDefinition* p)
{
  return (p == G4MuonMinus::MuonMinus() || p == G4MuonPlus::MuonPlus());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuPairProductionModel::Initialise(const G4ParticleDefinition*,
                                         const G4DataVector& cuts)
{
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  G4double fixedEnergy = sqrt(lowKinEnergy*highKinEnergy);

  partialSumSigma.clear();
  for (size_t i=0; i<numOfCouples; i++) {
    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
    const G4Material* material = couple->GetMaterial();
    G4DataVector* dv = ComputePartialSumSigma(material, fixedEnergy,
                             G4std::min(cuts[i], 0.25*highKinEnergy));
    partialSumSigma.push_back(dv);
  }
  if(!samplingTablesAreFilled) MakeSamplingTables();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuPairProductionModel::ComputeDEDX(const G4Material* material,
                                              const G4ParticleDefinition* p,
                                                    G4double kineticEnergy,
                                                    G4double cutEnergy)
{
  G4double dedx = 0.0;
  if(minPairEnergy >= cutEnergy) return dedx;
  G4double cut = cutEnergy;
  if(kineticEnergy <= cutEnergy) cut = kineticEnergy;

  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector();

  //  loop for elements in the material
  for (size_t i=0; i<material->GetNumberOfElements(); i++) {

    G4double Z = (*theElementVector)[i]->GetZ();

    G4double loss = ComputMuPairLoss(Z, kineticEnergy, cut);

    dedx += loss*theAtomicNumDensityVector[i];
  }
  if(dedx < 0.) dedx = 0.;
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuPairProductionModel::ComputMuPairLoss(G4double Z,
                                                   G4double tkin, G4double cutEnergy)
{
  static const
  G4double xgi[] ={ 0.0199,0.1017,0.2372,0.4083,0.5917,0.7628,0.8983,0.9801};
  static const
  G4double wgi[] ={ 0.0506,0.1112,0.1569,0.1813,0.1813,0.1569,0.1112,0.0506};
  static const G4double ak1=6.9;
  static const G4double ak2=1.0;
  static const G4double sqrte = sqrt(exp(1.));
  static const G4double aaa = log(minPairEnergy);
  G4double z13 = pow(Z,0.333333333);

  G4double loss = 0.0 ;

  G4double particleMass = (G4MuonPlus::MuonPlus())->GetPDGMass();
  G4double tmax = tkin + particleMass*(1.-0.75*sqrte*z13);

  //  G4cout << "###DEDX tkin= " << tkin << " tmax= " << tmax << " tmin= " << minPairEnergy << G4endl;

  G4double cut = cutEnergy;
  if(tmax <= cutEnergy) cut = tmax;
  if(cut <= minPairEnergy) return loss;

  // calculate the rectricted loss
  // numerical integration in log(PairEnergy)
  G4double bbb = log(cut) ;
  G4int    kkk = (G4int)((bbb-aaa)/ak1+ak2);
  if(kkk > 8) kkk = 8;
  G4double hhh = (bbb-aaa)/(G4double)kkk ;
  G4double x = aaa;

  //  G4cout << "###DEDX tkin= " << tkin << " cut= " << cut << " kkk= " << kkk << G4endl;

  for (G4int l=0 ; l<kkk; l++)
  {

    for (G4int ll=0; ll<8; ll++)
    {
      G4double ep = exp(x+xgi[ll]*hhh);
      //  G4cout << "ep= " << ep << G4endl;
      loss += wgi[ll]*ep*ep*ComputeDMicroscopicCrossSection(tkin, Z, ep);
    }
    x += hhh;
  }
  loss *= hhh ;
  //  cout << "### tmax= " << tmax << " hhh= " << hhh << " loss= " << loss << endl;
  if (loss < 0.) loss = 0.;
  return loss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuPairProductionModel::ComputeMicroscopicCrossSection(
                                           G4double tkin,
                                           G4double Z,
                                           G4double A,
                                           G4double cut)

{
  static const G4double ak1=6.9 ;
  static const G4double ak2=1.0 ;
  static const G4double sqrte = sqrt(exp(1.)) ;
  static const G4double
  xgi[]={ 0.0199,0.1017,0.2372,0.4083,0.5917,0.7628,0.8983,0.9801 };
  static const G4double
  wgi[]={ 0.0506,0.1112,0.1569,0.1813,0.1813,0.1569,0.1112,0.0506 };
  G4double z13 = pow(Z,0.333333333);

  G4double cross = 0. ;

  G4double particleMass = (G4MuonPlus::MuonPlus())->GetPDGMass();
  G4double tmax = tkin + particleMass*(1.-0.75*sqrte*z13);

  if(tmax <= cut) return cross;

  G4double aaa = log(cut);
  G4double bbb = log(tmax);
  G4int kkk = (G4int)((bbb-aaa)/ak1 + ak2);
  if(kkk > 8) kkk = 8;
  G4double hhh = (bbb-aaa)/float(kkk);
  G4double x = aaa;

  //  G4cout << "###Cross tkin= " << tkin << " cut= " << cut << " kkk= " << kkk << G4endl;

  for(G4int l=0; l<kkk; l++)
  {

    for(G4int i=0; i<8; i++)
    {
      G4double ep = exp(x + xgi[i]*hhh);

      cross += ep*wgi[i]*ComputeDMicroscopicCrossSection(tkin, Z, ep);
    }
    aaa += hhh;
  }

  cross *=hhh;
  if(cross < 0.0) cross = 0.0;

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuPairProductionModel::ComputeDMicroscopicCrossSection(
                                           G4double tkin,
                                           G4double Z,
                                           G4double pairEnergy)
 // Calculates the  differential (D) microscopic cross section
 // using the cross section formula of R.P. Kokoulin (18/01/98)
{

  static const G4double
  xgi[] ={ 0.0199,0.1017,0.2372,0.4083,0.5917,0.7628,0.8983,0.9801 };

  static const G4double
  wgi[] ={ 0.0506,0.1112,0.1569,0.1813,0.1813,0.1569,0.1112,0.0506 };

  G4double cross = 0.;

  G4double particleMass = (G4MuonPlus::MuonPlus())->GetPDGMass();
  G4double totalEnergy  = tkin + particleMass;
  G4double energyLoss   = totalEnergy - pairEnergy;
  G4double a = 6.*particleMass*particleMass/(totalEnergy*energyLoss) ;
  G4double b = 4.*electron_mass_c2/pairEnergy;

  G4double tmn = (b+2.*a*(1.-b))/(1.+(1.-a)*sqrt(1.-b));

  if(tmn <= 0.) return cross;

  //  G4cout << "a= " << a << " b= " << b << " tmn= " << tmn << G4endl;


  tmn = log(tmn);

  //  G4cout << "a= " << a << " b= " << b << " tmn= " << tmn << G4endl;


  // Gaussian integration in ln(1-ro) ( with 8 points)
  for (G4int i=0; i<7; i++)
  {
    G4double ro = 1.-exp(tmn*xgi[i]) ;

    cross += wgi[i]*(1.-ro)*ComputeDDMicroscopicCrossSection(tkin,Z,pairEnergy,ro);
    //    cout << "ro= " << ro << " cross= " << cross << endl;
  }

  cross *= -tmn ;

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuPairProductionModel::ComputeDDMicroscopicCrossSection(
                                           G4double tkin,
                                           G4double Z,
                                           G4double pairEnergy,
                                           G4double asymmetry)
 // Calculates the  differential (D) microscopic cross section
 // using the cross section formula of R.P. Kokoulin (18/01/98)
{
  static const G4double sqrte = sqrt(exp(1.)) ;

  G4double bbbtf= 183. ;
  G4double bbbh = 202.4 ;
  G4double g1tf = 1.95e-5 ;
  G4double g2tf = 5.3e-5 ;
  G4double g1h  = 4.4e-5 ;
  G4double g2h  = 4.8e-5 ;

  G4double particleMass = (G4MuonPlus::MuonPlus())->GetPDGMass();
  G4double totalEnergy  = tkin + particleMass;
  G4double energyLoss   = totalEnergy - pairEnergy;
  G4double massratio    = particleMass/electron_mass_c2 ;
  G4double massratio2   = massratio*massratio ;

  G4double z13 = pow(Z,0.333333333);
  G4double z23 = z13*z13 ;

  G4double c3 = 3.*sqrte*particleMass/4. ;

  G4double DDCrossSection = 0. ;

  if(energyLoss <= c3*z13) return DDCrossSection ;

  G4double c7 = 4.*electron_mass_c2 ;
  G4double c8 = 6.*particleMass*particleMass ;
  G4double alf = c7/pairEnergy ;
  G4double a3 = 1. - alf ;

  if(a3 <= 0.) return DDCrossSection ;

 // zeta calculation
  G4double bbb,g1,g2,zeta1,zeta2,zeta,z2 ;
  if( Z < 1.5 )
  {
    bbb = bbbh ;
    g1  = g1h ;
    g2  = g2h ;
  }
  else
  {
    bbb = bbbtf ;
    g1  = g1tf ;
    g2  = g2tf ;
  }
  zeta1 = 0.073 * log(totalEnergy/(particleMass+g1*z23*totalEnergy))-0.26 ;
  if( zeta1 > 0.)
  {
    zeta2 = 0.058*log(totalEnergy/(particleMass+g2*z13*totalEnergy))-0.14 ;
    zeta  = zeta1/zeta2 ;
  }
  else
  {
    zeta = 0. ;
  }

  z2 = Z*(Z+zeta) ;

  G4double screen0 = 2.*electron_mass_c2*sqrte*bbb/(z13*pairEnergy) ;
  G4double a0 = totalEnergy*energyLoss ;
  G4double a1 = pairEnergy*pairEnergy/a0 ;
  G4double bet = 0.5*a1 ;
  G4double xi0 = 0.25*massratio2*a1 ;
  G4double del = c8/a0 ;

  G4double romin = 0. ;
  G4double romax = (1.-del)*sqrt(1.-c7/pairEnergy) ;

  if((asymmetry < romin) || (asymmetry > romax)) return DDCrossSection ;

  G4double a4 = 1.-asymmetry ;
  G4double a5 = a4*(2.-a4) ;
  G4double a6 = 1.-a5 ;
  G4double a7 = 1.+a6 ;
  G4double a9 = 3.+a6 ;
  G4double xi = xi0*a5 ;
  G4double xii = 1./xi ;
  G4double xi1 = 1.+xi ;
  G4double screen = screen0*xi1/a5 ;

  G4double yeu = 5.-a6+4.*bet*a7 ;
  G4double yed = 2.*(1.+3.*bet)*log(3.+xii)-a6-a1*(2.-a6) ;
  G4double yel = 1.+yeu/yed ;
  G4double ale=log(bbb/z13*sqrt(xi1*yel)/(1.+screen*yel)) ;
  G4double cre = 0.5*log(1.+2.25/(massratio2*z23)*xi1*yel) ;
  G4double be ;
  if(xi <= 1.e3)
    be = ((2.+a6)*(1.+bet)+xi*a9)*log(1.+xii)+(a5-bet)/xi1-a9;
  else
    be = (3.-a6+a1*a7)/(2.+xi) ;
  G4double fe = (ale-cre)*be ;
  if( fe < 0.)
    fe = 0. ;

  G4double ymu = 4.+a6 +3.*bet*a7 ;
  G4double ymd = a7*(1.5+a1)*log(3.+xi)+1.-1.5*a6 ;
  G4double ym1 = 1.+ymu/ymd ;
  G4double alm_crm = log(bbb*massratio/(1.5*z23*(1.+screen*ym1))) ;
  G4double a10,bm ;
  if( xi >= 1.e-3)
  {
    a10 = (1.+a1)*a5 ;
    bm  = (a7*(1.+1.5*bet)-a10*xii)*log(xi1)+xi*(a5-bet)/xi1+a10 ;
  }
  else
    bm = (5.-a6+bet*a9)*(xi/2.) ;
  G4double fm = alm_crm*bm ;
  if( fm < 0.)
    fm = 0. ;

  DDCrossSection = (fe+fm/massratio2) ;

  DDCrossSection *= 4.*fine_structure_const*fine_structure_const
                   *classic_electr_radius*classic_electr_radius/(3.*pi) ;

  DDCrossSection *= z2*energyLoss/(totalEnergy*pairEnergy) ;


  return DDCrossSection ;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuPairProductionModel::CrossSection(const G4Material* material,
                                               const G4ParticleDefinition* p,
                                                     G4double kineticEnergy,
                                                     G4double cutEnergy,
                                                     G4double maxEnergy)
{
  G4double cross = 0.0;

  G4double tmax = G4std::min(maxEnergy, kineticEnergy);
  if(cutEnergy >= tmax) return cross;

  const G4ElementVector* theElementVector = material->GetElementVector() ;
  const G4double* theAtomNumDensityVector = material->GetAtomicNumDensityVector();

  for (size_t i=0; i<material->GetNumberOfElements(); i++) {

    G4double Z = (*theElementVector)[i]->GetZ();
    G4double A = (*theElementVector)[i]->GetA()/(g/mole) ;

    G4double cr = ComputeMicroscopicCrossSection(kineticEnergy, Z, A, cutEnergy);

    if(maxEnergy < kineticEnergy) {
      cr -= ComputeMicroscopicCrossSection(kineticEnergy, Z, A, maxEnergy);
    }
    cross += theAtomNumDensityVector[i] * cr;
  }

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DataVector* G4MuPairProductionModel::ComputePartialSumSigma(
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

void G4MuPairProductionModel::MakeSamplingTables()
{
  static const G4double sqrte = sqrt(exp(1.)) ;
  G4double particleMass = (G4MuonPlus::MuonPlus())->GetPDGMass();

  for (G4int iz=0; iz<nzdat; iz++)
  {
    G4double atomicNumber = zdat[iz];
    G4double z13 = exp(log(atomicNumber)/3.) ;

    for (G4int it=0; it<ntdat; it++)
    {
      G4double kineticEnergy = tdat[it];
      G4double maxPairEnergy = kineticEnergy+particleMass*(1.-0.75*sqrte*z13) ;

      G4double CrossSection = 0.0 ;

      G4double ymin = -5. ;
      G4double ymax = 0. ;
      G4double dy = (ymax-ymin)/NBIN ;

      G4double y = ymin - 0.5*dy ;
      G4double yy = ymin - dy ;
      G4double x = exp(y);
      G4double fac = exp(dy);
      G4double dx = exp(yy)*(fac - 1.0);

      if(maxPairEnergy > minPairEnergy) {
        G4double c = log(maxPairEnergy/minPairEnergy) ;

        for (G4int i=0 ; i<NBIN; i++)
        {
          y += dy ;
          x *= fac;
          dx*= fac;
          G4double ep = minPairEnergy*exp(c*x) ;
          CrossSection += ep*dx*ComputeDMicroscopicCrossSection(
                                      kineticEnergy, atomicNumber, ep);
          ya[i]=y ;
          proba[iz][it][i] = CrossSection ;

        }
      } else {

        for (G4int i=0 ; i<NBIN; i++)
	  {
          y += dy ;
          ya[i]=y ;
          proba[iz][it][i] = 0.0 ;
	  }
      }

      ya[NBIN]=0. ;

      proba[iz][it][NBIN] = CrossSection ;

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

G4DynamicParticle* G4MuPairProductionModel::SampleSecondary(
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
                                   G4double,
                                   G4double)
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4std::vector<G4DynamicParticle*>* G4MuPairProductionModel::SampleSecondaries(
                             const G4MaterialCutsCouple* couple,
                             const G4DynamicParticle* aDynamicParticle,
                                   G4double minEnergy,
                                   G4double maxEnergy)
{
   static const G4double esq = sqrt(exp(1.));
   G4double kineticEnergy = aDynamicParticle->GetKineticEnergy();
   G4double particleMass  = aDynamicParticle->GetDefinition()->GetPDGMass();
   G4ParticleMomentum ParticleDirection =
                                      aDynamicParticle->GetMomentumDirection();

   // select randomly one element constituing the material
   const G4Element* anElement = SelectRandomAtom(couple);

   // limits of the energy sampling
   G4double totalEnergy = kineticEnergy + particleMass ;
   //G4double TotalMomentum = sqrt(KineticEnergy*(TotalEnergy+particleMass)) ;
   G4double Z3 = anElement->GetIonisation()->GetZ3() ;
   G4double maxPairEnergy = totalEnergy-0.75*esq*particleMass*Z3 ;
   if(maxPairEnergy > maxEnergy) maxPairEnergy = maxEnergy;

   // check against insufficient energy
   if(minEnergy >= maxPairEnergy) return 0;

   // sample e-e+ energy, pair energy first
   G4double PairEnergy,x,yc,y ;
   // G4int iZ,iT;
   G4int iy ;

   // select sampling table ;
   G4double lnZ = log(anElement->GetZ()) ;
   G4double delmin = 1.e10 ;
   G4double del ;
   G4int izz = 0;
   G4int itt = 0;
   G4int NBINminus1 = NBIN-1;
   for (G4int iz=0; iz<nzdat; iz++)
   {
     del = abs(lnZ-log(zdat[iz])) ;
     if(del<delmin)
     {
        delmin=del ;
        izz=iz ;
     }
   }
   delmin = 1.e10 ;
   for (G4int it=0; it<ntdat; it++)
   {
     del = abs(log(kineticEnergy)-log(tdat[it])) ;
     if(del<delmin)
     {
       delmin=del;
       itt=it ;
     }
   }

   if( minEnergy <= minPairEnergy)
     iy = 0 ;
   else
   {
     G4double xc = log(minEnergy/minPairEnergy)/log(maxPairEnergy/minPairEnergy) ;
     yc = log(xc) ;

     iy = -1 ;
     do {
         iy += 1 ;
        } while ((ya[iy] < yc )&&(iy < NBINminus1)) ;
   }

   G4double norm = proba[izz][itt][iy] ;

   G4double r = norm+G4UniformRand()*(1.-norm) ;

   iy -= 1 ;
   do {
        iy += 1 ;
      } while ((proba[izz][itt][iy] < r)&&(iy < NBINminus1)) ;

   //sampling is uniformly in y in the bin
   if( iy < NBIN )
     y = ya[iy] + G4UniformRand() * ( ya[iy+1] - ya[iy]) ;
   else
     y = ya[iy] ;

   x = exp(y) ;

   PairEnergy = minPairEnergy*exp(x*log(maxPairEnergy/minPairEnergy)) ;

  // sample r=(E+-E-)/PairEnergy  ( uniformly .....)
   G4double rmax = (1.-6.*particleMass*particleMass/(totalEnergy*
                                               (totalEnergy-PairEnergy)))
                                       *sqrt(1.-minPairEnergy/PairEnergy) ;
   r = rmax * (-1.+2.*G4UniformRand()) ;

  // compute energies from PairEnergy,r
   G4double ElectronEnergy=(1.-r)*PairEnergy/2. ;
   G4double PositronEnergy=(1.+r)*PairEnergy/2. ;

   //  angles of the emitted particles ( Z - axis along the parent particle)
   //      (mean theta for the moment)
   G4double Teta = electron_mass_c2/totalEnergy ;

   G4double Phi  = twopi * G4UniformRand() ;
   G4double dirx = sin(Teta)*cos(Phi) , diry = sin(Teta)*sin(Phi) ,
            dirz = cos(Teta) ;

   G4double ElectronMomentum , PositronMomentum ;
   //G4double finalPx,finalPy,finalPz ;
   G4double ElectKineEnergy = ElectronEnergy - electron_mass_c2 ;

   ElectronMomentum = sqrt(ElectKineEnergy*(ElectronEnergy+electron_mass_c2));
   G4ThreeVector ElectDirection ( dirx, diry, dirz );
   ElectDirection.rotateUz(ParticleDirection);

   // create G4DynamicParticle object for the particle1
   G4DynamicParticle* aParticle1= new G4DynamicParticle();
   aParticle1->SetDefinition(G4Electron::Electron());
   aParticle1->SetMomentumDirection(ElectDirection);
   aParticle1->SetKineticEnergy(ElectKineEnergy);


   G4double PositKineEnergy = PositronEnergy - electron_mass_c2 ;
   PositronMomentum = sqrt(PositKineEnergy*(PositronEnergy+electron_mass_c2));

   G4ThreeVector PositDirection ( -dirx, -diry, dirz );
   PositDirection.rotateUz(ParticleDirection);

   // create G4DynamicParticle object for the particle2
   G4DynamicParticle* aParticle2= new G4DynamicParticle();
   aParticle2->SetDefinition(G4Positron::Positron());
   aParticle2->SetMomentumDirection(PositDirection);
   aParticle2->SetKineticEnergy(PositKineEnergy);


  G4std::vector<G4DynamicParticle*>* vdp = new G4std::vector<G4DynamicParticle*>;
  vdp->push_back(aParticle1);
  vdp->push_back(aParticle2);

  return vdp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4Element* G4MuPairProductionModel::SelectRandomAtom(
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


