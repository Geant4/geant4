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
// 24.11.08 V. Grichine - first implementation
//

#include "G4ComponentGGNuclNuclXsc.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadTmpUtil.hh"
#include "G4HadronNucleonXsc.hh"


G4ComponentGGNuclNuclXsc::G4ComponentGGNuclNuclXsc() 
 : G4VComponentCrossSection("Glauber-Gribov nucleus nucleus"),
//   fUpperLimit(100000*GeV),
   fLowerLimit(0.1*MeV),
   fRadiusConst(1.08*fermi),  // 1.1, 1.3 ?
   fTotalXsc(0.0), fElasticXsc(0.0), fInelasticXsc(0.0), fProductionXsc(0.0),
   fDiffractionXsc(0.0),
    cacheDP(G4Proton::Proton(),G4ParticleMomentum(1.,0,0),0),
    dProton(G4Proton::Proton(),G4ParticleMomentum(1.,0,0),0),
    dNeutron(G4Neutron::Neutron(),G4ParticleMomentum(1.,0,0),0)
// , fHadronNucleonXsc(0.0)
{
  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  hnXsc = new G4HadronNucleonXsc();
}


G4ComponentGGNuclNuclXsc::~G4ComponentGGNuclNuclXsc()
{
  delete hnXsc;
}

////////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::GetTotalIsotopeCrossSection(const G4ParticleDefinition* aParticle,
				       G4double kinEnergy,
				       G4int Z, G4int A)
{
  cacheDP.SetDefinition(aParticle);
  cacheDP.SetKineticEnergy(kinEnergy);
  fInelasticXsc = GetZandACrossSection(&cacheDP, Z, A);
  return fTotalXsc;
}

//////////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::GetTotalElementCrossSection(const G4ParticleDefinition* aParticle,
				       G4double kinEnergy, 
				       G4int Z, G4double A)
{
  cacheDP.SetDefinition(aParticle);
  cacheDP.SetKineticEnergy(kinEnergy);
  fInelasticXsc = GetZandACrossSection(&cacheDP, Z, G4int(A));
  return fTotalXsc;
}

////////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::GetInelasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4int A)
{
  cacheDP.SetDefinition(aParticle);
  cacheDP.SetKineticEnergy(kinEnergy);
  fInelasticXsc = GetZandACrossSection(&cacheDP, Z, A);
  return fInelasticXsc;
}

/////////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::GetInelasticElementCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4double A)
{
  cacheDP.SetDefinition(aParticle);
  cacheDP.SetKineticEnergy(kinEnergy);
  fInelasticXsc = GetZandACrossSection(&cacheDP, Z, G4int(A));
  return fInelasticXsc;
}

//////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::GetElasticElementCrossSection(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4double A)
{
  cacheDP.SetDefinition(aParticle);
  cacheDP.SetKineticEnergy(kinEnergy);
  fInelasticXsc = GetZandACrossSection(&cacheDP, Z, G4int(A));
  return fElasticXsc;
}

///////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::GetElasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4int A)
{
  cacheDP.SetDefinition(aParticle);
  cacheDP.SetKineticEnergy(kinEnergy);
  fInelasticXsc = GetZandACrossSection(&cacheDP, Z, A);
  return fElasticXsc;
}

////////////////////////////////////////////////////////////////
 
G4double G4ComponentGGNuclNuclXsc::ComputeQuasiElasticRatio(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4int A)
{
  cacheDP.SetDefinition(aParticle);
  cacheDP.SetKineticEnergy(kinEnergy);
  fInelasticXsc = GetZandACrossSection(&cacheDP, Z, A);
  G4double ratio = 0.;

  if(fInelasticXsc > 0.)
  {
    ratio = (fInelasticXsc - fProductionXsc)/fInelasticXsc;
    if(ratio < 0.) ratio = 0.;
  }
  return ratio;
}
 
//////////////////////////////////////////////////////////////////////

void
G4ComponentGGNuclNuclXsc::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4ComponentGGNuclNuclXsc calculates total, inelastic and\n"
          << "elastic cross sections for nucleus-nucleus collisions using\n"
          << "the Glauber model with Gribov corrections.  It is valid for\n"
          << "all incident energies above 100 keV./n";
}

/////////////////////////////////////////////////////////////////////

G4bool 
G4ComponentGGNuclNuclXsc::IsElementApplicable(const G4DynamicParticle* aDP, 
					      G4int Z, const G4Material*)
{
  G4bool applicable = false;
  G4double kineticEnergy = aDP->GetKineticEnergy();

  if (kineticEnergy >= fLowerLimit && Z > 1) applicable = true;
  return applicable;
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculates total and inelastic Xsc, derives elastic as total - inelastic 
// accordong to Glauber model with Gribov correction calculated in the dipole 
// approximation on light cone. Gaussian density helps to calculate rest 
// integrals of the model. [1] B.Z. Kopeliovich, nucl-th/0306044 


G4double G4ComponentGGNuclNuclXsc::
GetElementCrossSection(const G4DynamicParticle* aParticle, G4int Z,
		       const G4Material*)
{
  G4int A = G4lrint(G4NistManager::Instance()->GetAtomicMassAmu(Z));
  return GetZandACrossSection(aParticle, Z, A);
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculates total and inelastic Xsc, derives elastic as total - inelastic 
// accordong to Glauber model with Gribov correction calculated in the dipole 
// approximation on light cone. Gaussian density of point-like nucleons helps 
// to calculate rest integrals of the model. [1] B.Z. Kopeliovich,
// nucl-th/0306044 + simplification above


G4double G4ComponentGGNuclNuclXsc::
GetZandACrossSection(const G4DynamicParticle* aParticle,
                     G4int tZ, G4int tA)
{
  G4double xsection;
  G4double sigma;
  G4double cofInelastic = 2.4;
  G4double cofTotal = 2.0;
  G4double nucleusSquare;
  G4double cB;
  G4double ratio;

  G4double pZ = aParticle->GetDefinition()->GetPDGCharge();
  G4double pA = aParticle->GetDefinition()->GetBaryonNumber();

  G4double pTkin = aParticle->GetKineticEnergy();  
  pTkin /= pA;

  G4double pN = pA - pZ;
  if( pN < 0. ) pN = 0.;

  G4double tN = tA - tZ;
  if( tN < 0. ) tN = 0.;

  G4double tR = GetNucleusRadius( G4double(tZ),G4double(tA) );  
  G4double pR = GetNucleusRadius(pZ,pA); 

  cB = GetCoulombBarier(aParticle, G4double(tZ), G4double(tA), pR, tR);

  if ( cB > 0. ) 
  {
    dProton.SetKineticEnergy(pTkin);
    dNeutron.SetKineticEnergy(pTkin);

    sigma = (pZ*tZ+pN*tN)*hnXsc->GetHadronNucleonXscNS(&dProton, theProton);

    G4double ppInXsc = hnXsc->GetInelasticHadronNucleonXsc();

    sigma += (pZ*tN+pN*tZ)*hnXsc->GetHadronNucleonXscNS(&dNeutron, theProton);

    G4double npInXsc = hnXsc->GetInelasticHadronNucleonXsc();

    // G4cout<<"ppInXsc = "<<ppInXsc/millibarn<<"; npInXsc = "<<npInXsc/millibarn<<G4endl;
    // G4cout<<"npTotXsc = "<<hnXsc->GetTotalHadronNucleonXsc()/millibarn<<"; npElXsc = "
    //                      <<hnXsc->GetElasticHadronNucleonXsc()/millibarn<<G4endl;

    nucleusSquare = cofTotal*pi*( pR*pR + tR*tR );   // basically 2piRR

    ratio      = sigma/nucleusSquare;
    xsection   = nucleusSquare*G4Log( 1. + ratio );
    fTotalXsc  = xsection;
    fTotalXsc *= cB;

    fInelasticXsc = nucleusSquare*G4Log( 1. + cofInelastic*ratio )/cofInelastic;

    fInelasticXsc *= cB;
    fElasticXsc   = fTotalXsc - fInelasticXsc;

    // if (fElasticXsc < DBL_MIN) fElasticXsc = DBL_MIN;
    /*    
    G4double difratio = ratio/(1.+ratio);

    fDiffractionXsc = 0.5*nucleusSquare*( difratio - G4Log( 1. + difratio ) );
    */
    // production to be checked !!! edit MK xsc

    //sigma = (pZ*tZ+pN*tN)*GetHadronNucleonXscMK(theProton, pTkin, theProton) +
    //      (pZ*tN+pN*tZ)*GetHadronNucleonXscMK(theProton, pTkin, theNeutron);

    sigma = (pZ*tZ+pN*tN)*ppInXsc + (pZ*tN+pN*tZ)*npInXsc;
 
    ratio = sigma/nucleusSquare;
    fProductionXsc = nucleusSquare*G4Log( 1. + cofInelastic*ratio )/cofInelastic;

    if (fElasticXsc < 0.) fElasticXsc = 0.;
  }
  else
  {
    fInelasticXsc  = 0.;
    fTotalXsc      = 0.;
    fElasticXsc    = 0.;
    fProductionXsc = 0.;
  }
  return fInelasticXsc;   // xsection; 
}

///////////////////////////////////////////////////////////////////////////////
//
//

G4double G4ComponentGGNuclNuclXsc::
GetCoulombBarier(const G4DynamicParticle* aParticle, G4double tZ, G4double tA,
                 G4double pR, G4double tR)
{
  G4double ratio;
  G4double pZ = aParticle->GetDefinition()->GetPDGCharge();

  G4double pTkin = aParticle->GetKineticEnergy();
  // G4double pPlab = aParticle->GetTotalMomentum();
  G4double pM    = aParticle->GetDefinition()->GetPDGMass();
  // G4double tM    = tZ*proton_mass_c2 + (tA-tZ)*neutron_mass_c2; // ~ 1% accuracy
  G4double tM    = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass( G4int(tZ), G4int(tA) );
  G4double pElab = pTkin + pM;
  G4double totEcm  = std::sqrt(pM*pM + tM*tM + 2.*pElab*tM);
  // G4double pPcm  = pPlab*tM/totEcm;
  // G4double pTcm  = std::sqrt(pM*pM + pPcm*pPcm) - pM;
  G4double totTcm  = totEcm - pM -tM;

  G4double bC    = fine_structure_const*hbarc*pZ*tZ;
           bC   /= pR + tR;
           bC   /= 2.;  // 4., 2. parametrisation cof ??? vmg

	   // G4cout<<"pTkin = "<<pTkin/GeV<<"; pPlab = "
	   // <<pPlab/GeV<<"; bC = "<<bC/GeV<<"; pTcm = "<<pTcm/GeV<<G4endl;

  if( totTcm <= bC ) ratio = 0.;
  else             ratio = 1. - bC/totTcm;

  // if(ratio < DBL_MIN) ratio = DBL_MIN;
  if( ratio < 0.) ratio = 0.;

  // G4cout <<"ratio = "<<ratio<<G4endl;
  return ratio;
}


//////////////////////////////////////////////////////////////////////////
//
// Return single-diffraction/inelastic cross-section ratio

G4double G4ComponentGGNuclNuclXsc::
GetRatioSD(const G4DynamicParticle* aParticle, G4double tA, G4double tZ)
{
  G4double sigma, cofInelastic = 2.4, cofTotal = 2.0, nucleusSquare, ratio;

  G4double pZ = aParticle->GetDefinition()->GetPDGCharge();
  G4double pA = aParticle->GetDefinition()->GetBaryonNumber();

  G4double pTkin = aParticle->GetKineticEnergy();  
  pTkin /= pA;

  G4double pN = pA - pZ;
  if( pN < 0. ) pN = 0.;

  G4double tN = tA - tZ;
  if( tN < 0. ) tN = 0.;

  G4double tR = GetNucleusRadius(tZ,tA);  
  G4double pR = GetNucleusRadius(pZ,pA); 

  sigma = (pZ*tZ+pN*tN)*GetHadronNucleonXscNS(theProton, pTkin, theProton) +
          (pZ*tN+pN*tZ)*GetHadronNucleonXscNS(theProton, pTkin, theNeutron);

  nucleusSquare = cofTotal*pi*( pR*pR + tR*tR );   // basically 2piRR
  ratio = sigma/nucleusSquare;
  fInelasticXsc = nucleusSquare*G4Log(1. + cofInelastic*ratio)/cofInelastic;
  G4double difratio = ratio/(1.+ratio);

  fDiffractionXsc = 0.5*nucleusSquare*( difratio - G4Log( 1. + difratio ) );

  if (fInelasticXsc > 0.) ratio = fDiffractionXsc/fInelasticXsc;
  else                    ratio = 0.;

  return ratio; 
}

//////////////////////////////////////////////////////////////////////////
//
// Return quasi-elastic/inelastic cross-section ratio

G4double G4ComponentGGNuclNuclXsc::
GetRatioQE(const G4DynamicParticle* aParticle, G4double tA, G4double tZ)
{
  G4double sigma, cofInelastic = 2.4, cofTotal = 2.0, nucleusSquare, ratio;

  G4double pZ = aParticle->GetDefinition()->GetPDGCharge();
  G4double pA = aParticle->GetDefinition()->GetBaryonNumber();

  G4double pTkin = aParticle->GetKineticEnergy();  
  pTkin /= pA;

  G4double pN = pA - pZ;
  if( pN < 0. ) pN = 0.;

  G4double tN = tA - tZ;
  if( tN < 0. ) tN = 0.;

  G4double tR = GetNucleusRadius(tZ,tA);  
  G4double pR = GetNucleusRadius(pZ,pA); 

  sigma = (pZ*tZ+pN*tN)*GetHadronNucleonXscNS(theProton, pTkin, theProton) +
          (pZ*tN+pN*tZ)*GetHadronNucleonXscNS(theProton, pTkin, theNeutron);

  nucleusSquare = cofTotal*pi*( pR*pR + tR*tR );   // basically 2piRR
  ratio = sigma/nucleusSquare;
  fInelasticXsc = nucleusSquare*G4Log(1. + cofInelastic*ratio)/cofInelastic;

  //  sigma = GetHNinelasticXsc(aParticle, tA, tZ);
  ratio = sigma/nucleusSquare;
  fProductionXsc = nucleusSquare*G4Log(1. + cofInelastic*ratio)/cofInelastic;

  if (fInelasticXsc > fProductionXsc) ratio = (fInelasticXsc-fProductionXsc)/fInelasticXsc;
  else                                ratio = 0.;
  if ( ratio < 0. )                   ratio = 0.;

  return ratio; 
}

///////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to differnt parametrisations:
// [2] E. Levin, hep-ph/9710546
// [3] U. Dersch, et al, hep-ex/9910052
// [4] M.J. Longo, et al, Phys.Rev.Lett. 33 (1974) 725 

G4double 
G4ComponentGGNuclNuclXsc::GetHadronNucleonXsc(const G4DynamicParticle* aParticle, 
                                              const G4Element* anElement)
{
  G4int At = G4lrint(anElement->GetN());  // number of nucleons 
  G4int Zt = G4lrint(anElement->GetZ());  // number of protons
  return GetHadronNucleonXsc(aParticle, At, Zt);
}




///////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to differnt parametrisations:
// [2] E. Levin, hep-ph/9710546
// [3] U. Dersch, et al, hep-ex/9910052
// [4] M.J. Longo, et al, Phys.Rev.Lett. 33 (1974) 725 

G4double 
G4ComponentGGNuclNuclXsc::GetHadronNucleonXsc(const G4DynamicParticle* aParticle, 
                                                   G4int At, G4int Zt)
{
  G4double xsection = 0.;

  G4double targ_mass = G4ParticleTable::GetParticleTable()->
  GetIonTable()->GetIonMass(Zt, At);
  targ_mass = 0.939*GeV;  // ~mean neutron and proton ???

  G4double proj_mass = aParticle->GetMass();
  G4double proj_momentum = aParticle->GetMomentum().mag();
  G4double sMand = CalcMandelstamS ( proj_mass , targ_mass , proj_momentum );

  sMand /= GeV*GeV;  // in GeV for parametrisation
  proj_momentum /= GeV;
  const G4ParticleDefinition* pParticle = aParticle->GetDefinition();

  if(pParticle == theNeutron) // as proton ??? 
  {
    xsection = G4double(At)*(21.70*G4Pow::GetInstance()->powA(sMand,0.0808) + 56.08*G4Pow::GetInstance()->powA(sMand,-0.4525));
  } 
  else if(pParticle == theProton) 
  {
    xsection = G4double(At)*(21.70*G4Pow::GetInstance()->powA(sMand,0.0808) + 56.08*G4Pow::GetInstance()->powA(sMand,-0.4525));
  } 
 
  xsection *= millibarn;
  return xsection;
}



///////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to PDG parametrisation (2005):
// http://pdg.lbl.gov/2006/reviews/hadronicrpp.pdf
//  At = number of nucleons,  Zt = number of protons 

G4double 
G4ComponentGGNuclNuclXsc::GetHadronNucleonXscPDG(const G4ParticleDefinition* pParticle, 
                                                 G4double sMand, 
                                                 const G4ParticleDefinition* tParticle)
{
  G4double xsection = 0.;
  // G4bool pORn = (tParticle == theProton || nucleon == theNeutron  );  
  G4bool proton = (tParticle == theProton);
  G4bool neutron = (tParticle == theNeutron);

  // General PDG fit constants

  G4double s0   = 5.38*5.38; // in Gev^2
  G4double eta1 = 0.458;
  G4double eta2 = 0.458;
  G4double B    = 0.308;

  // const G4ParticleDefinition* pParticle = aParticle->GetDefinition();
  
  if(pParticle == theNeutron) // proton-neutron fit 
  {
    if ( proton )
    {
      xsection = ( 35.80 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                  + 40.15*G4Pow::GetInstance()->powA(sMand,-eta1) - 30.*G4Pow::GetInstance()->powA(sMand,-eta2));
    }
    if ( neutron )
    {
      xsection = (35.45 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
	     + 42.53*G4Pow::GetInstance()->powA(sMand,-eta1) - 33.34*G4Pow::GetInstance()->powA(sMand,-eta2)); // pp for nn
    }
  } 
  else if(pParticle == theProton) 
  {
    if ( proton )
    {
      xsection = (35.45 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                + 42.53*G4Pow::GetInstance()->powA(sMand,-eta1) - 33.34*G4Pow::GetInstance()->powA(sMand,-eta2));

    }
    if ( neutron )
    {
      xsection = (35.80 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                  + 40.15*G4Pow::GetInstance()->powA(sMand,-eta1) - 30.*G4Pow::GetInstance()->powA(sMand,-eta2));
    }
  } 
  xsection *= millibarn; // parametrised in mb
  return xsection;
}


///////////////////////////////////////////////////////////////////////////////
//
// Returns nucleon-nucleon cross-section based on N. Starkov parametrisation of
// data from mainly http://wwwppds.ihep.su:8001/c5-6A.html database
// projectile nucleon is pParticle with pTkin shooting target nucleon tParticle

G4double 
G4ComponentGGNuclNuclXsc::GetHadronNucleonXscNS(const G4ParticleDefinition* pParticle, 
                                                 G4double pTkin, 
                                                 const G4ParticleDefinition* tParticle)
{
  G4double xsection(0);
  // G4double Delta;   DHW 19 May 2011: variable set but not used
  G4double A0, B0;
  G4double hpXscv(0);
  G4double hnXscv(0);

  G4double targ_mass = tParticle->GetPDGMass();
  G4double proj_mass = pParticle->GetPDGMass(); 

  G4double proj_energy   = proj_mass + pTkin; 
  G4double proj_momentum = std::sqrt(pTkin*(pTkin+2*proj_mass));

  G4double sMand = CalcMandelstamS ( proj_mass , targ_mass , proj_momentum );

  sMand         /= GeV*GeV;  // in GeV for parametrisation
  proj_momentum /= GeV;
  proj_energy   /= GeV;
  proj_mass     /= GeV;

  // General PDG fit constants

  //  G4double s0   = 5.38*5.38; // in Gev^2
  //  G4double eta1 = 0.458;
  //  G4double eta2 = 0.458;
  //  G4double B    = 0.308;
  
  if( proj_momentum >= 373.)
  {
      return GetHadronNucleonXscPDG(pParticle,sMand,tParticle);
  }
  else if( proj_momentum >= 10. ) // high energy: pp = nn = np
    // if( proj_momentum >= 2.)
  {
    //  Delta = 1.;    // DHW 19 May 2011: variable set but not used
    //  if (proj_energy < 40.) Delta = 0.916+0.0021*proj_energy;

    //AR-12Aug2016  if (proj_momentum >= 10.) {
      B0 = 7.5;
      A0 = 100. - B0*G4Log(3.0e7);

      xsection = A0 + B0*G4Log(proj_energy) - 11
                + 103*G4Pow::GetInstance()->powA(2*0.93827*proj_energy + proj_mass*proj_mass+
                  0.93827*0.93827,-0.165);        //  mb
    //AR-12Aug2016  }
  }
  else // low energy pp = nn != np
  {
      if(pParticle == tParticle) // pp or nn      // nn to be pp
      {
        if( proj_momentum < 0.73 )
        {
          hnXscv = 23 + 50*( G4Pow::GetInstance()->powA( G4Log(0.73/proj_momentum), 3.5 ) );
        }
        else if( proj_momentum < 1.05  )
        {
          hnXscv = 23 + 40*(G4Log(proj_momentum/0.73))*
                         (G4Log(proj_momentum/0.73));
        }
        else  // if( proj_momentum < 10.  )
        {
          hnXscv = 39.0 + 
              75*(proj_momentum - 1.2)/(G4Pow::GetInstance()->powA(proj_momentum,3.0) + 0.15);
        }
        xsection = hnXscv;
      }
      else  // pn to be np
      {
        if( proj_momentum < 0.8 )
        {
          hpXscv = 33+30*G4Pow::GetInstance()->powA(G4Log(proj_momentum/1.3),4.0);
        }      
        else if( proj_momentum < 1.4 )
        {
          hpXscv = 33+30*G4Pow::GetInstance()->powA(G4Log(proj_momentum/0.95),2.0);
        }
        else    // if( proj_momentum < 10.  )
        {
          hpXscv = 33.3+
              20.8*(G4Pow::GetInstance()->powA(proj_momentum,2.0)-1.35)/
                 (G4Pow::GetInstance()->powA(proj_momentum,2.50)+0.95);
        }
        xsection = hpXscv;
      }
  }
  xsection *= millibarn; // parametrised in mb
  return xsection;
}

/////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon inelastic cross-section based on FTF-parametrisation 

G4double 
G4ComponentGGNuclNuclXsc::GetHNinelasticXscVU(const G4DynamicParticle* aParticle, 
                                              G4int At, G4int Zt)
{
  G4int PDGcode = aParticle->GetDefinition()->GetPDGEncoding();
  G4int absPDGcode = std::abs(PDGcode);
  G4double Elab = aParticle->GetTotalEnergy();              
                          // (s - 2*0.88*GeV*GeV)/(2*0.939*GeV)/GeV;
  G4double Plab = aParticle->GetMomentum().mag();            
                          // std::sqrt(Elab * Elab - 0.88);

  Elab /= GeV;
  Plab /= GeV;

  G4double LogPlab    = G4Log( Plab );
  G4double sqrLogPlab = LogPlab * LogPlab;

  //G4cout<<"Plab = "<<Plab<<G4endl;

  G4double NumberOfTargetProtons  = Zt; 
  G4double NumberOfTargetNucleons = At;
  G4double NumberOfTargetNeutrons = NumberOfTargetNucleons - NumberOfTargetProtons;

  if(NumberOfTargetNeutrons < 0.) NumberOfTargetNeutrons = 0.;

  G4double Xtotal = 0., Xelastic = 0., Xinelastic =0.;

  if( absPDGcode > 1000 )  //------Projectile is baryon --------
  {
       G4double XtotPP = 48.0 +  0. *G4Pow::GetInstance()->powA(Plab, 0.  ) +
                         0.522*sqrLogPlab - 4.51*LogPlab;

       G4double XtotPN = 47.3 +  0. *G4Pow::GetInstance()->powA(Plab, 0.  ) +
                         0.513*sqrLogPlab - 4.27*LogPlab;

       G4double XelPP  = 11.9 + 26.9*G4Pow::GetInstance()->powA(Plab,-1.21) +
                         0.169*sqrLogPlab - 1.85*LogPlab;

       G4double XelPN  = 11.9 + 26.9*G4Pow::GetInstance()->powA(Plab,-1.21) +
                         0.169*sqrLogPlab - 1.85*LogPlab;

       Xtotal          = ( NumberOfTargetProtons  * XtotPP +
                           NumberOfTargetNeutrons * XtotPN  );

       Xelastic        = ( NumberOfTargetProtons  * XelPP  +
                           NumberOfTargetNeutrons * XelPN   );
  }

  Xinelastic = Xtotal - Xelastic;
  if(Xinelastic < 0.) Xinelastic = 0.;

  return Xinelastic*= millibarn;
}

///////////////////////////////////////////////////////////////////////////////
//
//

G4double 
G4ComponentGGNuclNuclXsc::GetNucleusRadius(const G4DynamicParticle* , 
                                           const G4Element* anElement)
{
  G4double At = anElement->GetN();
  G4double oneThird = 1.0/3.0;
  G4double cubicrAt = G4Pow::GetInstance()->powA (At, oneThird); 

  G4double R;  // = fRadiusConst*cubicrAt;
  R = fRadiusConst*cubicrAt;

  G4double meanA  = 21.;
  G4double tauA1  = 40.; 
  G4double tauA2  = 10.; 
  G4double tauA3  = 5.; 

  G4double a1 = 0.85;
  G4double b1 = 1. - a1;

  G4double b2 = 0.3;
  G4double b3 = 4.;

  if (At > 20.)   // 20.
  {
    R *= ( a1 + b1*G4Exp( -(At - meanA)/tauA1) ); 
  }
  else if (At > 3.5)
  {
    R *= ( 1.0 + b2*( 1. - G4Exp( (At - meanA)/tauA2) ) ); 
  }
  else 
  {
    R *= ( 1.0 + b3*( 1. - G4Exp( (At - meanA)/tauA3) ) ); 
  }

  return R;
}

///////////////////////////////////////////////////////////////////////////////
//
//

G4double 
G4ComponentGGNuclNuclXsc::GetNucleusRadius(G4double Zt, G4double At)
{
  G4double R;
  R = GetNucleusRadiusDE(Zt,At);
  // R = GetNucleusRadiusRMS(Zt,At);

  return R;
}

///////////////////////////////////////////////////////////////////

G4double 
G4ComponentGGNuclNuclXsc::GetNucleusRadiusGG(G4double At)
{
  G4double oneThird = 1.0/3.0;
  G4double cubicrAt = G4Pow::GetInstance()->powA (At, oneThird); 

  G4double R;  // = fRadiusConst*cubicrAt;  
  R = fRadiusConst*cubicrAt;

  G4double meanA = 20.;
  G4double tauA  = 20.;

  if ( At > 20.)   // 20.
  {
    R *= ( 0.8 + 0.2*G4Exp( -(At - meanA)/tauA) ); 
  }
  else
  {
    R *= ( 1.0 + 0.1*( 1. - G4Exp( (At - meanA)/tauA) ) ); 
  }

  return R;
}

/////////////////////////////////////////////////////////////////////////////
//
//

G4double 
G4ComponentGGNuclNuclXsc::GetNucleusRadiusDE(G4double Z, G4double A)
{
  // algorithm from diffuse-elastic

  G4double R, r0, a11, a12, a13, a2, a3;

  a11 = 1.26;  // 1.08, 1.16
  a12 = 1.;  // 1.08, 1.16
  a13 = 1.12;  // 1.08, 1.16
  a2 = 1.1;
  a3 = 1.;

  // Special rms radii for light nucleii

  if (A < 50.)
  {
    if     (std::abs(A-1.) < 0.5)                         return 0.89*fermi; // p
    else if(std::abs(A-2.) < 0.5)                         return 2.13*fermi; // d
    else if(std::abs(Z-1.) < 0.5 && std::abs(A-3.) < 0.5) return 1.80*fermi; // t

    else if(std::abs(Z-2.) < 0.5 && std::abs(A-3.) < 0.5) return 1.96*fermi; // He3
    else if(std::abs(Z-2.) < 0.5 && std::abs(A-4.) < 0.5) return 1.68*fermi; // He4

    else if(std::abs(Z-3.) < 0.5)                         return 2.40*fermi; // Li7
    else if(std::abs(Z-4.) < 0.5)                         return 2.51*fermi; // Be9

    else if( 10.  < A && A <= 16. ) r0  = a11*( 1 - G4Pow::GetInstance()->powA(A, -2./3.) )*fermi;   // 1.08*fermi;
    else if( 15.  < A && A <= 20. ) r0  = a12*( 1 - G4Pow::GetInstance()->powA(A, -2./3.) )*fermi;
    else if( 20.  < A && A <= 30. ) r0  = a13*( 1 - G4Pow::GetInstance()->powA(A, -2./3.) )*fermi;
    else                            r0  = a2*fermi;

    R = r0*G4Pow::GetInstance()->powA( A, 1./3. );
  }
  else
  {
    r0 = a3*fermi;

    R  = r0*G4Pow::GetInstance()->powA(A, 0.27);
  }
  return R;
}


/////////////////////////////////////////////////////////////////////////////
//
// RMS radii from e-A scattering data

G4double 
G4ComponentGGNuclNuclXsc::GetNucleusRadiusRMS(G4double Z, G4double A)
{

  if     (std::abs(A-1.) < 0.5)                         return 0.89*fermi; // p
  else if(std::abs(A-2.) < 0.5)                         return 2.13*fermi; // d
  else if(std::abs(Z-1.) < 0.5 && std::abs(A-3.) < 0.5) return 1.80*fermi; // t

  else if(std::abs(Z-2.) < 0.5 && std::abs(A-3.) < 0.5) return 1.96*fermi; // He3
  else if(std::abs(Z-2.) < 0.5 && std::abs(A-4.) < 0.5) return 1.68*fermi; // He4

  else if(std::abs(Z-3.) < 0.5)                         return 2.40*fermi; // Li7
  else if(std::abs(Z-4.) < 0.5)                         return 2.51*fermi; // Be9

  else                               return 1.24*G4Pow::GetInstance()->powA(A, 0.28 )*fermi; // A > 9
}


///////////////////////////////////////////////////////////////////////////////
//
//

G4double G4ComponentGGNuclNuclXsc::CalculateEcmValue(const G4double mp, 
                                                     const G4double mt, 
                                                     const G4double Plab)
{
  G4double Elab = std::sqrt ( mp * mp + Plab * Plab );
  G4double Ecm  = std::sqrt ( mp * mp + mt * mt + 2 * Elab * mt );
  // G4double Pcm  = Plab * mt / Ecm;
  // G4double KEcm = std::sqrt ( Pcm * Pcm + mp * mp ) - mp;

  return Ecm ; // KEcm;
}


///////////////////////////////////////////////////////////////////////////////
//
//

G4double G4ComponentGGNuclNuclXsc::CalcMandelstamS(const G4double mp, 
                                                   const G4double mt, 
                                                   const G4double Plab)
{
  G4double Elab = std::sqrt ( mp * mp + Plab * Plab );
  G4double sMand  = mp*mp + mt*mt + 2*Elab*mt ;

  return sMand;
}

//
//
///////////////////////////////////////////////////////////////////////////////
