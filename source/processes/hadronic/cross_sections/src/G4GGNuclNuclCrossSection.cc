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

#include "G4GGNuclNuclCrossSection.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadTmpUtil.hh"


G4GGNuclNuclCrossSection::G4GGNuclNuclCrossSection() 
 : G4VCrossSectionDataSet("Glauber-Gribov nucleus nucleus"),
   fUpperLimit(100000*GeV), fLowerLimit(0.1*MeV),
   fRadiusConst(1.08*fermi),  // 1.1, 1.3 ?
   fTotalXsc(0.0), fElasticXsc(0.0), fInelasticXsc(0.0), fProductionXsc(0.0),
   fDiffractionXsc(0.0), fHadronNucleonXsc(0.0)
{
  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
}


G4GGNuclNuclCrossSection::~G4GGNuclNuclCrossSection()
{}


G4bool 
G4GGNuclNuclCrossSection::IsApplicable(const G4DynamicParticle* aDP, 
				       const G4Element*  anElement)
{
  G4int Z = G4lrint(anElement->GetZ());
  G4int N = G4lrint(anElement->GetN());
  return IsIsoApplicable(aDP, Z, N);
} 

///////////////////////////////////////////////////////////////////////////////
//
//

G4bool 
G4GGNuclNuclCrossSection::IsIsoApplicable(const G4DynamicParticle* aDP, 
					  G4int Z, G4int)
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


G4double G4GGNuclNuclCrossSection::
GetCrossSection(const G4DynamicParticle* aParticle, const G4Element* anElement,
                G4double T)
{
  G4int Z = G4lrint(anElement->GetZ());
  G4int N = G4lrint(anElement->GetN());
  return GetZandACrossSection(aParticle, Z, N, T);
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculates total and inelastic Xsc, derives elastic as total - inelastic 
// accordong to Glauber model with Gribov correction calculated in the dipole 
// approximation on light cone. Gaussian density of point-like nucleons helps 
// to calculate rest integrals of the model. [1] B.Z. Kopeliovich,
// nucl-th/0306044 + simplification above


G4double G4GGNuclNuclCrossSection::
GetZandACrossSection(const G4DynamicParticle* aParticle,
                     G4int tZ, G4int tA, G4double)
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

  G4double tR = GetNucleusRadius(tA);  
  G4double pR = GetNucleusRadius(pA); 

  cB = GetCoulombBarier(aParticle, G4double(tZ), G4double(tA), pR, tR);
  if (cB > 0.) {

    sigma = (pZ*tZ+pN*tN)*GetHadronNucleonXscNS(theProton, pTkin, theProton) +
          (pZ*tN+pN*tZ)*GetHadronNucleonXscNS(theProton, pTkin, theNeutron);

    nucleusSquare = cofTotal*pi*( pR*pR + tR*tR );   // basically 2piRR

    ratio = sigma/nucleusSquare;
    xsection =  nucleusSquare*std::log( 1. + ratio );
    fTotalXsc = xsection;
    fTotalXsc *= cB;

    fInelasticXsc = nucleusSquare*std::log( 1. + cofInelastic*ratio )/cofInelastic;

    fInelasticXsc *= cB;
    fElasticXsc   = fTotalXsc - fInelasticXsc;

    // if (fElasticXsc < DBL_MIN) fElasticXsc = DBL_MIN;
    /*    
    G4double difratio = ratio/(1.+ratio);

    fDiffractionXsc = 0.5*nucleusSquare*( difratio - std::log( 1. + difratio ) );
    */
    // production to be checked !!! edit MK xsc

    //sigma = (pZ*tZ+pN*tN)*GetHadronNucleonXscMK(theProton, pTkin, theProton) +
    //      (pZ*tN+pN*tZ)*GetHadronNucleonXscMK(theProton, pTkin, theNeutron);

    sigma = (pZ*tZ+pN*tN)*GetHadronNucleonXscNS(theProton, pTkin, theProton) +
          (pZ*tN+pN*tZ)*GetHadronNucleonXscNS(theProton, pTkin, theNeutron);
 
    ratio = sigma/nucleusSquare;
    fProductionXsc = nucleusSquare*std::log( 1. + cofInelastic*ratio )/cofInelastic;

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

G4double G4GGNuclNuclCrossSection::
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

G4double G4GGNuclNuclCrossSection::
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

  G4double tR = GetNucleusRadius(tA);  
  G4double pR = GetNucleusRadius(pA); 

  sigma = (pZ*tZ+pN*tN)*GetHadronNucleonXscNS(theProton, pTkin, theProton) +
          (pZ*tN+pN*tZ)*GetHadronNucleonXscNS(theProton, pTkin, theNeutron);

  nucleusSquare = cofTotal*pi*( pR*pR + tR*tR );   // basically 2piRR
  ratio = sigma/nucleusSquare;
  fInelasticXsc = nucleusSquare*std::log(1. + cofInelastic*ratio)/cofInelastic;
  G4double difratio = ratio/(1.+ratio);

  fDiffractionXsc = 0.5*nucleusSquare*( difratio - std::log( 1. + difratio ) );

  if (fInelasticXsc > 0.) ratio = fDiffractionXsc/fInelasticXsc;
  else                    ratio = 0.;

  return ratio; 
}

//////////////////////////////////////////////////////////////////////////
//
// Return quasi-elastic/inelastic cross-section ratio

G4double G4GGNuclNuclCrossSection::
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

  G4double tR = GetNucleusRadius(tA);  
  G4double pR = GetNucleusRadius(pA); 

  sigma = (pZ*tZ+pN*tN)*GetHadronNucleonXscNS(theProton, pTkin, theProton) +
          (pZ*tN+pN*tZ)*GetHadronNucleonXscNS(theProton, pTkin, theNeutron);

  nucleusSquare = cofTotal*pi*( pR*pR + tR*tR );   // basically 2piRR
  ratio = sigma/nucleusSquare;
  fInelasticXsc = nucleusSquare*std::log(1. + cofInelastic*ratio)/cofInelastic;

  //  sigma = GetHNinelasticXsc(aParticle, tA, tZ);
  ratio = sigma/nucleusSquare;
  fProductionXsc = nucleusSquare*std::log(1. + cofInelastic*ratio)/cofInelastic;

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
G4GGNuclNuclCrossSection::GetHadronNucleonXsc(const G4DynamicParticle* aParticle, 
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
G4GGNuclNuclCrossSection::GetHadronNucleonXsc(const G4DynamicParticle* aParticle, 
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
    xsection = G4double(At)*(21.70*std::pow(sMand,0.0808) + 56.08*std::pow(sMand,-0.4525));
  } 
  else if(pParticle == theProton) 
  {
    xsection = G4double(At)*(21.70*std::pow(sMand,0.0808) + 56.08*std::pow(sMand,-0.4525));
  } 
 
  xsection *= millibarn;
  return xsection;
}


///////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to PDG parametrisation (2005):
// http://pdg.lbl.gov/2006/reviews/hadronicrpp.pdf

G4double 
G4GGNuclNuclCrossSection::GetHadronNucleonXscPDG(const G4DynamicParticle* aParticle, 
                                                  const G4Element* anElement)
{
  G4int At = G4lrint(anElement->GetN());  // number of nucleons 
  G4int Zt = G4lrint(anElement->GetZ());  // number of protons
  return GetHadronNucleonXscPDG( aParticle, At, Zt );
}


///////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to PDG parametrisation (2005):
// http://pdg.lbl.gov/2006/reviews/hadronicrpp.pdf
//  At = number of nucleons,  Zt = number of protons 

G4double 
G4GGNuclNuclCrossSection::GetHadronNucleonXscPDG(const G4DynamicParticle* aParticle, 
                                                 G4int At, G4int Zt)
{
  G4double xsection = 0.;

  G4double Nt = At-Zt;              // number of neutrons
  if (Nt < 0.) Nt = 0.;  

  G4double targ_mass = G4ParticleTable::GetParticleTable()->
  GetIonTable()->GetIonMass(Zt, At);
  targ_mass = 0.939*GeV;  // ~mean neutron and proton ???

  G4double proj_mass     = aParticle->GetMass(); 
  G4double proj_momentum = aParticle->GetMomentum().mag();
  G4double sMand = CalcMandelstamS ( proj_mass , targ_mass , proj_momentum );
  sMand         /= GeV*GeV;  // in GeV for parametrisation

  // General PDG fit constants

  G4double s0   = 5.38*5.38; // in Gev^2
  G4double eta1 = 0.458;
  G4double eta2 = 0.458;
  G4double B    = 0.308;

  const G4ParticleDefinition* pParticle = aParticle->GetDefinition();
  
  if(pParticle == theNeutron) // proton-neutron fit 
  {
    xsection = G4double(Zt)*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                  + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2));
    xsection += Nt*(35.45 + B*std::pow(std::log(sMand/s0),2.) 
	     + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2)); // pp for nn
  } 
  else if(pParticle == theProton) 
  {
    xsection = G4double(Zt)*(35.45 + B*std::pow(std::log(sMand/s0),2.) 
                + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2));

    xsection += Nt*(35.80 + B*std::pow(std::log(sMand/s0),2.) 
                  + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2));
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
G4GGNuclNuclCrossSection::GetHadronNucleonXscNS(G4ParticleDefinition* pParticle, 
                                                 G4double pTkin, 
                                                 G4ParticleDefinition* tParticle)
{
  G4double xsection(0), Delta, A0, B0;
  G4double hpXsc(0);
  G4double hnXsc(0);

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
  
  if( proj_momentum >= 10. ) // high energy: pp = nn = np
    // if( proj_momentum >= 2.)
  {
    Delta = 1.;
    if (proj_energy < 40.) Delta = 0.916+0.0021*proj_energy;

    if (proj_momentum >= 10.) {
      B0 = 7.5;
      A0 = 100. - B0*std::log(3.0e7);

      xsection = A0 + B0*std::log(proj_energy) - 11
                + 103*std::pow(2*0.93827*proj_energy + proj_mass*proj_mass+
                  0.93827*0.93827,-0.165);        //  mb
    }
  }
  else // low energy pp = nn != np
  {
      if(pParticle == tParticle) // pp or nn      // nn to be pp
      {
        if( proj_momentum < 0.73 )
        {
          hnXsc = 23 + 50*( std::pow( std::log(0.73/proj_momentum), 3.5 ) );
        }
        else if( proj_momentum < 1.05  )
        {
          hnXsc = 23 + 40*(std::log(proj_momentum/0.73))*
                         (std::log(proj_momentum/0.73));
        }
        else  // if( proj_momentum < 10.  )
        {
          hnXsc = 39.0 + 
              75*(proj_momentum - 1.2)/(std::pow(proj_momentum,3.0) + 0.15);
        }
        xsection = hnXsc;
      }
      else  // pn to be np
      {
        if( proj_momentum < 0.8 )
        {
          hpXsc = 33+30*std::pow(std::log(proj_momentum/1.3),4.0);
        }      
        else if( proj_momentum < 1.4 )
        {
          hpXsc = 33+30*std::pow(std::log(proj_momentum/0.95),2.0);
        }
        else    // if( proj_momentum < 10.  )
        {
          hpXsc = 33.3+
              20.8*(std::pow(proj_momentum,2.0)-1.35)/
                 (std::pow(proj_momentum,2.50)+0.95);
        }
        xsection = hpXsc;
      }
  }
  xsection *= millibarn; // parametrised in mb
  return xsection;
}

/////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon inelastic cross-section based on FTF-parametrisation 

G4double 
G4GGNuclNuclCrossSection::GetHNinelasticXscVU(const G4DynamicParticle* aParticle, 
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

  G4double LogPlab    = std::log( Plab );
  G4double sqrLogPlab = LogPlab * LogPlab;

  //G4cout<<"Plab = "<<Plab<<G4endl;

  G4double NumberOfTargetProtons  = Zt; 
  G4double NumberOfTargetNucleons = At;
  G4double NumberOfTargetNeutrons = NumberOfTargetNucleons - NumberOfTargetProtons;

  if(NumberOfTargetNeutrons < 0.) NumberOfTargetNeutrons = 0.;

  G4double Xtotal = 0., Xelastic = 0., Xinelastic =0.;

  if( absPDGcode > 1000 )  //------Projectile is baryon --------
  {
       G4double XtotPP = 48.0 +  0. *std::pow(Plab, 0.  ) +
                         0.522*sqrLogPlab - 4.51*LogPlab;

       G4double XtotPN = 47.3 +  0. *std::pow(Plab, 0.  ) +
                         0.513*sqrLogPlab - 4.27*LogPlab;

       G4double XelPP  = 11.9 + 26.9*std::pow(Plab,-1.21) +
                         0.169*sqrLogPlab - 1.85*LogPlab;

       G4double XelPN  = 11.9 + 26.9*std::pow(Plab,-1.21) +
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
G4GGNuclNuclCrossSection::GetNucleusRadius(const G4DynamicParticle* , 
                                           const G4Element* anElement)
{
  G4double At = anElement->GetN();
  G4double oneThird = 1.0/3.0;
  G4double cubicrAt = std::pow (At, oneThird); 

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
    R *= ( a1 + b1*std::exp( -(At - meanA)/tauA1) ); 
  }
  else if (At > 3.5)
  {
    R *= ( 1.0 + b2*( 1. - std::exp( (At - meanA)/tauA2) ) ); 
  }
  else 
  {
    R *= ( 1.0 + b3*( 1. - std::exp( (At - meanA)/tauA3) ) ); 
  }

  return R;
}

///////////////////////////////////////////////////////////////////////////////
//
//

G4double 
G4GGNuclNuclCrossSection::GetNucleusRadius(G4double At)
{
  G4double R;
  R = GetNucleusRadiusDE(At);

  return R;
}

///////////////////////////////////////////////////////////////////

G4double 
G4GGNuclNuclCrossSection::GetNucleusRadiusGG(G4double At)
{
  G4double oneThird = 1.0/3.0;
  G4double cubicrAt = std::pow (At, oneThird); 

  G4double R;  // = fRadiusConst*cubicrAt;  
  R = fRadiusConst*cubicrAt;

  G4double meanA = 20.;
  G4double tauA  = 20.;

  if ( At > 20.)   // 20.
  {
    R *= ( 0.8 + 0.2*std::exp( -(At - meanA)/tauA) ); 
  }
  else
  {
    R *= ( 1.0 + 0.1*( 1. - std::exp( (At - meanA)/tauA) ) ); 
  }

  return R;
}


G4double 
G4GGNuclNuclCrossSection::GetNucleusRadiusDE(G4double A)
{
  // algorithm from diffuse-elastic

  G4double R, r0, a11, a12, a13, a2, a3;

  a11 = 1.26;  // 1.08, 1.16
  a12 = 1.;  // 1.08, 1.16
  a13 = 1.12;  // 1.08, 1.16
  a2 = 1.1;
  a3 = 1.;


  if (A < 50.)
  {
    if( 10  < A && A <= 15. )     r0  = a11*( 1 - std::pow(A, -2./3.) )*fermi;   // 1.08*fermi;
    else if( 15  < A && A <= 20 ) r0  = a12*( 1 - std::pow(A, -2./3.) )*fermi;
    else if( 20  < A && A <= 30 ) r0  = a13*( 1 - std::pow(A, -2./3.) )*fermi;
    else                          r0  = a2*fermi;

    R = r0*std::pow( A, 1./3. );
  }
  else
  {
    r0 = a3*fermi;

    R  = r0*std::pow(A, 0.27);
  }
  return R;
}


///////////////////////////////////////////////////////////////////////////////
//
//

G4double G4GGNuclNuclCrossSection::CalculateEcmValue(const G4double mp, 
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

G4double G4GGNuclNuclCrossSection::CalcMandelstamS(const G4double mp, 
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
