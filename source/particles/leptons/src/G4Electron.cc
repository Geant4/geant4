// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Electron.cc,v 1.2 1999-12-15 14:51:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo
// **********************************************************************
//  Added particle definitions, H.Kurashige, 19 April 1996
//  Added SetCuts implementation, L.Urban, 30 May 1996
//  Revised, G.Cosmo, 6 June 1996
//  Code uses operators (+=, *=, ++, -> etc.) correctly, P. Urban, 26/6/96
//  Add ElectronDefinition() H.Kurashige 4 July 1996
// ----------------------------------------------------------------------

#include "g4std/fstream"
#include "g4std/iomanip"
    
#include "G4Electron.hh"
// ######################################################################
// ###                         ELECTRON                               ###
// ######################################################################

G4Electron::G4Electron(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable )
 : G4VLepton( aName,mass,width,charge,iSpin,iParity,
	      iConjugation,iIsospin,iIsospin3,gParity,pType,
              lepton,baryon,encoding,stable,lifetime,decaytable )
{
}

// ......................................................................
// ...                 static member definitions                      ...
// ......................................................................
//     
//    Arguments for constructor are as follows
//               name             mass          width         charge
//             2*spin           parity  C-conjugation
//          2*Isospin       2*Isospin3       G-parity
//               type    lepton number  baryon number   PDG encoding
//             stable         lifetime    decay table 

G4Electron G4Electron::theElectron(
		 "e-",  0.51099906*MeV,       0.0*MeV,    -1.*eplus, 
		    1,               0,             0,          
		    0,               0,             0,             
	     "lepton",               1,             0,          11,
		 true,            -1.0,          NULL
);

G4Electron* G4Electron::ElectronDefinition(){return &theElectron;}
// initialization for static cut values
G4double   G4Electron::theElectronLengthCut = -1.0;
G4double*  G4Electron::theElectronKineticEnergyCuts = NULL;

// **********************************************************************
// ************************* ComputeLoss ********************************
// **********************************************************************
G4double G4Electron::ComputeLoss(G4double AtomicNumber,
                                 G4double KineticEnergy) const
{
  static G4double Z;  
  static G4double taul, ionpot, ionpotlog;
  const  G4double cbr1=0.02, cbr2=-5.7e-5, cbr3=1., cbr4=0.072;
  const  G4double Tlow=10.*keV, Thigh=1.*GeV;

  static G4double bremfactor= 0.1 ;

  //  calculate dE/dx for electrons
  if( abs(AtomicNumber-Z)>0.1 )
  {
    Z = AtomicNumber;
    taul = Tlow/GetPDGMass();
    ionpot = 1.6e-5*MeV*exp(0.9*log(Z))/GetPDGMass();
    ionpotlog = log(ionpot);
  } 


  G4double tau = KineticEnergy/GetPDGMass();
  G4double dEdx;

  if(tau<taul) {
    G4double t1 = taul+1.;
    G4double t2 = taul+2.;
    G4double tsq = taul*taul;
    G4double beta2 = taul*t2/(t1*t1);
    G4double f = 1.-beta2+log(tsq/2.)
                  +(0.5+0.25*tsq+(1.+2.*taul)*log(0.5))/(t1*t1);
    dEdx = (log(2.*taul+4.)-2.*ionpotlog+f)/beta2;
    dEdx = twopi_mc2_rcl2*Z*dEdx;
    G4double clow = dEdx*sqrt(taul);
    dEdx = clow/sqrt(KineticEnergy/GetPDGMass());
  } else {
    G4double t1 = tau+1.;
    G4double t2 = tau+2.;
    G4double tsq = tau*tau;
    G4double beta2 = tau*t2/(t1*t1);
    G4double f = 1.-beta2+log(tsq/2.)
                   +(0.5+0.25*tsq+(1.+2.*tau)*log(0.5))/(t1*t1);
    dEdx = (log(2.*tau+4.)-2.*ionpotlog+f)/beta2;
    dEdx = twopi_mc2_rcl2*Z*dEdx;

    // loss from bremsstrahlung follows
    G4double cbrem = (cbr1+cbr2*Z)
                       *(cbr3+cbr4*log(KineticEnergy/Thigh));
    cbrem = Z*(Z+1.)*cbrem*tau/beta2;

    cbrem *= bremfactor ;

    dEdx += twopi_mc2_rcl2*cbrem;
  }

  return dEdx;
}

// **********************************************************************
// *********************** BuildRangeVector *****************************
// **********************************************************************

void G4Electron::BuildRangeVector(const G4Material* aMaterial,
				  const G4LossTable* aLossTable,
				  G4double       maxEnergy,
				  G4double       aMass,
                                  G4PhysicsLogVector* rangeVector)
{
  //  create range vector for a material
  const G4double tlim = 10.*keV;
  const G4int maxnbint = 100;

  const G4ElementVector* elementVector = aMaterial->GetElementVector();
  const G4double* atomicNumDensityVector = aMaterial->GetAtomicNumDensityVector();
  G4int NumEl = aMaterial->GetNumberOfElements();

  // calculate parameters of the low energy part first
  G4int i;
  G4double loss=0.;
  for (i=0; i<NumEl; i++)
  {
    G4bool isOut;
    G4int IndEl = (*elementVector)(i)->GetIndex();
    loss += atomicNumDensityVector[i]*
           (*aLossTable)[IndEl]->GetValue(tlim,isOut);
  }
  G4double taulim = tlim/aMass;
  G4double clim = sqrt(taulim)*loss;
  G4double taumax = maxEnergy/aMass;

  // now the range vector can be filled

  for ( i=0; i<TotBin; i++)
  {
    G4double LowEdgeEnergy = rangeVector->GetLowEdgeEnergy(i);
    G4double tau = LowEdgeEnergy/aMass;

    if ( tau <= taulim ) {
      G4double Value = 2.*aMass*tau*sqrt(tau)/(3.*clim);
      rangeVector->PutValue(i,Value);
    } else {
      G4double rangelim = 2.*aMass*taulim*sqrt(taulim)/(3.*clim);
      G4double ltaulow = log(taulim);
      G4double ltauhigh = log(tau);
      G4double ltaumax = log(taumax);
      G4int    nbin = G4int(maxnbint*(ltauhigh-ltaulow)/(ltaumax-ltaulow));
      if( nbin < 1 ) nbin = 1;
      G4double Value = RangeLogSimpson(elementVector, atomicNumDensityVector,
                                       aLossTable,    aMass,
				       ltaulow,       ltauhigh,
                                       nbin,          NumEl)    + rangelim;
      rangeVector->PutValue(i,Value);
    }
  }
} 






