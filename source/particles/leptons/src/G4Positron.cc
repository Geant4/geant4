// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Positron.cc,v 1.3 2000-02-27 06:23:41 kurasige Exp $
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
//  Added SetCuts implementation, L.Urban, 12 June 1996
//  Code uses operators (+=, *=, ++, -> etc.) correctly, P. Urban, 26/6/96
//  Add PositronDefinition(), H.Kurashige 4 July 1996
// ----------------------------------------------------------------------

#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4Positron.hh"
    
// ######################################################################
// ###                         POSITRON                               ###
// ######################################################################

G4Positron::G4Positron(
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
  SetParticleSubType("e");
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
G4Positron G4Positron::thePositron(
		 "e+",  0.51099906*MeV,       0.0*MeV,    +1.*eplus, 
		    1,               0,             0,          
		    0,               0,             0,             
	     "lepton",              -1,             0,          -11,
		 true,            -1.0,          NULL
);

G4Positron* G4Positron::PositronDefinition() {return &thePositron;}
// initialization for static cut values
G4double   G4Positron::thePositronLengthCut = -1.0;
G4double*  G4Positron::thePositronKineticEnergyCuts = NULL;

// **********************************************************************
// ************************* ComputeLoss ********************************
// **********************************************************************

G4double G4Positron::ComputeLoss(G4double AtomicNumber,
                                 G4double KineticEnergy) const
{
  static G4double Z;  
  static G4double taul, ionpot, ionpotlog;
  const  G4double cbr1=0.02, cbr2=-5.7e-5, cbr3=1., cbr4=0.072;
  const  G4double Tlow=10.*keV, Thigh=1.*GeV;
  static G4double bremfactor = 0.1 ;

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

  if(tau<taul)
  {
    G4double t1 = taul+1.;
    G4double t2 = taul+2.;
    G4double tsq = taul*taul;
    G4double beta2 = taul*t2/(t1*t1);
    G4double     f = 2.*log(taul)
                     -(6.*taul+1.5*tsq-taul*(1.-tsq/3.)/t2-tsq*(0.5-tsq/12.)/
                       (t2*t2))/(t1*t1);
    dEdx = (log(2.*taul+4.)-2.*ionpotlog+f)/beta2;
    dEdx = twopi_mc2_rcl2*Z*dEdx;
    G4double clow = dEdx*sqrt(taul);
    dEdx = clow/sqrt(KineticEnergy/GetPDGMass());
  } else {
    G4double t1 = tau+1.;
    G4double t2 = tau+2.;
    G4double tsq = tau*tau;
    G4double beta2 = tau*t2/(t1*t1);
    G4double f = 2.*log(tau)
                 - (6.*tau+1.5*tsq-tau*(1.-tsq/3.)/t2-tsq*(0.5-tsq/12.)/
                     (t2*t2))/(t1*t1);
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
void G4Positron::BuildRangeVector(const G4Material* aMaterial,
				  const G4LossTable* aLossTable,
				  G4double       maxEnergy,
				  G4double       aMass,
                                  G4PhysicsLogVector* rangeVector)
{
  G4Electron* pElectron = G4Electron::ElectronDefinition();
  pElectron->BuildRangeVector(aMaterial,
			      aLossTable,
			      maxEnergy,
			      aMass,
                              rangeVector);
}
