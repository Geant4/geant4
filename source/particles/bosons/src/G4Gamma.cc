// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Gamma.cc,v 1.1 1999-01-07 16:09:48 gunter Exp $
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
//  Added SetCuts implementation , L.Urban , 12 June 1996
//  Code uses operators (+=, *=, ++, -> etc.) correctly, P. Urban, 26/6/96
//  Add GammaDefinition(), H.Kurashige, 4 July 1996

// ----------------------------------------------------------------------

#include "G4ios.hh"
#include <fstream.h>
#include <iomanip.h>
 
#include "G4Gamma.hh"

// ######################################################################
// ###                            GAMMA                               ###
// ######################################################################

G4Gamma::G4Gamma(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable )
 : G4VBoson( aName,mass,width,charge,iSpin,iParity,
             iConjugation,iIsospin,iIsospin3,gParity,pType,
             lepton,baryon,encoding,stable,lifetime,decaytable )
   //-- members initialisation for SetCuts ------------------------------
{
   const G4double DefaultLowestEnergy  = 1.0*keV;
   const G4double DefaultHighestEnergy = 1.0*GeV;
   const G4int DefaultTotBin = 100;

   //-- members initialisation for SetCuts ------------------------------
   LowestEnergy  = DefaultLowestEnergy;
   HighestEnergy = DefaultHighestEnergy;
   TotBin = DefaultTotBin;

   // Anti-particle of gamma is gamma itself
   SetAntiPDGEncoding(encoding);
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

G4Gamma G4Gamma::theGamma(
	      "gamma",          0.0*MeV,       0.0*MeV,         0.0, 
		    2,              -1,            -1,          
		    0,               0,             0,             
	      "gamma",               0,             0,          22,
		 true,             0.0,          NULL
);
G4Gamma*  G4Gamma::GammaDefinition() {return &theGamma;}
// initialization for static cut values
G4double   G4Gamma::theGammaLengthCut = -1.0;
G4double*  G4Gamma::theGammaKineticEnergyCuts = NULL;

// ***********************************************************************
// ******************* BuildAbsorptionLengthVector ***********************
// ***********************************************************************

void G4Gamma::BuildAbsorptionLengthVector(
                            const G4Material* aMaterial,
			    const G4CrossSectionTable* aCrossSectionTable,
			    G4double       ,     
			    G4double       ,
			    G4RangeVector* absorptionLengthVector )
{
  // fill the absorption length vector for this material
  // absorption length is defined here as
  //
  //    absorption length = 5./ macroscopic absorption cross section
  //
  const G4ElementVector* elementVector = aMaterial->GetElementVector();
  const G4double* atomicNumDensityVector = aMaterial->GetAtomicNumDensityVector();
  //  fill absorption length vector
  if (absorptionLengthVector == NULL) {
    G4String aErrorMessage("Error in G4Gamma::BuildAbsorptionLengthVector()");
    aErrorMessage += " NULL pointer is found in absorptionLengthVector \n";
    G4Exception((const char*)aErrorMessage);
  }
  G4int NumEl = aMaterial->GetNumberOfElements();
  G4double absorptionLengthMax = 0.0;
  for (G4int ibin=0; ibin<TotBin; ibin++)
  {
    G4double lowEdgeEnergy = absorptionLengthVector->GetLowEdgeEnergy(ibin);

    G4double SIGMA = 0. ;
    for (G4int iel=0; iel<NumEl; iel++)
    {
      G4bool isOut;
      G4int IndEl = (*elementVector)(iel)->GetIndex();
      SIGMA +=  atomicNumDensityVector[iel]*
                (*aCrossSectionTable)[IndEl]->GetValue(lowEdgeEnergy,isOut);
    }
    //  absorption length=5./SIGMA
    absorptionLengthVector->PutValue(ibin, 5./SIGMA);
    if (absorptionLengthMax < 5./SIGMA ) absorptionLengthMax = 5./SIGMA;
  }
  if ( theCutInMaxInteractionLength >= absorptionLengthMax ) {
      G4cout << "******** SetCuts for " << GetParticleName(); 
      G4cout << " ********************" << endl;
      G4cout << "The maximal meaningful cut is ";
      G4cout << absorptionLengthMax/mm << " mm." << endl;
      G4cout << "All the " << GetParticleName() << "will be killed !" << endl;
      G4cout << "in the material " << aMaterial->GetName() << "." << endl;
  }
}

// ***********************************************************************
// ********************** ComputeCrossSection ****************************
// ***********************************************************************

G4double G4Gamma::ComputeCrossSection(G4double AtomicNumber,
                                      G4double KineticEnergy) const
{
  //  Compute the "absorption" cross section of the photon "absorption"
  //  cross section means here the sum of the cross sections of the
  //  pair production, Compton scattering and photoelectric processes
  static G4double Z;  
  const  G4double t1keV = 1.*keV;
  const  G4double t200keV = 200.*keV;
  const  G4double t100MeV = 100.*MeV;

  static G4double s200keV, s1keV;
  static G4double tmin, tlow; 
  static G4double smin, slow;
  static G4double cmin, clow, chigh;
  //  compute Z dependent quantities in the case of a new AtomicNumber
  if(abs(AtomicNumber-Z)>0.1)
  {
    Z = AtomicNumber;
    G4double Zsquare = Z*Z;
    G4double Zlog = log(Z);
    G4double Zlogsquare = Zlog*Zlog;

    s200keV = (0.2651-0.1501*Zlog+0.02283*Zlogsquare)*Zsquare;
    tmin = (0.552+218.5/Z+557.17/Zsquare)*MeV;
    smin = (0.01239+0.005585*Zlog-0.000923*Zlogsquare)*exp(1.5*Zlog);
    cmin=log(s200keV/smin)/(log(tmin/t200keV)*log(tmin/t200keV));

    tlow = 0.2*exp(-7.355/sqrt(Z))*MeV;
    slow = s200keV*exp(0.042*Z*log(t200keV/tlow)*log(t200keV/tlow));
    s1keV = 300.*Zsquare;
    clow =log(s1keV/slow)/log(tlow/t1keV);

    chigh=(7.55e-5-0.0542e-5*Z)*Zsquare*Z/log(t100MeV/tmin);
  }

  //  calculate the cross section (using an approximate empirical formula)
  G4double s;
  if ( KineticEnergy<tlow ) {
    if(KineticEnergy<t1keV) s = slow*exp(clow*log(tlow/t1keV));
    else                    s = slow*exp(clow*log(tlow/KineticEnergy));
  } else if ( KineticEnergy<t200keV ) {
    s = s200keV
         * exp(0.042*Z*log(t200keV/KineticEnergy)*log(t200keV/KineticEnergy));
  } else if( KineticEnergy<tmin ){
    s = smin
         * exp(cmin*log(tmin/KineticEnergy)*log(tmin/KineticEnergy));
  } else {
    s = smin + chigh*log(KineticEnergy/tmin);
  }
  return s * barn;

}

// **********************************************************************
// ******************** ConvertCutToKineticEnergy ***********************
// **********************************************************************

G4double G4Gamma::ConvertCutToKineticEnergy(
                   G4RangeVector* absorptionLengthVector) const
{
  const G4double epsilon=0.01;
  const G4int NBIN=200;
  
  //  find max. absorption length and the corresponding energy (rmax,Tmax)
  G4double fac=log(HighestEnergy/LowestEnergy)/NBIN;
  fac=exp(fac);
  G4double T=LowestEnergy/fac;
  G4double Tmax = T;
  G4double rmax = -1.e10*mm;
  G4bool isOut;
  for (G4int ibin=0; ibin<NBIN; ibin++)
  {
    T=fac*T;
    G4double r=absorptionLengthVector->GetValue(T,isOut);
    if ( r>rmax ) 
    {
       Tmax=T;
       rmax=r;
    }
  }
  G4double T1 = LowestEnergy;
  G4double r1 = absorptionLengthVector->GetValue(T1,isOut);
  if ( theCutInMaxInteractionLength <= r1 ) {
     return T1;
  } else if ( theCutInMaxInteractionLength >= rmax ) {
      G4cout << "******** ConvertCutToKineticEnergy for " << GetParticleName(); 
      G4cout << " ********************" << endl;
      G4cout << "The cut energy is set " << DBL_MAX/GeV << "GeV " <<endl; 
      return  DBL_MAX;
  } else {
    G4double T2 = Tmax;
    G4double T3 = sqrt(T1*T2);
    G4double r3 = absorptionLengthVector->GetValue(T3,isOut);
    while ( abs(1.-r3/theCutInMaxInteractionLength)>epsilon ) {
      if ( theCutInMaxInteractionLength <= r3 ) 
        T2 = T3;
      else 
        T1 = T3;
    T3 = sqrt(T1*T2);
    r3 = absorptionLengthVector->GetValue(T3,isOut);
    }
    return T3;
  }
}
