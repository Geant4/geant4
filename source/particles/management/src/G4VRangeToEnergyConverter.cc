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
// $Id: G4VRangeToEnergyConverter.cc,v 1.2 2002-12-16 11:15:46 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file/  History:
//    5 Oct. 2002, H.Kuirashige : Structure created based on object model
// --------------------------------------------------------------

#include "G4VRangeToEnergyConverter.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4PhysicsLogVector.hh"

#include "G4ios.hh"
#include "g4std/iomanip"
#include "g4std/strstream"

// energy range
G4double  G4VRangeToEnergyConverter::LowestEnergy = 0.99e-3*MeV;
G4double  G4VRangeToEnergyConverter::HighestEnergy = 100.0e6*MeV;

G4VRangeToEnergyConverter::G4VRangeToEnergyConverter():
  theParticle(0), theLossTable(0), NumberOfElements(0), TotBin(200),
  verboseLevel(1)
{
}

G4VRangeToEnergyConverter::G4VRangeToEnergyConverter(const G4VRangeToEnergyConverter& right)
{
  *this = right;
}

G4VRangeToEnergyConverter & G4VRangeToEnergyConverter::operator=(const G4VRangeToEnergyConverter &right)
{
  if (this == &right) return *this;
  if (theLossTable) delete theLossTable;

  NumberOfElements = right.NumberOfElements;
  TotBin = right.TotBin;
  theParticle = right.theParticle;
  verboseLevel = right.verboseLevel;
  
  // create the loss table
  theLossTable = new G4LossTable();
  theLossTable->reserve(G4Element::GetNumberOfElements());  
  // fill the loss table
  for (size_t j=0; j<size_t(NumberOfElements); j++){
    G4LossVector* aVector= new
            G4LossVector(LowestEnergy, HighestEnergy, TotBin);
    for (size_t i=0; i<size_t(TotBin); i++) {
      G4double Value = (*((*right.theLossTable)[j]))[i];
      aVector->PutValue(i,Value);
    }
    theLossTable->insert(aVector);
  }
  return *this;
}


G4VRangeToEnergyConverter::~G4VRangeToEnergyConverter()
{ 
  if (theLossTable) delete  theLossTable;
}

G4int G4VRangeToEnergyConverter::operator==(const G4VRangeToEnergyConverter &right) const
{
  return this == &right;
}

G4int G4VRangeToEnergyConverter::operator!=(const G4VRangeToEnergyConverter &right) const
{
  return this != &right;
}


// **********************************************************************
// ************************* Convert  ***********************************
// **********************************************************************
G4double G4VRangeToEnergyConverter::Convert(G4double rangeCut, 
					    const G4Material* material) 
{
  //???????????? G4double Charge = theParticle->GetPDGCharge();
  G4double Mass   = theParticle->GetPDGMass();
  G4double theKineticEnergyCuts = 0.;
 
  // Build the energy loss table
  if (theLossTable ==0) BuildLossTable();
  
  // Build range vector for every material, convert cut into energy-cut,
  // fill theKineticEnergyCuts and delete the range vector
  //???????????? G4double tune = 0.025*mm*g/cm3 ,lowen = 30.*keV ; 

  G4int idx = material->GetIndex(); 
  G4double density = material->GetDensity() ;
  if(density > 0.) {
    G4RangeVector* rangeVector = new G4RangeVector(LowestEnergy, HighestEnergy, TotBin);
    BuildRangeVector(material, HighestEnergy, Mass, rangeVector);
    theKineticEnergyCuts = ConvertCutToKineticEnergy(rangeVector, rangeCut, idx);
    if(theKineticEnergyCuts < LowestEnergy) {
      theKineticEnergyCuts = LowestEnergy ;
    }
    delete rangeVector;
  }
  return theKineticEnergyCuts;
}

// **********************************************************************
// ************************ SetEnergyRange  *****************************
// **********************************************************************
void G4VRangeToEnergyConverter::SetEnergyRange(G4double lowedge, 
					       G4double highedge)
{
  // check LowestEnergy/ HighestEnergy 
  if ( (lowedge<0.0)||(highedge<=lowedge) ){
    G4cerr << "Error in G4VRangeToEnergyConverter::SetEnergyRange";
    G4cerr << " :  illegal energy range" << "(" << lowedge/GeV;
    G4cerr << "," << highedge/GeV << ") [GeV]" << G4endl;
  } else {
    LowestEnergy = lowedge;
    HighestEnergy = highedge;
  }
}


G4double G4VRangeToEnergyConverter::GetLowEdgeEnergy()
{
  return LowestEnergy;
}
    

G4double G4VRangeToEnergyConverter::GetHighEdgeEnergy()
{
  return HighestEnergy;
}

// **********************************************************************
// ************************ RangeLinSimpson *****************************
// **********************************************************************
G4double G4VRangeToEnergyConverter::RangeLinSimpson(
                                     const G4ElementVector* elementVector,
                                     const G4double* atomicNumDensityVector,
                                     G4double aMass,   
                                     G4double taulow, G4double tauhigh)
{
  // Simpson numerical integration, linear binning
  G4double dtau = (tauhigh-taulow)/TotBin;
  G4double Value=0.;
  for (size_t i=0; i<=size_t(TotBin); i++){
    G4double taui=taulow+dtau*i;
    G4double ti=aMass*taui;
    G4double lossi=0.;
    size_t nEl = elementVector->size();
    for (size_t j=0; j<nEl; j++) {
      G4bool isOut;
      G4int IndEl = (*elementVector)[j]->GetIndex();
      lossi += atomicNumDensityVector[j]*
              (*theLossTable)[IndEl]->GetValue(ti,isOut);
   }
    if ( i==0 ) {
      Value += 0.5/lossi;
    } else {
      if ( i<size_t(TotBin) ) Value += 1./lossi;
      else            Value += 0.5/lossi;
    }
  }
  Value *= aMass*dtau;

  return Value;
}


// **********************************************************************
// ************************ RangeLogSimpson *****************************
// **********************************************************************
G4double G4VRangeToEnergyConverter::RangeLogSimpson(
                                     const G4ElementVector* elementVector,
                                     const G4double* atomicNumDensityVector,
                                     G4double aMass,   
                                     G4double ltaulow, G4double ltauhigh)
{
  // Simpson numerical integration, logarithmic binning
  G4double ltt = ltauhigh-ltaulow;
  G4double dltau = ltt/TotBin;
  G4double Value = 0.;
  for (size_t i=0; i<=size_t(TotBin); i++){
    G4double ui = ltaulow+dltau*i;
    G4double taui = exp(ui);
    G4double ti = aMass*taui;
    G4double lossi = 0.;
    size_t nEl = elementVector->size();
    for (size_t j=0; j<nEl; j++) {
      G4bool isOut;
      G4int IndEl = (*elementVector)[j]->GetIndex();
      lossi += atomicNumDensityVector[j]*
              (*theLossTable)[IndEl]->GetValue(ti,isOut);
    }
    if ( i==0 ) {
      Value +=  0.5*taui/lossi;
    } else {
      if ( i<size_t(TotBin) ) Value += taui/lossi;
      else Value +=  0.5*taui/lossi;
    }
  }
  Value *= aMass*dltau;

  return Value;
}

// **********************************************************************
// ************************ BuildLossTable ******************************
// **********************************************************************
//   create Energy Loss Table for charged particles 
//   (cross section tabel for neutral )
void G4VRangeToEnergyConverter::BuildLossTable()
{
   //  Build dE/dx tables for elements
  if (size_t(NumberOfElements) != G4Element::GetNumberOfElements()) {
    if (theLossTable!=0) delete theLossTable;
    theLossTable =0; 
    NumberOfElements = 0;
  }

  if (NumberOfElements ==0) {
    NumberOfElements = G4Element::GetNumberOfElements();
    theLossTable = new G4LossTable();
    theLossTable->reserve(G4Element::GetNumberOfElements());
#ifdef G4VERBOSE
    if (GetVerboseLevel()>2) {
      G4cout << "G4VRangeToEnergyConverter::BuildLossTable() ";
      G4cout << "Create theLossTable[" << theLossTable << "]";
      G4cout << " NumberOfElements=" << NumberOfElements <<G4endl;
    }
#endif
  }

  // fill the loss table
  for (size_t j=0; j<size_t(NumberOfElements); j++){
    G4double Value;
    G4LossVector* aVector= new
            G4LossVector(LowestEnergy, HighestEnergy, TotBin);
    for (size_t i=0; i<size_t(TotBin); i++) {
      Value = ComputeLoss(  (*G4Element::GetElementTable())[j]->GetZ(),
                            aVector->GetLowEdgeEnergy(i)
                          );
      aVector->PutValue(i,Value);
    }
    theLossTable->insert(aVector);
  }
}

// **********************************************************************
// ************************** ComputeLoss *******************************
// **********************************************************************
G4double G4VRangeToEnergyConverter::ComputeLoss(G4double AtomicNumber,
						G4double KineticEnergy) const
{
  //  calculate dE/dx

  static G4double Z;  
  static G4double ionpot, tau0, taum, taul, ca, cba, cc;

  G4double  z2Particle = theParticle->GetPDGCharge()/eplus;
  z2Particle *=  z2Particle;
  if (z2Particle < 0.1) return 0.0;

  if( abs(AtomicNumber-Z)>0.1 ){
    // recalculate constants
    Z = AtomicNumber;
    G4double Z13 = exp(log(Z)/3.);
    tau0 = 0.1*Z13*MeV/proton_mass_c2;
    taum = 0.035*Z13*MeV/proton_mass_c2;
    taul = 2.*MeV/proton_mass_c2;
    ionpot = 1.6e-5*MeV*exp(0.9*log(Z));
    cc = (taul+1.)*(taul+1.)*log(2.*electron_mass_c2*taul*(taul+2.)/ionpot)/(taul*(taul+2.))-1.;
    cc = 2.*twopi_mc2_rcl2*Z*cc*sqrt(taul);
    ca = cc/((1.-0.5*sqrt(tau0/taum))*tau0);
    cba = -0.5/sqrt(taum);
  }

  G4double tau = KineticEnergy/theParticle->GetPDGMass();
  G4double dEdx;
  if ( tau <= tau0 ) {
    dEdx = ca*(sqrt(tau)+cba*tau);
  } else {
    if( tau <= taul ) {
      dEdx = cc/sqrt(tau);
    } else {
      dEdx = (tau+1.)*(tau+1.)*
	     log(2.*electron_mass_c2*tau*(tau+2.)/ionpot)/(tau*(tau+2.))-1.;
      dEdx = 2.*twopi_mc2_rcl2*Z*dEdx;
    }
  }
  return dEdx*z2Particle ;
}

// **********************************************************************
// ************************ BuildRangeVector ****************************
// **********************************************************************
void G4VRangeToEnergyConverter::BuildRangeVector(
                                  const G4Material* aMaterial,
                                  G4double       maxEnergy,
                                  G4double       aMass,
                                  G4RangeVector* rangeVector)
{
  //  create range vector for a material
  const G4double tlim=2.*MeV, t1=0.1*MeV, t2=0.025*MeV; 
  const G4int  maxnbint=100;
 
  const G4ElementVector* elementVector = aMaterial->GetElementVector();
  const G4double* atomicNumDensityVector = aMaterial->GetAtomicNumDensityVector();

  G4int NumEl = aMaterial->GetNumberOfElements();

  // calculate parameters of the low energy part first
  G4double loss1=0.;
  G4double loss2=0.;
  size_t i;
  for (i=0; i<size_t(NumEl); i++) {
    G4bool isOut;
    G4int IndEl = (*elementVector)[i]->GetIndex();
    loss1 += atomicNumDensityVector[i]*
            (*theLossTable)[IndEl]->GetValue(t1,isOut);
    loss2 += atomicNumDensityVector[i]*
            (*theLossTable)[IndEl]->GetValue(t2,isOut);
  }
  G4double tau1 = t1/proton_mass_c2;
  G4double sqtau1 = sqrt(tau1);
  G4double ca = (4.*loss2-loss1)/sqtau1;
  G4double cb = (2.*loss1-4.*loss2)/tau1;
  G4double cba = cb/ca;
  G4double taulim = tlim/proton_mass_c2;
  G4double taumax = maxEnergy/aMass;
  G4double ltaumax = log(taumax);

  // now we can fill the range vector....
  G4double  rmax = 0.0;
  for (i=0; i<size_t(TotBin); i++) {
    G4double  LowEdgeEnergy = rangeVector->GetLowEdgeEnergy(i);
    G4double  tau = LowEdgeEnergy/aMass;
    G4double  Value;
 
    if ( tau <= tau1 ){
      Value =2.*aMass*log(1.+cba*sqrt(tau))/cb;
    } else {
      Value = 2.*aMass*log(1.+cba*sqtau1)/cb;
      if ( tau <= taulim ) {
        G4int nbin = (G4int)(maxnbint*(tau-tau1)/(taulim-tau1));
        if ( nbin<1 ) nbin = 1;
        Value += RangeLinSimpson(elementVector,atomicNumDensityVector, 
                                 aMass,
                                 tau1, tau);
      } else {
        Value += RangeLinSimpson(elementVector,atomicNumDensityVector,
                                   aMass,
				   tau1, taulim);
        G4double ltaulow  = log(taulim);
        G4double ltauhigh = log(tau);
        G4int nbin = (G4int)(maxnbint*(ltauhigh-ltaulow)/(ltaumax-ltaulow));
        if ( nbin<1 ) nbin = 1;
        Value += RangeLogSimpson(elementVector,atomicNumDensityVector,
				 aMass,
				 ltaulow, ltauhigh);
      }
    }
    rangeVector->PutValue(i,Value); 
    if (rmax < Value) rmax = Value;
  }
}

// **********************************************************************
// ****************** ConvertCutToKineticEnergy *************************
// **********************************************************************
G4double G4VRangeToEnergyConverter::ConvertCutToKineticEnergy(
				    G4RangeVector* rangeVector,
				    G4double       theCutInLength, 
				    size_t         materialIndex
				                              ) const
{
  const G4double epsilon=0.01;

  //  find max. range and the corresponding energy (rmax,Tmax)
  G4double rmax= -1.e10*mm;
  G4double Tmax= HighestEnergy;
  G4double fac = exp( log(HighestEnergy/LowestEnergy)/TotBin );
  G4double T=LowestEnergy/fac;
  G4bool isOut;

  for (size_t ibin=0; ibin<size_t(TotBin); ibin++) {
    T *= fac;
    G4double r=rangeVector->GetValue(T,isOut);
    if ( r>rmax )    {
       Tmax=T;
       rmax=r;
    }
  }

  // check cut in length is smaller than range max
  if ( theCutInLength >= rmax )  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4VRangeToEnergyConverter::ConvertCutToKineticEnergy ";
      G4cout << "  for " << theParticle->GetParticleName() << G4endl;
      G4cout << "The cut in range [" << theCutInLength/mm << " (mm)]  ";
      G4cout << " is too big  " << G4endl; 
      G4cout << "The cut in energy is set" << DBL_MAX/GeV << "GeV " <<G4endl; 
    }
#endif
    return  DBL_MAX;
  }
  
  // convert range to energy
  G4double T1 = LowestEnergy;
  G4double r1 = rangeVector->GetValue(T1,isOut);
  if ( theCutInLength <= r1 ) return T1;

  G4double T2 = Tmax ;
  G4double T3 = sqrt(T1*T2);
  G4double r3 = rangeVector->GetValue(T3,isOut);
  while ( abs(1.-r3/theCutInLength)>epsilon ) {
    if ( theCutInLength <= r3 ) {
      T2 = T3;
    } else {
      T1 = T3;
    }
    T3 = sqrt(T1*T2);
    r3 = rangeVector->GetValue(T3,isOut);
  }
  return T3;
}

