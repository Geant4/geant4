// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleWithCuts.cc,v 1.6 1999-12-15 14:51:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
//      Hisaya Kurashige, 21 Oct 1996
//      Hisaya Kurashige, 04 Jan 1997
// --------------------------------------------------------------
//   New Physics scheme           8 Jan. 1997  H.Kurahige
//   The cut in kinetic energy is set to zero for vacuum,
//   bug in ConvertCutToKineticEnergy is corrected . L.Urban 04/04/97
//   remove sprintf              10 Nov. 1997 H.Kurashige
//   remove BuildPhysicTabel()   06 June 1998 H.Kurashige
//   bug in CalcEnergyCuts is corrected , 22 June 1998 L.Urban
//   modify CalcEnergyCuts 09 Nov. 1998, L.Urban
// ------------------------------------------------------------
#include "globals.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ios.hh"

#include "g4std/strstream"


G4double  G4ParticleWithCuts::LowestEnergy = 0.99e-3*MeV;
G4double  G4ParticleWithCuts::HighestEnergy = 100.0e6*MeV;


G4ParticleWithCuts::G4ParticleWithCuts(
		const G4String&  aName,  
                G4double         mass,     
                G4double         width,
                G4double         charge,   
                G4int            iSpin,
                G4int            iParity,
                G4int            iConjugation,
                G4int            iIsospin,   
                G4int            iIsospinZ, 
                G4int            gParity,
                const G4String&  pType,
                G4int            lepton,
                G4int            baryon,
                G4int            encoding,
                G4bool           stable,
                G4double         lifetime,
                G4DecayTable     *decaytable,
		G4bool           shortlived)
	: G4ParticleDefinition(aName, mass, width, charge, iSpin, iParity,
                               iConjugation, iIsospin, iIsospinZ, gParity,
                               pType, lepton, baryon, encoding, stable,
                               lifetime, decaytable, shortlived), 
         NumberOfElements(0),
	 theLossTable(0),
   //-- members initialisation for SetCuts ------------------------------
         theCutInMaxInteractionLength(-1.0),
         theKineticEnergyCuts(0)
{  		   
   // -- set ApplyCutsFlag in default   ----------
         SetApplyCutsFlag(false);

   //-- default values for SetCuts ------------------------------
   //    Lowest/Highest energy is defined in MeV
         TotBin = 200;
}

G4ParticleWithCuts::~G4ParticleWithCuts()
{ 
  if (theKineticEnergyCuts) delete [] theKineticEnergyCuts;
  if (theLossTable) delete  theLossTable;
}

// **********************************************************************
// ************************ RangeLinSimpson *****************************
// **********************************************************************

G4double G4ParticleWithCuts::RangeLinSimpson(
				     const G4ElementVector* elementVector,
                                     const G4double* atomicNumDensityVector,
				     const G4LossTable* aLossTable,
	                             G4double aMass,   
			             G4double taulow, G4double tauhigh,
                                     G4int nbin, G4int NumEl)
{
  // Simpson numerical integration, linear binning
  G4double dtau = (tauhigh-taulow)/nbin;
  G4double Value=0.;
  for (G4int i=0; i<=nbin; i++)
  {
    G4double taui=taulow+dtau*i;
    G4double ti=aMass*taui;
    G4double lossi=0.;
    for (G4int j=0; j<NumEl; j++)
    {
      G4bool isOut;
      G4int IndEl = (*elementVector)(j)->GetIndex();
      lossi += atomicNumDensityVector[j]*
              (*aLossTable)[IndEl]->GetValue(ti,isOut);
    }
    if ( i==0 )
      Value += 0.5/lossi;
    else {
      if ( i<nbin ) Value += 1./lossi;
      else          Value += 0.5/lossi;
    }
  }
  Value *= aMass*dtau;

  return Value;
}

// **********************************************************************
// ************************ RangeLogSimpson *****************************
// **********************************************************************

G4double G4ParticleWithCuts::RangeLogSimpson(
				     const G4ElementVector* elementVector,
                                     const G4double* atomicNumDensityVector,
				     const G4LossTable* aLossTable,
                                     G4double aMass,   
				     G4double ltaulow, G4double ltauhigh,
                                     G4int nbin, G4int NumEl)
{
  // Simpson numerical integration, logarithmic binning
  G4double ltt = ltauhigh-ltaulow;
  G4double dltau = ltt/nbin;

  G4double Value = 0.;
  for (G4int i=0; i<=nbin; i++)
  {
    G4double ui = ltaulow+dltau*i;
    G4double taui = exp(ui);
    G4double ti = aMass*taui;
    G4double lossi = 0.;
    
    for (G4int j=0; j<NumEl; j++)
    {
      G4bool isOut;
      G4int IndEl = (*elementVector)(j)->GetIndex();
      lossi += atomicNumDensityVector[j]*
              (*aLossTable)[IndEl]->GetValue(ti,isOut);
    }
    if ( i==0 )
      Value +=  0.5*taui/lossi;
    else {
      if ( i<nbin ) Value += taui/lossi;
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
void G4ParticleWithCuts::BuildLossTable()
{
   //  Build dE/dx tables for elements
  if (NumberOfElements ==0) 
  {
    NumberOfElements = G4Element::GetNumberOfElements();
    theLossTable = new
               G4LossTable(G4Element::GetNumberOfElements());
#ifdef G4VERBOSE
    if (GetVerboseLevel()>2) {
      G4cout << "G4ParticleWithCuts::BuildLossTable() ";
      G4cout << "Create theLossTable[" << theLossTable << "]";
      G4cout << " NumberOfElements=" << NumberOfElements <<G4endl;
    }
#endif
  } else {
    if (NumberOfElements != G4Element::GetNumberOfElements()){
      char errMsg[1024];
      G4std::ostrstream errOs(errMsg,1024);
      errOs << "Error in G4ParticlWithCuts::BuildLossTable()";
      errOs << "[" <<  this->GetParticleName() << "] ";
      errOs << "  :  inconsistent G4Element::GetNumberOfElements =" << G4Element::GetNumberOfElements();
      errOs << " previous value" << NumberOfElements << '\0';
      G4Exception(errMsg);
    }
  }

  // fill the loss table
  for (G4int J=0; J<NumberOfElements; J++)
  {
    G4double Value;
    G4LossVector* aVector= new
            G4LossVector(LowestEnergy, HighestEnergy, TotBin);
    for (G4int i=0; i<TotBin; i++)
    {
      Value = ComputeLoss(
		           (*G4Element::GetElementTable())[J]->GetZ(),
			    aVector->GetLowEdgeEnergy(i)
	         	  );
     aVector->PutValue(i,Value);
    }
    theLossTable->insert(aVector);
  }
}
// **********************************************************************
// ****************** ConvertCutToKineticEnergy *************************
// **********************************************************************

G4double G4ParticleWithCuts::ConvertCutToKineticEnergy(G4RangeVector* rangeVector) const
{
  const G4double epsilon=0.01;
  const G4int NBIN=200;

  //  find max. range and the corresponding energy (rmax,Tmax)
  G4double Tmax=HighestEnergy;
  G4double rmax=-1.e10*mm;

  G4double fac=log(HighestEnergy/LowestEnergy)/NBIN;
  fac=exp(fac);
  G4double T=LowestEnergy/fac;
  G4bool isOut;
  for (G4int ibin=0; ibin<NBIN; ibin++)
  {
    T=fac*T;
    G4double r=rangeVector->GetValue(T,isOut);
    if ( r>rmax )
    {
       Tmax=T;
       rmax=r;
    }
  }

  G4double T1 = LowestEnergy;
  G4double r1 = rangeVector->GetValue(T1,isOut);
  if ( theCutInMaxInteractionLength <= r1 )
    return T1;

  if ( theCutInMaxInteractionLength >= rmax )
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "Error in G4ParticleWithCuts::ConvertCutToKineticEnergy" <<G4endl;
      G4cout << "******** ConvertCutToKineticEnergy for " << GetParticleName(); 
      G4cout << " ********************" << G4endl;
      G4cout << "The cut energy is set " << DBL_MAX/GeV << "GeV " <<G4endl; 
    }
#endif
    return  DBL_MAX;
  } else {
    G4double T2 = Tmax ;
    G4double T3 = sqrt(T1*T2);
    G4double r3 = rangeVector->GetValue(T3,isOut);
    while ( abs(1.-r3/theCutInMaxInteractionLength)>epsilon )
    {
      if ( theCutInMaxInteractionLength <= r3 ) {
	T2 = T3;
      } else {
	T1 = T3;
      }
      T3 = sqrt(T1*T2);
      r3 = rangeVector->GetValue(T3,isOut);
    }
    return T3;
  }
}

// **********************************************************************
// **************************** SetCuts *********************************
// **********************************************************************

void  G4ParticleWithCuts::CalcEnergyCuts(G4double aCut) 
{
  char errMsg[1024];
  G4std::ostrstream errOs(errMsg,1024);
  // check LowestEnergy/ HighestEnergy/TotBin 
  if (TotBin<1) {
    errOs <<  "Error in G4ParticlWithCuts::G4ParticlWithCuts" ;
    errOs << "[" << this->GetParticleName() << "]";
    errOs << " :  not defined or illegal TotBin [" << TotBin << "]" << '\0';
    G4Exception(errMsg);
  }
  if ( (LowestEnergy<0.0)||(HighestEnergy<=LowestEnergy) ){
    errOs << "Error in G4ParticlWithCuts::G4ParticlWithCuts";
    errOs << "[" << this->GetParticleName() << "]";
    errOs << " :  illegal energy range" << "(" << LowestEnergy/GeV;
    errOs << "," << HighestEnergy/GeV << ") [GeV]" << '\0';
    G4Exception(errMsg);
  }

  // Set cut in stopping range
  theCutInMaxInteractionLength = aCut;

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();

  // Create the vector of cuts in energy
  // corresponding to the stopping range cut
  if(theKineticEnergyCuts) delete [] theKineticEnergyCuts;
  theKineticEnergyCuts = new G4double [materialTable->length()];
                    
  G4double Charge = this->GetPDGCharge() ;

  G4bool useProtonCut =
           ((GetParticleName() != "gamma" ) &&
            (GetParticleName() != "e-"    ) &&
            (GetParticleName() != "e+"    ) &&
            (GetParticleName() != "mu-"   ) &&
            (GetParticleName() != "mu+"   ) &&
            (GetParticleName() != "proton" ) &&
            (GetParticleName() != "anti_proton" ) && 
            (Charge != 0.)                             );

  static G4ParticleDefinition* theProton =0;   

  // check if the proton exists or not 
  if ((useProtonCut) && (theProton ==0)) {
    theProton =   G4ParticleTable::GetParticleTable()->FindParticle("proton");
    if (theProton ==0) {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0) {
	    G4cout << " G4ParticleWithCuts::CalcEnergyCuts  ";
	    G4cout << " proton is not defined !!" << G4endl;
      }
#endif
      useProtonCut = false;
    }
  } 

  if (useProtonCut) {
    // check if cuts for the proton are defined or not
    if (theProton->GetEnergyCuts()==0) {
      errOs << " G4ParticleWithCuts::CalcEnergyCuts  ";
      errOs << "   proton energy cut is not defined !!" << '\0';    
      G4Exception(errMsg);  
    }
    //  check if the cut in range is same as one fro the proton
    useProtonCut =  ( abs(aCut-theProton->GetLengthCuts())<1.*nanometer );  
  }	

  if (useProtonCut) {
    // use energy cuts for Proton
#ifdef G4VERBOSE
    if (GetVerboseLevel()>2) {
      G4cout << " G4ParticleWithCuts: [" << GetParticleName() <<"]";
      G4cout << " uses Proton Cut " << G4endl;
    }
#endif                                                      
    G4double ChargeSquare = Charge*Charge/(eplus*eplus) ;
    G4double massRatio = proton_mass_c2/(this->GetPDGMass()) ;
    
    for (G4int J=0; J<materialTable->length(); J +=1) {
      G4double protonEnergyCut = (theProton->GetEnergyCuts())[J] ;           
      // cut energy is rescaled by using charge and mass ratio
      theKineticEnergyCuts[J] = ChargeSquare*protonEnergyCut/massRatio ; 
      if(theKineticEnergyCuts[J] < LowestEnergy) {
        theKineticEnergyCuts[J] = LowestEnergy ;
      }
    }
  } else {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>2) {
      G4cout << " G4ParticleWithCuts: [" << GetParticleName() <<"]";
      G4cout << " calcurate by using its own loss table  " << G4endl;
    }
#endif
    // Build the energy loss table
    BuildLossTable();
    
    // Build range vector for every material, convert cut into energy-cut,
    // fill theKineticEnergyCuts and delete the range vector
    G4double tune = 0.025*mm*g/cm3 ,lowen = 30.*keV ; 
    G4double density ;
    for (G4int J=0; J<materialTable->length(); J +=1){
      G4RangeVector* rangeVector = new
	G4RangeVector(LowestEnergy, HighestEnergy, TotBin);
      G4Material* aMaterial = (*materialTable)[J];
      density = aMaterial->GetDensity() ;
      if(density == 0.) {
	theKineticEnergyCuts[J] = 0. ;
      } else {
	this->BuildRangeVector(aMaterial, this->theLossTable, 
			       this->HighestEnergy, this->GetPDGMass(),
			       rangeVector);
	theKineticEnergyCuts[J] = ConvertCutToKineticEnergy(rangeVector);
	
	if(    ((GetParticleName()=="e-")||(GetParticleName()=="e+"))
	   && (theKineticEnergyCuts[J] < lowen) ) {
	  theKineticEnergyCuts[J] /= (1.+tune/(aCut*density)) ;
	}
	if(theKineticEnergyCuts[J] < LowestEnergy) {
        theKineticEnergyCuts[J] = LowestEnergy ;
      }
      }
    delete rangeVector;
    }
     
    // Delete energy loss table
    theLossTable->clearAndDestroy();
  }
}

// **********************************************************************
// ************************** ComputeLoss *******************************
// **********************************************************************

G4double G4ParticleWithCuts::ComputeLoss(G4double AtomicNumber,
                                          G4double KineticEnergy) const
{
  //  calculate dE/dx

  static G4double Z;  
  static G4double ionpot, tau0, taum, taul, ca, cba, cc;

  G4double  z2Particle = GetPDGCharge()/eplus;
  z2Particle *=  z2Particle;
  if (z2Particle < 0.1) return 0.0;

  if( abs(AtomicNumber-Z)>0.1 )
  {
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

  G4double tau = KineticEnergy/GetPDGMass();
  G4double dEdx;
  if ( tau <= tau0 )
    dEdx = ca*(sqrt(tau)+cba*tau);
  else
  {
    if( tau <= taul )
      dEdx = cc/sqrt(tau);
    else
    {
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

void G4ParticleWithCuts::BuildRangeVector(
                                  const G4Material* aMaterial,
				  const G4LossTable*   aLossTable,
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
  if (rangeVector == 0) {
    char errMsg[1024];
    G4std::ostrstream errOs(errMsg,1024);
    errOs << "Error in G4ParticleWithCuts::BuildRangeVector()";
    errOs << "[" << this->GetParticleName() << "] ";
    errOs << " :  0 pointer is found in absorptionLengthVector" << '\0';
    G4Exception(errMsg);
  }

  // calculate parameters of the low energy part first
  G4double loss1=0.;
  G4double loss2=0.;
  G4int i;
  for (i=0; i<NumEl; i++)
  {
    G4bool isOut;
    G4int IndEl = (*elementVector)(i)->GetIndex();
    loss1 += atomicNumDensityVector[i]*
            (*aLossTable)[IndEl]->GetValue(t1,isOut);
    loss2 += atomicNumDensityVector[i]*
            (*aLossTable)[IndEl]->GetValue(t2,isOut);
  }
  G4double tau1 = t1/proton_mass_c2;
  G4double sqtau1 = sqrt(tau1);
  G4double ca = (4.*loss2-loss1)/sqtau1;
  G4double cb = (2.*loss1-4.*loss2)/tau1;
  G4double cba = cb/ca;
  G4double taulim = tlim/proton_mass_c2;
  G4double taumax = maxEnergy/aMass;
  G4double ltaulim = log(taulim);
  G4double ltaumax = log(taumax);

  // now we can fill the range vector....
  G4double  rmax = 0.0;
  for (i=0; i<TotBin; i++)
  {
    G4double  LowEdgeEnergy = rangeVector->GetLowEdgeEnergy(i);
    G4double  tau = LowEdgeEnergy/aMass;
    G4double  Value;
 
    if ( tau <= tau1 ){
      Value =2.*aMass*log(1.+cba*sqrt(tau))/cb;
    } else {
      Value = 2.*aMass*log(1.+cba*sqtau1)/cb;
      if ( tau <= taulim )
      {
        G4int nbin = (G4int)(maxnbint*(tau-tau1)/(taulim-tau1));
        if ( nbin<1 ) nbin = 1;
        Value += RangeLinSimpson(elementVector,atomicNumDensityVector, 
				 aLossTable, aMass,
                                 tau1, tau,
                                 nbin, NumEl);
      } else {
        Value += RangeLinSimpson(elementVector,atomicNumDensityVector,
				   aLossTable, aMass,
                                   tau1, taulim,
                                   maxnbint, NumEl);
        G4double ltaulow  = log(taulim);
        G4double ltauhigh = log(tau);
        G4int nbin = (G4int)(maxnbint*(ltauhigh-ltaulow)/(ltaumax-ltaulow));
        if ( nbin<1 ) nbin = 1;
        Value += RangeLogSimpson(elementVector,atomicNumDensityVector,
				aLossTable, aMass,
				ltaulow, ltauhigh,
                                nbin, NumEl);
      }
    }
    rangeVector->PutValue(i,Value); 
    if (rmax < Value) rmax = Value;
  }
  if ( theCutInMaxInteractionLength >= rmax) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "Error in G4ParticleWithCuts::BuildRangeVector()" << G4endl;
      G4cout << " SetCuts for " << GetParticleName() << G4endl;
      G4cout << "The maximal meaningful cut is " << rmax/mm << " mm." << G4endl;
      G4cout << "All the " << GetParticleName() << "will be killed !" << G4endl;
      G4cout << "in the material " << aMaterial->GetName() << "." << G4endl;
    }
#endif
  }
} 























