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
//
// $Id: G4VRangeToEnergyConverter.cc 93090 2015-10-05 13:14:43Z gcosmo $
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file/  History:
//    5 Oct. 2002, H.Kuirashige : Structure created based on object model
// --------------------------------------------------------------

#include "G4VRangeToEnergyConverter.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4PhysicsLogVector.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

G4double  G4VRangeToEnergyConverter::LowestEnergy = 0.99e-3*MeV;
G4double  G4VRangeToEnergyConverter::HighestEnergy = 100.0e6*MeV;
G4double  G4VRangeToEnergyConverter::MaxEnergyCut = 10.0*GeV;

G4VRangeToEnergyConverter::G4VRangeToEnergyConverter():
  theParticle(0), theLossTable(0), NumberOfElements(0), TotBin(300),
  verboseLevel(1)
{
  fMaxEnergyCut = 0.;
}

G4VRangeToEnergyConverter::G4VRangeToEnergyConverter(const G4VRangeToEnergyConverter& right) :  theParticle(right.theParticle), theLossTable(0), TotBin(right.TotBin)
{
  fMaxEnergyCut = 0.;
  if (theLossTable) {
    theLossTable->clearAndDestroy();
    delete theLossTable;
    theLossTable=0;
  }
  *this = right;
}

G4VRangeToEnergyConverter & G4VRangeToEnergyConverter::operator=(const G4VRangeToEnergyConverter &right)
{
  if (this == &right) return *this;
  if (theLossTable) {
    theLossTable->clearAndDestroy();
    delete theLossTable;
    theLossTable=0;
 }

  fMaxEnergyCut = right.fMaxEnergyCut;
  NumberOfElements = right.NumberOfElements;
  theParticle = right.theParticle;
  verboseLevel = right.verboseLevel;
  
  // create the loss table
  theLossTable = new G4LossTable();
  theLossTable->reserve(G4Element::GetNumberOfElements());  
  // fill the loss table
  for (size_t j=0; j<size_t(NumberOfElements); j++){
    G4LossVector* aVector = new G4LossVector(LowestEnergy, MaxEnergyCut, TotBin);
    for (size_t i=0; i<=size_t(TotBin); i++) {
      G4double Value = (*((*right.theLossTable)[j]))[i];
      aVector->PutValue(i,Value);
    }
    theLossTable->insert(aVector);
  }

  // clean up range vector store
  for (size_t idx=0; idx<fRangeVectorStore.size(); idx++){
    delete fRangeVectorStore.at(idx);
  }
  fRangeVectorStore.clear();

  // copy range vector store
  for (size_t j=0; j<((right.fRangeVectorStore).size()); j++){
    G4RangeVector* vector = (right.fRangeVectorStore).at(j);
    G4RangeVector* rangeVector = 0; 
    if (vector !=0 ) {
      rangeVector = new G4RangeVector(LowestEnergy, MaxEnergyCut, TotBin);
      fMaxEnergyCut = MaxEnergyCut;   
      for (size_t i=0; i<=size_t(TotBin); i++) {
	G4double Value = (*vector)[i];
	rangeVector->PutValue(i,Value);
      }
    }
    fRangeVectorStore.push_back(rangeVector);
  }
  return *this;
}


G4VRangeToEnergyConverter::~G4VRangeToEnergyConverter()
{ 
  Reset();
  // Comment out Reset() for MT application  

/////MA  // delete loss table without deleteing vectors  
/////MA  if (theLossTable) {  
/////MA    delete theLossTable;
/////MA  }
/////MA  theLossTable=0;
/////MA  NumberOfElements=0;
/////MA  
/////MA  //clear RangeVectorStore without deleteing vectors
/////MA  fRangeVectorStore.clear();

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
#ifdef G4VERBOSE
    if (GetVerboseLevel()>3) {
      G4cout << "G4VRangeToEnergyConverter::Convert() ";
      G4cout << "Convert for " << material->GetName() 
	     << " with Range Cut " << rangeCut/mm << "[mm]" << G4endl;
    }
#endif

  G4double theKineticEnergyCuts = 0.;

  if (fMaxEnergyCut != MaxEnergyCut) {
    fMaxEnergyCut = MaxEnergyCut;      
    // clear loss table and renge vectors
    Reset();
  }
 
  // Build the energy loss table
  BuildLossTable();
  
  // Build range vector for every material, convert cut into energy-cut,
  // fill theKineticEnergyCuts and delete the range vector
  static const G4double tune = 0.025*mm*g/cm3 ,lowen = 30.*keV ; 

  // check density
  G4double density = material->GetDensity() ;
  if(density <= 0.) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4VRangeToEnergyConverter::Convert() ";
      G4cout << material->GetName() << "has zero density "
	     << "( " << density << ")" << G4endl;
    }
#endif
    return 0.;
  }
 
   // initialize RangeVectorStore
  const G4MaterialTable* table = G4Material::GetMaterialTable();
  G4int ext_size = table->size() - fRangeVectorStore.size();
  for (int i=0; i<ext_size; i++) fRangeVectorStore.push_back(0);
  
  // Build Range Vector
  G4int idx = material->GetIndex(); 
  G4RangeVector* rangeVector = fRangeVectorStore.at(idx);
  if (rangeVector == 0) {
    rangeVector = new G4RangeVector(LowestEnergy, MaxEnergyCut, TotBin); 
    BuildRangeVector(material, rangeVector);
    fRangeVectorStore.at(idx) = rangeVector;
  }

  // Convert Range Cut ro Kinetic Energy Cut 
  theKineticEnergyCuts = ConvertCutToKineticEnergy(rangeVector, rangeCut, idx);
  
  if( ((theParticle->GetParticleName()=="e-")||(theParticle->GetParticleName()=="e+"))
      && (theKineticEnergyCuts < lowen) ) {
    //  corr. should be switched on smoothly   
    theKineticEnergyCuts /= (1.+(1.-theKineticEnergyCuts/lowen)*
			     tune/(rangeCut*density)); 
  }
  
  if(theKineticEnergyCuts < LowestEnergy) {
    theKineticEnergyCuts = LowestEnergy ;
  } else if(theKineticEnergyCuts > MaxEnergyCut) {
    theKineticEnergyCuts = MaxEnergyCut;
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
#ifdef G4VERBOSE
    G4cerr << "Error in G4VRangeToEnergyConverter::SetEnergyRange";
    G4cerr << " :  illegal energy range" << "(" << lowedge/GeV;
    G4cerr << "," << highedge/GeV << ") [GeV]" << G4endl;
#endif
    G4Exception( "G4VRangeToEnergyConverter::SetEnergyRange()",
		 "ProcCuts101",
		 JustWarning, "Illegal energy range ");
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
// ******************* Get/SetMaxEnergyCut  *****************************
// **********************************************************************
G4double G4VRangeToEnergyConverter::GetMaxEnergyCut()
{
  return MaxEnergyCut;
}

void G4VRangeToEnergyConverter::SetMaxEnergyCut(G4double value)
{
  MaxEnergyCut = value;
}

// **********************************************************************
// ************************ Reset  **************************************
// **********************************************************************
void G4VRangeToEnergyConverter::Reset()
{
  // delete loss table
  if (theLossTable) {  
    theLossTable->clearAndDestroy();
    delete theLossTable;
  }
  theLossTable=0;
  NumberOfElements=0;
  
  //clear RangeVectorStore
  for (size_t idx=0; idx<fRangeVectorStore.size(); idx++){
    delete fRangeVectorStore.at(idx);
  }
  fRangeVectorStore.clear();
} 


// **********************************************************************
// ************************ BuildLossTable ******************************
// **********************************************************************
//   create Energy Loss Table for charged particles 
//   (cross section tabel for neutral )
void G4VRangeToEnergyConverter::BuildLossTable()
{
  if (size_t(NumberOfElements) == G4Element::GetNumberOfElements()) return;
  
  // clear Loss table and Range vectors
  Reset();

  //  Build dE/dx tables for elements
  NumberOfElements = G4Element::GetNumberOfElements();
  theLossTable = new G4LossTable();
  theLossTable->reserve(G4Element::GetNumberOfElements());
#ifdef G4VERBOSE
  if (GetVerboseLevel()>3) {
    G4cout << "G4VRangeToEnergyConverter::BuildLossTable() ";
    G4cout << "Create theLossTable[" << theLossTable << "]";
    G4cout << " NumberOfElements=" << NumberOfElements <<G4endl;
  }
#endif
  
  
  // fill the loss table
  for (size_t j=0; j<size_t(NumberOfElements); j++){
    G4double Value;
    G4LossVector* aVector= 0;
    aVector= new G4LossVector(LowestEnergy, MaxEnergyCut, TotBin);
    for (size_t i=0; i<=size_t(TotBin); i++) {
      Value = ComputeLoss(  (*G4Element::GetElementTable())[j]->GetZ(),
			    aVector->Energy(i)
			    );
      aVector->PutValue(i,Value);
    }
    theLossTable->insert(aVector);
  }
}

// **********************************************************************
// ************************ BuildRangeVector ****************************
// **********************************************************************
void G4VRangeToEnergyConverter::BuildRangeVector(const G4Material* aMaterial,
					     G4PhysicsLogVector* rangeVector)
{
  //  create range vector for a material
  const G4ElementVector* elementVector = aMaterial->GetElementVector();
  const G4double* atomicNumDensityVector = aMaterial->GetAtomicNumDensityVector();
  G4int NumEl = aMaterial->GetNumberOfElements();

  // calculate parameters of the low energy part first
  size_t i;
  std::vector<G4double> lossV;
  for ( size_t ib=0; ib<=size_t(TotBin); ib++) {
    G4double loss=0.;
    for (i=0; i<size_t(NumEl); i++) {
      G4int IndEl = (*elementVector)[i]->GetIndex();
      loss += atomicNumDensityVector[i]*
	        (*((*theLossTable)[IndEl]))[ib];
    }
    lossV.push_back(loss);
  }
   
  // Integrate with Simpson formula with logarithmic binning
  G4double dltau = 1.0;
  if (LowestEnergy>0.) {
      G4double ltt =std::log(MaxEnergyCut/LowestEnergy);
      dltau = ltt/TotBin;
  }

  G4double s0 = 0.;
  G4double Value;
  for ( i=0; i<=size_t(TotBin); i++) {
    G4double t = rangeVector->GetLowEdgeEnergy(i);
    G4double q = t/lossV[i];
    if (i==0) s0 += 0.5*q;
    else s0 += q;
    
    if (i==0) {
       Value = (s0 + 0.5*q)*dltau ;
    } else {
      Value = (s0 - 0.5*q)*dltau ;
    }
    rangeVector->PutValue(i,Value);
  }
} 

// **********************************************************************
// ****************** ConvertCutToKineticEnergy *************************
// **********************************************************************
G4double G4VRangeToEnergyConverter::ConvertCutToKineticEnergy(
				    G4RangeVector* rangeVector,
				    G4double       theCutInLength, 
#ifdef G4VERBOSE
				    size_t         materialIndex
#else
                                    size_t
#endif
				                              ) const
{
  const G4double epsilon=0.01;

  //  find max. range and the corresponding energy (rmax,Tmax)
  G4double rmax= -1.e10*mm;

  G4double T1 = LowestEnergy;
  G4double r1 =(*rangeVector)[0] ;

  G4double T2 = MaxEnergyCut;

  // check theCutInLength < r1 
  if ( theCutInLength <= r1 ) {  return T1; }

  // scan range vector to find nearest bin 
  // ( suppose that r(Ti) > r(Tj) if Ti >Tj )
  for (size_t ibin=0; ibin<=size_t(TotBin); ibin++) {
    G4double T=rangeVector->GetLowEdgeEnergy(ibin);
    G4double r=(*rangeVector)[ibin];
    if ( r>rmax )   rmax=r;
    if (r <theCutInLength ) {
      T1 = T;
      r1 = r;
    } else if (r >theCutInLength ) {
      T2 = T;
      break;
    }
  }

  // check cut in length is smaller than range max
  if ( theCutInLength >= rmax )  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>2) {
      G4cout << "G4VRangeToEnergyConverter::ConvertCutToKineticEnergy ";
      G4cout << "  for " << theParticle->GetParticleName() << G4endl;
      G4cout << "The cut in range [" << theCutInLength/mm << " (mm)]  ";
      G4cout << " is too big  " ;
      G4cout << " for material  idx=" << materialIndex <<G4endl; 
    }
#endif
    return  MaxEnergyCut;
  }
  
  // convert range to energy
  G4double T3 = std::sqrt(T1*T2);
  G4double r3 = rangeVector->Value(T3);
  const size_t MAX_LOOP = 1000; 
  for (size_t loop_count=0; loop_count<MAX_LOOP; ++loop_count){
    if (std::fabs(1.-r3/theCutInLength)<epsilon ) break;
    if ( theCutInLength <= r3 ) {
      T2 = T3;
    } else {
      T1 = T3;
    }
    T3 = std::sqrt(T1*T2);
    r3 = rangeVector->Value(T3);
  }

  return T3;
}

