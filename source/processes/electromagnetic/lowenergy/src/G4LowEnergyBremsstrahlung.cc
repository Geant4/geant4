// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyBremsstrahlung.cc,v 1.24 2000-04-19 13:30:27 lefebure Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      ------------ G4LowEnergyBremsstrahlung: low energy modifications --------
//                   by Alessandra Forti, March 1999
//
// **************************************************************
// 
// 18.04.2000 V.L.
// - First implementation of continuous energy loss.
// 17.02.2000 Veronique Lefebure
//  - correct bug : the gamma energy was not deposited when the gamma was 
//    not produced when its energy was < CutForLowEnergySecondaryPhotons
//
// Added Livermore data table construction methods A. Forti
// Modified BuildMeanFreePath to read new data tables A. Forti
// Modified PostStepDoIt to insert sampling with with EEDL data A. Forti
// Added SelectRandomAtom A. Forti
// Added map of the elements A. Forti
// --------------------------------------------------------------

#include "G4LowEnergyBremsstrahlung.hh"
#include "G4EnergyLossTables.hh"
#include "G4UnitsTable.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
// constructor
 
G4LowEnergyBremsstrahlung::G4LowEnergyBremsstrahlung(const G4String& processName)
  : G4eLowEnergyLoss(processName),      // initialization
    theCrossSectionTable(0),
    theMeanFreePathTable(0),
    ATable(0),
    BTable(0),
    ZNumVec(0),
    lowEnergyCut(0.1*eV),
    CutForLowEnergySecondaryPhotons(0.)
{ 
    LowestKineticEnergy  = GetLowerBoundEloss();
    HighestKineticEnergy = GetUpperBoundEloss();
    TotBin = GetNbinEloss();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
// destructor
 
G4LowEnergyBremsstrahlung::~G4LowEnergyBremsstrahlung()
{
     if (theMeanFreePathTable) {
        theMeanFreePathTable->clearAndDestroy();
        delete theMeanFreePathTable;
     }
     if (theCrossSectionTable) {

        delete theCrossSectionTable;
     }

     if(ZNumVec){
       
       ZNumVec->clear();
       delete ZNumVec;
     }
   
     if (ATable) {

        delete ATable;
     }

     if (BTable) {

        delete BTable;
     }

   if (&PartialSumSigma) {

      PartialSumSigma.clearAndDestroy();
   }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// SET CUT FOR LOW ENERGY SECONDARY PHOTONS A. FORTI
void G4LowEnergyBremsstrahlung::SetCutForLowEnSecPhotons(G4double cut){

  CutForLowEnergySecondaryPhotons = cut;
}

  // METHOD BELOW  FROM STANDARD E_M PROCESSES CODE
void G4LowEnergyBremsstrahlung::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{

 
  BuildZVec();

  // energy sampling formula coefficient
  BuildATable();
  BuildBTable();
  BuildCrossSectionTable() ;
  
  BuildLossTable(aParticleType) ;
  if (&aParticleType==G4Electron::Electron()){

    RecorderOfElectronProcess[CounterOfElectronProcess] = (*this).theLossTable ;
    CounterOfElectronProcess++;
    PrintInfoDefinition();  
  }
  else{

    RecorderOfPositronProcess[CounterOfPositronProcess] = (*this).theLossTable ;
    CounterOfPositronProcess++;
   }

  BuildMeanFreePathTable() ;
 
  BuildDEDXTable(aParticleType) ;
 
 
 


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  // CONSTRUCT THE CROSS SECTION TABLE FOR THE ELEMENTS MAPPED IN ZNUMVEC. 
void G4LowEnergyBremsstrahlung::BuildCrossSectionTable(){
 
  if (theCrossSectionTable) {
    
    delete theCrossSectionTable; 
  }

  theCrossSectionTable = new G4SecondLevel();
  G4int dataNum = 2;
 
  for(G4int TableInd = 0; TableInd < ZNumVec->entries(); TableInd++){
    
    G4int AtomInd = (G4int) (*ZNumVec)[TableInd];
    
    G4FirstLevel* oneAtomCS = util.BuildFirstLevelTables(AtomInd, dataNum, "brem/br-cs-");
    
    theCrossSectionTable->insert(oneAtomCS);
    
  }//end for on atoms
}

// CONSTRUCT THE TABLE OF THE FIRST PARAMETER OF THE SAMPLING FORMULA 
void G4LowEnergyBremsstrahlung::BuildATable(){

  if (ATable) {
    
    delete ATable; 
  }
  G4int dataNum = 2;
  ATable = util.BuildSecondLevelTables(0,dataNum,"brem/br-co-a");

}

// CONSTRUCT THE TABLE OF THE PARAMETERS OF THE FORMULA OF THE 
// SECOND PARAMETER OF THE SAMPLING FORMULA
void G4LowEnergyBremsstrahlung::BuildBTable(){

  if (BTable) {
    
    delete BTable; 
  }
  G4int dataNum = 2;
  BTable = util.BuildFirstLevelTables(0, dataNum, "brem/br-co-b");

}

// Vector mapping the existing elements in the material table
// needed at initialization time to load only the necessary data
void G4LowEnergyBremsstrahlung::BuildZVec(){

  const G4MaterialTable* theMaterialTable=G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();

  if(ZNumVec){

    ZNumVec->clear();
    delete ZNumVec;
  }

  ZNumVec = new G4Data(); 
  for (G4int J=0 ; J < numOfMaterials; J++){ 
 
    const G4Material* material= (*theMaterialTable)[J];        
    const G4ElementVector* theElementVector = material->GetElementVector();
    const G4int NumberOfElements = material->GetNumberOfElements() ;

    for (G4int iel=0; iel<NumberOfElements; iel++ ){

      G4double Zel = (*theElementVector)(iel)->GetZ();

      if(ZNumVec->contains(Zel) == FALSE){

	ZNumVec->insert(Zel);
      }
      else{
	
	continue;
      }
    }
  }
}

void G4LowEnergyBremsstrahlung::BuildLossTable(const G4ParticleDefinition& aParticleType)
{
   //  Build table for energy loss due to soft brems
   //  the tables are built for *MATERIALS*
 
  
  //  create table
  
  if (theLossTable) { 
      theLossTable->clearAndDestroy();
      delete theLossTable;
  }
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  const G4int numOfMaterials = theMaterialTable->length();
  theLossTable = new G4PhysicsTable(numOfMaterials);
  
  //  loop for materials
  
  for (G4int J=0; J<numOfMaterials; J++){
    
      // create physics vector and fill it
      G4PhysicsLogVector* aVector = new G4PhysicsLogVector(LowestKineticEnergy,
							   HighestKineticEnergy,
							   TotBin);
      // get material parameters needed for the energy loss calculation
      const G4Material* material= (*theMaterialTable)[J];

      const G4double Tcut = G4Gamma::Gamma()->GetCutsInEnergy()[material->GetIndex()] ;
      G4cout<<"*** LE Bremsstrahlung using Gamma Tcut = "<<Tcut
            <<" for material "<< material->GetName()
	    <<G4endl;
      const G4ElementVector* theElementVector = material->GetElementVector();
      const G4int NumberOfElements = material->GetNumberOfElements() ;
      const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector();
      
      // now comes the loop for the kinetic energy values
      for (G4int i = 0 ; i < TotBin ; i++){

           const G4double LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;
           G4double ionloss = 0.;          
	   // loop for elements in the material
           for (G4int iel=0; iel<NumberOfElements; iel++ ){
                const G4double Z = (*theElementVector)(iel)->GetZ();
                ionloss += GetEnergyLossWithCut(Z,LowEdgeEnergy,Tcut)*
                           theAtomicNumDensityVector[iel] ;
	   }	
           aVector->PutValue(i,ionloss) ;
      }
      theLossTable->insert(aVector);
    }
}

      
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

//
// METHOD BELOW  FROM STANDARD E_M PROCESSES CODE MODIFIED TO USE 
// LIVERMORE DATA (using log-log interpolation as reported in stepanek paper) 
//
void G4LowEnergyBremsstrahlung::BuildMeanFreePathTable()

// Build  mean free path tables for the gamma emission by e- or e+.
// tables are Build for MATERIALS. see GENERAL part of processes in GEANT4
  // manual
{
   G4double FixedEnergy = (LowestKineticEnergy + HighestKineticEnergy)/2.;

   //create table
   if (theMeanFreePathTable) {
       theMeanFreePathTable->clearAndDestroy();
       delete theMeanFreePathTable;
   }

   G4double NumbOfMaterials = G4Material::GetNumberOfMaterials();
   const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
   G4Material* material;
   G4double* CutInKineticEnergy = G4Gamma::Gamma()->GetCutsInEnergy() ;

   PartialSumSigma.resize(NumbOfMaterials);

   G4double LowEdgeEnergy , Value;
   theMeanFreePathTable = new G4PhysicsTable(NumbOfMaterials);
   G4PhysicsLogVector* ptrVector;

   for ( G4int J=0 ; J < NumbOfMaterials; J++ ){ 
     
     //create physics vector then fill it ....
     ptrVector = new G4PhysicsLogVector(LowestKineticEnergy, HighestKineticEnergy,
					TotBin ) ;
     
     material= (*theMaterialTable)(J);
     const G4ElementVector* theElementVector = material->GetElementVector();
     const G4double* theAtomNumDensityVector = material->GetAtomicNumDensityVector();   
     const G4double Threshold = CutInKineticEnergy[J] ;
        
     for ( G4int i = 0 ; i < TotBin ; i++ ){
       
       LowEdgeEnergy = ptrVector->GetLowEdgeEnergy( i ) ;
       const G4double BigPath= DBL_MAX;
       G4double SIGMA = 0 ;
       
       for ( G4int k=0 ; k < material->GetNumberOfElements() ; k++ ){ 
	 
	 G4int AtomIndex = (G4int) (*theElementVector)(k)->GetZ();
         G4double interCrsSec = GetCrossSectionWithCut(AtomIndex, LowEdgeEnergy,Threshold);
	 SIGMA += theAtomNumDensityVector[k]*interCrsSec;
       }       
       
       Value = SIGMA<=0.0 ? BigPath : 1./SIGMA ;
       ptrVector->PutValue( i , Value ) ;

     }
     
     theMeanFreePathTable->insert( ptrVector );
     
     // Compute the PartialSumSigma table at a given fixed energy
     ComputePartialSumSigma(FixedEnergy, material,Threshold) ;       
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//
// METHOD BELOW  FROM STANDARD E_M PROCESSES CODE MODIFIED TO USE 
// LIVERMORE DATA (using log-log interpolation as reported in stepanek paper)
//
void G4LowEnergyBremsstrahlung::ComputePartialSumSigma(const G4double KineticEnergy,
						       const G4Material* aMaterial,
						       const G4double Threshold)

// Build the table of cross section per element. The table is built for MATERIALS.
// This table is used by DoIt to select randomly an element in the material. 

{
   G4int Imate = aMaterial->GetIndex();
   G4int NbOfElements = aMaterial->GetNumberOfElements();
   const G4ElementVector* theElementVector = aMaterial->GetElementVector(); 
   const G4double* theAtomNumDensityVector = aMaterial->GetAtomicNumDensityVector();

   PartialSumSigma(Imate) = new G4ValVector(NbOfElements);

   G4double SIGMA = 0. ;

   for ( G4int Ielem=0 ; Ielem < NbOfElements ; Ielem++ ){

     G4int AtomIndex = (G4int) (*theElementVector)(Ielem)->GetZ();
     
     G4double interCrsSec = GetCrossSectionWithCut(AtomIndex,KineticEnergy,Threshold);
     
     SIGMA += theAtomNumDensityVector[Ielem]*interCrsSec;
	 
     PartialSumSigma(Imate)->insert(SIGMA);
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4LowEnergyBremsstrahlung::PostStepDoIt(const G4Track& trackData,
							   const G4Step& stepData){

  // This parametrization is derived from : 
  // Migdal corrections (dielectric suppression). 
  // Migdal: Phys Rev 103:1811 (1956); Messel & Crawford: Pergamon Press (1970)
  //     
  
  
  aParticleChange.Initialize(trackData);
  
  G4Material* aMaterial=trackData.GetMaterial() ;
  
  
  const G4DynamicParticle* aDynamicParticle=trackData.GetDynamicParticle();
  G4double charge = aDynamicParticle->GetDefinition()->GetPDGCharge();   
  
  G4double ElectKinEn = aDynamicParticle->GetKineticEnergy();

  if(ElectKinEn <= LowestKineticEnergy){
    
    aParticleChange.SetStatusChange(fStopAndKill);
    aParticleChange.SetEnergyChange(0.);
    aParticleChange.SetLocalEnergyDeposit(ElectKinEn);
    
    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

  }

  G4ParticleMomentum ElectDirection = aDynamicParticle->GetMomentumDirection();
  
  // Gamma production cut in this material
  G4double GammaEnergyCut = (G4Gamma::GetCutsInEnergy())[aMaterial->GetIndex()];
  
  
  // check against insufficient energy
  if (ElectKinEn < GammaEnergyCut){    
      aParticleChange.SetEnergyChange(ElectKinEn);  
      aParticleChange.SetLocalEnergyDeposit(0.);
      return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
  }
  
  // select randomly one element constituing the material  
  G4Element* anElement = SelectRandomAtom(aMaterial);
  
  // limits of the energy sampling
  G4double TotalEnergy = ElectKinEn + electron_mass_c2;
  G4double TotalEnergysquare = TotalEnergy*TotalEnergy ;
  
  //
  // The emitted gamma energy is from EEDL data fitted with A/E+B function.
  // Original formula A/E+B+C*E and sampling methods are reported by  J. Stepanek 
  // formula has been modified by A. Forti and S. Giani. 

  // 
  //  sample the energy of the emitted gamma 
  //  
  G4double p1 = 0, p2 = 0;
  G4double coeffA = 0, coeffB = 0;
  G4int AtomicNum = (G4int) anElement->GetZ();
  coeffA = ComputeA(AtomicNum, ElectKinEn);
  coeffB = ComputeB(AtomicNum, ElectKinEn);
  
  //const G4double minEn = lowEnergyCut; 
  const G4double minEn = GammaEnergyCut;


  p1 = coeffA*log(ElectKinEn/minEn);
  p2 = coeffB*(ElectKinEn - minEn); 
  
  G4double IntegrProb = p1+p2;  
  G4double R1 = G4UniformRand()*IntegrProb;
  
  G4double GammaEnergy;

  if(R1 <= p1){ 
    
    G4double R2 = G4UniformRand();
    GammaEnergy = ElectKinEn*pow((minEn/ElectKinEn),R2);
    /// stepanek does: GammaEnergy = exp(R2*log(ElectKinEn/minEn)+log(ElectKinEn));
  }
  else if ((p1 < R1) && (R1 <= p1+p2)){
    
    G4double R2 = G4UniformRand();
    GammaEnergy = ElectKinEn - R2*(ElectKinEn - minEn);
  }

/*
 G4double R1 = minEn + G4UniformRand()*(ElectKinEn- minEn); 
 G4double Max = coeffA/minEn + coeffB;
 G4double R2 = G4UniformRand()*Max;
 while (coeffA/R1 + coeffB < R2){
        R1 = minEn + G4UniformRand()*(ElectKinEn- minEn);
        R2 = G4UniformRand()*Max;
 }
 G4double GammaEnergy = R1;
*/
 
  //**********************//
  // Angular distribution //
  //**********************//
  
  //  angles of the emitted gamma. ( Z - axis along the parent particle)
  //  universal distribution suggested by L. Urban (Geant3 manual (1993) Phys211),
  //  derived from Tsai distribution (Rev Mod Phys 49,421(1977))

  if(GammaEnergy < minEn){
     G4cerr<<"Problem with bremsstrahlung gamma energy sampling: Energy<cut:"
           <<GammaEnergy<<" < "<<minEn
	   <<G4endl;
  }
   
    G4double u;
    const G4double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;
    
    if (9./(9.+d) > G4UniformRand()) u = - log(G4UniformRand()*G4UniformRand())/a1 ;
    else                          u = - log(G4UniformRand()*G4UniformRand())/a2 ;
    
    G4double Teta = u*electron_mass_c2/TotalEnergy ;
    G4double Phi  = twopi * G4UniformRand() ;
    G4double dirx = sin(Teta)*cos(Phi) , diry = sin(Teta)*sin(Phi) , dirz = cos(Teta) ;
    
    G4ThreeVector GammaDirection ( dirx, diry, dirz);
    
    GammaDirection.rotateUz(ElectDirection);   
  
    //
    // Update the incident particle 
    //
    
    G4double NewKinEnergy = ElectKinEn - GammaEnergy;      
    
    //      
    ///final state electron:
    //
    if (NewKinEnergy > 0.){
      
      aParticleChange.SetMomentumChange( ElectDirection );
      aParticleChange.SetEnergyChange( NewKinEnergy );

    } 
    else{
      
      aParticleChange.SetEnergyChange( 0. );
      if (charge<0.){
	
	aParticleChange.SetStatusChange(fStopAndKill);
      }
      else{
	
	aParticleChange.SetStatusChange(fStopButAlive);
      }    
    }
    //      
    ///emitted photon:
    //
    if(GammaEnergy <  GammaEnergyCut){

	 aParticleChange.SetLocalEnergyDeposit(GammaEnergy); 
    }
    else{

	// create G4DynamicParticle object for the Gamma 
	G4DynamicParticle* aGamma= new G4DynamicParticle (G4Gamma::Gamma(),
							  GammaDirection, GammaEnergy);
	
	aParticleChange.SetNumberOfSecondaries(1);
	aParticleChange.AddSecondary(aGamma); 
	aParticleChange.SetLocalEnergyDeposit(0.);
    }
  
    
#ifdef G4VERBOSE
  if(verboseLevel > 15){

    G4cout<<"LE Bremsstrahlung PostStepDoIt"<<G4endl;
  }
#endif
  return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}

G4double G4LowEnergyBremsstrahlung::GetEnergyLossWithCut(const G4double AtomicNumber,
                                                         const G4double KineticEnergy,
                                                         const G4double Tcut){
  const G4double minEn = lowEnergyCut ;
  if(minEn == 0.) G4cerr<<"Minimum Gamma energy should be finite"<<G4endl;

  //  shortcut ..........................
  if(Tcut <= minEn) return 0. ;

  G4double CrossSection = GetCrossSection(AtomicNumber,KineticEnergy) ;
  //  shortcut ..........................
  if(CrossSection <= 0.) return 0. ;
  
  G4double loss = 0.;
   // 
   //  energy spectrum of the emitted gamma 
   //  
   G4double MeanTinc;
   MeanTinc = KineticEnergy;
   const G4double MeanCS = GetCrossSection(AtomicNumber,MeanTinc);
   const G4double coeffA = ComputeA(AtomicNumber, MeanTinc);
   const G4double coeffB = ComputeB(AtomicNumber, MeanTinc);
   //
   //integration of  T*dSigma/dT between Tmin = minEn and Tcut
   //
   G4double Tmax;
   //
   //integration of  T*dSigma/dT between Tmin = minEn and Tcut
   //
   Tmax = Tcut;
   if(Tmax>MeanTinc) Tmax = MeanTinc;
   G4double SmallLoss = 0.;
   SmallLoss = 0.5*coeffB*(Tmax*Tmax - minEn*minEn) + coeffA*(Tmax-minEn);
   if(SmallLoss < 0.) G4cerr<<"Problem with integration of gamma spectrum: SmallLoss = "<<SmallLoss<<G4endl;
   //
   //integration of dSigma/dT between Tmin = minEn and KineticEnergy 
   //
   Tmax = MeanTinc;
   G4double norm = coeffB*(Tmax-minEn) + coeffA*log(Tmax/minEn);
   if(norm <= 0.) G4cerr<<"Problem with integration of gamma spectrum: norm = "<<norm<<G4endl;
  
   SmallLoss *= MeanCS/norm ;
   loss+=SmallLoss;
  return loss ;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4LowEnergyBremsstrahlung::GetCrossSection(const G4double AtomicNumber,
                                                    const G4double KineticEnergy){
						    
         const G4FirstLevel* oneAtomCS 
	       = (*theCrossSectionTable)[ZNumVec->index(AtomicNumber)];
	 
         return util.DataLogInterpolation(KineticEnergy, 
	                             (*(*oneAtomCS)[0]), 
                                     (*(*oneAtomCS)[1]) )*barn;
						    
}
G4double G4LowEnergyBremsstrahlung::GetCrossSectionWithCut(const G4double AtomicNumber,
					    	           const G4double KineticEnergy,
							   const G4double Tcut){
    if(KineticEnergy<=Tcut) return 0.;
    G4double Tmin = Tcut;
    if(Tcut<lowEnergyCut) Tmin = lowEnergyCut;
    G4double Tmax = KineticEnergy;
    
    G4double CrossSection = GetCrossSection(AtomicNumber,KineticEnergy) ;
    if(CrossSection <= 0.) return 0.;

    const G4double coeffA = ComputeA(AtomicNumber, KineticEnergy);
    const G4double coeffB = ComputeB(AtomicNumber, KineticEnergy);
    
    G4double fraction = coeffB*(Tmax-Tmin) + coeffA*log(Tmax/Tmin);
    if(fraction <= 0.) G4cerr<<"Problem with integration of gamma spectrum: fraction = "<<fraction<<G4endl;
    G4double norm = coeffB*(Tmax-lowEnergyCut) + coeffA*log(Tmax/lowEnergyCut);
    if(norm <= 0.) G4cerr<<"Problem with integration of gamma spectrum: norm = "<<norm<<G4endl;
    fraction /= norm;
    
    return CrossSection*fraction;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// METHOD BELOW  FROM STANDARD E_M PROCESSES CODE MODIFIED TO USE 
// LIVERMORE DATA (using log-log interpolation as reported in stepanek paper)
G4Element* G4LowEnergyBremsstrahlung::SelectRandomAtom(G4Material* aMaterial) const
{


  const G4int Index = aMaterial->GetIndex();
  const G4int NumberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();

  G4double rval = G4UniformRand()*((*PartialSumSigma(Index))(NumberOfElements-1));
  for ( G4int i=0; i < NumberOfElements; i++ )
    if (rval <= (*PartialSumSigma(Index))(i)) return ((*theElementVector)(i));
  return (*theElementVector)(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyBremsstrahlung::PrintInfoDefinition()
{
  G4String comments = "Total cross sections from EEDL database,";
           comments += "Gamma energy sampled from a parametrised formula.";
           comments += "First implementation of the continuous dE/dx part.";  
           comments += "\n At present it can be used for electrons only ";
           comments += " in the energy range [250eV,100GeV]";
                     
	   G4cout << G4endl << GetProcessName() << ":  " << comments<<G4endl;

}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....














