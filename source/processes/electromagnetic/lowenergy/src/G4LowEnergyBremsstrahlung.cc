// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyBremsstrahlung.cc,v 1.19 2000-01-26 09:50:00 lefebure Exp $
// $Id: G4LowEnergyBremsstrahlung.cc,v 1.19 2000-01-26 09:50:00 lefebure Exp $
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
// Added Livermore data table construction methods A. Forti
// Modified BuildMeanFreePath to read new data tables A. Forti
// Modified PostStepDoIt to insert sampling with with EEDL data A. Forti
// Added SelectRandomAtom A. Forti
// Added map of the elements A. Forti
// --------------------------------------------------------------

// This Class Header
#include "G4LowEnergyBremsstrahlung.hh"

// Collaborating Class Headers
#include "G4EnergyLossTables.hh"
#include "G4UnitsTable.hh"
#include "G4Gamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
// constructor
 
G4LowEnergyBremsstrahlung::G4LowEnergyBremsstrahlung(const G4String& processName)
  : G4eEnergyLoss(processName),      // initialization
    theCrossSectionTable(0),
    theMeanFreePathTable(0),
    ATable(0),
    BTable(0),
    ZNumVec(0),
    LowestKineticEnergy (250.*eV),
    HighestKineticEnergy(100.*GeV),
    lowEnergyCut(0.1*eV),
    CutForLowEnergySecondaryPhotons(0.),
    TotBin(200)
{ 


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

void G4LowEnergyBremsstrahlung::SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins)
{
  LowestKineticEnergy = lowE;  HighestKineticEnergy = highE; TotBin = nBins;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// SET CUT FOR LOW ENERGY SECONDARY PHOTONS A. FORTI
void G4LowEnergyBremsstrahlung::SetCutForLowEnSecPhotons(G4double cut){

  CutForLowEnergySecondaryPhotons = cut;
}

  // METHOD BELOW  FROM STANDARD E_M PROCESSES CODE
void G4LowEnergyBremsstrahlung::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{

    BuildLossTable(aParticleType) ;
 
  if (&aParticleType==G4Electron::Electron()){

    RecorderOfElectronProcess[CounterOfElectronProcess] = (*this).theLossTable ;
    CounterOfElectronProcess++;
  }
  else{

    RecorderOfPositronProcess[CounterOfPositronProcess] = (*this).theLossTable ;
    CounterOfPositronProcess++;
   }

  BuildZVec();

  BuildCrossSectionTable() ;
  BuildMeanFreePathTable() ;

  BuildDEDXTable  (aParticleType) ;

  // energy sampling formula coefficient
  BuildATable();
  BuildBTable();

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

//  Build table for energy loss due to soft brems
//  tables are built for *MATERIALS* already in the standard processes
// to be changed when the new energy loss will be calculated.
//
// // METHOD BELOW  FROM STANDARD E_M PROCESSES LEFT BUT AT THE MOMENT NOT USED
//
void G4LowEnergyBremsstrahlung::BuildLossTable(const G4ParticleDefinition& aParticleType)
{
  G4double KineticEnergy,TotalEnergy,bremloss,Z,x,
           losslim,loss,rate,natom,Cut;

  const G4double MinKinEnergy = 1.*keV;
  const G4double MinCut = 1.*keV;
  const G4double Thigh = 100.*GeV;
  const G4double Cuthigh = 50.*GeV;
  const G4double Factorhigh = 36./(1450.*GeV);
  const G4double coef1 = -0.5, coef2 = 2./9.;

  ParticleMass = aParticleType.GetPDGMass() ;
  G4double* GammaCutInKineticEnergy = G4Gamma::Gamma()->GetEnergyCuts();
  
  //  create table
  
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length() ;

  if (theLossTable) { theLossTable->clearAndDestroy();
                         delete theLossTable;
                    }
                       
  theLossTable = new G4PhysicsTable(numOfMaterials);

//  loop for materials

  for (G4int J=0; J<numOfMaterials; J++)
    {
     // create physics vector and fill it

     G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                               LowestKineticEnergy,HighestKineticEnergy,TotBin);

     // get elements in the material
     const G4Material* material = (*theMaterialTable)[J];
 
     const G4ElementVector* theElementVector = material->GetElementVector();
     const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector();
     const G4int NumberOfElements = material->GetNumberOfElements();

       //  loop for the kinetic energy values
       for (G4int i=0; i<TotBin; i++){

          KineticEnergy = aVector->GetLowEdgeEnergy(i) ;
          TotalEnergy = KineticEnergy+ParticleMass ;

          Cut = GammaCutInKineticEnergy[J] ;
          if (Cut < MinCut) Cut = MinCut ;
          if (Cut > KineticEnergy) Cut = KineticEnergy ;

          bremloss = 0.;

          if (KineticEnergy>MinKinEnergy)
            {
             if (Cut > KineticEnergy) Cut = KineticEnergy ;

             //  loop for elements in the material
             for (G4int iel=0; iel<NumberOfElements; iel++)
               {
                Z=(*theElementVector)(iel)->GetZ();
                natom = theAtomicNumDensityVector[iel] ;
                if (KineticEnergy <= Thigh)
                  {
                   //loss for MinKineticEnergy<KineticEnergy<=100 GeV
                   x=log(TotalEnergy/ParticleMass);
                   loss = ComputeBremLoss(Z,natom,KineticEnergy,Cut,x) ;
                   if (&aParticleType==G4Positron::Positron())
                      loss *= ComputePositronCorrFactorLoss(Z,KineticEnergy,Cut) ;   
                  }
                else
                  {
                   // extrapolation for KineticEnergy>100 GeV
                   x=log(Thigh/ParticleMass) ; 
                   if (Cut<Thigh)
                     {
                      losslim = ComputeBremLoss(Z,natom,Thigh,Cut,x) ;
                      if (&aParticleType==G4Positron::Positron())
                         loss *= ComputePositronCorrFactorLoss(Z,Thigh,Cut) ;   
                      rate = Cut/TotalEnergy ;
                      loss = losslim*(1.+coef1*rate+coef2*rate*rate) ;
                      rate = Cut/Thigh ;
                      loss /= (1.+coef1*rate+coef2*rate*rate) ;
                     }
                   else
                     {
                      losslim = ComputeBremLoss(Z,natom,Thigh,Cuthigh,x) ;  
                      if (&aParticleType==G4Positron::Positron())
                         loss *= ComputePositronCorrFactorLoss(Z,Thigh,Cuthigh) ;   
                      rate = Cut/TotalEnergy ;
                      loss = losslim*(1.+coef1*rate+coef2*rate*rate) ;
                      loss *= Factorhigh*Cut ;
                     }

                  }
                bremloss += natom*loss;
               }

            }

           // now compute the correction due to the LPM effect
           const G4double MigdalConstant = classic_electr_radius*
                                           electron_Compton_length*
                                           electron_Compton_length/pi ;

           const G4double LPMconstant = fine_structure_const*electron_mass_c2*
                                electron_mass_c2/(8.*pi*hbarc) ;
	   
           const G4double kmin = 1.*eV ;
           const G4double klim = 1.*keV ;

           G4double LPMEnergy = LPMconstant*(material->GetRadlen()) ;
           G4double TotalEnergysquare = TotalEnergy*TotalEnergy ;
           G4double LPMGammaEnergyLimit = TotalEnergysquare/LPMEnergy ;

           if(LPMGammaEnergyLimit > klim)
           {
             G4double kmax = G4std::min(Cut,LPMGammaEnergyLimit) ;

             G4double floss = 0. ;
             G4int nmax = 1000 ;
             G4int nn ;
             G4double vmin=log(kmin);
             G4double vmax=log(Cut) ;
             nn = int(nmax*(vmax-vmin)/(log(HighestKineticEnergy)-vmin)) ;
             G4double u,uu,s2lpm,sp,fac,c,v,dv,w ;
             dv = (vmax-vmin)/nn ;
             v = vmin-dv ;
             for(G4int n=0; n<=nn; n++)
             {
               v += dv ;
               u = exp(v) ;
               uu = u*u ;
               if(u<=kmax)
               {
                 s2lpm=LPMEnergy*u/TotalEnergysquare ;
                 sp=uu/(uu+MigdalConstant*TotalEnergysquare*
                           (material->GetElectronDensity())) ;
                 w=s2lpm*(1.+1./sp) ;
                 fac=0.5*(sqrt(w*w+4.*s2lpm)-w)/sp;
                 if(fac>1.)
                 fac=1. ;
               }
               else
               {
                 fac=1. ;
               }

               fac *= uu*u ;

               if((n==0)||(n==nn))
                 c=0.5;
               else
                 c=1.;

               fac *= c ;
               floss += fac ;
             }

             floss *=dv*3./(Cut*Cut*Cut-kmin*kmin*kmin) ;
             if(floss > 1.) floss = 1. ; 

             // correct the loss
             bremloss *= floss ;
          }
  
          if(bremloss < 0.) bremloss = 0. ;
          aVector->PutValue(i,bremloss);  
        }

       theLossTable->insert(aVector);
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//
// METHOD BELOW  FROM STANDARD E_M PROCESSES CODE LEFT BUT AT THE MOMENT NOT USED
//
G4double G4LowEnergyBremsstrahlung::ComputeXYPolynomial(G4double x,  G4double y,
                                             G4int xSize, G4int ySize,
                                             const G4double coeff[])
{
  // Computes the polynomial (1 y y^2 ...) * matrix * (1 x x^2 ...) .
  // xSize and ySize are the dimensions of the matrix,
  // coeff containts the elements, stored row-wise.    

  G4double* a= new G4double[xSize];
  G4int i, j;

  for (i=0; i<xSize; i++) a[i]= 0.0;  

  G4int index= 0; G4double yy= 1.0;
  for (j=0; j<ySize; j++)
    { for (i=0; i<xSize; i++) a[i]+= coeff[index++]*yy;     
      yy*= y;
    }
  
  G4double r= a[0]; G4double xx= x;
  for (i=1; i<xSize; i++) { r+= a[i]*xx; xx*= x;}  
  
  delete[] a;
  return r;
}                                             

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//
// METHOD BELOW  FROM STANDARD E_M PROCESSES CODE LEFT BUT AT THE MOMENT NOT USED
//
G4double G4LowEnergyBremsstrahlung::ComputeBremLoss(G4double Z,G4double natom,
                         G4double T,G4double Cut,G4double x)

// compute loss due to soft brems 
// 'Migdal' version , this is the default in GEANT3 

{
  const G4double beta=0.99,ksi=2.51,ve=0.00004 ;
  const G4double corrfac = classic_electr_radius*electron_Compton_length*electron_Compton_length/pi  ;

  static const G4double
  CMbarn[]= {
    -0.960613e-1, 0.631029e-1,-0.142819e-1, 0.150437e-2,-0.733286e-4, 0.131404e-5,
     0.859343e-1,-0.529023e-1, 0.131899e-1,-0.159201e-2, 0.926958e-4,-0.208439e-5,
    -0.684096e+1, 0.370364e+1,-0.786752e0,  0.822670e-1,-0.424710e-2, 0.867980e-4,
    -0.200856e+1, 0.129573e+1,-0.306533e0,  0.343682e-1,-0.185931e-2, 0.392432e-4, 
     0.127538e+1,-0.515705e0,  0.820644e-1,-0.641997e-2, 0.245913e-3,-0.365789e-5,
     0.115792e0, -0.463143e-1, 0.725442e-2,-0.556266e-3, 0.208049e-4,-0.300895e-6};

  static const G4double
  CPbarn[]= {
    -0.960613e-1, 0.631029e-1,-0.142819e-1, 0.150437e-2,-0.733286e-4, 0.131404e-5,
     0.859343e-1,-0.529023e-1, 0.131899e-1,-0.159201e-2, 0.926958e-4,-0.208439e-5,
    -0.271082e-1, 0.173949e-1,-0.452531e-2, 0.569405e-3,-0.344856e-4, 0.803964e-6,
     0.419855e-2,-0.277188e-2, 0.737658e-3,-0.939463e-4, 0.569748e-5,-0.131737e-6,
    -0.318752e-3, 0.215144e-3,-0.579787e-4, 0.737972e-5,-0.441485e-6, 0.994726e-8,
     0.938233e-5,-0.651642e-5, 0.177303e-5,-0.224680e-6, 0.132080e-7,-0.288593e-9};

  static const G4double
  CCMbarn[]= {
    -0.245667e-3, 0.833406e-4,-0.129217e-4, 0.915099e-6,-0.247179e-7,
     0.147696e-3,-0.498793e-4, 0.402375e-5, 0.989281e-7,-0.133378e-7,
    -0.737702e-2, 0.333057e-2,-0.553141e-3, 0.402464e-4,-0.107977e-5,
    -0.641533e-2, 0.290113e-2,-0.477641e-3, 0.342008e-4,-0.900582e-6,
     0.574303e-5, 0.908521e-4,-0.256900e-4, 0.239921e-5,-0.741271e-7};

  static const G4double
  CCPbarn[]= {
    -0.245667e-3, 0.833406e-4,-0.129217e-4, 0.915099e-6,-0.247179e-7,
     0.147696e-3,-0.498793e-4, 0.402375e-5, 0.989281e-7,-0.133378e-7,
    -0.341260e-4, 0.971711e-5,-0.172031e-6,-0.119455e-6, 0.704166e-8,
     0.341740e-5,-0.775867e-6,-0.653231e-7, 0.225605e-7,-0.114860e-8,
    -0.119391e-6, 0.194885e-7, 0.588959e-8,-0.127589e-8, 0.608247e-10};

  G4double CM[36],CP[36],CCM[25],CCP[25]; //Set the unit: barn
   
  for (G4int i=0; i<36; i++)    { CM[i] = CMbarn[i]*barn;
                                  CP[i] = CPbarn[i]*barn;
                                }
  for (G4int ii=0; ii<25; ii++) { CCM[ii] = CCMbarn[ii]*barn;
                                  CCP[ii] = CCPbarn[ii]*barn;
                                }
  //  -----------------------------------------------------------

  G4double TotalEnergy = T + electron_mass_c2;
  G4double y=log(Cut/(ve*TotalEnergy));
  
  G4double loss;
  
  if (y <= 0.) loss = ComputeXYPolynomial(x, y, 6, 6, CM)
                     + Z * ComputeXYPolynomial(x, y, 5, 5, CCM);
  else         loss = ComputeXYPolynomial(x, y, 6, 6, CP)
                     + Z * ComputeXYPolynomial(x, y, 5, 5, CCP);

  G4double rate = TotalEnergy/Cut ;
  G4double corr = 1./(1.+corrfac*natom*rate*rate) ;

  G4double factor = pow(Cut*corr/T,beta);
  factor *= Z*(Z+ksi)*TotalEnergy*TotalEnergy/(TotalEnergy+electron_mass_c2) ;

  loss   *= factor ;

  return loss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//
// METHOD BELOW  FROM STANDARD E_M PROCESSES LEFT BUT AT THE MOMENT NOT USED
//

G4double G4LowEnergyBremsstrahlung::ComputePositronCorrFactorLoss(
                            G4double Z,G4double KineticEnergy,G4double GammaCut)

//calculates the correction factor for the energy loss due to bremsstrahlung for positrons
//the same correction is in the (discrete) bremsstrahlung 

{
  static const G4double K = 132.9416*eV ;
  static const G4double a1=4.15e-1, a3=2.10e-3, a5=54.0e-5 ;

  G4double x   = log(KineticEnergy/(K*Z*Z)), x2 = x*x, x3 = x2*x;
  G4double eta = 0.5+atan(a1*x+a3*x3+a5*x3*x2)/pi;
  G4double e0  = GammaCut/KineticEnergy;
  
  G4double factor(0.);
  if (e0!=1.0) { factor=log(1.-e0)/eta; factor=exp(factor);}  
  factor = eta*(1.-factor)/e0;

  return factor;
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
   if (theMeanFreePathTable) {theMeanFreePathTable->clearAndDestroy();
                              delete theMeanFreePathTable;
                             }

   G4double NumbOfMaterials = G4Material::GetNumberOfMaterials();
   const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
   G4Material* material;

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
     
     for ( G4int i = 0 ; i < TotBin ; i++ ){
       
       LowEdgeEnergy = ptrVector->GetLowEdgeEnergy( i ) ;
       const G4double BigPath= DBL_MAX;
       G4double SIGMA = 0 ;
       
       for ( G4int k=0 ; k < material->GetNumberOfElements() ; k++ ){ 
	 
	 G4int AtomIndex = (G4int) (*theElementVector)(k)->GetZ();
	 const G4FirstLevel* oneAtomCS
	   = (*theCrossSectionTable)[ZNumVec->index(AtomIndex)];
	 
	 G4double interCrsSec = util.DataLogInterpolation(LowEdgeEnergy, 
							  (*(*oneAtomCS)[0]), 
							  (*(*oneAtomCS)[1]))*barn;

	 SIGMA += theAtomNumDensityVector[k]*interCrsSec;
       }       
       
       Value = SIGMA<=0.0 ? BigPath : 1./SIGMA ;
       ptrVector->PutValue( i , Value ) ;

     }
     
     theMeanFreePathTable->insert( ptrVector );
     
     // Compute the PartialSumSigma table at a given fixed energy
     ComputePartialSumSigma(FixedEnergy, material) ;       
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//
// METHOD BELOW  FROM STANDARD E_M PROCESSES CODE MODIFIED TO USE 
// LIVERMORE DATA (using log-log interpolation as reported in stepanek paper)
//
void G4LowEnergyBremsstrahlung::ComputePartialSumSigma(G4double KineticEnergy,
						       const G4Material* aMaterial)

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
     const G4FirstLevel* oneAtomCS
       = (*theCrossSectionTable)[ZNumVec->index(AtomIndex)];
     
     G4double interCrsSec = util.DataLogInterpolation(KineticEnergy, 
						      (*(*oneAtomCS)[0]), 
						      (*(*oneAtomCS)[1]))*barn;
     
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
  // MIGDAL constant and LPM effect LEFT FROM STANDARD PROCESS
  //     
  
  const G4double MigdalConstant = classic_electr_radius
    *electron_Compton_length
    *electron_Compton_length/pi;


  const G4double LPMconstant = fine_structure_const*electron_mass_c2*
    electron_mass_c2/(8.*pi*hbarc) ;
  
  aParticleChange.Initialize(trackData);
  
  G4Material* aMaterial=trackData.GetMaterial() ;
  
  G4double LPMEnergy = LPMconstant*(aMaterial->GetRadlen()) ;
  
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
    
    aParticleChange.SetMomentumChange( ElectDirection );
    aParticleChange.SetEnergyChange( ElectKinEn );
    aParticleChange.SetLocalEnergyDeposit (0.); 
    aParticleChange.SetNumberOfSecondaries(0);
    
    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
  }
  
  // select randomly one element constituing the material  
  G4Element* anElement = SelectRandomAtom(aMaterial);
  
  // limits of the energy sampling
  G4double TotalEnergy = ElectKinEn + electron_mass_c2;
  G4double TotalEnergysquare = TotalEnergy*TotalEnergy ;
  G4double LPMGammaEnergyLimit = TotalEnergysquare/LPMEnergy ;
  
  //
  // The emitted gamma energy is from EEDL data fitted with A/E+B function.
  // Original formula A/E+B+C*E and sampling methods are reported by  J. Stepanek 
  // formula has been modified by A. Forti and S. Giani. 

  // 
  //  sample the energy of the emitted gamma for electron kinetic energy
  //  
  G4double p1 = 0, p2 = 0;
  G4double coeffA = 0, coeffB = 0;
  G4int AtomicNum = (G4int) anElement->GetZ();
  coeffA = ComputeA(AtomicNum, ElectKinEn);
  coeffB = ComputeB(AtomicNum, ElectKinEn);
  
  p1 = coeffA*log(ElectKinEn/lowEnergyCut);
  p2 = coeffB*(ElectKinEn - lowEnergyCut); 
  
  G4double IntegrProb = p1+p2;  
  G4double R1 = G4UniformRand()*IntegrProb;
  
  G4double GammaEnergy;

  if(R1 <= p1){ 
    
    G4double R2 = G4UniformRand();
    GammaEnergy = ElectKinEn*pow((lowEnergyCut/ElectKinEn),R2);   
  }
  else if ((p1 < R1) && (R1 <= p1+p2)){
    
    G4double R2 = G4UniformRand();
    GammaEnergy = ElectKinEn - R2*(ElectKinEn - lowEnergyCut);
  }
  
  // now comes the supression due to the LPM effect (gamma production suppression
  // due to the multiple scattering of the electron) SEE ABOVE

  if(GammaEnergy < LPMGammaEnergyLimit){
    
    G4double S2LPM = LPMEnergy*GammaEnergy/TotalEnergysquare ;
    G4double Spol  = GammaEnergy*GammaEnergy/(GammaEnergy*GammaEnergy +
					      MigdalConstant*(aMaterial->GetElectronDensity())*
					      TotalEnergysquare) ;
    G4double w = S2LPM*(1.+1./Spol) ;
    G4double Supr = 0.5*(sqrt(w*w+4.*S2LPM)-w)/Spol ;
    
    if (G4UniformRand() > Supr )
      GammaEnergy = 0. ;
  }
  
  //protection: DO NOT PRODUCE a gamma with energy 0. !
  if (GammaEnergy <= 0.){
 
    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
  }
  
  //**********************//
  // Angular distribution //
  //**********************//
  
  //  angles of the emitted gamma. ( Z - axis along the parent particle)
  //  universal distribution suggested by L. Urban (Geant3 manual (1993) Phys211),
  //  derived from Tsai distribution (Rev Mod Phys 49,421(1977))

  if(GammaEnergy > CutForLowEnergySecondaryPhotons){
    
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
    
    if (NewKinEnergy > 0.){
      
      aParticleChange.SetMomentumChange( ElectDirection );
      aParticleChange.SetEnergyChange( NewKinEnergy );

      if(GammaEnergy <  GammaEnergyCut){

	 aParticleChange.SetLocalEnergyDeposit(GammaEnergy); 
      }
      else{

	// create G4DynamicParticle object for the Gamma 
	G4DynamicParticle* aGamma= new G4DynamicParticle (G4Gamma::Gamma(),
							  GammaDirection, GammaEnergy);
	
	aParticleChange.SetNumberOfSecondaries(1);
	aParticleChange.AddSecondary(aGamma); 
	aParticleChange.SetLocalEnergyDeposit(0);
      }
    } 
    else{
      
      aParticleChange.SetEnergyChange( 0. );
      aParticleChange.SetLocalEnergyDeposit (0.);
      if (charge<0.){
	
	aParticleChange.SetStatusChange(fStopAndKill);
      }
      else{
	
	aParticleChange.SetStatusChange(fStopButAlive);
      }    
    }   
  }
  else{

    aParticleChange.SetNumberOfSecondaries(0);
  }
    
#ifdef G4VERBOSE
  if(verboseLevel > 15){

    G4cout<<"LE Bremsstrahlung PostStepDoIt"<<G4endl;
  }
#endif
  return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
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
  G4String comments = "Total cross sections from EEDL database";
           comments += "Good description from 1 eV to 100 GeV.\n";
           comments += "Gamma energy sampled from a parametrised formula.";
                     
	   G4cout << G4endl << GetProcessName() << ":  " << comments<<G4endl;

}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....














