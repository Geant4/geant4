// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4hLowEnergyIonisation.cc,v 1.1 1999-07-20 17:36:50 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4hLowEnergyIonisation physics process -------
//                by Vladimir Ivanchenko, 14 July 1999 
// **************************************************************
// It is the first implementation of the NEW IONISATION PROCESS.
// It calculates the ionisation of charged hadrons.
// **************************************************************
// --------------------------------------------------------------


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
#include "G4hLowEnergyIonisation.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4hLowEnergyIonisation::G4hLowEnergyIonisation(const G4String& processName)
   : G4hIonisation(processName),
     theMeanFreePathTable(NULL),
     LowestKineticEnergy(0.10*keV),
     HighestKineticEnergy(100.*TeV),
     ZieglerLowEnergy(1.*keV),
     ZieglerHighEnergy(2.*MeV),
     DEDXtable("Ziegler1977H"),
     TotBin(200),
     theProton (G4Proton::Proton()),
     theAntiProton (G4AntiProton::AntiProton()),
     theElectron ( G4Electron::Electron() ),
     twoln10(2.*log(10.)),
     Factor(twopi_mc2_rcl2),
     bg2lim(0.0169), 
     taulim(8.4146e-3),
     RateMass(electron_mass_c2/proton_mass_c2),
     ProtonMassAMU(1.007276),
     ZieglerFactor(eV*cm2*1.0e-15) 
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
     
G4hLowEnergyIonisation::~G4hLowEnergyIonisation() 
{
     if (theMeanFreePathTable) {
        theMeanFreePathTable->clearAndDestroy();
        delete theMeanFreePathTable;
     }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::SetDEDXTableName(const G4String& dedxTable)
{
  if(dedxTable == "Ziegler1977H") { 
     return;

  } else if(dedxTable == "Ziegler1977He") {
     DEDXtable = "Ziegler1977He";

  } else {
    cout << "G4hLowEnergyIonisation Warning: There is no table with the name ="
         << dedxTable; 
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::BuildLossTable(const G4ParticleDefinition& aParticleType)
{
  // cuts for  electron ....................
  DeltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;

  G4double LowEdgeEnergy , ionloss, ionlossBB;
  G4double Eexc, tau, paramA, paramB;
  static const G4MaterialTable* theMaterialTable=
                                   G4Material::GetMaterialTable();

  //  create table

    cout << "\n ### Construct new tables for proton energy loss" 
         << "\n     ZieglerHighEnergy = " << ZieglerHighEnergy 
         << " MeV" << " \n";

  G4int numOfMaterials = theMaterialTable->length();

  if ( theLossTable) {
     theLossTable->clearAndDestroy();
     delete theLossTable;
  }
  theLossTable = new G4PhysicsTable(numOfMaterials);

  //  loop for materials

  for (G4int J=0; J<numOfMaterials; J++)
  {

    // create physics vector and fill it

    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                    LowestKineticEnergy, HighestKineticEnergy, TotBin);
 
    // get material parameters needed for the energy loss calculation
   
    const G4Material* material= (*theMaterialTable)[J];

    // get  electron cut in kin. energy for the material

    DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[J] ;

    // define the boundary energies between using Stopping Power Tables and
    // Geant4 methods which includes delta-ray production

    G4double tauHighEnergy = ZieglerHighEnergy ;
    G4double tauLowEnergy = ZieglerLowEnergy ;

    // define constants A and B for this material  

    tau    = tauLowEnergy/keV ;                     // tau in keV
    paramA = GetZieglerLoss(material, tau)/sqrt(tau) ; 
    tau    = tauHighEnergy/keV ;                    // tau in keV
    ionloss = GetZieglerLoss(material, tau) ; 
    tau    = tauHighEnergy/proton_mass_c2 ;         // tau is relative energy
    ionlossBB = GetBetheBlochLoss(material, tau) ; 
    paramB =  ionloss/ionlossBB - 1.0 ; 

    cout << "\n ### " << material->GetName() << ": loss = " 
         << ionloss << "; paramB = " << paramB 
         << "; tauHighEnergy = " << tauHighEnergy << " MeV \n";

    // now comes the loop for the kinetic energy values

    for (G4int i = 0 ; i < TotBin ; i++)
    {
      LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;

      if ( LowEdgeEnergy < tauHighEnergy ) {
      //  low energy part , parametrized energy loss formulae
     
        tau = LowEdgeEnergy/keV ;                        // tau in keV

        if ( LowEdgeEnergy < tauLowEnergy ) {

        // The model of free electron gas
          ionloss = GetFreeElectronGasLoss(paramA, tau) ;

        } else {
	  // Ziegler parametrisation
          ionloss = GetZieglerLoss(material, tau) ; 
	}
      } else {

	// high energy part , Bethe-Bloch formula
        tau = LowEdgeEnergy/proton_mass_c2 ;           // tau is relative energy   
        ionloss = GetBetheBlochLoss(material, tau) ; 
        ionloss *= (1.0 + paramB*tauHighEnergy/LowEdgeEnergy) ;
      }

      // now put the loss into the vector
      //IV  cout << "  E = " << LowEdgeEnergy/MeV << " MeV; ion = " << ionloss << endl ;

      aVector->PutValue(i,ionloss) ;
    }

    // insert vector for this material into the table
    theLossTable->insert(aVector) ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetZieglerLoss(const G4Material* material, 
                                                   G4double tau)
{
    G4double ionloss, ion;

    // get elements in the actual material,
    const G4ElementVector* theElementVector=
                   material->GetElementVector() ;
    const G4double* theAtomicNumDensityVector=
                   material->GetAtomicNumDensityVector() ;
    const G4int NumberOfElements=
                   material->GetNumberOfElements() ;

    ionloss = 0.0 ;

    //  loop for the elements in the material
    for (G4int iel=0; iel<NumberOfElements; iel++)
    {
      const G4Element* element = (*theElementVector)(iel) ;
      G4double A1 = ProtonMassAMU ;
      G4double Z2 = element->GetZ() ;
      G4double A2 = element->GetA() ;
      G4int iz = int(Z2) ;
      if( iz <= 0 ) iz = 1 ;
      if( iz > 92 ) iz = 92 ; 
    
      // Electronic Stopping Power 
      // Choose the parametrisation
 
      if(DEDXtable == "Ziegler1977H") { 
        ion = GetStoppingPower1977H(iz, tau/ProtonMassAMU) ; 
        ion *= theAtomicNumDensityVector[iel]*ZieglerFactor ;

      } else if(DEDXtable == "Ziegler1977He") {
	// This must be modified!!!
        ion = GetStoppingPower1977He(iz, tau) ; 
        ion *= theAtomicNumDensityVector[iel]*ZieglerFactor ;
      }
      // Nuclear Stopping Power
      G4double ionn = GetStoppingPower1977n(1.0, Z2, A1, A2, tau) ;
        ion += ionn*theAtomicNumDensityVector[iel]*ZieglerFactor ;

      ionloss += ion ;
    }

    // Correction due to delta-electrons energy loss 
    // Bethe-Bloch formulae is used
    G4double tauU = tau * keV / proton_mass_c2 ;
    ionloss -= GetDeltaRaysEnergy(material, tauU) ;

  if ( ionloss <= 0.) ionloss = 0. ;

  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetBetheBlochLoss(const G4Material* material, 
                                                   G4double tau)
{
  G4double ionloss ;
  G4double taul = material->GetIonisation()->GetTaul() ;

  if ( tau < taul ) {

  //  low energy part , parametrized L.Urban energy loss formulae

    const G4ElementVector* theElementVector=
                   material->GetElementVector() ;
    const G4double* theAtomicNumDensityVector=
                   material->GetAtomicNumDensityVector() ;
    const G4int NumberOfElements=
                   material->GetNumberOfElements() ;

    ionloss = 0. ;
        
    //  loop for the elements in the material
    for (G4int iel=0; iel<NumberOfElements; iel++)
      {
      const G4Element* element = (*theElementVector)(iel) ;
      ionloss += GetUrbanModel(element, tau) * theAtomicNumDensityVector[iel] ;
      }

  } else {
    // Standard Bethe-Bloch formulae

    // some local variables 

    G4double gamma,bg2,beta2,Tmax,rcut,x,delta,sh ;
    G4double ElectronDensity = material->GetElectronDensity();
    G4double Eexc = material->GetIonisation()->GetMeanExcitationEnergy();
    G4double Eexc2 = Eexc*Eexc ;
    G4double Cden = material->GetIonisation()->GetCdensity();
    G4double Mden = material->GetIonisation()->GetMdensity();
    G4double Aden = material->GetIonisation()->GetAdensity();
    G4double X0den = material->GetIonisation()->GetX0density();
    G4double X1den = material->GetIonisation()->GetX1density();
    G4double* ShellCorrectionVector;
    ShellCorrectionVector = material->GetIonisation()->
                                          GetShellCorrectionVector();

    gamma = tau + 1.0 ;
    bg2 = tau*(tau+2.0) ;
    beta2 = bg2/(gamma*gamma) ;
    Tmax = 2.*electron_mass_c2*bg2/(1.+2.*gamma*RateMass+RateMass*RateMass) ;

    if ( DeltaCutInKineticEnergyNow < Tmax)
      rcut = DeltaCutInKineticEnergyNow/Tmax ;
    else
      rcut = 1.;

    ionloss = log(2.*electron_mass_c2*bg2*Tmax/Eexc2)+log(rcut)-(1.+rcut)*beta2 ;

    // density correction 

    x = log(bg2)/twoln10 ;
    if ( x < X0den )
        delta = 0. ;
    else 
        {
          delta = twoln10*x - Cden ;
          if ( x < X1den )
            delta += Aden*pow((X1den-x),Mden) ;
        } 

    // shell correction 
         
    if ( bg2 > bg2lim ) {
        sh = 0. ;      
        x = 1. ;
        for (G4int k=0; k<=2; k++) {
            x *= bg2 ;
            sh += ShellCorrectionVector[k]/x;
          }
    } else {
        sh = 0. ;      
        x = 1. ;
        for (G4int k=0; k<=2; k++) {
            x *= bg2lim ;
            sh += ShellCorrectionVector[k]/x;
          }
        sh *= log(tau/taul)/log(taulim/taul) ;     
        }

    // now you can compute the total ionization loss

  //IV  cout << "ionloss = " << ionloss << "; delta = " << delta << "; sh = " << sh << endl;
  //IV cout << "Factor = " << Factor << "; Edensity = " << ElectronDensity << "; beta2 = " << beta2 << endl;

    ionloss -= delta + sh ;
    ionloss *= Factor*ElectronDensity/beta2 ;
  }
  if ( ionloss <= 0.) ionloss = 0. ;

  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetDeltaRaysEnergy(const G4Material* material, 
                                                   G4double tau)
{
    G4double ionloss = 0.0 ;

    // some local variables 

    G4double gamma,bg2,beta2,Tmax,x ;
    G4double ElectronDensity = material->GetElectronDensity();

    gamma = tau + 1.0 ;
    bg2 = tau*(tau+2.0) ;
    beta2 = bg2/(gamma*gamma) ;
    Tmax = 2.*electron_mass_c2*bg2/(1.+2.*gamma*RateMass+RateMass*RateMass) ;

    if ( DeltaCutInKineticEnergyNow < Tmax) {
      x = DeltaCutInKineticEnergyNow / Tmax ;
      ionloss = (x - log(x) - 1.0 ) * Factor * ElectronDensity / beta2 ;
    }
    return ionloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetUrbanModel(const G4Element* element, 
                                                     G4double tau)
{
  // Parametrisation of low energy protons energy loss by L.Urban 
  // for the class G4hIonisation

  G4double ionloss = 0. ;

  if ( tau < element->GetIonisation()->GetTau0()) {  
    ionloss = ( element->GetIonisation()->GetAlow()*sqrt(tau)
              + element->GetIonisation()->GetBlow()*tau) ;
  }
  else {
    ionloss = element->GetIonisation()->GetClow()/sqrt(tau) ;
  }

  return ionloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetStoppingPower1977H(G4int iz, G4double E)
{
  G4double ionloss ;
  G4int i = iz-1 ;  // index of atom

  // The data and the fit from: 
  // H.H.Andersen & J.F.Ziegler Hydrogen Stopping Powers and
  // Ranges in All Elements, Vol.3, Pergamon Press, 1977
 
  static G4double A[92][12] = {
1.262,1.440,  242.6,12000.,0.115900,0.0005099,54360.0,-5.0520,2.0490,-0.30440,0.019660,-0.0004659,
1.229,1.397,  484.5,5873.0,0.052250,0.0010200,24510.0,-2.1580,0.8278,-0.11720,0.007259,-0.000166,
1.411,1.600,  725.6,3013.0,0.045780,0.0015300,21470.0,-0.5831,0.5620,-0.11830,0.009298,-0.0002498,
2.248,2.590,  966.0, 153.8,0.034750,0.0020390,16300.0, 0.2779,0.1745,-0.05684,0.005155,-0.0001488,
2.474,2.815, 1206.0,1060.0,0.028550,0.0025490,13450.0,-2.4450,1.2830,-0.22050,0.015600,-0.0003930,
2.631,2.989, 1445.0, 957.2,0.028190,0.0030590,13220.0,-4.3800,2.0440,-0.32830,0.022210,-0.0005417,
2.954,3.350, 1683.0,1900.0,0.025130,0.0035690,11790.0,-5.0540,2.3250,-0.37130,0.025060,-0.0006109,
2.652,3.000, 1920.0,2000.0,0.022300,0.0040790,10460.0,-6.7340,3.0190,-0.47480,0.031710,-0.0007669,
2.085,2.352, 2157.0,2634.0,0.018160,0.0045890,8517.0, -5.5710,2.4490,-0.37810,0.024830,-0.0005919,
1.951,2.199, 2393.0,2699.0,0.015680,0.0050990,7353.0, -4.4080,1.8790,-0.28140,0.017960,-0.0004168,
2.542,2.869, 2628.0,1854.0,0.014720,0.0056090,6905.0, -4.9590,2.0730,-0.30540,0.019210,-0.0004403,
3.792,4.293, 2862.0,1009.0,0.013970,0.0061180,6551.0, -5.5100,2.2660,-0.32950,0.020470,-0.0004637,
4.154,4.739, 2766.0, 164.5,0.020230,0.0066280,6309.0, -6.0610,2.4600,-0.35350,0.021730,-0.0004871,
4.150,4.700, 3329.0, 550.0,0.013210,0.0071380,6194.0, -6.2940,2.5380,-0.36280,0.022200,-0.0004956,
3.232,3.647, 3561.0,1560.0,0.012670,0.0076480,5942.0, -6.5270,2.6160,-0.37210,0.022670,-0.0005040,
3.447,3.891, 3792.0,1219.0,0.012110,0.0081580,5678.0, -6.7610,2.6940,-0.38140,0.023140,-0.0005125,
5.047,5.714, 4023.0, 878.6,0.011780,0.0086680,5524.0, -6.9940,2.7730,-0.39070,0.023610,-0.0005209,
5.731,6.500, 4253.0, 530.0,0.011230,0.0091780,5268.0, -7.2270,2.8510,-0.40000,0.024070,-0.0005294,
5.151,5.833, 4482.0, 545.7,0.011290,0.0096870,5295.0, -7.4400,2.9230,-0.40940,0.024620,-0.0005411,
5.521,6.252, 4710.0, 553.3,0.011120,0.0102000,5214.0, -7.6530,2.9950,-0.41870,0.025160,-0.0005529,
5.201,5.884, 4938.0, 560.9,0.009995,0.0107100,4688.0, -8.0120,3.1230,-0.43500,0.026050,-0.0005707,
4.862,5.496, 5165.0, 568.5,0.009474,0.0112200,4443.0, -8.3710,3.2510,-0.45130,0.026940,-0.0005886,
4.480,5.055, 5391.0, 952.3,0.009117,0.0117300,4276.0, -8.7310,3.3790,-0.46760,0.027830,-0.0006064,
3.983,4.489, 5616.0,1336.0,0.008413,0.0122400,3946.0, -9.0900,3.5070,-0.48380,0.028720,-0.0006243,
3.469,3.907, 5725.0,1461.0,0.008829,0.0127500,3785.0, -9.4490,3.6350,-0.50010,0.029610,-0.0006421,
3.519,3.963, 6065.0,1243.0,0.007782,0.0132600,3650.0, -9.8090,3.7630,-0.51640,0.030500,-0.0006600,
3.140,3.535, 6288.0,1372.0,0.007361,0.0137700,3453.0,-10.1700,3.8910,-0.53270,0.031390,-0.0006779,
3.553,4.004, 6205.0, 555.1,0.008763,0.0142800,3297.0,-10.5300,4.0190,-0.54900,0.032290,-0.0006957,
3.696,4.175, 4673.0, 387.8,0.021880,0.0147900,3174.0,-11.1800,4.2520,-0.57910,0.033990,-0.0007314,
4.210,4.750, 6953.0, 295.2,0.006809,0.0153000,3194.0,-11.5700,4.3940,-0.59800,0.035060,-0.0007537,
5.041,5.697, 7173.0, 202.6,0.006725,0.0158100,3154.0,-11.9500,4.5370,-0.61690,0.036130,-0.0007759,
5.554,6.300, 6496.0, 110.0,0.009689,0.0163200,3097.0,-12.3400,4.6800,-0.63580,0.037210,-0.0007981,
5.323,6.012, 7611.0, 292.5,0.006447,0.0168300,3024.0,-12.7200,4.8230,-0.65470,0.038280,-0.0008203,
5.847,6.656, 7395.0, 117.5,0.007684,0.0173400,3006.0,-13.1100,4.9650,-0.67350,0.039350,-0.0008425,
5.611,6.335, 8046.0, 365.2,0.006244,0.0178500,2928.0,-13.4000,5.0830,-0.69060,0.040420,-0.0008675,
6.411,7.250, 8262.0, 220.0,0.006087,0.0183600,2855.0,-13.6900,5.2000,-0.70760,0.041500,-0.0008925,
5.694,6.429, 8478.0, 292.9,0.006087,0.0188600,2855.0,-13.9200,5.2660,-0.71400,0.041730,-0.0008943,
6.339,7.159, 8693.0, 330.3,0.006003,0.0193700,2815.0,-14.1400,5.3310,-0.72050,0.041960,-0.0008962,
6.407,7.234, 8907.0, 367.8,0.005889,0.0198800,2762.0,-14.3600,5.3970,-0.72690,0.042190,-0.0008980,
6.734,7.603, 9120.0, 405.2,0.005765,0.0203900,2704.0,-14.5900,5.4630,-0.73330,0.042420,-0.0008998,
6.902,7.791, 9333.0, 442.7,0.005587,0.0209000,2621.0,-16.2200,6.0940,-0.82250,0.047910,-0.0010240,
6.425,7.248, 9545.0, 480.2,0.005367,0.0214100,2517.0,-17.8500,6.7250,-0.91160,0.053390,-0.0011480,
6.799,7.671, 9756.0, 517.6,0.005315,0.0219200,2493.0,-17.9600,6.7520,-0.91350,0.053410,-0.001147,
6.108,6.887, 9966.0, 555.1,0.005151,0.0224300,2416.0,-18.0700,6.7790,-0.91540,0.053420,-0.0011450,
5.924,6.677,10180.0, 592.5,0.004919,0.0229400,2307.0,-18.1800,6.8060,-0.91730,0.053430,-0.0011430,
5.238,5.900,10380.0, 630.0,0.004758,0.0234500,2231.0,-18.2800,6.8330,-0.91920,0.053450,-0.0011420,
5.623,6.354, 7160.0, 337.6,0.013940,0.0239600,2193.0,-18.3900,6.8600,-0.92110,0.053460,-0.0011400,
5.814,6.554,10800.0, 355.5,0.004626,0.0244700,2170.0,-18.6200,6.9150,-0.92430,0.053400,-0.0011340,
6.230,7.024,11010.0, 370.9,0.004540,0.0249800,2129.0,-18.8500,6.9690,-0.92750,0.053350,-0.0011270,
6.410,7.227,11210.0, 386.4,0.004474,0.0254900,2099.0,-19.0700,7.0240,-0.93080,0.053290,-0.0011210,
7.500,8.480, 8608.0, 348.0,0.009074,0.0260000,2069.0,-19.5700,7.2250,-0.96030,0.055180,-0.0011650,
6.979,7.871,11620.0, 392.4,0.004402,0.0265100,2065.0,-20.0700,7.4260,-0.98990,0.057070,-0.0012090,
7.725,8.716,11830.0, 394.8,0.004376,0.0270200,2052.0,-20.5600,7.6270,-1.01900,0.058960,-0.0012540,
8.231,9.289,12030.0, 397.3,0.004384,0.0275300,2056.0,-21.0600,7.8280,-1.04900,0.060850,-0.0012980,
7.287,8.218,12230.0, 399.7,0.004447,0.0280400,2086.0,-20.4000,7.5400,-1.00400,0.057820,-0.0012240,
7.899,8.911,12430.0, 402.1,0.004511,0.0285500,2116.0,-19.7400,7.2520,-0.95880,0.054790,-0.0011510,
8.041,9.071,12630.0, 404.5,0.004540,0.0290600,2129.0,-19.0800,6.9640,-0.91360,0.051760,-0.0010770,
7.489,8.444,12830.0, 406.9,0.004420,0.0295700,2073.0,-18.4300,6.6770,-0.86840,0.048720,-0.0010030,
7.291,8.219,13030.0, 409.3,0.004298,0.0300800,2016.0,-17.7700,6.3890,-0.82330,0.045690,-0.0009292,
7.098,8.000,13230.0, 411.8,0.004182,0.0305900,1962.0,-17.1100,6.1010,-0.77810,0.042660,-0.0008553,
6.910,7.786,13430.0, 414.2,0.004058,0.0311000,1903.0,-16.4500,5.8130,-0.73300,0.039630,-0.0007815,
6.728,7.580,13620.0, 416.6,0.003976,0.0316100,1865.0,-15.7900,5.5260,-0.68780,0.036600,-0.0007077,
6.551,7.380,13820.0, 419.0,0.003877,0.0321200,1819.0,-15.1300,5.2380,-0.64260,0.033570,-0.0006339,
6.739,7.592,14020.0, 421.4,0.003863,0.0326300,1812.0,-14.4700,4.9500,-0.59750,0.030530,-0.0005601,
6.212,6.996,14210.0, 423.9,0.003725,0.0331400,1747.0,-14.5600,4.9840,-0.60220,0.030820,-0.0005668,
5.517,6.210,14400.0, 426.3,0.003632,0.0336500,1703.0,-14.6500,5.0180,-0.60690,0.031110,-0.0005734,
5.219,5.874,14600.0, 428.7,0.003498,0.0341600,1640.0,-14.7400,5.0510,-0.61170,0.031410,-0.0005801,
5.071,5.706,14790.0, 433.0,0.003405,0.0346700,1597.0,-14.8300,5.0850,-0.61640,0.031700,-0.0005867,
4.926,5.542,14980.0, 433.5,0.003342,0.0351800,1567.0,-14.9100,5.1190,-0.62110,0.031990,-0.0005933,
4.787,5.386,15170.0, 435.9,0.003292,0.0356900,1544.0,-15.0000,5.1530,-0.62580,0.032280,-0.0006000,
4.893,5.505,15360.0, 438.4,0.003243,0.0362000,1521.0,-15.0900,5.1860,-0.63050,0.032570,-0.0006066,
5.028,5.657,15550.0, 440.8,0.003195,0.0367100,1499.0,-15.1800,5.2200,-0.63530,0.032860,-0.0006133,
4.738,5.329,15740.0, 443.2,0.003186,0.0372200,1494.0,-15.2700,5.2540,-0.64000,0.033150,-0.0006199,
4.574,5.144,15930.0, 442.4,0.003144,0.0377300,1475.0,-15.6700,5.3920,-0.65770,0.034180,-0.0006426,
5.200,5.851,16120.0, 441.6,0.003122,0.0382400,1464.0,-16.0700,5.5290,-0.67550,0.035210,-0.0006654,
5.070,5.704,16300.0, 440.9,0.003082,0.0387500,1446.0,-16.4700,5.6670,-0.69320,0.036240,-0.0006881,
4.945,5.563,16490.0, 440.1,0.002965,0.0392600,1390.0,-16.8800,5.8040,-0.71100,0.037270,-0.0007109,
4.476,5.034,16670.0, 439.3,0.002871,0.0397700,1347.0,-17.2800,5.9420,-0.72870,0.038300,-0.0007336,
4.856,5.460,18320.0, 438.5,0.002542,0.0402800,1354.0,-17.0200,5.8460,-0.71490,0.037400,-0.0007114,
4.308,4.843,17040.0, 487.8,0.002882,0.0407900,1352.0,-17.8400,6.1830,-0.76590,0.040760,-0.0007925,
4.723,5.311,17220.0, 537.0,0.002913,0.0413000,1366.0,-18.6600,6.5200,-0.81690,0.044110,-0.0008737,
5.319,5.982,17400.0, 586.3,0.002871,0.0418100,1347.0,-19.4800,6.8570,-0.86780,0.047470,-0.0009548,
5.956,6.700,17800.0, 677.0,0.002660,0.0423200,1336.0,-19.5500,6.8710,-0.86860,0.047480,-0.0009544,
6.158,6.928,17770.0, 586.3,0.002812,0.0428300,1319.0,-19.6200,6.8840,-0.86940,0.047480,-0.0009540,
6.204,6.979,17950.0, 586.3,0.002776,0.0433400,1302.0,-19.6900,6.8980,-0.87020,0.047490,-0.0009536,
6.181,6.954,18120.0, 586.3,0.002748,0.0438500,1289.0,-19.7600,6.9120,-0.87100,0.047490,-0.0009532,
6.949,7.820,18300.0, 586.3,0.002737,0.0443600,1284.0,-19.8300,6.9260,-0.87180,0.047500,-0.0009528,
7.506,8.448,18480.0, 586.3,0.002727,0.0448700,1279.0,-19.9000,6.9400,-0.87260,0.047510,-0.0009524,
7.649,8.609,18660.0, 586.3,0.002697,0.0453800,1265.0,-19.9700,6.9530,-0.87330,0.047510,-0.0009520,
7.710,8.679,18830.0, 586.3,0.002641,0.0458900,1239.0,-20.0400,6.9670,-0.87410,0.047520,-0.0009516,
7.407,8.336,19010.0, 586.3,0.002603,0.0464000,1221.0,-20.1100,6.9810,-0.87490,0.047520,-0.0009512,
7.290,8.204,19180.0, 586.3,0.002573,0.0469100,1207.0,-20.1800,6.9950,-0.87570,0.047530,-0.0009508
  };

  if ( E < 10.0 ) {
     ionloss = A[i][0] * sqrt(E) ;

  } else if ( E < 1000.0 ) {
     G4double Slow  = A[i][1] * pow(E, 0.45) ;
     G4double Shigh = log( 1.0 + A[i][3]/E + A[i][4]*E ) * A[i][2]/E ;
     ionloss = Slow*Shigh / (Slow + Shigh) ; 

  } else {
     G4double le = log(E) ;
     G4double gam = 1.0 + E * ProtonMassAMU * keV /proton_mass_c2 ;
     G4double beta2 = 1.0 - 1.0/ (gam*gam) ;
     ionloss = ( log(A[i][6]*beta2/(1.0 - beta2)) - beta2 -
               A[i][7] - A[i][8]*le - A[i][9]*le*le - A[i][10]*le*le*le -
               A[i][11]*le*le*le*le ) * A[i][5]/beta2 ;
  }

  if ( ionloss <= 0.) ionloss = 0. ;

  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetStoppingPower1977He(G4int iz, G4double E)
{
  G4double ionloss ;
  G4int i = iz-1 ;  // index of atom

  // The He4 data and the fit from: 
  // J.F.Ziegler, Helium Stopping Powers and
  // Ranges in All Elemental Matter, Vol.4, Pergamon Press, 1977
 
  static G4double A[92][9] = {
0.9661,0.4126,  6.92, 8.831,  2.582, 2.371, 0.5462,   -0.07932,-0.006853,
2.027, 0.2931, 26.34, 6.66,   0.3409,2.809, 0.4847,   -0.08756,-0.007281,
1.42,  0.49,   12.25,32.,     9.161, 3.095, 0.4434,   -0.09259,-0.007459,
2.206, 0.51,   15.32, 0.25,   8.995, 3.28,  0.4188,   -0.09564,-0.007604,
3.691, 0.4128, 18.48,50.72,   9.,    3.426, 0.4,      -0.09796,-0.007715,
4.232, 0.3877, 22.99,35.,     7.993, 3.588, 0.3921,   -0.09935,-0.007804,
2.51,  0.4752, 38.26,13.02,   1.892, 3.759, 0.4094,   -0.09646,-0.007661,
1.766, 0.5261, 37.11,15.24,   2.804, 3.782, 0.3734,   -0.1011, -0.007874,
1.533, 0.531,  40.44,18.41,   2.718, 3.816, 0.3504,   -0.1046, -0.008074,
1.183, 0.55,   39.83,17.49,   4.001, 3.863, 0.3342,   -0.1072, -0.008231,
9.894, 0.3081, 23.65, 0.384, 92.93,  3.898, 0.3191,   -0.1086, -0.008271,
4.3,   0.47,   34.3,  3.3,   12.74,  3.961, 0.314,    -0.1091, -0.008297,
2.5,   0.625,  45.7,  0.1,    4.359, 4.024, 0.3113,   -0.1093, -0.008306,
2.1,   0.65,   49.34, 1.788,  4.133, 4.077, 0.3074,   -0.1089, -0.008219,
1.729, 0.6562, 53.41, 2.405,  3.845, 4.124, 0.3023,   -0.1094, -0.00824,
1.402, 0.6791, 58.98, 3.528,  3.211, 4.164, 0.2964,   -0.1101, -0.008267,
1.117, 0.7044, 69.69, 3.705,  2.156, 4.21,  0.2936,   -0.1103, -0.00827,
0.9172,0.724,  79.44, 3.648,  1.646, 4.261, 0.2994,   -0.1085, -0.008145,
8.554, 0.3817, 83.61,11.84,   1.875, 4.3,   0.2903,   -0.1103, -0.008259,
6.297, 0.4622, 65.39,10.14,   5.036, 4.334, 0.2897,   -0.1102, -0.008245,
5.307, 0.4918, 61.74,12.4,    6.665, 4.327, 0.2707,   -0.1127, -0.00837,
4.71,  0.5087, 65.28, 8.806,  5.948, 4.34,  0.2618,   -0.1138, -0.00842,
6.151, 0.4524, 83.,  18.31,   2.71,  4.361, 0.2559,   -0.1145, -0.008447,
6.57,  0.4322, 84.76,15.53,   2.779, 4.349, 0.24,     -0.1166, -0.00855,
5.738, 0.4492, 84.61,14.18,   3.101, 4.362, 0.2327,   -0.1174, -0.008588,
5.013, 0.4707, 85.58,16.55,   3.211, 4.375, 0.2253,   -0.1185, -0.008648,
4.32,  0.4947, 76.14,10.85,   5.441, 4.362, 0.2069,   -0.1214, -0.008815,
4.652, 0.4571, 80.73,22.,     4.952, 4.346, 0.1857,   -0.1249, -0.009021,
3.114, 0.5236, 76.67, 7.62,   6.385, 4.355, 0.18,     -0.1255, -0.009045,
3.114, 0.5236, 76.67, 7.62,   7.502, 4.389, 0.1806,   -0.1253, -0.009028,
3.114, 0.5236, 76.67, 7.62,   8.514, 4.407, 0.1759,   -0.1258, -0.009054,
5.746, 0.4662, 79.24, 1.185,  7.993, 4.419, 0.1694,   -0.1267, -0.009094,
2.792, 0.6346,106.1,  0.2986, 2.331, 4.412, 0.1545,   -0.1289, -0.009202,
4.667, 0.5095,124.3,  2.102,  1.667, 4.419, 0.1448,   -0.1303, -0.009269,
2.44,  0.6346,105.,   0.83,   2.851, 4.436, 0.1443,   -0.1299, -0.009229,
1.491, 0.7118,120.6,  1.101,  1.877, 4.478, 0.1608,   -0.1262, -0.008983,
11.72, 0.3826,102.8,  9.231,  4.371, 4.489, 0.1517,   -0.1278, -0.009078,
7.126, 0.4804,119.3,  5.784,  2.454, 4.514, 0.1551,   -0.1268, -0.009005,
11.61, 0.3955,146.7,  7.031,  1.423, 4.533, 0.1568,   -0.1261, -0.008945,
10.99, 0.41,  163.9,  7.1,    1.052, 4.548, 0.1572,   -0.1256, -0.008901,
9.241, 0.4275,163.1,  7.954,  1.102, 4.553, 0.1544,   -0.1255, -0.008883,
9.276, 0.418, 157.1,  8.038,  1.29,  4.548, 0.1485,   -0.1259, -0.008889,
3.999, 0.6152, 97.6,  1.297,  5.792, 4.489, 0.1128,   -0.1309, -0.009107,
4.306, 0.5658, 97.99, 5.514,  5.754, 4.402, 0.06656,  -0.1375, -0.009421,
3.615, 0.6197, 86.26, 0.333,  8.689, 4.292, 0.01012,  -0.1459, -0.009835,
5.8,   0.49,  147.2,  6.903,  1.289, 4.187,-0.04539,  -0.1542, -0.01025,
5.6,   0.49,  130.,  10.,     2.844, 4.577, 0.13,     -0.1285, -0.009067,
3.55,  0.6068,124.7,  1.112,  3.119, 4.583, 0.1253,   -0.1291, -0.009084,
3.6,   0.62,  105.8,  0.1692, 6.026, 4.58,  0.1174,   -0.1301, -0.009129,
5.4,   0.53,  103.1,  3.931,  7.767, 4.581, 0.111,    -0.1309, -0.009161,
3.97,  0.6459,131.8,  0.2233, 2.723, 4.582, 0.1046,   -0.1317, -0.009193,
3.65,  0.64,  126.8,  0.6834, 3.411, 4.6,   0.1052,   -0.1315, -0.009178,
3.118, 0.6519,164.9,  1.208,  1.51,  4.614, 0.1043,   -0.1315, -0.009175,
2.031, 0.7181,153.1,  1.362,  1.958, 4.619, 0.09769,  -0.1325, -0.009231,
14.4,  0.3923,152.5,  8.354,  2.597, 4.671, 0.1136,   -0.1298, -0.009078,
10.99, 0.4599,138.4,  4.811,  3.726, 4.706, 0.1206,   -0.1287, -0.009009,
16.6,  0.3773,224.1,  6.28,   0.9121,4.732, 0.1244,   -0.128,  -0.008968,
10.54, 0.4533,159.3,  4.832,  2.529, 4.722, 0.1156,   -0.1292, -0.00903,
10.33, 0.4502,162.,   5.132,  2.444, 4.71,  0.106,    -0.1305, -0.0091, 
10.15, 0.4471,165.6,  5.378,  2.328, 4.698, 0.09647,  -0.1319, -0.009169,
9.976, 0.4439,168.,   5.721,  2.258, 4.681, 0.08536,  -0.1335, -0.009252,
9.804, 0.4408,176.2,  5.675,  1.997, 4.676, 0.07819,  -0.1345, -0.009302,
14.22, 0.363, 228.4,  7.024,  1.016, 4.663, 0.06867,  -0.1358, -0.009373,
9.952, 0.4318,233.5,  5.065,  0.9244,4.676, 0.06861,  -0.1357, -0.009363,
9.272, 0.4345,210.,   4.911,  1.258, 4.649, 0.05362,  -0.1379, -0.00948,
10.13, 0.4146,225.7,  5.525,  1.055, 4.634, 0.04335,  -0.1394, -0.009558,
8.949, 0.4304,213.3,  5.071,  1.221, 4.603, 0.02679,  -0.1418, -0.00969,
11.94, 0.3783,247.2,  6.655,  0.849, 4.584, 0.01494,  -0.1436, -0.009783,
8.472, 0.4405,195.5,  4.051,  1.604, 4.576, 0.007043, -0.1447, -0.009841,
8.301, 0.4399,203.7,  3.667,  1.459, 4.571, 0.0007046,-0.1456, -0.009886,
6.567, 0.4858,193.,   2.65,   1.66,  4.566,-0.005626, -0.1464, -0.00993,
5.951, 0.5016,196.1,  2.662,  1.589, 4.561,-0.01197,  -0.1473, -0.009975,
7.495, 0.4523,251.4,  3.433,  0.8619,4.572,-0.012,    -0.1472, -0.009965,
6.335, 0.4825,255.1,  2.834,  0.8228,4.569,-0.01755,  -0.148,  -0.01,   
4.314, 0.5558,214.8,  2.354,  1.263, 4.573,-0.01992,  -0.1482, -0.01001,
4.02,  0.5681,219.9,  2.402,  1.191, 4.57, -0.02547,  -0.149,  -0.01005,
3.836, 0.5765,210.2,  2.742,  1.305, 4.528,-0.04613,  -0.1521, -0.01022,
4.68,  0.5247,244.7,  2.749,  0.8962,4.494,-0.0637,   -0.1548, -0.01037,
3.223, 0.5883,232.7,  2.954,  1.05,  4.564,-0.027,    -0.1471, -0.009852,
2.892, 0.6204,208.6,  2.415,  1.416, 4.546,-0.04963,  -0.1523, -0.01022,
4.728, 0.5522,217.,   3.091,  1.386, 4.594,-0.03339,  -0.1496, -0.01006,
6.18,  0.52,  170.,   4.,     3.224, 4.608,-0.02886,  -0.1485, -0.00999,
9.,    0.47,  198.,   3.8,    2.032, 4.624,-0.02639,  -0.1481, -0.009971,
2.324, 0.6997,216.,   1.599,  1.399, 4.636,-0.02422,  -0.1477, -0.009939,
1.961, 0.7286,223.,   1.621,  1.296, 4.648,-0.02172,  -0.1471, -0.009903,
1.75,  0.7427,350.1,  0.9789, 0.5507,4.662,-0.1192,   -0.1752, -0.01196,
10.31, 0.4613,261.2,  4.738,  0.9899,4.69, -0.009867, -0.1449, -0.009771,
7.962, 0.519, 235.7,  4.347,  1.313, 4.715,-0.002113, -0.1435, -0.009689,
6.227, 0.5645,231.9,  3.961,  1.379, 4.729, 0.001392, -0.1428, -0.009644,
5.246, 0.5947,228.6,  4.027,  1.423, 4.729,-0.0005983,-0.143,  -0.009647,
5.408, 0.5811,235.7,  3.961,  1.358, 4.738, 0.001075, -0.1425, -0.009618,
5.218, 0.5828,245.,   3.838,  1.25,  4.751, 0.004244, -0.1419, -0.009576
  };

  if ( E < 1.0 ) {
     G4double Slow  = A[i][0]  ;
     G4double E1 = 1.0/1000.0 ;
     G4double Shigh = log( 1.0 + A[i][3]/E1 + A[i][4]*E ) * A[i][2]/E1 ;
     ionloss = sqrt(E) * Slow*Shigh / (Slow + Shigh) ; 

  } else if ( E < 10000.0 ) {
     G4double Slow  = A[i][0] * pow(E, A[i][1]) ;
     G4double E1 = E/1000.0 ;
     G4double Shigh = log( 1.0 + A[i][3]/E1 + A[i][4]*E ) * A[i][2]/E1 ;
     ionloss = Slow*Shigh / (Slow + Shigh) ; 

  } else {
     G4double le = log(1.0/E) ;
     ionloss = exp( A[i][5] + A[i][6]*le + A[i][7]*le*le + A[i][8]*le*le*le) ;
  }

  if ( ionloss <= 0.) ionloss = 0. ;

  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetStoppingPower1977n(G4double Z1, G4double Z2, 
                                          G4double M1, G4double M2, G4double E)
{
  G4double ionloss ;

  // The fit of nuclear stopping from: 
  // J.F.Ziegler, Helium Stopping Powers and
  // Ranges in All Elemental Matter, Vol.4, Pergamon Press, 1977

  G4double rm = (M1 + M2) * sqrt( pow(Z1, 2.0/3.0) + pow(Z2, 2.0/3.0) ) ;

  G4double er = 32.53 * M2 * E / ( Z1 * Z2 * rm ) ;  // reduced energy

  if ( er < 0.01 ) {
     ionloss = sqrt(er) * 1.593 ; 

  } else if ( E < 10.0 ) {
     ionloss = 1.7 * sqrt(er) * log(er + exp(1.0)) / 
               (1.0 + 6.8 * er + 3.4 * pow(er, 1.5)) ; 

  } else {
     ionloss = log(0.47 * er) * 0.5 / er  ;
  }

  ionloss *= 8.462 * Z1 * Z2 * M1 / rm ; // Return to [ev/(10^15 atoms/cm^2]
  if ( ionloss <= 0.) ionloss = 0. ;

  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::PrintInfoDefinition()
{
  G4String comments = "  Knock-on electron cross sections . ";
           comments += "\n         Good description above the mean excitation energy.\n";
           comments += "         delta ray energy sampled from  differential Xsection.";

  G4cout << endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << G4BestUnit(LowestKineticEnergy,
                                                  "Energy")
         << " to " << G4BestUnit(HighestKineticEnergy,"Energy")
         << " in " << TotBin << " bins."
         << "\n        Low energy approximation implemented in the function " << DEDXtable
         << "\n  19.07 from " << G4BestUnit(ZieglerLowEnergy,"Energy")
         << " to " << G4BestUnit(ZieglerHighEnergy,"Energy") << " \n" << endl ;
}






