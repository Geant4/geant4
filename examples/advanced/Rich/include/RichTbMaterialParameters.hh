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
// Rich advanced example for Geant4
// RichTbMaterialParameters.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbMaterialParameters_h
#define RichTbMaterialParameters_h 1

#include "globals.hh"
#include "RichTbGeometryParameters.hh"
#include "RichTbAnalysisManager.hh"
#include "RichTbRunConfig.hh"
#include "AerogelTypeSpec.hh"

extern void InitializeRichTbMaterial();
extern void HistoRichTbMaterialProperties(RichTbRunConfig* RConfig);


extern std::vector<G4double> InitializeHpdQE(G4int);
extern std::vector<G4double> InitializeHpdWaveL(G4int);
extern std::vector<G4double> InitN2RefIndex(G4double, G4double);
extern std::vector<G4double> InitN2RefPhotW();
extern std::vector<G4double> InitAgelPhotW();
extern std::vector<G4double> InitializePhotonMomentumVector();
extern std::vector<G4int> getDeadPixelList(G4int ihpdNum , G4int IsectNum);
extern std::vector<G4double>GetAerogelRScatLength(AerogelType);
extern G4double GetCurrentBulkTrans(G4double currentMatRefIndex,
				       G4double currentNeighbourRefIndex, 
				       G4double MaxTotMeasuredTransmission);

// extern void GetGlassD263FilterTrans();
//
// the following value should be calculated in terms of the
// fundamental constants in the future. It is the conversion
// factor between the wavelength of a photon in nanometers
// and its energy in eV.
static const G4double PhotMomWaveConv=1243.125;

// Limits of Photon Energy  and number of bins for the
// Photon energy range.
//   static const G4double PhotonMinEnergy=1.5*eV;
static const G4double PhotonMinEnergy=1.3*eV;
static const G4double PhotonMaxEnergy=7.3*eV;
   static const G4int NumPhotWaveLengthBins = 1000;
   static const G4int NumPhotonRichMirrorReflWaveLengthBins=63;
   static const G4int NumAerogelRefIndexPhotonEnergyBins=37;
//   static const G4int NumFilterGlassD263WaveLengthBins=10;
// Defintion of STP pressure and temp

//static const G4double Pressure_STP=1.013*bar;
//static const G4double Temperature_STP=273.*kelvin
static const G4double  GasPressure_STP=STP_Pressure;
static const G4double GasTemperature_STP=STP_Temperature;
// Ref Index of nitrogen using sellmeir parametrization
static const G4double SellN2E1=13.414;
static const G4double SellN2E2=23.215;
static const G4double SellN2F1=921.28;
static const G4double SellN2F2=3569.60;
static const G4double GasMolWeightN2=28.02;      //unit is grams
static const G4double GasRhoN2atSTP=0.00125053;  //unit is gramPercm3
//
//
//Mirror reflectivity
// In the following, the bins at 100 nm, and 1000nm are
// defined just for convenience of interpolation. 
// They are not measured points.
static const G4double PhotonWavelengthRefl[]=
{100.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0,
 290.0, 300.0, 310.0, 320.0, 330.0, 340.0, 350.0, 360.0, 370.0, 380.0,
 390.0, 400.0, 410.0, 420.0, 430.0, 440.0, 450.0, 460.0, 470.0, 480.0,
 490.0, 500.0, 510.0, 520.0, 530.0, 540.0, 550.0, 560.0, 570.0, 580.0,
 590.0, 600.0, 610.0, 620.0, 630.0, 640.0, 650.0, 660.0, 670.0, 680.0,
 690.0, 700.0, 710.0, 720.0, 730.0, 740.0, 750.0, 760.0, 770.0, 780.0,
 790.0, 800.0, 1000.0 };

static const G4double RichTbMirrorReflectivity[]=
 {0.0, 0.9106, 0.9232, 0.9285, 0.9314, 0.9323, 0.9312, 0.9287, 0.9264,
 0.9234, 0.9195, 0.9156, 0.9109, 0.9066, 0.9022, 0.8981, 0.8925, 0.8883,
 0.8836, 0.8796, 0.8756, 0.8727, 0.8697, 0.8672, 0.8653, 0.8636, 0.8624,
 0.8612, 0.8608, 0.8601, 0.8601, 0.8601, 0.8600, 0.8603, 0.8603, 0.8604,
 0.8605, 0.8608, 0.8609, 0.8608, 0.8608, 0.8606, 0.8604, 0.8600, 0.8598,
 0.8591, 0.8581, 0.8573, 0.8563, 0.8549, 0.8535, 0.8517, 0.8497, 0.8475,
 0.8447, 0.8417, 0.8382, 0.8388, 0.8296, 0.8258, 0.8204, 0.8172, 0.8172 };

static const G4double RichTbMirrorEfficiency[]=
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
   0.0, 0.0, 0.0 };

// Transmission in quartz

static const  G4int NumPhotonRichTbGasQuartzWSurfaceWaveLengthBins=10;
static const  G4double RichGasQuartzWSurfacePhotMom[]=
   {1.0*eV,2.0*eV, 3.0*eV,4.0*eV,5.0*eV,6.0*eV,7.0*eV,8.0*eV,
    9.0*eV,10.0*eV};
static const  G4double RichTbGasQuartzWSurfaceReflectivity[]=
    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
static const  G4double RichTbGasQuartzWSurfaceEfficiency[]=
    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

// Transmission in the filters
static const  G4int NumPhotonRichTbFilterSurfaceWaveLengthBins=10;
static const  G4double RichTbFilterSurfacePhotMom[]=
   {1.0*eV,2.0*eV, 3.0*eV,4.0*eV,5.0*eV,6.0*eV,7.0*eV,8.0*eV,
    9.0*eV,10.0*eV};
static const  G4double RichTbFilterSurfaceReflectivity[]=
    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
static const  G4double RichTbFilterSurfaceEfficiency[]=
    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
// In the G4Example only 1 type of filter is used.
// In the LHCb implementation several types of filters
// were used. The program was originally setup for upto
// 6 types of filters.
// Now for the Glass D263 filter
static const G4int NumPhotBinGlassD263Trans=601; 
// static const G4double RefIndexGlassD263=1.50;
// artifically low value for the ref index of filter used for the
// G4example. If one uses the original value of 1.5, the cherenkov
// photons created in the filter by the 9 GeV/c pions will have cherenkov
// angle larger than the critical angle and hence would cause
// an (almost) infinite loop of total internal reflections. These 
// photons eventually die from photon absorption, 
// but they use a lot of cpu time.
// In the lhcb implementation the G4OpboundaryProcess is modified to
// to take care of this. But for the G4example this complication is
// avoided by using an artificially low value for the ref idex of the filter. 
// SE.  15-11-2002.
static const G4double RefIndexGlassD263=1.30;
// for initlializing the arrays for the filter tramsmission
// the following value is used. The arrays are later resized 
// to the appropiate value.
static const G4int NumFilterTransBins= NumPhotBinGlassD263Trans;
//
//Now for the different type of Aerogel Tiles
//In the G4example the same Aerogel type is simply repeated 5 times.
//In the LHCb implementation 5 different types are used.
static const G4double AerogelTypeAClarity
 =0.00719*micrometer*micrometer*micrometer*micrometer/cm;
static const G4double AerogelTypeATotTrans=0.9368;
static const G4double AerogelTypeANominalRefIndex=1.03066;

static const G4double AerogelTypeBClarity= AerogelTypeAClarity;
static const G4double AerogelTypeBTotTrans=AerogelTypeATotTrans;
static const G4double AerogelTypeBNominalRefIndex=AerogelTypeANominalRefIndex;
//
static const G4double AerogelTypeCClarity=AerogelTypeAClarity;
static const G4double AerogelTypeCTotTrans=AerogelTypeATotTrans;
static const G4double AerogelTypeCNominalRefIndex=AerogelTypeANominalRefIndex;
static const G4double AerogelTypeDClarity=AerogelTypeAClarity;
static const G4double AerogelTypeDTotTrans=AerogelTypeATotTrans;
static const G4double AerogelTypeDNominalRefIndex=AerogelTypeANominalRefIndex;
static const G4double AerogelTypeEClarity=AerogelTypeAClarity;
static const G4double AerogelTypeETotTrans=AerogelTypeATotTrans;
static const G4double AerogelTypeENominalRefIndex=AerogelTypeANominalRefIndex;
static const G4double AerogelReferenceRefIndWavelength[]=
{400.0*nanometer,400.0*nanometer,400.0*nanometer,400.0*nanometer,
 400.0*nanometer};
static const G4double StdAerogelNominalRefIndex=1.034;
// for test try tdr aerogel
//static const G4double StdAerogelNominalRefIndex=1.03123653;
// static const char* StdAerogelRefIndFile =
// "/afs/cern.ch/user/s/seaso/mycmt/RichTb/v5/inputData/aerogelRefIndex.txt";
// The following is the nominal ref index at STP for nitrogen.
// the aerogel is kept inside the nitrogen;
static const G4double NitrogenNominalRefIndex=1.000298;
//
// Quantum efficiency of the photocathodes
// the following line already in geometryparameters.hh
//static const G4int NumberOfHpds=4;

static const G4int NumHpdTot= NumberOfHpds;
static const G4int NumQEbins=41;
// for now all HPDs have the same wavelength bins.
// for the G4Example the QE is multiplied by 1.08 to account for
// the fresnel losses at the HPD Input Quartz window.
// In the G4 example, the correction to
// account for the QE reduction in the periphery of the HPD
// is not done, for simplicity. 
// In the LHCb implementation this taken care of more accurately.
static const G4double HpdQEReductionFactor=1.08;
static const G4double HpdQEWaveL[]=
{200.0,210.0,220.0,230.0,240.0,250.0,260.0,270.0,280.0,290.0,300.0,
 310.0,320.0,330.0,340.0,350.0,360.0,370.0,380.0,390.0,400.0,410.0,
 420.0,430.0,440.0,450.0,460.0,470.0,480.0,490.0,500.0,510.0,520.0,
 530.0,540.0,550.0,560.0,570.0,580.0,590.0,600.0};
static const G4double Hpd0QEPerCent[]=
{0.0,0.0,2.4264,9.0559,14.5723,18.9039,24.2476,27.9504,29.0910,
 31.5069,28.2676,26.8652,26.0060,25.2934,24.5034,23.4996,22.6084,
 21.0618,19.6225,18.2942,16.4564,14.3433,12.5074,10.6036,8.7953,
 7.2646,6.1266,5.0693,3.9949,3.0248,2.2300,1.5261,0.8610,0.4803,
 0.2870,0.1320,0.0500,0.0135,0.0,0.0,0.0};
static const G4double Hpd1QEPerCent[]= 
{0.0,1.0810,3.1196,6.4212,10.4983,14.3031,17.9612,20.6541,22.1570,
 25.5633,25.3696,25.2293,25.1988,24.9773,24.8055,24.4088,24.1165,
 23.3368,22.4608,21.8726,20.4717,18.9725,17.5489,15.8493,14.1541,
 12.5689,11.1047,9.7609,8.4922,7.3264,6.2243,5.0602,3.5248,2.3012,
 1.5403,1.0813,0.6823,0.4553,0.2524,0.1227,0.0583};
static const G4double Hpd2QEPerCent[]=
{0.0,0.0,1.3648,4.1624,6.8814,9.5704,12.5728,15.4286,16.7373,19.9520,
 20.0918,20.3359,20.6562,20.7517,20.6236,20.1992,19.8387,18.6558,17.5653,
 16.6784,15.2398,13.6456,12.1764,10.5722,9.0127,7.6529,6.4820,5.4517,
 4.4408,3.5450,2.7400,1.9825,1.3073,0.7549,0.4712,0.3031,0.1675,
 0.0803,0.0391,0.0169,0.0094};
static const G4double Hpd3QEPerCent[]=
{0.0,11.3501,12.1319,13.2398,16.2100,18.7656,21.7482,24.3610,25.0792,
 28.6009,28.2349,28.4114,28.5484,28.6614,28.8812,28.7106,28.1548,28.0591,
 25.5305,21.8913,19.6042,18.3709,17.3953,17.1329,16.2153,14.8937,14.6245,
 13.0102,13.8088,13.2000,12.0881,10.7617,8.1836,5.9729,4.6335,3.6822,
 2.9047,2.2591,1.7516,1.2502,0.8219};
// only the linear term is used for now
// static const G4double HpdDemagConst[]={ }
static const G4double HpdDemagLinearTerm[]=
{2.46433268,2.28343310,2.29306727,2.37845184};
// In the G4example the quadratic term is neglected and hence set to zero.
static const G4double HpdDemagQuadraticTerm[]=
{0.0,0.0,0.0,0.0};
static const G4double HpdDemagErrorPercent=0.0;
// in the G4example  a PSF factor is used, just an example.
static const G4double PadHpdPSFsigma=100.0*micrometer;
//quartz Transmission for 10 mm thickness.
static const G4double QuTransDataThickness =10.0*mm;
static const G4int NumPhotbinQuartzTrans=6;
static const G4double QuartzTransWL[]=
{180.0 ,185.0, 190.0, 200.0, 220.0,1000.0};
static const G4double QuartzTransmis[]=
{0.80, 0.90, 0.95 ,0.98 , 1.0,1.0};
static const G4double PhCathodeNominalTransmission=0.52;

//The dead Pixel List is set to zero for the G4Example.
//Hence all pixels are set to be active.
// The following is just a maximum number of deadpixels for array sizes.
static const G4int MaxNumDeadPixelPerHpdSect=50;

// Now for the back scattering
static const     G4double backscaprob=0.18;
// Now for the approximate and adhoc way of evaluating the effect of
// backscattering.
static const G4double NsigmaInPedCut=4.0;
static const G4double SignalToNoiseInData=10.0;

#endif 











