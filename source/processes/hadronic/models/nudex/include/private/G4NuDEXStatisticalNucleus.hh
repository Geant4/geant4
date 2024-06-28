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
// -------------------------------------------------------------------
//
//      Author:        E.Mendoza
// 
//      Creation date: May 2024
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//  NuDEX code (https://doi.org/10.1016/j.nima.2022.167894)
// 


#ifndef NUDEXSTATISTICALNUCLEUS_HH
#define NUDEXSTATISTICALNUCLEUS_HH 1

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "G4NuDEXRandom.hh"


class G4NuDEXLevelDensity;
class G4NuDEXInternalConversion;
class G4NuDEXPSF;


//This define remains:
//#define GENERATEEXPLICITLYALLLEVELSCHEME 1

//Class to obtain the level density for each excitation energy, spin, and parity
//All energies in MeV, all times in s
//Some of the class methods could be functions out of the class

struct Level{
  G4double Energy;
  G4int spinx2;
  G4bool parity; //true/false --> positive,negative
  unsigned int seed;
  G4int KnownLevelID;
  G4int NLevels;
  G4double Width;
};



//multipolarity of a transition is ...,-2,-1,0,1,2,... --> ...,M2,M1,Unk,E1,E2,...

struct KnownLevel{
  G4int id;
  G4double Energy;
  G4int spinx2;
  G4bool parity; //true/false --> positive,negative
  G4double T12; //half life - seconds
  G4int Ndecays;
  G4double* decayFraction;
  std::string* decayMode;
  G4int NGammas;
  G4int *FinalLevelID,*multipolarity;
  G4double *Eg,*cumulPtot,*Pg,*Pe,*Icc;
};



G4int ComparisonLevels(const void* va, const void* vb);
void CopyLevel(Level* a,Level* b);
void CopyLevel(KnownLevel* a,Level* b);


class G4NuDEXStatisticalNucleus{

public:
  G4NuDEXStatisticalNucleus(G4int Z,G4int A);
  ~G4NuDEXStatisticalNucleus();

public:
  //Initialize everything. All the required files should be in dirname.
  //some of the data could also be in inputfname
  G4int Init(const char* dirname,const char* inputfname=0);

  //If InitialLevel==-1 then we start from the thermal capture level
  //If ExcitationEnergy>0 then is the excitation energy of the nucleus
  //If ExcitationEnergy<0 then is a capture reaction of a neutron with energy -ExcitationEnergy (MeV)
  G4int GenerateCascade(G4int InitialLevel,G4double ExcitationEnergy,std::vector<char>& pType,std::vector<double>& pEnergy,std::vector<double>& pTime);

  G4int GetClosestLevel(G4double Energy,G4int spinx2,G4bool parity); //if spinx2<0, then retrieves the closest level of any spin and parity
  G4double GetLevelEnergy(G4int i_level);
  void GetSnAndI0(G4double &sn,G4double &i0){sn=Sn; i0=I0;}
  Level* GetLevel(G4int i_level);
  void ChangeLevelSpinParityAndBR(G4int i_level,G4int newspinx2,G4bool newParity,G4int nlevels,G4double width,unsigned int seed=0); //if nlevels or width are negative they don't change. If seed (to generate the BR) is 0 it does not change.
  void ChangeThermalCaptureLevelBR(G4double LevelEnergy,G4double absoluteIntensity);

  void SetSomeInitalParameters(G4int LDtype=-1,G4int PSFFlag=-1,G4double MaxSpin=-1,G4int minlevelsperband=-1,G4double BandWidth_MeV=0,G4double maxExcEnergy=0,G4int BrOption=-1,G4int sampleGammaWidths=-1,unsigned int aseed1=0,unsigned int aseed2=0,unsigned int aseed3=0);
  void SetInitialParameters02(G4int knownLevelsFlag=-1,G4int electronConversionFlag=-1,G4double primGamNormFactor=-1,G4double primGamEcut=-1,G4double ecrit=-1);
  void SetBandWidth(G4double bandWidth){ if(bandWidth==0){bandWidth=-1;} BandWidth=bandWidth;} //So it is not re-written with the lib-params.
  void SetBrOption(G4int BrOption){BROpt=BrOption;}
  void SetRandom1Seed(unsigned int seed){theRandom1->SetSeed(seed); Rand1seedProvided=true;}
  void SetRandom2Seed(unsigned int seed){theRandom2->SetSeed(seed); Rand2seedProvided=true;}
  void SetRandom3Seed(unsigned int seed){theRandom3->SetSeed(seed); Rand3seedProvided=true;}
  
  G4NuDEXRandom* GetRandom3(){return theRandom3;}
  G4bool HasBeenInitialized(){return hasBeenInitialized;}


  //-------------------------------------------------------
  //Print:
  void PrintAll(std::ostream &out);
  
  void PrintParameters(std::ostream &out);
  void PrintKnownLevels(std::ostream &out);
  void PrintLevelDensity(std::ostream &out);
  void PrintLevelScheme(std::ostream &out);
  void PrintThermalPrimaryTransitions(std::ostream &out);
  void PrintPSF(std::ostream &out);
  void PrintICC(std::ostream &out);
  void PrintTotalCumulBR(G4int i_level,std::ostream &out);
  void PrintBR(G4int i_level,G4double MaxExcEneToPrint_MeV,std::ostream &out);
  void PrintInput01(std::ostream &out);
  //----------------
  void PrintKnownLevelsInDEGENformat(std::ostream &out);
  void PrintLevelSchemeInDEGENformat(const char* fname,G4int MaxLevelID=-1);
  //-------------------------------------------------------

  
private:
  //-------------------------------------------------------
  //Used by Init():
  //Read different data from files (do it in this order). If returnval<0 --> error reading file or nucleus not present in the file:
  G4int ReadSpecialInputFile(const char* fname);
  G4int ReadGeneralStatNuclParameters(const char* fname);
  G4double ReadEcrit(const char* fname);
  G4double ReadKnownLevels(const char* fname);
  void CreateLevelScheme();
  G4int InsertHighEnergyKnownLevels();
  void ComputeKnownLevelsMissingBR();
  void MakeSomeParameterChecks01();
  //-------------------------------------------------------
  G4double TakeTargetNucleiI0(const char* fname,G4int& check);
  void CreateThermalCaptureLevel(unsigned int seed=0); //If seed (to generate the BR) is 0 it does not change.
  void GenerateThermalCaptureLevelBR(const char* dirname);
  //-------------------------------------------------------

  //-------------------------------------------------------
  //cascade generation:
  G4double ComputeDecayIntensities(G4int i_level,G4double* cumulativeBR=0,G4double randnumber=-1,G4double TotGR=-1,G4bool AllowE1=false);
  G4int SampleFinalLevel(G4int i_level,G4int& multipolarity,G4double &icc_fac,G4int nTransition);
  G4int GetMultipolarity(Level* theInitialLevel,Level* theFinalLevel);
  //-------------------------------------------------------


private:
  //-------------------------------------------------------
  //Used to create the unknown Levels:
  G4int GenerateLevelsInBigRange(G4double Emin,G4double Emax,G4int spinx2,G4bool parity,Level* someLevels,G4int MaxNLevelsToFill); //salen sin ordenar
  G4int GenerateLevelsInSmallRange(G4double Emin,G4double Emax,G4int spinx2,G4bool parity,Level* someLevels,G4int MaxNLevelsToFill); //salen sin ordenar
  G4int GenerateWignerLevels(G4double Emin,G4double Emax,G4int spinx2,G4bool parity,Level* someLevels,G4int MaxNLevelsToFill); //salen ordenados
  G4int GenerateBandLevels(G4int bandmin,G4int bandmax,G4int spinx2,G4bool parity,Level* someLevels,G4int MaxNLevelsToFill);
  G4int GenerateAllUnknownLevels(Level* someLevels,G4int MaxNLevelsToFill); //salen ordenados
  G4int CreateBandsFromLevels(G4int thisNLevels,Level* someLevels,G4int spinx2,G4bool parity); 
  G4int EstimateNumberOfLevelsToFill(); //to estimate the length of "theLevels" vector
  //-------------------------------------------------------


private:

  //General info:
  G4int A_Int,Z_Int;
  G4double Sn,D0,I0; //I0 es el del nucleo A-1 (el que captura)
  G4bool hasBeenInitialized;
  std::string theLibDir;

  G4NuDEXRandom* theRandom1;  //To generate the unknown level scheme
  G4NuDEXRandom* theRandom2;  //To calculate the Gamma-rho values (i.e. to generate the branching ratios)
  G4NuDEXRandom* theRandom3;  //To generate the cascades
  unsigned int seed1,seed2,seed3;
  G4bool Rand1seedProvided,Rand2seedProvided,Rand3seedProvided;

  //--------------------------------------------------------------------------
  //Parameters which will define how the level scheme will be created:
  G4double Ecrit; //Energy between the known and unknown levels
  G4double MaxExcEnergy,BandWidth;
  G4int maxspinx2,NBands,MinLevelsPerBand; //maximum spin (x2) to consider, number of bands used to "rebin" the stat. part
  G4int LevelDensityType; //if negative or cero, use the default one.
  G4int PSFflag; // use IAEA PSF-data (PSFflag==0), use RIPL-3 data (PSFflag==1)
  G4double E_unk_min,E_unk_max; //min and max energy where the statistical part will be generated
  G4double Emin_bands,Emax_bands; //limites de energia para calcular las bandas de niveles
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //Level scheme:
  Level* theLevels; //known+unknown levels
  KnownLevel* theKnownLevels; // known levels
  G4int NKnownLevels,NUnknownLevels,NLevels,KnownLevelsVectorSize;
  Level theThermalCaptureLevel;
  G4int NLevelsBelowThermalCaptureLevel; //excluding the last one
  G4int KnownLevelsFlag;
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //Branching ratios:
  G4int BROpt,SampleGammaWidths;
  G4double* TotalGammaRho;
  G4double* theThermalCaptureLevelCumulBR;
  G4double** TotalCumulBR; //all BR
  G4double PrimaryGammasIntensityNormFactor;
  G4double PrimaryGammasEcut; //This variable can be used to avoid generating transitions close to the "Primary Gammas" region
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //LD,ICC, PSF:
  G4int ElectronConversionFlag;
  G4NuDEXLevelDensity* theLD;
  G4NuDEXInternalConversion* theICC;
  G4NuDEXPSF* thePSF;
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //for internal use, when generating the cascades:
  G4int theSampledLevel,theSampledMultipolarity;
  //--------------------------------------------------------------------------
};

//***************************************************************************************************************
//***************************************************************************************************************




#endif




