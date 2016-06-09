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
// $Id: G4ElasticHadrNucleusHE.hh,v 1.25 2006/06/29 20:09:01 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// G4ElasticHadrNucleusHe.hh

//  The generator of high energy hadron-nucleus elastic scattering
//  The hadron kinetic energy T > 1 GeV
//  N.  Starkov 2003.
//
//  Modifications:
//  19.05.04 Variant for G4 6.1: The 'ApplyYourself' was changed
//  14.11.05 The HE elastic scattering on proton is added (N.Starkov)
//  30.05.06 Version without use of elastic data (N.Starkov)
//

#ifndef G4ElasticHadrNucleusHE_h
#define G4ElasticHadrNucleusHE_h 1

#include <vector>

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4Ions.hh"
#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "Randomize.hh"
#include "G4Nucleus.hh"
#include "G4IonTable.hh"

#include "G4DiffElasticHadrNucleus.hh"
#include "G4IntegrHadrNucleus.hh"
#include "G4HadronicInteraction.hh"

#define   ONQ2     150      //  The number of steps on Q2
#define   ONE      5        //  The number of steps on E
#define   AreaNumb 6        //  The number of order steps on E
#define   ONQ2XE   ONQ2*ONE //  The dimension of a distr. func. array
#define   MaxN     10       //  The atomic number where the calculation
                            //  on the formula is changed on the integral
                            //  one

class ElasticData
{
public:

  G4String  hadrName;
  G4int     nuclAtomicNumber;
  G4double  TableE[ONE*AreaNumb];
  G4double  TableQ2[ONQ2];
  G4double  TableFQ2[ONQ2];
  G4double  TableCrossSec[ONQ2XE*AreaNumb];
  G4double  CurrentT;
  G4int     CurrentN, Nstep;
  G4double  RandMax, theMaxQ2;
  G4double  Weight, dQ2;
  G4double  theR1, R2, Pnucl, Aeff;
  NucleusParameters NucPar;     
    
  ElasticData() {;}

  ElasticData(G4String HadrName, G4int AtomWeight);

  virtual ~ElasticData(){;}

  ElasticData (const ElasticData &t)
     {
       G4int k;

       hadrName         = t.hadrName;
       nuclAtomicNumber = t.nuclAtomicNumber;

       for(k = 0; k<ONE*AreaNumb; k++) TableE[k] = t.TableE[k];

       for(k = 0; k<ONQ2; k++) TableQ2[k] = t.TableQ2[k];

       for(k = 0; k< ONQ2XE*AreaNumb; k++)
                TableCrossSec[k] = t.TableCrossSec[k];

       theR1    = t.theR1;
       R2    = t.R2;
       Aeff  = t.Aeff;
       Pnucl = t.Pnucl;
       theMaxQ2 = t.theMaxQ2;
       dQ2   = t.dQ2;
       NucPar = t.NucPar;

     }

  ElasticData & operator=(const ElasticData &t)
     {
       G4int k;
       if(this!=&t)
	 {
          hadrName         = t.hadrName;
          nuclAtomicNumber = t.nuclAtomicNumber;
          for(k = 0; k<ONE*AreaNumb; k++)
	              TableE[k] = t.TableE[k];

          for(k = 0; k<ONQ2; k++) TableQ2[k] = t.TableQ2[k];

          for(k = 0; k< ONQ2XE*AreaNumb; k++)
                      TableCrossSec[k] = t.TableCrossSec[k];

          theR1 = t.theR1;
          R2    = t.R2;
          Aeff  = t.Aeff;
          Pnucl = t.Pnucl;
          theMaxQ2 = t.theMaxQ2;
          NucPar = t.NucPar;
          dQ2    = t.dQ2;

	 }
       return *this;
     }

  void Clean()
     {
       G4int k;

       hadrName         = "Nothing";
       nuclAtomicNumber = 0;
       for(k = 0; k<ONE*AreaNumb; k++) TableE[k] = 0;

       for(k = 0; k<ONQ2; k++) TableQ2[k] = 0;

       for(k = 0; k< ONQ2XE*AreaNumb; k++) TableCrossSec[k] = 0;

          theR1    = 0;
          R2    = 0;
          Aeff  = 0;
          Pnucl = 0;
          theMaxQ2 = 0;
     }

  G4double GetQ2limit(G4double aR1);
  G4int    GetNumberE(G4double E);
};

//  ############################################################
class G4ElasticHadrNucleusHE : public G4DiffElasticHadrNucleus,
			       public G4HadronicInteraction
{
public:
  G4ElasticHadrNucleusHE(const G4ParticleDefinition * aHadron,
                                      G4Nucleus            * aNucleus);

  G4ElasticHadrNucleusHE();

  virtual ~G4ElasticHadrNucleusHE() {;}

  G4HadFinalState * ApplyYourself( const G4HadProjectile  &aTrack,
				   G4Nucleus        &aNucleus);

  G4double RandomElastic0();

  G4double RandomElastic1( const G4DynamicParticle *   aHadron,
			   const ElasticData       *   aData);

  G4double SampleT(const G4ParticleDefinition* p,
		   G4double pTotLabMomentum, G4int Z, G4int N);
  G4double SampleT1(const G4ParticleDefinition* p,
		    G4double pTotLabMomentum, G4int Z, G4int N);

  G4bool GetHadronNucleusData(G4DynamicParticle * aParticle,
			      G4Nucleus         * aNucleus,
			      ElasticData & ElD );
private:
public:
  G4int   ReadOfData(G4ParticleDefinition * aParticle,
		     G4Nucleus            * aNucleus);

  G4double GetQ2limit(G4double aR1);

  void  CreationArray(const G4DynamicParticle * aHadron,
		      G4Nucleus         * aNucleus);

  void  ArrayForHeavy(const G4DynamicParticle * aHadron,
		      G4Nucleus         * aNucleus);

  void  ArrayForLight(const G4DynamicParticle * aHadron,
                                   G4Nucleus         * aNucleus);

  G4double InterPol(G4double X1, G4double X2, G4double X3,
		    G4double Y1, G4double Y2, G4double Y3, 
		    G4double X);
  G4double HadronNucleusQ2(G4DynamicParticle * aHadron,
			   G4Nucleus         & aNucleus);
  G4double HadronNucleusQ2_2(G4DynamicParticle * aHadron,
			     G4int A, G4int kk, 
			     ElasticData * pElD);
  G4double GetLightFq2(G4int N, G4double Q, G4int Step,
		       NucleusParameters * NP);

  G4int    GetBinom(G4int m, G4int n);
//  ======================================================
  
  G4DynamicParticle  aHad;

  G4float  GetFt(G4double T);
  G4float  GetDistrFun(G4double Q2);
  G4double GetQ2(G4double Ran);
  G4double GetQ2_2(G4int  N, G4double * Q, 
		   G4double * F, G4double R);
  G4double HadronProtonQ2(G4DynamicParticle * aHadron);
                                                                          
  G4double  Weight;
  G4int     HadrCode;
  G4String  HadronName;
  //  G4double  RR1, R2, Pnucl, Aeff;

  std::vector<ElasticData> SetOfElasticData;
  ElasticData       ElD;
  NucleusParameters NucPar;

  G4IonTable                * MyIonTable;
  G4DiffElasticHadrNucleus    aDiffElHadNcls;

//     G4HadFinalState  FinState; 

  G4int     Nstep,         // The number of steps on Q2
            iKindWork,     // 
            iContr,        //
            iPoE;          // The number of steps on E
  G4int     iTypeWork, CurrentN;
  G4double  aNucleon;
//              ,* pTableCrSec,   //  The array of distr. func.
//                              //  at all energies
//             * pTableE;       //  The array of E values
  G4double  iQ2[ONQ2],     //  The array of Q2 values
            pTableCrSec[ONQ2XE*AreaNumb],
            pTableE[ONE*AreaNumb],
            iIntgr[ONQ2],  //  The array of distr. func.   
                                  //    at one energy
            Factorials1[250]; // The array for factorials
  G4double  dEbeg1, dEend1, dQ2, maxQ2, RandMax;
};     //   The end of the class description
//  ######################################################
#endif

