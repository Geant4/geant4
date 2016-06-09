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
// $Id: G4ElasticHadrNucleusHE.hh,v 1.15 2005/06/10 13:28:07 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// G4ElasticHadrNucleusHe.hh

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
     G4double  TableCrossSec[ONQ2XE*AreaNumb];
    
     ElasticData() {;}
    ~ElasticData(){;}

     ElasticData (const ElasticData &t)
     {
       G4int k;

          hadrName         = t.hadrName;
          nuclAtomicNumber = t.nuclAtomicNumber;
          for(k = 0; k<ONE*AreaNumb; k++)
	       TableE[k] = t.TableE[k];

          for(k = 0; k<ONQ2; k++)
	       TableQ2[k] = t.TableQ2[k];

          for(k = 0; k< ONQ2XE*AreaNumb; k++)
	    TableCrossSec[k] = t.TableCrossSec[k];
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

          for(k = 0; k<ONQ2; k++)
	       TableQ2[k] = t.TableQ2[k];

          for(k = 0; k< ONQ2XE*AreaNumb; k++)
	    TableCrossSec[k] = t.TableCrossSec[k];
         }
     return *this;
     }
   };

   class G4ElasticHadrNucleusHE : public G4DiffElasticHadrNucleus,
                                  public G4HadronicInteraction
   {
 public:
         G4ElasticHadrNucleusHE(const G4ParticleDefinition * aHadron,
                                      G4Nucleus            * aNucleus);

         G4ElasticHadrNucleusHE();

        ~G4ElasticHadrNucleusHE() {;}

         G4HadFinalState * ApplyYourself( const G4HadProjectile  &aTrack,
                                                G4Nucleus        &aNucleus);

         G4double RandomElastic0();

         G4double RandomElastic1( const G4DynamicParticle *   aHadron,
                                  const ElasticData       *   aData);

 private:
         G4int   ReadOfData(G4ParticleDefinition * aParticle,
                            G4Nucleus            * aNucleus);

         G4String GetHadronName(const G4DynamicParticle * aHadron);

         G4double GetQ2limit(G4double R1);

         void  CreationArray(const G4DynamicParticle * aHadron,
                                   G4Nucleus         * aNucleus);

         void  ArrayForHeavy(const G4DynamicParticle * aHadron,
                                   G4Nucleus         * aNucleus);

         void  ArrayForLight(const G4DynamicParticle * aHadron,
                                   G4Nucleus         * aNucleus);

         G4double InterPol(G4double X1, G4double X2, G4double X3,
                           G4double Y1, G4double Y2, G4double Y3, 
                           G4double X);
//    +++++++++++++++++++++++++++++++++++++++++++++++++++++
   G4double Factorial1(G4int N)
     {
        G4double  Res;
                  Res = 1;
              if(N == 0) return Res;

          if(N < 100) for(G4int M = 1; M<=N; M++)  Res = Res*M;         

           else  Res = 2.50662827*std::exp(static_cast<double>(-N-1))*
                         std::pow(static_cast<double>(N+1),N+0.5)*
                         (1+1/12/(N+1)+1/288/(N+1)/(N+1)-
                         139/51840/(N+1)/(N+1)/(N+1)-
                         571/2488320/(N+1)/(N+1)/(N+1)/(N+1));
              return Res;
     }
//     ++++++++++++++++++++++++++++++++++++++++++++++++++
         std::vector<ElasticData> SetOfElasticData;

         G4IonTable                * MyIonTable;
         G4DiffElasticHadrNucleus    aDiffElHadNcls;
//         G4HadFinalState  FinState; 

         G4int     Nstep,         // The number of steps on Q2
                   iKindWork,     // 
                   iContr,        //
                   iPoE;          // The number of steps on E
         G4int     iTypeWork;
         G4double  aNucleon;  //, 
//                 * pTableCrSec,   //  The array of distr. func.
//                                  //  at all energies
//                 * pTableE;       //  The array of E values
         G4double  iQ2[ONQ2],     //  The array of Q2 values
                   pTableCrSec[ONQ2XE*AreaNumb],
                   pTableE[ONE*AreaNumb],
                   iIntgr[ONQ2],  //  The array of distr. func.   
                                  //    at one energy
                   Factorials1[250]; // The array for factorials
         G4double  dEbeg1, dEend1, dQ2, maxQ2;

  };     //   The end of the class description

#endif
