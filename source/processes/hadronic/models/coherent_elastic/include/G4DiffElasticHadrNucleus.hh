#ifndef G4DiffElasticHadrNucleus_h
#define G4DiffElasticHadrNucleus_h 1

#include "G4Nucleus.hh"
#include "G4DynamicParticle.hh"
#include "G4HadronValues.hh"

class G4DiffElasticHadrNucleus: public G4HadronValues
 {
    
  public:
    G4DiffElasticHadrNucleus() : G4HadronValues() {;}
   ~G4DiffElasticHadrNucleus()   {;}

    G4double HadrNuclDifferCrSec(const
                                 G4DynamicParticle*  aHadron,
                                 G4Nucleus *         aNucleus, 
                                 G4double            Q2);
    
  private:

    G4double binom(G4int N, G4int M)
      {
       
          G4double  Fact1 = 1;
       if ((N>1) & (N>=M)) 
          {
           G4int ii;
           for( ii = 1; ii<=M; ii++)
              Fact1 = Fact1*(N-ii+1)/ii;
          }
        return Fact1;
       }

   void Factors()
    {

     for(G4int ii = 0; ii<50; ii++)
        {
         G4double  Sum1 = 0;
         G4double  Fac1 = 1;

           for(G4int ll = 0; ll<=ii; ii++)
             {
              Fac1 = binom(ii,ll);

////              Fac1 := Fac1*ll/(ii-ll+1);
////              Sum1 := Sum1+Fac1*(ii-ll+1)/(ii-ll);
                G4double Fac3 = 1;
                G4double Sum3 = 1;
                 for(G4int mm = 1; mm<=ii-ll; mm++)
                   {
                     Fac3 = binom(ii-ll,mm);
                     Sum3 = Sum3 + 1/Fac3;
//                   Fac3 := Fac3*mm/(ii-ll-mm+1);
//                   Sum3 := Sum3 + Fac3;
                   }  //  mm
                  Sum1 = Sum1 + Sum3/Fac1;
             }      //  ll

             Mnoj[ii] = Sum1;
         }           //  ii
             Mnoj[0] = 1;
  }   //   Factors

    G4double Mnoj[50];

};

#endif
