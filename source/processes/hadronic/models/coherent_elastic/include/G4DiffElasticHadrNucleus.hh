#ifndef G4DiffElasticHadrNucleus_h
#define G4DiffElasticHadrNucleus_h 1

#include "G4Nucleus.hh"
#include "G4DynamicParticle.hh"
#include "G4HadronValues.hh"

class G4DiffElasticHadrNucleus: public G4HadronValues
 {
    
  public:
    G4DiffElasticHadrNucleus() : G4HadronValues() {Factors();}

   ~G4DiffElasticHadrNucleus()   {;}

    G4double HadrNuclDifferCrSec(
                          const  G4DynamicParticle *  aHadron,
                                 G4Nucleus         *  aNucleus, 
                                 G4double            Q2);
    
    void GetNucleusParameters(G4Nucleus   * aNucleus);

  private:

// ++++++++++++++++++++++++++++++++++++++++++++++++++
    G4double binom(G4int N, G4int M)
      {       
          G4double  Fact1 = 1;
       if ((N>1) & (N>=M)) 
          {
              Fact1 = Factorials[N]/Factorials[M]/
                      Factorials[N-M];
          }
        return Fact1;
       }

// +++++++++++++++++++++++++++++++++++++++++++++++++++
   G4double Factorial(G4int N)
     {
        G4double  Res;
                  Res = 1;
              if(N == 0) return Res;

          if(N < 100) for(G4int M = 1; M<=N; M++)  Res = Res*M;         

          if(N >= 100)  Res = 2.50662827*exp(-N-1)*pow(N+1,N+0.5)*
                         (1+1/12/(N+1)+1/288/(N+1)/(N+1)-
                         139/51840/(N+1)/(N+1)/(N+1)-
                         571/2488320/(N+1)/(N+1)/(N+1)/(N+1));

              return Res;
            }

// +++++++++++++++++++++++++++++++++++++++++++++++++++
   void Factors()
    {
     G4int ii, ll, mm;
     G4double Sum1, Fac1, Fac3, Sum3;

          Factorials[0] = 1;

     for( ii = 1; ii<250; ii++)
        {
          if(ii >= 100) Mnoj[ii] = 3.03;  //  there is the saturation
          else
           {
           Sum1 = 0;
           Fac1 = 1;

               Factorials[ii] = Factorial(ii);

           for( ll = 0; ll<=ii; ll++)
             {
              Fac1 = binom(ii,ll);

                 Fac3 = 1;
                 Sum3 = 1;
                 for( mm = 1; mm<=ii-ll; mm++)
                   {
                     Fac3 = binom(ii-ll,mm);
                     Sum3 = Sum3 + 1/Fac3;
                   }  //  mm
                  Sum1 = Sum1 + Sum3/Fac1;
             }      //  ll
             Mnoj[ii] = Sum1;
           }       //  else
         }           //  ii
             Mnoj[0] = 1;

  }   //   Factors
//  --------------------------------------------------------

     G4double Mnoj[250], Factorials[250];

  public:
    
     G4double  R1, R2, Pnucl;
  };

#endif
