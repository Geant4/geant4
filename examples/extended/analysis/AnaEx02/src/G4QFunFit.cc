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
// $Id: G4QFunFit.cc,v 1.1 2006-01-11 15:43:52 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------- 

//#define debug
//#define tdebug
//#define edebug
//#define vdebug

#include "G4QFunFit.hh"

G4QFunFit::G4QFunFit(G4bool f, std::vector<G4double>* AI, std::vector<G4double>* DI,
                   		std::vector<G4double>* EX, G4VQFunction* usefun, G4int Np,
                     std::vector<G4double>* MI, std::vector<G4double>* MA)
:NS(Np),fulik(f),A(AI),DA(DI),AMN(MI),AMX(MA),EXDA(EX),fumfun(usefun),miF(false),maF(false)
{
  AM = 1.e37;                        // @@ Make a possibility to change by modifiers
  RP = 1.e-6;                        // @@ ...
  AMM= 1.e18;                        // @@ ...
  //#ifdef edebug
  if(AM<=0) G4cout<<"%Warning% G4QFunFit::: AM="<<AM<<" <= 0 (sqrt, 1./)"<<G4endl;
  //#endif
  NA=A->size();                      // a#of parameters
  NP=EXDA->size()/NS;                // Y,DY,Xi length  
  G4int I=0;                         // Prototype of the index
  if(!AMN)
  {
    AMN = new std::vector<G4double>; // create a new vector of minimum values
    for(I=0; I<NA; I++) AMN->push_back(-AM); // fill the minimum by the default values
    miF=true;
		}
  if(!AMX)
  {
    AMX = new std::vector<G4double>;          // create a new vector of minimum values
    for(I=0; I<NA; I++) AMX->push_back(AM);  // fill the minimum by the default values
    maF=true;
		}
  //#ifdef edebug
  G4int ND=DA->size();
  if(NA!=ND) G4cout<<"%Warning% G4QFunFit:::lengthOf A="<<NA<<", lengthOf DA="<<ND<<G4endl;
  //#endif  
  APR=AM/NA;                         // squared big number
		APS=std::sqrt(APR);                // big number
  AP=1./APR;                         // small number
  SIGMA = new std::vector<G4double>;
  PL0 = new std::vector<G4double>;
  DF = new std::vector<G4double>;
  for(I=0; I<NA; I++) PL0->push_back((*DA)[I]); // PL0==DA initialization
}

G4QFunFit::~G4QFunFit()
{
  if(miF) AMN->clear(); delete AMN;
  if(maF) AMX->clear(); delete AMX;
  SIGMA->clear(); delete SIGMA;
  PL0->clear(); delete PL0;
  DF->clear(); delete DF;
}

G4double G4QFunFit::FIT(G4int N1, G4int N2, G4int N3, G4double EPS, G4int IT)
{
  G4double S=0.;                     // Prototype of HI2/2
  //N1=2;                            // Step Reduction Limit @@ Make a Set/Get functions
  //N2=1;                            // Step Increase Limit @@ Make a Set/Get functions
  //#ifdef debug
  G4cout<<"G4QFunFit::FIT:S="<<S<<",N1="<<N1<<",N2="<<N2<<",N3="<<N3<<",EPS="<<EPS<<",IT="
        <<IT<<G4endl;
  //#endif
  if(EPS<0)
		{
    EPS=-EPS;                        // SIGMA  is always calculated
    G4cout<<"%Warning% G4QFunFit::FIT: Negative EPS="<<EPS<<" is made psitive"<<G4endl;
  }
  else if(!EPS)
		{
    EPS=0.01;                        // SIGMA  is always calculated
    G4cout<<"%Warning% G4QFunFit::FIT: Zero EPS is made standard: EPS=0.01"<<G4endl;
  }
  //#ifdef tdebug
  G4Timer* timer = new G4Timer();    // create a timer
  G4double tot_time=0.;              // total time calculation
  G4double ini_time=0.;              // initialization time
  G4int    iter=0;                   // number of iterations for the time counting
  timer->Start();                    // start the timer
  //#endif
  //#ifdef vdebug
  if(IT>0) // Title of program and definitions @@ Des not look necesary
		{
  		G4cout<<G4endl<<"G4QFunFit::FIT: Definitions: "<<G4endl;
  		G4cout<<"Function minimisation by subroutine G4QFunFit "
          <<"with following print-out"<<G4endl;
    G4cout<<"S = value of the objective function,"
          <<" EC = expected change in S during the next iteration"<<G4endl;
    G4cout<<"KAPPA = estimated distance to minimum,"
          <<" LAMBDA = step length modifier"<<G4endl<<G4endl;
  }
  //#endif
  G4int I=0;                         // Prototype of the index
  G4int FIXFLG=0;                    // initialize "#of fixed parameters" +
  G4int ENDFLG=0;                    // exit condition for the fit
  G4bool hasafix=false;              // INDFLG(2) just a flag of existing fixed parameters?
  G4int IFIX1=-1;                    // Prototype C++ for (0...) fixed parameter # (global)
  G4int NONFIX=-1;                   // Prototype C++ for (FI) counter of not fixed (?)
  G4int NN2=-1;                      // Instead of =0 (WHILE compensated)
  G4int NN3=-1;                      // Instead of =0 (WHILE compensated)
  R.clear();                         // PL-reduced diagonal elements of Z-matrix
  PL.clear();                        // Changable copy of the predefined PL0 Step Limiter
  SIGMA->clear();
  for(I=0; I<NA; I++)                // @@ vector functions can be used for initialization
  {
    R.push_back(0.);                 // Reset (initialize) working vector
    PL.push_back((*PL0)[I]);         // Copy start values of PL0 to variable PL
    SIGMA->push_back(0.);            // Initialize errors of parameters
  }
  G4double AIMAX=0.;                 // Fake prototype of the value of the fixed parameter
  G4int     IMAX=-1;                 // Prototype of the fixed parameter to fill by AIMAX
  G4double SP=0;                     // ProtOf SmallNumber inRespect to S (SR=1.E-6*ABS(S))
  G4double GT=0.;                    // Prototype of expectation of the HI^2 change on step
  G4bool first=false;                // First skip of the initial modification in WHILE
  while(!ENDFLG)                     // 3
  {
    //#ifdef debug
    G4cout<<"G4QFunFit::FIT:WHILE, N3="<<NN3<<",N="<<NA<<",I="<<IMAX<<",A="<<AIMAX<<G4endl;
    for(I=0; I<NA; I++) G4cout<<",A["<<I<<"]="<<(*A)[I];
    G4cout<<G4endl;
    for(I=0; I<NA; I++) G4cout<<",D["<<I<<"]="<<(*DA)[I];
    G4cout<<G4endl;
    //#endif
    if(first)
    {
      for(I=0; I<NA; I++) (*A)[I]+=(*DA)[I]; // Modify the parameter values
      if (IMAX>=0) (*A)[IMAX]=AIMAX; // Fix a parameter if the fixed parameter exixts
    }
    first=true;                      // Never skip any more
    fumfun->SetA(A);                 // Fill set of parameters
    G4double OLDS=S;
    NN2++;
    NN3++;
    G4double T1=1.;
    for(G4int NN1=0; NN1<N1; NN1++)  // Step Reduction LOOP 4
    {
      S=0.;
      G4int N=0;                     // counter of #of free parameters (starts from 0!)
      G.clear();                     // for reset
      PL.clear();                    // for reset
      for(I=0; I<NA; I++)
      {
        G.push_back(0.);             // Reset
        if((*PL0)[I]>0.)             // Skip pre-fixed parameters if negative
        {
          N++;                       // Count free parameters
										if(PL[I]>0.) (*PL0)[I]=PL[I]; // @@ ? PL0 is changed !
        }
      }
      G4int NN=N*(N+1)/2;            // New triangle matrix length for free parameters
      //#ifdef debug
      G4cout<<"G4QFunFit::FIT: N="<<N<<", NN="<<NN<<", S="<<S<<G4endl;
      //#endif
      Z0.clear();                    // Is calculated in SGZ
      if(NN) for(I=0; I<NN; I++) Z0.push_back(0); // Initialize Z0
      fixedpr=-1;                    // INDFLG(1)=0 == No fixed parameters
      // ----------------------- CALCULATE OBJECTIVE FUNCTION ----------------------
      S=SGZ();
      //#ifdef tdebug
      timer->Stop();                 // stop the timer
      G4double calc_time=timer->GetUserElapsed();
      G4double peri_time=0.;
      if(iter)
						{
        tot_time+=calc_time;
        peri_time=tot_time/iter;
      }
      else ini_time=calc_time;
      iter++;
      G4cout<<"G4QFunFit::FIT: it="<<iter<<", CalcTime="<<calc_time<<", init="<<ini_time
            <<", perIter="<<peri_time<<G4endl;
      timer->Start();                // start the timer
      //#endif
      SP=RP*std::fabs(S);            // Small number in respect to S (SR=1.E-6*ABS(S))
      //#ifdef debug
      G4cout<<"G4QFunFit::FIT: after SGZ, S="<<S<<", SP="<<SP<<", NS="<<NS<<G4endl;
      //#endif
      if(NN3>0&&NN1<=N1)             // Step reduction for not first iteration (19)
						{
        G4double W=S-OLDS-GT;
        G4double T=W+W;              // change in HI^2 with correction for the expectation
        //#ifdef debug
        G4cout<<"G4QFunFit::FIT: fpr="<<fixedpr<<",OS="<<OLDS<<",T="<<T<<",G="<<GT<<G4endl;
        for(I=0; I<NA; I++) G4cout<<",A["<<I<<"]="<<(*A)[I];
        G4cout<<G4endl;
        for(I=0; I<NA; I++) G4cout<<",D["<<I<<"]="<<(*DA)[I];
        G4cout<<G4endl;
        //#endif
        if(fixedpr>=0 || std::fabs(S-OLDS)>SP || -GT>SP && 0.59*T>=-GT)
        {
          if(fixedpr==-1)
										{
            if(!T) G4cout<<"G4QFunFit::FIT: T=0, fixedpr=-1"<<G4endl;
            T=-GT/T;
            if(T<0.25) T=0.25;
            //#ifdef debug
            G4cout<<"G4QFunFit::FIT: Precalc(no fixed pars) T="<<T<<", GT="<<GT<<G4endl;
          //#endif
          }
          else T=0.25;
          GT*=T;                     // Correction of the expectation
          T1*=T;                     // accumulate resulting step reduction (T1 compare 1)
          NN2=-1;                    // New life with the step increase strategy
          for(I=0; I<NA; I++) if(PL[I]>0.) // The NewReducedStep after unsuccessful attempt
          {
            (*A)[I]-=(*DA)[I];       // Recover the old values of parameters
            PL[I]*=T;                // Use new strategy with reduced steps
            (*DA)[I]*=T;             // Use new reduced step sizes
            (*A)[I]+=(*DA)[I];       // make the new reduced step
          }
          //#ifdef debug
          G4cout<<"G4QFunFit::FIT: After correction T="<<T<<G4endl;
          for(I=0; I<NA; I++) G4cout<<",A["<<I<<"]="<<(*A)[I];
          G4cout<<G4endl;
          for(I=0; I<NA; I++) G4cout<<",D["<<I<<"]="<<(*DA)[I];
          G4cout<<G4endl;
          //#endif
          fumfun->SetA(A);           // Fill set of parameters
          NN1++;                     // Increment a counter of reduced steps
        }
        else break;                  // Go out from the 4 Loop
      }
      else break;                    // Go out from the 4 Loop
    }                                // 4: End of the long (N1) Loop (step reduction limit)
    G4double ALAMBD=1.;              // Prototype of the Step Modification Factor
    G4double AKAPPA=0.;              // Prototype of the Estim Distance to Minimum
    G4bool lb77=true;                // Flag of emulation "go to 77"
    G4bool lb19=false;               // Prototype to enter in the "go to 19" while loop
    if(fixedpr>=0) ENDFLG=-4;        // Negative value is an algorithmic argument (?)
    else while(!lb19)                // C++ LOOP to emulate back go to 19
    {
      lb77=false;                    // kill the short cut
      lb19=true;                     // To exit in the end of the loop
      G4int K1=0;                    // pointer for Z0 matrix (first in new line)
      G4int K2=0;                    // pointer for Z matrix (increment after each equiv.)
      G4int I1=0;                    // line pointer for reduced (by PL0) Z0 matrix
      G4int I2=0;                    // line pointer for already reduced (by PL) Z matrix
      Z.clear();                     // Reduced matrix for only free parameters
      for(I=0; I<NA; I++)            // Count real number of free parameters to save space
      {
        if((*PL0)[I]>0.)             // Skip to use Z0
								{
          if(!PL[I])
          {
            PL[I]=(*PL0)[I];         // Recover the boundary fixed parameter
            //#ifdef debug
            G4cout<<"G4QFunFit::FIT: par # "<<I<<" is made free once more"<<G4endl;
            //#endif
          }
          if(PL[I]<=0.) K1+=I1+1;    // Parameter is lockally fixed, 1 line forward in Z0
          else
										{
            if((*A)[I]>=(*AMX)[I]&&G[I]<0.|| (*A)[I]<=(*AMN)[I]&&G[I]>0.)// par<minVpar>max
												{
              PL[I]=0.;              // fix the parameter temporary (for this step)
              K1+=I1;
            }
            else
												{
              for(G4int J=0; J<=I; J++) if ((*PL0)[J]>0.) // check global fixations
														{
                if (PL[J]>0.)        // check lockal fixations
																{
                  Z.push_back(Z0[K1]); // When it's checked -> get rid of K2
                  if(J==I&&(Z0[K1]<0.||K1+K1!=(I1+1)*(I1+2)-2||K2+K2!=(I2+1)*(I1+2)-2))
  		                G4cout<<"%Warning% G4QFunFit::FIT: "<<K1+1<<"="<<(I1+1)*(I+2)/2<<", Z("
                       <<Z.size()-1<<"="<<K2<<"="<<(I2+1)*(I2+2)/2-1<<")="<<Z0[K1]<<G4endl;
                  K2++;
                }
                K1++;
              }
              I2++;
            }
          } // 29
          I1++;
        }
      } // End of I<NA LOOP
      // ------------ INVERT Z  ???? for what ????? ----------------
      I1=0;                         // Start from 0 to calculate #of PL-reduced values
      G4int L=0;
      R.clear();
      JR.clear();
      for(I=0; I<NA; I++)
      {
        if(PL[I]>0)
						  {
          R.push_back(Z[L]);         // Extract diagonal elements of Z-matrix (PL reduced)
          JR.push_back(I);           // Remember full line pointers for reduced counting
          //#ifdef debug
          G4cout<<"G4QFunFit::FIT:R["<<I<<"]=Z["<<L<<"]="<<Z[L]<<",JR["<<I1<<"]="<<JR[I1]
                <<",PL="<<PL[I]<<G4endl;
          //#endif
          I1++;                      // Increment free parameters
          L+=I1+1;                   // Update diagonal pointer [0,2,5,9,14,...] (OK)
        } // End of I<NA LOOP (PL>0)
        else R.push_back(0.);        // Skip
      }
      G4int N=I1;                    // A#of free parameters (PL<=0 reduced) (@@ Why not N)
      //#ifdef debug
      G4cout<<"G4QFunFit::FIT: NA="<<NA<<", free="<<N;
      for(I=0; I<NA; I++) G4cout<<",A["<<I<<"]="<<(*A)[I];
      //for(I=0; I<NA; I++) G4cout<<",A["<<I<<"]="<<fumfun->GetA(I);
      G4cout<<G4endl;
      //#endif
#ifdef mdebug
      // Remember the direct matrix
      std::vector<G4double>* z = new std::vector<G4double>;
      G4int nz=Z.size();             // Length of the Z matrix
      for(I=0; I<nz; I++) z->push_back(Z[I]);
      G4cout<<"G4QFunFit::FIT: nz="<<nz<<", z0="<<(*z)[0]<<", z1="<<(*z)[1]<<G4endl;
#endif
      MCONV(N);                      // Z is the input, Z is the output (Z is kept in Z0)
#ifdef mdebug
      // Check that the result is a unit matrix
      G4int na=NA;
      G4cout<<"G4QFunFit::FIT:"<<nz<<"="<<NA*(NA+1)/2<<",Z0="<<Z[0]<<",Z1="<<Z[1]<<G4endl;
      while(nz<na*(na+1)/2)
      {
        G4cout<<"G4QFunFit::FIT: na="<<na<<", nz="<<nz<<"="<<na*(na+1)/2<<G4endl;
        na--;
      }
      G4cout<<"G4QFunFit::FIT: final na="<<na<<", nz="<<nz<<"="<<na*(na+1)/2<<G4endl;
      G4int ner=0;
      G4double er=0.;
      G4double ers=0.;
      for(I=0; I<na; I++)
      {
        for(G4int J=0; J<na; J++)
        {
          G4double sum=0.;
          for(G4int K=0; K<na; K++)
          {
            G4int I0=I;
            G4int K0=K;
            if(I0<K0)                // K0<=I0
												{
              K0=I;
              I0=K;
            }
            G4int J1=J;
            G4int K1=K;
            if(J1<K1)                // K1<=J1
												{
              K1=J;
              J1=K;
            }
            G4int IK=I0*(I0+1)/2+K0;
            G4int JK=J1*(J1+1)/2+K1;
            sum+=(*z)[IK]*Z[JK];
          }
										G4cout<<"G4QFunFit::FIT: sum["<<I<<","<<J<<"]="<<sum<<G4endl;
          if(I==J) er=std::fabs(1.-sum);
          else     er=std::fabs(sum);
          if(er>0.1)                 // 10% check
										{
            ner++;
            ers+=er;
          }
        }
        G4cout<<"G4QFunFit::FIT: *************** End of line ********************"<<G4endl;
      }
      if(ner)
						{
        G4cout<<"G4QFunFit::FIT: NErrors="<<ner<<", SE="<<ers<<G4endl;
        G4Exception("Matrix is not converted");
      }
      delete z;                      // ************** End Of Check Print *****************
#endif
      //#ifdef debug
      G4cout<<"G4QFunFit::FIT: fixedpr="<<fixedpr<<G4endl;
      //#endif
      if(fixedpr>=0)                 // INDFLG(1) Is any parameter fixed
						{
        fixedpr=-1;                  // Reinitialize a flag of fixed parameters (par #)
        hasafix=true;                // INDFLG(2) Put a flag of existing fixed parameter
        FIXFLG++;                    // 49 Increment #of fixed parameters
        //#ifdef debug
        G4cout<<"G4QFunFit::FIT: fixed prarameter is found FIXFLG="<<FIXFLG<<G4endl;
        //#endif
        NONFIX=0;
        lb19=false;
      }
      if(lb19)
      {
        G4double AKAP=0.;            // Fake C++ prototype
        G4double RR=0.;              // Fake C++ prototype
        //#ifdef debug
        G4cout<<"G4QFunFit::FIT:CalcTheorStep";
        for(I=0; I<NA; I++) G4cout<<",PL["<<I<<"]="<<PL[I];
        G4cout<<G4endl;
        //#endif
        // ---------------------------- CALCULATE THEORETICAL STEP TO MINIMUM ----------
        I1=0;
        for(I=0; I<NA; I++)
						  {
          RR=0.;
          if(PL[I]>0)
								  {
										  G4int L1=0;
            for(L=0; L<NA; L++)
										  {
              if(PL[L]>0)
								      {
                G4int K=0;
                if(I1>L1) K=(I1+1)*I1/2+L1;
														  else      K=(L1+1)*L1/2+I1;
                //#ifdef debug
                G4cout<<"G4QFunFit::FIT: I1="<<I1<<",G[L="<<L<<"]="<<G[L]<<",Z[K="<<K<<"]="
                      <<Z[K]<<G4endl;
                //#endif
														  if(G[L] && Z[K]) RR-=G[L]*Z[K];
                L1++;
              } 
            }
            I1++;
            //#ifdef debug
            if(!RR) G4cout<<"G4QFunFit::FIT: RR=0 for I="<<I<<G4endl;
            //#endif
          } // End PL>0 if
          //#ifdef debug
          G4cout<<"G4QFunFit::FIT:D["<<I<<"]="<<RR<<G4endl;
          //#endif
          (*DA)[I]=RR;
        } // End of Step Calculation LOOP
        //#ifdef debug
        G4cout<<"G4QFunFit::FIT: Check of parameters on boundary"<<G4endl;
        //#endif
        // ---------------------------- CHECK OF PARAMETERS ON BOUNDARY ------------
        G4double AFIX=0.;            // Prototype of max AKAP=ABS(DA(I)/SIGI)
        G4int    IFIX=-1;            // initialize the fixed parameter# (local used for PL)
        I1=0;                        // Line #
        L=0;                         // Diagonal element pointer
        for(I=0;I<NA;I++)if(PL[I]>0) // LOOP over free parameters
						  {
          G4double SIGI=1.;          // Default value
          if(Z[L]!=1.)
          {
            SIGI=std::sqrt(std::fabs(Z[L])); // always positive
            R[I]*=Z[L];              // diag_el_of_Covar * diag_el_of_Converted
          }
          (*SIGMA)[I]=SIGI;          // Error = sqrt(|diagonal_element|)
          //#ifdef debug
          G4cout<<"G4QFunFit::FIT:A["<<I<<"]="<<(*A)[I]<<",max="<<(*AMX)[I]<<",min="
                <<(*AMN)[I]<<", DA="<<(*DA)[I]<<",R="<<R[I]<<", L="<<L<<G4endl;
          //#endif
          if((*A)[I]>=(*AMX)[I]&&(*DA)[I]>0. || (*A)[I]<=(*AMN)[I]&&(*DA)[I]<0.) // OnBound
								  {
            if(!SIGI) G4cout<<"%Warning% G4QFunFit::FIT: SIGI=0"<<G4endl;
            AKAP=std::fabs((*DA)[I]/SIGI); // ProposedStep in ErrorUnits (to fix the worst)
            //#ifdef debug
            G4cout<<"G4QFunFit::FIT: AKAP="<<AKAP<<" V AFIX="<<AFIX<<G4endl;
            //#endif
            if(AKAP>AFIX)
										  {
              AFIX=AKAP;             // finallyfix a parameter with maximum AFIX
              IFIX=I;                // Prepare to fix only parameter wich comes beyond
              IFIX1=I;               // Remember for the future use
              //#ifdef debug
              G4cout<<"G4QFunFit::FIT: >>>> Parameter I="<<I<<" can be fixed <<<<"<<G4endl;
              //#endif
            }
								  } // 46
          I1++;                     // Line number in the ReducedMatrix Z
          L+=I1+1;                  // OK -> Pointer to diagonal elements of ReducedMatrix
        } // End if PL, end LOOP
        //#ifdef debug
        G4cout<<"G4QFunFit::FIT: Is IFIX="<<IFIX<<" >= 0 ?"<<G4endl;
        //#endif
        if(IFIX>=0)
						  {
          //#ifdef debug
          G4cout<<"G4QFunFit::FIT: par# "<<IFIX<<" is fixed, AFIX="<<AFIX<<G4endl;
          //#endif
          PL[IFIX]=-1.;
          FIXFLG++;
          NONFIX=0;
          lb19=false;
        } // 50
        if(lb19)                     // 49
        {
          //#ifdef debug
          G4cout<<"G4QFunFit::FIT: Calculate StepModification & EstimDistToMinHi2"<<G4endl;
          //#endif
          ALAMBD=1.;                 // Step modification factor: Min[AL=ABS(BI/DA(I))]
          AKAPPA=0.;                 // EstimDistance to min: Max[AKAP=ABS(DA[I]/SIGMA[I])]
          IMAX=-1;                   // Parameter # for max AL=max(A-AMN/AMX-A V PL)/DA
          for(I=0; I<NA; I++) if(PL[I]>0)
						    {
            G4double BM=(*AMX)[I]-(*A)[I];
            G4double ABI=(*A)[I]+PL[I];
            G4double ABM=(*AMX)[I];
            //#ifdef debug
            G4cout<<"G4QFunFit::FIT: DA["<<I<<"]="<<(*DA)[I]<<", BM="<<BM<<G4endl;
            //#endif
            if((*DA)[I]<=0.)
										  {
              BM=(*A)[I]-(*AMN)[I];
              ABI=(*A)[I]-PL[I];
              ABM=(*AMN)[I];
            }
            G4double BI=PL[I];
            //#ifdef debug
            G4cout<<"G4QFunFit::FIT: BI="<<BI<<" V BM="<<BM<<G4endl;
            //#endif
            if(BI>BM)
										  {
              BI=BM;
              ABI=ABM;
            }
            //#ifdef debug
            G4cout<<"G4QFunFit::FIT: ===> DA["<<I<<"]="<<(*DA)[I]<<", BI="<<BI<<G4endl;
            //#endif
            if(std::fabs((*DA)[I]) > BI)
										  {
              if(!(*DA)[I]) G4cout<<"%Warning% G4QFunFit::FIT: DA["<<I<<"]=0"<<G4endl;
              G4double AL=std::fabs(BI/(*DA)[I]);
              if(ALAMBD>AL)
												  {
                IMAX=I;
                AIMAX=ABI;
                ALAMBD=AL;
              }
            }
            if(!(*SIGMA)[I]) G4cout<<"%Warning% G4QFunFit::FIT: SIGMA["<<I<<"]=0"<<G4endl;
            AKAP=std::fabs((*DA)[I]/(*SIGMA)[I]);
            //#ifdef debug
            G4cout<<"G4QFunFit::FIT: AKAP="<<AKAP<<", AKAPPA="<<AKAPPA<<", I="<<I<<G4endl;
            //#endif
            if(AKAP>AKAPPA) AKAPPA=AKAP;
          } // End Pl if & end of I<NA LOOP
          //#ifdef debug
          G4cout<<"G4QFunFit::FIT: NewCorrectedStep AMM="<<AMM<<",ALAMBD="<<ALAMBD<<G4endl;
          //#endif
          // ---------------------------- CALCULATE NEW CORRECTED STEP --------------
          RR=0.;
          G4double AMB=AMM;          // Start with the big value
          if(ALAMBD>0.) AMB=0.25/ALAMBD;
          for(I=0; I<NA; I++)
						    {
            G4double RS=PL[I];       // Incrementable copy
            if(RS>0.)
										  {
              //#ifdef debug
              G4cout<<"G4QFunFit::FIT:"<<I<<",NN2="<<NN2<<" V N2="<<N2<<", R="<<(*DA)[I]/RS
                    <<" V AMB="<<AMB<<", DA["<<I<<"]="<<(*DA)[I]<<G4endl;
              //#endif
              if(NN2>N2 && std::fabs((*DA)[I]/RS)>=AMB)
												  {
                RS+=RS;
                PL[I]=RS+RS;         // PL[I]=4*PL[I]
                T1=4.;               // Multiplication factor
                //#ifdef debug
                G4cout<<"G4QFunFit::FIT:PL["<<I<<"] is multiplied ByFactorOf "<<T1<<G4endl;
                //#endif
              }
              //(*DA)[I]*=ALAMBD;    // @@ Modification !!
              RR+=(*DA)[I]*G[I];
              //#ifdef debug
              G4cout<<"G4QFunFit::FIT:Modified DA["<<I<<"]="<<(*DA)[I]<<",RR="<<RR<<G4endl;
              //#endif
            }
          } // End I<NA LOOP
          GT=RR;                     // Calculate expectation of HI^2 change on the step
          //---------------- CHECK IF MINIMUM ATTAINED AND SET EXIT MODE -----------------
          //#ifdef debug
          G4cout<<"G4QFunFit::FIT:GT="<<GT<<",SP="<<SP<<",T1="<<T1<<",AL="<<ALAMBD<<G4endl;
          //#endif
          if(-GT<=SP && T1<1. && ALAMBD<1.) ENDFLG=-1; // No further decrease of HI2
          //#ifdef debug
          G4cout<<"G4QFunFit::FIT: ENDFLG="<<ENDFLG<<", AKAPPA="<<AKAPPA<<" V EPS="<<EPS
                <<",NONFIX="<<NONFIX<<",FIXFLG="<<FIXFLG<<G4endl;
          //#endif
          if(ENDFLG>=0&&FIXFLG && (AKAPPA<EPS&&!ENDFLG || NONFIX>FIXFLG)) // go to 78
										{
            ENDFLG=0;
            //ENDFLG=1;
            FIXFLG=0;            // Makes all parameters free (no fixed parameters)
            IFIX1=-1;            // fixed parameter # (-1 = Reset) (global)
            //#ifdef debug
            G4cout<<"G4QFunFit::FIT: all parameters are free"<<G4endl;
            //#endif
            for(I=0; I<NA; I++) PL[I]=(*PL0)[I]; // Recover PL from scratch (PL0)
            hasafix=false;       // Reset a flag of existing fixed parameter
            lb19=false;
          }
						  } // label 19 skip 2
      } // label 19 skip 1
    } // while(lb19=false)(lable 19)
    //#ifdef debug
    G4cout<<"G4QFunFit::FIT:**WHILE, ENDFLG="<<ENDFLG<<", AKAPPA="<<AKAPPA<<" V EPS="<<EPS
										<<",FIXFLG="<<FIXFLG<<", lb77="<<lb77<<G4endl;
    //#endif
    if(ENDFLG>=0 || lb77)                 // NOT go to 77 (including far away short cut)
				{
      if(AKAPPA<EPS && !FIXFLG) ENDFLG=1; // go to 77 Done(ENDFLG=1)SmallStep(allPars EXIT)
      else if(AKAPPA<EPS && ENDFLG>=0 && IFIX1>0) // go to 76
						{
        NONFIX++;                         // 76 Increment bad tryes NONFIX
        ENDFLG=0;                         // (REPETE the WHILE)
      }
    }
    // 77
    if(!ENDFLG && NN3>=N3) ENDFLG=-3;      // OK, but the iterations are exhosted (EXIT)
    if(ENDFLG>0 && hasafix>0) ENDFLG=-2;  // OK, but atLeastOneParameter is fixed (EXIT)
    MONITO(S,NN3,IT,ENDFLG,GT,AKAPPA,ALAMBD);
    //if(ENDFLG==1 && FIXFLG) 
  } // while(!ENDFLAG) (3)
  //#ifdef tdebug
  delete timer;
  //#endif
  //#ifdef debug
  G4cout<<"G4QFunFit::FIT: return ENDFLAG="<<ENDFLG<<G4endl;
  //#endif
  if(ENDFLG<0) S=-S;                 // if S<0 ENDFLG can be read by GetENDFLG()
  return S;
}

// Inverts the positive definite triang-packed symmetric matrix Z by the square-root method
void G4QFunFit::MCONV(G4int N)
{
  if(N<1)                            // N here is length of the PL reduced cov matrix Z
		{
    G4cout<<"%Warning% G4QFunFit::MCONV: N="<<N<<" <1, matrix is not changed"<<G4endl;
    return;
  }
#ifdef debug
  else G4cout<<"G4QFunFit::MCONV: called, N="<<N<<G4endl;
#endif
  G4int NM1=N-1;                     // Pointer to the last row
  G4int I;                           // Prototype of the LOOP index
  G4int K;                           // Prototype of the LOOP index
  G4int L;                           // Prototype of the LOOP index
  G4int LL;                          // Prototype
  G4int LI;                          // Prototype
  G4int NL;                          // Prototype
  for(I=0; I<N; I++)
		{
#ifdef debug
    G4cout<<"G4QFunFit::MCONV: JR="<<JR.size()<<", Z="<<Z.size()<<", R="<<R.size()<<G4endl;
#endif
    G4int IR=JR[I];                  // Recover "not PL-reduced pointer" (I is PL reduced)
    G4int NI=(I+1)*I/2;              // pointer to the first element in the line
    G4int II=NI+I;                   // Diagonal element of this line
				//#ifdef debug
    G4cout<<"G4QFunFit::MCONV: JR["<<I<<"]="<<JR[I]<<", Z["<<II<<"]="<<Z[II]<<", R["
          <<IR<<"]="<<R[IR]<<G4endl;
				//#endif
    if(Z[II] < RP*std::fabs(R[IR]) || Z[II] < AP)
				{
						//#ifdef debug
      G4cout<<"G4QFunFit::MCONV:Small Z["<<II<<"]="<<Z[II]<<",R["<<IR<<"]="<<R[IR]<<G4endl;
						//#endif
      MEXIT(IR,N,II);                // Example of go to 19
      return;
    }
    // Relatively small diagonal element = V = abs small element => skip additions
    if(Z[II]<=0.) G4cout<<"%Warning% G4QFunFit::MCONV: Z["<<II<<"] <= 0, AP="<<AP<<G4endl;
    G4double ZII=1./std::sqrt(Z[II]);// Convert diagonal element for normalization
    Z[II]=ZII;                       // Replace the diagonal element
    if(II>NI)for(NL=II-1;NL>=NI;NL--)// Norm all except for the diagonal element
				{
      if(Z[NL]) Z[NL]*=ZII;         // Normalize (#0) elements of the row
#ifdef debug
      G4cout<<"G4QFunFit::MCONV: Z["<<NL<<"]="<<Z[NL]<<" or APS="<<APS<<G4endl;
#endif
      if(std::fabs(Z[NL]) >= APS)    // Z(NL) is too big (can't estimate errors)
						{
								//#ifdef debug
        G4cout<<"G4QFunFit::MCONV: Z["<<NL<<"]="<<Z[NL]<<", NI="<<NI<<G4endl;
								//#endif
        // JR is pointers in a full line for reduced counting
        IR=JR[NL-NI-1];              // IR=JR(NL-NI) ? @@ CHECK !
        II=NL;                       // for the proper print
        MEXIT(IR,N,NL);              // (go to 19) II=NL just for the print
        return;        
      }
    }
#ifdef debug
    G4cout<<"G4QFunFit::MCONV: I="<<I<<" or NM1="<<NM1<<G4endl;
#endif
    if(NM1>I)                        // not for the last in the loop
    {
      for(K=NM1; K>I; K--)           // through the rest of rows above I
						{
      		G4int NK=(K+1)*K/2;          // Pointer to the begin of the row
        NL=NK;                       // Modified dublicate for the L-LOOP
        G4int KK=NK+I;               // I-th member in K-th line (K>I not yet normed)
        G4double D=Z[KK]*ZII;        // norm the rest of the line, but remember it later
#ifdef debug
        G4cout<<"G4QFunFit::MCONV: K="<<K<<", Z["<<KK<<"]="<<Z[KK]<<", D="<<D<<G4endl;
#endif
								if(D)
								{
          G4double C=D*ZII;          // overnorming for second L>I element
#ifdef debug
          G4cout<<"G4QFunFit::MCONV: I="<<I<<", K="<<K<<", C="<<C<<G4endl;
#endif
          for(L=K; L>I; L--)         // this K >=L >I part includes diagonal element (K,K)
										{
            LL=NK+L;                 // L-th member in K-th line
            LI=NL+I;                 // I-th member in L-th line
#ifdef debug
            G4cout<<"G4QFunFit::MCONV: L="<<L<<", Z["<<LL<<"]="<<Z[LL]<<", Z["<<LI<<"]="
                  <<Z[LI]<<G4endl;
#endif
            Z[LL]-=Z[LI]*C;          // Calculate a diagonal product
            NL-=L;                   // Jump back to previous line (@@ ?)
          }
          // ------------ L=I is skiped ---------- it is taken into account az Z(LL) itself
          if(I>0) for(L=I; L>=0; L--)// Not for the firs, previous part
										{
            LL=NK+L;                 // L-th member in the K-th line
            LI=NI+L;                 // L-th element in the I-th line (already normed)
            Z[LL]-=Z[LI]*D;          // Continue calculating the diagonal product
          }
          // As a result <L-th in K-th line> = <it> * <L-th line>
          Z[KK]=-C;
        }
        else Z[KK]=0.;
      }
    }
  }
  for(I=0; I<N; I++) for(K=I; K<N; K++)
		{
    NL=(K+1)*K/2;
    G4int KI=NL+I;
    G4double E=0.;
#ifdef debug
    G4cout<<"G4QFunFit::MCONV: I="<<I<<", NL="<<NL<<", KI="<<KI<<G4endl;
#endif
    for(L=K; L<N; L++)
				{
      LI=NL+I;
      G4int LK=NL+K;
#ifdef debug
      G4cout<<"G4QFunFit::MCONV:"<<L<<",["<<LI<<"]="<<Z[LI]<<",["<<LK<<"]="<<Z[LK]<<G4endl;
#endif
      E+=Z[LI]*Z[LK];
      NL+=L+1;
    }
    Z[KI]=E;
  }
  return;
}

void G4QFunFit::MEXIT(G4int IR, G4int N, G4int II) // Emergency exit action for MCONV
{
  PL[IR]=-2;                         // here we really need to know an address in full row
  //#ifdef debug
  G4cout<<"G4QFunFit::MEXIT(MCONV): par#"<<IR<<" is fixed N="<<N<<", Z["<<II<<"]="<<Z[II]
        <<",R="<<R[IR]<<G4endl;
  //#endif
  R[IR]=0.;
  fixedpr=IR;                        // filled by a number in a full row !
  return;
}

G4double G4QFunFit::SGZ()        // S=ObjectiveFunction, G=GradientOf S, Z=CovarianceMatrix
{
  G4double S=0.;                     // Hi2/2
  G4int K=NS;                        // A number of the Experimental Points
  G4int K2=0;                        // Initial pointer to the current Point 
  G4int I=0;                         // Prototype of the index
  G4int J=0;                         // Prototype of the index
  std::vector<G4double>* X = new std::vector<G4double>; // vector of arguments
  for(G4int L1=0; L1<K; L1++)        // LOOP over points
		{
    G4int K1=K2;                     // Assistant pointer for the Likelihood
    G4int NX=NP-2;                   // A#of arguments
#ifdef debug
    G4cout<<"G4QFunFit::SGZ: L1="<<L1<<", K2="<<K2<<", NX="<<NX<<",IF2="<<fulik<<G4endl;
#endif
    if(fulik)                        // Likihood: Likelehood does not have y, dy (?)
				{
      NX=NP;                         // All line is only arguments (no y, dy)
      K1-=2;                         // To shift the copy-pointer by the y, dy places
    }
    G4int K12=K1+2;                  // Shift to arguments
    for(I=0; I<NX; I++)              // Copy arguments
				{
      if(L1) (*X)[I]=(*EXDA)[K12+I]; // Just a copy @@ Copy directly to fumfun->SetX(I,*)
      else X->push_back((*EXDA)[K12+I]); // First time initialize the empyt X-vector (@@X?)
    }
#ifdef debug
    G4cout<<"G4QFunFit::SGZ: K1="<<K1<<", FLIK="<<fulik;
				for(J=0; J<NX; J++) G4cout<<", X["<<J<<"]="<<(*X)[J]; // (@@X?)
				for(J=0; J<NA; J++) G4cout<<", P["<<J<<"]="<<(*PL0)[J];
    G4cout<<G4endl;
#endif
    fumfun->SetX(X);                 // Fill set of arguments
    G4double Y=fumfun->Y(DF,PL0);    // Transfer the DF-pointer to avector of derivitives
#ifdef debug
    G4cout<<"G4QFunFit::SGZ: FLIK="<<fulik<<", Y="<<Y;
    for(J=0; J<NX; J++) G4cout<<", DF["<<J<<"]="<<(*DF)[J];
    G4cout<<G4endl;
#endif
    G4double RR=0.;                  // Prototype of shift/error
    G4double SIG=Y;                  // Prototype    
    if(fulik)                      // In case of not Likelihood (INDFLG(3)#0)
				{
      if(Y<=0)                       // Negative Likelihood - abandon attempt
						{
        fixedpr=0;                   // -VE or Y=0 in maximum liklehood method
        S=1.e127;                    // (value ?)
      }
      S-=std::log(Y);                // Collect the likelyhood function
      Y=-Y;                          // Y is not used in this function any more!
      SIG=Y;                         // Can error be negative?
      RR=1.;                         // Shift/error ratio is formal
    }
    else                             // Minimization of HI^2 (INDFLG[2]=INDFLG(3)=0)
				{
      SIG=(*EXDA)[K2+1];             // Error of Y (second value in the data line)
      Y-=(*EXDA)[K2];                // Shift exY from theorY (K1 in original,K1=K2 forHi2)
      if(!SIG)G4cout<<"%Warning% G4QFunFit::SGZ:***SIG=0***, K1="<<K1<<", K2="<<K2<<G4endl;
      RR=Y/SIG;                      // shift/error ratio (can be negative)
      S+=RR*RR/2;                    // S=HI^2/2 (!)
#ifdef debug
      G4cout<<"G4QFunFit::SGZ:****> S="<<S<<", dY="<<Y<<", SIG="<<SIG<<", RR="<<RR<<G4endl;
#endif
    }
    G4int N=-1;                      // Counter of the active parameters (N-tot=N_max+1)
    for(J=0; J<NA; J++) if((*PL0)[J]>0) // LOOP only over active parameters
				{
      N++;                           // Counted as an active parameter
      if((*DF)[J])                   // At this point a derivitive for the parameter J != 0
						{
        (*DF)[N]=(*DF)[J]/SIG;       // Norm & compress DF. Forget about NA after this LOOP
        if(RR) G[J]+=(*DF)[N]*RR;    // G-matrix=dhHi2/da_j=Sum_i[(dy/da(j))*dy_i/sig_i^2]
      }
      else (*DF)[N]=0.;              // Active parameter does not change the function Y
    }
    N++;                             // C++ correction (N_tot)
    G4int L=0;                       // Pointer for the Z0-matrix (triangular, condenced)
    if(N) for(I=0; I<N; I++)         // LOOP over the active parameters
				{
      G4double RS=(*DF)[I];          // Get the derivitive for the parameter
      if(RS)
						{
								for(J=0; J<=I; J++)
        {
          RR=(*DF)[J];
          if(RR) Z0[L]+=RS*RR;
          L++;
        }
      }
      else L+=I+1;                   // (+1 ?) Skip all row
    }
    K2+=NP;
  }
  //#ifdef debug
  G4cout<<"G4QFunFit::SGZ:=>S="<<S<<", Z0="<<Z0[0]<<", Z2="<<Z0[2]<<", Z5="<<Z0[5]<<G4endl;
  //#endif
  delete X;                          // (@@X?)
  return S;
}

void G4QFunFit::ERRORF()             // Additional calculation of HI2 & errors
{
  //#ifdef vdebug
  G4cout<<"G4QFunFit::ERRORF: Errors for each DataPoint using FinalParameterValues"<<G4endl
        <<"POINT    FITTED    Y STANDARD  CONTRIBUTION  X CO-ORDINATES"<<G4endl
        <<"  #     Y  VALUE    DEVIATION   TO CHISQ/2   OF DATA POINTS"<<G4endl;
  //#endif
  std::vector<G4double>* X = new std::vector<G4double>; // vector of arguments
  G4int K1=-NP;
  ERROR.clear();                     // Reset the error vector (@@)
  G4int I=0;
  for (G4int J=0; J<NS; J++)         // LOOP over points
		{
    K1+=NP;                          //  Data line pointer
    G4double FUN=(*EXDA)[K1];        //  Measured function value
    G4double SIG=(*EXDA)[K1+1];      //  Function error value
    if(!SIG) G4cout<<"%Warning% G4QFunFit::ERRORF: SIG=0 for K1="<<K1<<G4endl;
    G4int NX=NP-2;                   // Length of the argument vector (@@ struct)
    G4int K12=K1+2;                  // Shift to arguments
    for(I=0; I<NX; I++)              // Copy arguments in the X vector
				{
      if(J) (*X)[I]=(*EXDA)[K12+I];   // Just a copy @@ Copy directly to fumfun->SetX(I,*)
      else X->push_back((*EXDA)[K12+I]); // First time initialize the empyt X-vector (@@X?)
    }
    fumfun->SetX(X);                 // Fill set of arguments
    G4double Y=fumfun->Y(DF,PL0);    // Transfer the DF-pointer to avector of derivitives
    G4double ER=SCAL();              // Calculates ER = variance of theoretical function Y
    if(ER<0.) G4cout<<"%Warning% G4QFunFit::ERRORF: ER="<<ER<<" < 0"<<G4endl;
    ER=std::sqrt(ER);                // Diagonal elements are sigma=dx^2
    ERROR.push_back(ER);             // Make a vector of errors (@@ For what?)
    G4double XI=(FUN-Y)/SIG;         // Intermediate acceleration
    XI=XI*XI/2;                      // Hi^2/2
    //#ifdef vdebug
    G4cout<<" "<<J<<"  "<<Y<<"  "<<ER<<"  "<<XI;
    for(I=0; I<NX; I++) G4cout<<" "<<(*X)[I]; // (@@X?)
    G4cout<<G4endl;                  // @@ fit spaces
    //#endif
  }
  delete X;                          // (@@X?)
  return;
}

G4double G4QFunFit::SCAL()           // Calculates ER = variance of theoretical function Y
{
  G4double ER=0.;
  G4int N=0;                         // Counter of the free parameters
  G4int I=0;                         // Prototype of the index
  G4int K=0;                         // Prototype of the index
  for(I=0; I<NA; I++) if(PL[I]>0)    // Count free (not fixed) parameters
		{
    DF[N]=DF[I];
    N++;                             // @@ make the same trick in SGZ
  }   
  if(!N) return ER;
  for(I=0; I<N; I++) for(G4int J=0; J<N; J++)
		{
    if(I<J) K=(J+1)*J/2+I;           // Symmetric natrix solution for the triangular packed
    else    K=(I+1)*I/2+J;
    ER+=Z[K]*(*DF)[I]*(*DF)[J];
  }
  return ER;
}

void G4QFunFit::MONITO(G4double S, G4int NN3, G4int IT, G4int ENDFLG, G4double GT,
                    G4double AKAPPA, G4double ALAMBD) // IterationProgressOutput @@PutInFIT
{
  G4int NM=0;                        // C++ initialization
  G4int I=0;                         // Prototype of the index (& more)
  if(IT>=0 && ((IT && (NN3<=0 || NM>=0)) || ENDFLG))
		{
    //#ifdef vdebug
    G4cout<<"ITERATION #"<<NN3<<", S="<<S<<", EC="<<GT<<", KAPPA="<<AKAPPA<<", LAMBDA="
          <<ALAMBD<<G4endl<<"PARAMETER  PARAMETER   STANDARD   CORRELATION"
                  <<G4endl<<" NUMBER      VALUE     DEVIATION   FACTOR"<<G4endl;
    //#endif
    for(I=0; I<NA; I++)
				{
      if((*PL0)[I]>0)
						{
        if(PL[I]>0)
          G4cout<<"    "<<I<<"      "<<(*A)[I]<<"     "<<(*SIGMA)[I]<<"  "<<R[I]<<G4endl;
        else if(PL[I]>-2) G4cout<<"    "<<I<<"      "<<(*A)[I]<<"     "<<(*SIGMA)[I]<<" "
                                <<R[I]<<" ON BOUNDARY"<<G4endl;
        else G4cout<<"    "<<I<<"      "<<(*A)[I]<<"  => INFINITE ERROR ESTIMATED"<<G4endl;
      }
      else G4cout<<"    "<<I<<"      "<<(*A)[I]<<"  => THIS PARAMETER IS FIXED"<<G4endl;
    }
    NM=-IT;
  }
  else if(IT<0) NM=-IT;
  NM++;
  if(ENDFLG<0 && IT>0)
		{
				G4cout<<"%Warning% G4QFunFit: MINIMISATION IS TERMINATED: ";
    switch(ENDFLG)
    {
				case -1:
						G4cout<<"No further decrease of S=HI^2/2"<<G4endl;
						break;
				case -2:
						G4cout<<"Infinite error is estimated"<<G4endl;
						break;
				case -3:
						G4cout<<"Iteration limit is reached"<<G4endl;
						break;
				case -4:
						G4cout<<"Y<=0 is a logarithmic argument"<<G4endl;
						break;
				default:
						G4cout<<"Unknown error ENDFLG="<<ENDFLG<<G4endl;
						break;
    }
  }
  return;
}
