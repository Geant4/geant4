// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SimpleIntegration.cc,v 1.2 1999-11-16 17:31:11 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Implementation file for simple integration methods
//

#include "G4SimpleIntegration.hh"


G4int G4SimpleIntegration::fMaxDepth = 100 ;


G4SimpleIntegration::G4SimpleIntegration( function pFunction )
{
   fFunction = pFunction ;
}

G4SimpleIntegration::G4SimpleIntegration( function pFunction,
					  G4double pTolerance)
{
   fFunction = pFunction ;
   fTolerance = pTolerance ;
}


G4SimpleIntegration::~G4SimpleIntegration() 
{
   ;
}
       
       // Simple integration methods
       
G4double
G4SimpleIntegration::Trapezoidal(G4double xInitial,
                                 G4double xFinal,
			         G4int iterationNumber ) 
{
   G4int i ;
   G4double Step = (xFinal - xInitial)/iterationNumber ;
   G4double mean = (fFunction(xInitial) + fFunction(xFinal))*0.5 ;
   G4double x = xInitial ;
   for(i=1;i<iterationNumber;i++)
   {
      x += Step ;
      mean += fFunction(x) ;
   }
   return mean*Step ;
}

G4double 
G4SimpleIntegration::MidPoint(G4double xInitial,
                              G4double xFinal,
			      G4int iterationNumber ) 
{
   G4int i ;
   G4double Step = (xFinal - xInitial)/iterationNumber ;
   G4double x = xInitial + 0.5*Step;
   G4double mean = fFunction(x) ;
   for(i=1;i<iterationNumber;i++)
   {
      x += Step ;
      mean += fFunction(x) ;
   }
   return mean*Step ;   
}

G4double      
G4SimpleIntegration::Gauss(G4double xInitial,
                           G4double xFinal,
			   G4int iterationNumber ) 
{
   G4int i ;
   G4double x ;
   static G4double root = 1.0/sqrt(3.0) ;
   G4double Step = (xFinal - xInitial)/(2.0*iterationNumber) ;
   G4double delta = Step*root ;
   G4double mean = 0.0 ;
   for(i=0;i<iterationNumber;i++)
   {
      x = (2*i + 1)*Step ;
      mean += (fFunction(x+delta) + fFunction(x-delta)) ;
   }
   return mean*Step ;   
}

G4double    
G4SimpleIntegration::Simpson(G4double xInitial,
                             G4double xFinal,
			     G4int iterationNumber ) 
{
   G4int i ;
   G4double Step = (xFinal - xInitial)/iterationNumber ;
   G4double x = xInitial ;
   G4double xPlus = xInitial + 0.5*Step ;
   G4double mean = (fFunction(xInitial) + fFunction(xFinal))*0.5 ;
   G4double sum = fFunction(xPlus) ;
   for(i=1;i<iterationNumber;i++)
   {
      x     += Step ;
      xPlus += Step ;
      mean  += fFunction(x) ;
      sum   += fFunction(xPlus) ;
   }
   mean += 2.0*sum ;
   return mean*Step/3.0 ;   
}



       // Adaptive Gauss integration
       
G4double       
G4SimpleIntegration::AdaptGaussIntegration( G4double xInitial,
                                            G4double xFinal   ) 
{
   G4int depth = 0 ;
   G4double sum = 0.0 ;
   AdaptGauss(xInitial,xFinal,sum,depth) ;
   return sum ;
}
 			    

G4double
G4SimpleIntegration::Gauss( G4double xInitial,
                            G4double xFinal   ) 
{
   static G4double root = 1.0/sqrt(3.0) ;
   
   G4double xMean = (xInitial + xFinal)/2.0 ;
   G4double Step = (xFinal - xInitial)/2.0 ;
   G4double delta = Step*root ;
   G4double sum = (fFunction(xMean + delta) + fFunction(xMean - delta)) ;
   
   return sum*Step ;   
}


void    
G4SimpleIntegration::AdaptGauss( G4double xInitial,
                                 G4double xFinal,
				 G4double& sum,
			         G4int& depth      ) 
{
   if(depth >fMaxDepth)
   {
      G4Exception("Function varies too rapidly in G4SimpleIntegration::AdaptGauss") ;
   }
   G4double xMean = (xInitial + xFinal)/2.0 ;
   G4double leftHalf = Gauss(xInitial,xMean) ;
   G4double rightHalf = Gauss(xMean,xFinal) ;
   G4double full = Gauss(xInitial,xFinal) ;
   if(fabs(leftHalf+rightHalf-full) < fTolerance)
   {
      sum += full ;
   }
   else
   {
      depth++ ;
      AdaptGauss(xInitial,xMean,sum,depth) ;
      AdaptGauss(xMean,xFinal,sum,depth) ;
   }
}
