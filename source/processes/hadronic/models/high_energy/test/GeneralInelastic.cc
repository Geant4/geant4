#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HEAntiKaonZeroInelastic.hh"
#include "G4HEKaonZeroShortInelastic.hh"
#include "G4HEKaonZeroLongInelastic.hh"
#include "G4HEKaonMinusInelastic.hh"
#include "G4HEProtonInelastic.hh"
#include "G4HEAntiProtonInelastic.hh"
#include "G4HENeutronInelastic.hh"
#include "G4HEAntiNeutronInelastic.hh"
#include "G4HELambdaInelastic.hh"
#include "G4HEAntiLambdaInelastic.hh"
#include "G4HESigmaPlusInelastic.hh"
#include "G4HESigmaZeroInelastic.hh"
#include "G4HESigmaMinusInelastic.hh"
#include "G4HEAntiSigmaPlusInelastic.hh"
#include "G4HEAntiSigmaZeroInelastic.hh"
#include "G4HEAntiSigmaMinusInelastic.hh"
#include "G4HEXiZeroInelastic.hh"
#include "G4HEXiMinusInelastic.hh"
#include "G4HEAntiXiZeroInelastic.hh"
#include "G4HEAntiXiMinusInelastic.hh"
#include "G4HEOmegaMinusInelastic.hh"
#include "G4HEAntiOmegaMinusInelastic.hh"

#include "G4HEPlot.hh"

int main()
{

  G4HEPionPlusInelastic theProtonInelastic;


  G4double A = 207.;
  G4double Z = 82.; 
  G4double ppp = 100.;

  G4int maxpart = 1000;
  G4HEVector incidentParticle;
  G4HEVector* pv = new G4HEVector[maxpart];
  G4int iteration = 0;
  G4int i, iprocess, vecL; 
  G4double cost, sint, phi, xfey, pt, pTot;
  G4String aPlotFile = "General.plot";
  G4HEPlot* plot = new G4HEPlot[100];
  G4HEVector* pana = new G4HEVector[10];

  for(i=0; i<=25; i++)
    {
      plot[i].Init(20, 0.9, 0.01);
      plot[30+i].Init(20, -1., 0.1);
      plot[60+i].Init(20,  0., 0.1);
    }
  
    do 
     {
       for(i=0; i<maxpart; i++)
         {
           pv[i].setZero();
         }        
       iprocess = (G4int)(G4UniformRand()*26.);
       if(iprocess == 0)
         {
           G4HEPionPlusInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("PionPlus");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 1)
         {
           G4HEPionMinusInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("PionMinus");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 2)
         {
           G4HEKaonPlusInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("KaonPlus");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 3)
         {
           G4HEKaonZeroInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("KaonZero");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 4)
         {
           G4HEAntiKaonZeroInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("AntiKaonZero");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 5)
         {
           G4HEKaonZeroShortInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("KaonZeroShort");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 6)
         {
           G4HEKaonZeroLongInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("KaonZeroLong");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 7)
         {
           G4HEKaonMinusInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("KaonMinus");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 8)
         {
           G4HEProtonInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("Proton");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 9)
         {
           G4HEAntiProtonInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("AntiProton");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 10)
         {
           G4HENeutronInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("Neutron");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 11)
         {
           G4HEAntiNeutronInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("AntiNeutron");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 12)
         {
           G4HELambdaInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("Lambda");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 13)
         {
           G4HEAntiLambdaInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("AntiLambda");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 14)
         {
           G4HESigmaPlusInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("SigmaPlus");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 15)
         {
           G4HESigmaZeroInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("SigmaZero");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 16)
         {
           G4HESigmaMinusInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("SigmaMinus");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 17)
         {
           G4HEAntiSigmaPlusInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("AntiSigmaPlus");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 18)
         {
           G4HEAntiSigmaZeroInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("AntiSigmaZero");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 19)
         {
           G4HEAntiSigmaMinusInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("AntiSigmaMinus");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 20)
         {
           G4HEXiZeroInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("XiZero");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 21)
         {
           G4HEXiMinusInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("XiMinus");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 22)
         {
           G4HEAntiXiZeroInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("AntiXiZero");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 23)
         {
           G4HEAntiXiMinusInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("AntiXiMinus");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 24)
         {
           G4HEOmegaMinusInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("OmegaMinus");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else if(iprocess == 25)
         {
           G4HEAntiOmegaMinusInelastic anInteraction;
           anInteraction.SetMaxNumberOfSecondaries(maxpart);
           cost = -1. + 2.*G4UniformRand();
           sint = sqrt(1. -cost*cost);
           phi  = G4UniformRand()*M_2PI;           
           incidentParticle.setDefinition("AntiOmegaMinus");
           incidentParticle.setMomentumAndUpdate(ppp*sint*cos(phi),
                                                 ppp*sint*sin(phi),
                                                 ppp*cost          ); 
           if(iteration < 10) anInteraction.SetVerboseLevel(2);
           else               anInteraction.SetVerboseLevel(0);
           anInteraction.ApplyYourself( pv, incidentParticle, A, Z);
           vecL = anInteraction.GetNumberOfSecondaries();
         }
       else
         {
           cout << " no process number " << iprocess << endl;
           continue; 
         }  
       pana[1] = incidentParticle;
       pana[2].setDefinition("Proton");
       pana[2].setMomentumAndUpdate(0., 0., 0.);
       pana[3].Add(pana[1], pana[2]);
       pana[4].Lor(pana[1], pana[3]);
       pana[6].setZero();
       for(i=0; i<vecL; i++)
         {
           pana[6].Add(pana[6],pv[i]);
           pana[5].Lor(pv[i], pana[3]);
           cost = pv[i].CosAng(pana[4]);
           xfey = cost*pv[i].Length()/pana[4].Length();
           pt   = sqrt(1.-cost*cost)*pv[i].Length();
           plot[30+iprocess].Fill(xfey, 1.);
           plot[60+iprocess].Fill(pt, 1.);
         }
       pTot = pana[6].Length()/pana[1].Length();
       plot[iprocess].Fill(pTot, 1.); 

       iteration++; 

     }  while (iteration < 100000) ;

   for(i=0; i<=25; i++)
     {
        plot[i].DumpToFile(i, aPlotFile);
        plot[30+i].DumpToFile(30+i, aPlotFile);
        plot[60+i].DumpToFile(60+i, aPlotFile);
     } 

}


