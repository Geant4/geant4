// Test- Job simulates the experimental data given in
// 
// 
// 
// 
// 
// 
// 

#include "G4HEProtonInelastic.hh"
#include "G4HEPlot.hh"

int main()
{
//*******************************************************

  G4int NumberEvents = 1000000;

  G4String aFile = "ProtonAu.plot";

  G4HEProtonInelastic theInteraction;

  G4HEVector incidentParticle;
  incidentParticle.setDefinition("Proton");

//******************************************************* 
        
  G4double A, Z;
  A = 196.97;
  Z = 79.;
  G4int maxpart = 1000;
  G4HEVector* pv = new G4HEVector[maxpart];
  incidentParticle.setMomentumAndUpdate(0., 0., 200.); 
  theInteraction.SetMaxNumberOfSecondaries(maxpart);
  theInteraction.SetVerboseLevel(2);
  G4int i = 0;
  G4HEPlot* plot = new G4HEPlot[5]; 
  G4HEVector* pana = new G4HEVector[11];

  plot[ 1].Init(50, 0., 1.);
  plot[ 2].Init(50, 0., 0.1);
  plot[ 3].Init(50, 0., 0.1);
  plot[ 4].Init(50, 0., 0.1);
  G4int nevhmf = 0;
  G4int npi0hf = 0;

  pana[1] = incidentParticle;
  pana[2].setMass(0.9383);
  pana[2].setKineticEnergyAndUpdate(0.);
  pana[10].Add(pana[1], pana[2]);
  pana[3].Lor(pana[1], pana[10]);

  G4int iteration = 0;   
    do 
     {
       for(i=0; i<maxpart; i++)
	 {
           pv[i].setZero();
         }
       if(iteration == 5) theInteraction.SetVerboseLevel(0);
       theInteraction.ApplyYourself( pv, incidentParticle, A, Z);
       G4int vecLength = theInteraction.GetNumberOfSecondaries();
       G4int ngk = vecLength;
       G4int nchhmf = 0;
       G4int nhmf = 0; 
       
       for(i=0; i<ngk; i++)
         {
           if(pv[i].getName() == "SigmaZero")
             {
               pana[8].setZero();
               pana[9].setZero();
               pana[7] = pv[i];
               pana[7].Smul(pana[7],-1.);
               G4double ppp = pv[i].Length();
               pana[8].setDefinition("Lambda");
               pana[9].setDefinition("Gamma");
               G4double phmf = 0.074;       
               pana[8].setEnergy(sqrt(sqr(phmf)+sqr(pana[8].getMass())));
               pana[9].setEnergy(phmf);
               G4double cost = -1. + G4UniformRand()*2.;
               G4double sint = sqrt(1.-cost*cost);
               G4double phi = 6.28318531*G4UniformRand();
               pana[8].setMomentum(phmf*sint*sin(phi),
                                   phmf*sint*cos(phi),
                                   phmf*cost);
               pana[9].Smul(pana[8],-1.);
               pana[8].Lor(pana[8],pana[7]);
               pana[9].Lor(pana[9],pana[7]);
               pv[i] = pana[8];
               pv[ngk] = pana[9];
               ngk++;
	     }   
         }

       for(i=0; i<ngk; i++)
         {
           nhmf++;
           if(fabs(pv[i].getCharge()) > 0.5) nchhmf++;
         }
       if(nhmf <= 2) continue;
       nevhmf++;
       plot[1].Fill((G4double)nchhmf +0.5, 1.);       

       for(i=0; i<ngk; i++)
         {
           G4int ioff;
           if(pv[i].getName() == "PionZero") ioff = 1;
           else if(pv[i].getName() == "PionPlus") ioff = 2;
           else if(pv[i].getName() == "PionMinus") ioff = 3;
           else ioff = 0; 
           if(ioff > 0)
             {
               G4double tetahf = pv[i].Ang(pana[1]);
               if(tetahf < 1.e-6)
                 {
                   tetahf += fabs(theInteraction.normal())*1.e-6; 
                 }
               G4double pt2 = sqr(pv[i].getMomentum().x()) + sqr(pv[i].getMomentum().y());
               G4double wgthmf = 0.;
               if(pt2 > 1.e-3) wgthmf = 1./pt2;
               if(tetahf > 0.239111 && tetahf < 0.452040) plot[1+ioff].Fill(pt2,wgthmf);
             } 
         }               


       iteration++; 

     }  while (iteration < NumberEvents) ;

    cout << " number of useful events " << nevhmf << G4endl; 
   
    G4double wgthmf = 1./nevhmf;
    plot[1].Scale(wgthmf, plot[1]);
    wgthmf = 1800./(nevhmf*0.6*0.10*3.1415926*2.);
    wgthmf *= 1.3;    
    plot[2].Scale(wgthmf, plot[2]);    
    plot[3].Scale(wgthmf, plot[3]);    
    plot[4].Scale(wgthmf, plot[4]);    

    for (i=1; i<4; i++)
      {
        plot[i].DumpToFile( i, aFile);
      }

    delete pana;
    delete plot;
    delete pv; 
}


