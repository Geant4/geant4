// Test- Job simulates the experimental data given in
// M.R.Atayan et al., Z.Phys.C54, 247 (1992);
// The same comparison has been done with Geant3 MC
// H. Fesefeldt, August 1998.

#include "G4HEPionPlusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEPlot.hh"

int main()
{
//*******************************************************

  G4int NumberEvents = 100000;
  G4double incidentMomentum = 250.;
 
// PionPlus, KaonPlus at 250 GeV

  G4String aFile = "PionPlus.plot";

  G4HEPionPlusInelastic theInteraction;

  G4HEVector incidentParticle;
  incidentParticle.setDefinition("PionPlus");

//******************************************************* 
        
  G4double A, Z;
  A = 1.01;
  Z = 1.;
  G4int maxpart = 1000;
  G4HEVector* pv = new G4HEVector[maxpart];
  incidentParticle.setMomentumAndUpdate(0., 0., incidentMomentum); 
  theInteraction.SetMaxNumberOfSecondaries(maxpart);
  theInteraction.SetVerboseLevel(2);
  G4int i = 0;

  G4HEPlot* plot = new G4HEPlot[6]; 
  G4HEVector* pana = new G4HEVector[11];
  G4double* help7 = new G4double[maxpart];
  G4double* help8 = new G4double[maxpart];
  G4double* help9 = new G4double[maxpart];

  plot[ 1].Init(40, -1., 0.05);
  plot[ 2].Init(40, -1., 0.05);
  plot[ 3].Init(40, -1., 0.05);
  plot[ 4].Init(40, -1., 0.05);
  plot[ 5].Init(50,  0., 0.05);

  G4int nhmf = 0;
  G4int nevhmf = 0;

  pana[1] = incidentParticle;
  pana[2].setMass(0.9382);
  pana[2].setCharge(1.);
  pana[2].setKineticEnergyAndUpdate(0.);
  pana[3].Add(pana[1], pana[2]);
  pana[4].Lor(pana[1], pana[3]);

  G4int iteration = 0;   
    do 
     {
       for(i=0; i<maxpart; i++)
	 {
           pv[i].setZero();
           help7[i] = 0.;
           help8[i] = 0.;
           help9[i] = 0.;
         }
       if(iteration == 5) theInteraction.SetVerboseLevel(0);
       theInteraction.ApplyYourself( pv, incidentParticle, A, Z);
       G4int vecLength = theInteraction.GetNumberOfSecondaries();
       G4int ngk = vecLength;

       nhmf = 0;
       for(i=0; i<ngk; i++)
         {
           pana[5].Lor(pv[i],pana[3]);
           help7[i] = pana[5].getMomentum().z()/pana[4].getMomentum().z();
           help8[i] = pana[5].getEnergy();
           help9[i] = sqr(pana[5].getMomentum().x())+sqr(pana[5].getMomentum().y());
           if(fabs(pv[i].getCharge()) > 0.5) nhmf++;
         }     
       if(nhmf <= 2) continue;
       nevhmf++;  
       for(i=0; i<ngk; i++)
         {
           if(pv[i].getCharge() > 0.5) plot[1].Fill(help7[i],20.);
           else if(pv[i].getCharge() < -0.5) plot[2].Fill(help7[i],20.);
           else if(pv[i].getName() == "PionZero")
	     {
               plot[3].Fill(help7[i],20.);
               G4double wgthmf = 2.*help8[i]/pana[3].getMass();
               plot[4].Fill(help7[i],wgthmf*20.);
               if(help7[i] > 0.025) plot[5].Fill(help9[i],20.);
             }
         }                           

       iteration++; 

     }  while (iteration < NumberEvents) ;

    cout << " number of useful events " << nevhmf << G4endl; 
   
    G4double wgthmf = 20.94/nevhmf;
    if(incidentParticle.getName() == "KaonPlus") wgthmf = 17.72/nevhmf;
    for(i=1; i<=5; i++) 
     {
       plot[i].Scale(wgthmf, plot[i]);
       plot[i].DumpToFile(i, aFile);
     } 
    delete pana;
    delete plot;
    delete pv; 
    delete help7;
    delete help8;
    delete help9; 
}


