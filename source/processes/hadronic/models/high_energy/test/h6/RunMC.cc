// Test- Job simulates the experimental data given in
// M.Aquilar-Benitez et al., CERN-PPE/91-22.;
// The same comparison has been done with Geant3 MC
// H. Fesefeldt, August 1998.

#include "G4HEProtonInelastic.hh"
#include "G4HEPlot.hh"

int main()
{
//*******************************************************

  G4int NumberEvents = 100000;
  G4double incidentMomentum = 400.;

  G4String aFile = "Proton.plot";

  G4HEProtonInelastic theInteraction;

  G4HEVector incidentParticle;
  incidentParticle.setDefinition("Proton");

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

  G4HEPlot* plot = new G4HEPlot[52]; 
  G4HEVector* pana = new G4HEVector[11];
  G4double* help7 = new G4double[maxpart];
  G4double* help8 = new G4double[maxpart];
  G4double* help9 = new G4double[maxpart];

  plot[ 1].Init(30, 0., 1.);  
  plot[ 2].Init(20, 0., 0.05);
  plot[ 3].Init(20, 0., 0.05);
  plot[ 4].Init(20, 0., 0.05);
  plot[ 5].Init(20, 0., 0.05);
  plot[ 6].Init(20, 0., 0.05);
  plot[ 7].Init(20, 0., 0.05);
  plot[ 8].Init(20, 0., 0.05);
  plot[ 9].Init(20, 0., 0.05);
  plot[10].Init(20, 0., 0.05);
  plot[11].Init(20, 0., 0.05);
  plot[12].Init(20, 0., 0.05);
  plot[13].Init(20, 0., 0.05);
  plot[14].Init(20, 0., 0.05);
  plot[15].Init(20, 0., 0.05);
  plot[16].Init(50, 0., 0.1);
  plot[17].Init(50, 0., 0.1);
  plot[18].Init(50, 0., 0.1);
  plot[19].Init(50, 0., 0.1);
  plot[20].Init(50, 0., 0.1);
  plot[21].Init(50, 0., 0.1);
  plot[22].Init(50, 0., 0.1);
  plot[23].Init(20, 0., 0.05);
  plot[24].Init(20, 0., 0.05);
  plot[25].Init(20, 0., 0.05);
  plot[26].Init(20, 0., 0.05);
  plot[27].Init(20, 0., 0.05);
  plot[28].Init(20, 0., 0.05);
  plot[29].Init(20, 0., 0.05);
  plot[30].Init(20, 0., 0.0025);
  plot[31].Init(20, -10., 1.);
  plot[32].Init(25, 0., 0.2);
  plot[33].Init(25, 0., 0.2);
  plot[34].Init(25, 0., 0.2);
  plot[35].Init(25, 0., 0.2);
  plot[36].Init(25, 0., 0.2);
  plot[37].Init(25, 0., 0.2);
  plot[38].Init(25, 0., 0.2);
  plot[51].Init(50, 0., 0.004);

  G4int nevhmf = 0;

  pana[1] = incidentParticle;
  pana[2].setMass(0.9382);
  pana[2].setCharge(1.);
  pana[2].setKineticEnergyAndUpdate(0.);
  pana[3].Add(pana[1], pana[2]);
  pana[4].Lor(pana[1], pana[3]);
  pana[5].Lor(pana[2], pana[5]);

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
       if(ngk <=2) continue;
       pana[10].setZero();
       for(i=0; i<ngk; i++) pana[10].Add(pana[10], pv[i]);
       pana[5].Sub(pana[1],pana[10]);
       G4double dplab = pana[5].Length();
       dplab /= pana[1].Length();
       dplab = fabs(dplab);
       plot[51].Fill(dplab,1.);

       G4int nhmf = 0;
       G4int nchhmf = 0;
       G4int nchsum = 0;
       G4double xfast = -10.;
       G4int ifast = 0;

       for(i=0; i<ngk; i++)
         {
           pana[6].Lor(pv[i],pana[3]);
           help7[i] = pana[6].getMomentum().z()/pana[4].getMomentum().z();
           help8[i] = pana[6].getEnergy();
           G4double xup = pana[6].getEnergy() + pana[6].getMomentum().z();
           G4double xlo = pana[6].getEnergy() - pana[6].getMomentum().z();
           help9[i] = log(xup/xlo)/2.;
           if(fabs(pv[i].getCharge()) > 0.5) nchhmf++;
           if(pv[i].getCharge() > 0.5) nchsum++;
           if(pv[i].getCharge() < -0.5) nchsum--;
           nhmf++;
           if(help7[i] > xfast)
             {
               xfast = help7[i];
               ifast = i;
             } 
         }     
       if(nhmf <= 2) continue;
       G4double ran = G4UniformRand();
       if(nchhmf == 2 && ran < 0.983) continue;
       if(nchhmf == 4 && ran < 0.563) continue;
       if(nchhmf == 6 && ran < 0.201) continue;
       if(nchhmf == 8 && ran < 0.043) continue;
       G4double wgt = 1.;
       if(nchhmf == 4) wgt = 2.288;
       if(nchhmf == 6) wgt = 1.252;
       if(nchhmf == 8) wgt = 1.045;
       nevhmf++;
       plot[1].Fill((G4double)nchhmf+0.5, wgt);
       plot[31].Fill((G4double)nchhmf+0.5, wgt);
  
       for(i=0; i<ngk; i++)
         {
           G4int ioff = 0;
           if(pv[i].getName() == "PionZero") ioff = 1;
           else if(pv[i].getName() == "PionPlus") ioff = 2;
           else if(pv[i].getName() == "PionMinus") ioff = 3;
           else if(pv[i].getName() == "Proton") ioff = 4;
           else if(pv[i].getName() == "KaonPlus") ioff = 5;
           else if(pv[i].getName() == "KaonMinus") ioff = 6;
           else if(pv[i].getName() == "AntiProton") ioff = 7;
           if(ioff > 0)
             {
               G4double xhmf = help7[i];
               G4double yhmf = help9[i];
               if(help7[i] > 0.)
                 {
                   plot[1+ioff].Fill(xhmf, 20.*wgt);
                   if(ioff == 1) plot[30].Fill(xhmf, 400.*wgt);
                   G4double wgthmf = help8[i]/pana[4].getMomentum().z();
                   wgthmf *= 20.;
                   plot[8+ioff].Fill(xhmf, wgthmf*wgt);
                   G4double pt2hmf = sqr(pv[i].getMomentum().x()) + sqr(pv[i].getMomentum().y());
                   plot[15+ioff].Fill(pt2hmf, 10.*wgt);
                   wgthmf = pt2hmf*20.;
                   plot[22+ioff].Fill(xhmf, wgthmf*wgt);
                 }
               if(yhmf > 0.) plot[31+ioff].Fill( yhmf, 5.*wgt);
             }         
         }                           

       iteration++; 

     }  while (iteration < NumberEvents) ;

    cout << " number of useful events " << nevhmf << endl; 
   
    G4double wgthmf = 1./nevhmf;
    plot[1].Scale(wgthmf, plot[1]);
    wgthmf = 28.6/nevhmf;
    for(i=2; i<= 30; i++) plot[i].Scale(wgthmf, plot[i]);
    for(i=32; i<= 38; i++) plot[i].Scale(wgthmf, plot[i]);
    for(i=1; i<=7; i++) plot[22+i].Divide(1., 1., plot[22+i], plot[1+i]); 
    for(i=1; i<=38; i++) plot[i].DumpToFile(i, aFile);
    plot[51].DumpToFile(51, aFile);
    delete pana;
    delete plot;
    delete pv; 
    delete help7;
    delete help8;
    delete help9; 
}


