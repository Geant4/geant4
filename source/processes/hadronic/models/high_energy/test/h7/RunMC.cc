// Test- Job simulates the experimental data given in
// G. Donaldsen et al., Phys. Letters 73B, 375 (1978);
// F. Pauss et al., Z. f. Phys. C27, 211 (1985);
// R.P.Apsomon et al., Z. f. Phys. C52, 397 (1991);
// pi0 Production at high Energy 
// The same comparison has been done with Geant3 MC
// H. Fesefeldt, August 1998.

#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEKaonMinusInelastic.hh"
#include "G4HEPlot.hh"

int main()
{
//*******************************************************

  G4int NumberEvents = 100000;
  G4double incidentMomentum = 58.;
 
// PionPlus, PionMinus, KaonPlus, KaonMinus at 110 GeV
// PionMinus, KaonMinus at 58 GeV

  G4String aFile = "PionMinus58.plot";

  G4HEPionMinusInelastic theInteraction;

  G4HEVector incidentParticle;
  incidentParticle.setDefinition("PionMinus");

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

  G4HEPlot* plot = new G4HEPlot[30]; 
  G4HEVector* pana = new G4HEVector[11];
  G4double* help7 = new G4double[maxpart];
  G4double* help8 = new G4double[maxpart];

  plot[ 1].Init(30, 0., 1.);
  plot[ 2].Init(25, 0., 0.1);
  plot[ 3].Init(25, 0., 0.1);
  plot[ 4].Init(25, 0., 0.1);
  plot[ 5].Init(25, 0., 0.1);
  plot[ 6].Init(25, 0., 0.1);
  plot[ 7].Init(25, 0., 0.1);
  plot[ 8].Init(25, 0., 0.1);
  plot[ 9].Init(25, 0., 0.1);
  plot[10].Init(25, 0., 0.1);
  plot[11].Init(25, 0., 0.1);
  plot[12].Init(10, 0., 0.1);
  plot[13].Init(25, 0., 0.1);
  plot[14].Init(25, 0., 0.1);
  plot[15].Init(25, 0., 0.1);
  plot[16].Init(20, 0., 0.05);
  plot[17].Init(20, 0., 0.05);
  plot[18].Init(20, 0., 0.05);
  plot[19].Init(20, 0., 0.05);
  plot[20].Init(20, 0., 0.05);
  plot[21].Init(20, 0., 0.05);
  plot[22].Init(10, 0., 0.1);
  plot[23].Init(10, 0., 0.1);
  plot[24].Init(10, 0., 0.1);
  plot[25].Init(20, 0., 0.05);
  plot[26].Init(20, 0., 0.05);
  plot[27].Init(20, 0., 0.05);
  plot[28].Init(20, 0., 0.05);
  plot[29].Init(20, 0., 0.05);

  G4int nevhmf = 0;

  pana[1] = incidentParticle;
  pana[2].setMass(0.9383);
  pana[2].setCharge(1.);
  pana[2].setKineticEnergyAndUpdate(0.);
  pana[3].Add(pana[1], pana[2]);
  pana[4].Lor(pana[1], pana[3]);
  pana[5].Lor(pana[2], pana[3]);

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
       if(ngk <=2) continue;
       nevhmf++;

       pana[10].setZero();
       for(i=0; i<ngk; i++) pana[10].Add(pana[10], pv[i]);
       pana[6].Sub(pana[1],pana[10]);
       G4double dplab = pana[6].Length();
       dplab /= incidentMomentum;
       dplab = fabs(dplab);
       plot[29].Fill(dplab, 1.);

       for(i=0; i<ngk; i++)
         {
           help7[i] = pv[i].getMomentum().z()/pana[1].getMomentum().z();
           help8[i] = sqrt(sqr(pv[i].getMomentum().x()) + sqr(pv[i].getMomentum().y()));
           if(pv[i].getName() == "PionPlus")
             {
               G4double xhmf = help7[i];
               G4int ixhmf = (G4int)(xhmf/0.1) + 1;
               ixhmf = G4std::min(10,ixhmf);
               G4double pthmf = help8[i];
               if(pthmf < 0.05) pthmf = 0.025;
               G4double wgthmf = pv[i].getEnergy()/(pana[1].getMomentum().z()*pthmf);
               plot[1+ixhmf].Fill(help8[i],wgthmf);
               plot[12].Fill(xhmf,1.);
               if(ixhmf >=2 && ixhmf <= 7)
                 {
                   G4int ix1hmf = (ixhmf-2)/2 + 1;
                   plot[12+ix1hmf].Fill(help8[i],1.);
                 }            
             }
           G4int ioff = 0;
           if(pv[i].getName() == "PionZero") ioff = 1; 
           if(pv[i].getName() == "PionPlus") ioff = 2; 
           if(pv[i].getName() == "PionMinus") ioff = 3; 
           if(pv[i].getName() == "KaonPlus") ioff = 4; 
           if(pv[i].getName() == "KaonMinus") ioff = 5; 
           if(pv[i].getType() == "Baryon") ioff = 6;
           if(ioff > 0)
             {
               G4double xhmf = help7[i];
               G4double wgthmf = xhmf;
               plot[15+ioff].Fill(xhmf,wgthmf);
             }
           if(ioff > 0 && ioff <= 3) plot[21+ioff].Fill(help7[i],1.); 
         }               

       iteration++; 

     }  while (iteration < NumberEvents) ;

    cout << " number of useful events " << nevhmf << G4endl; 
   
    i = 0; 
    G4double sig = 0.;
    if(incidentParticle.getName() == "PionPlus") sig = 20.0;
    else if(incidentParticle.getName() == "PionMinus") sig = 20.0; 
    else if(incidentParticle.getName() == "KaonPlus") sig = 17.0; 
    else if(incidentParticle.getName() == "KaonMinus") sig = 17.0; 


    G4double wgthmf = 1./nevhmf;
    plot[ 1].Scale(wgthmf, plot[ 1]);
    wgthmf = sig/nevhmf;    
    plot[12].Scale(wgthmf, plot[12]);    
    plot[13].Scale(wgthmf, plot[13]);    
    plot[14].Scale(wgthmf, plot[14]);    
    plot[15].Scale(wgthmf, plot[15]);    
    plot[22].Scale(wgthmf, plot[22]);    
    plot[23].Scale(wgthmf, plot[23]);    
    plot[24].Scale(wgthmf, plot[24]);    
    plot[26].Scale(wgthmf, plot[26]);    
    plot[27].Scale(wgthmf, plot[27]);    
    plot[28].Scale(wgthmf, plot[28]);
    wgthmf = sig/(nevhmf*0.1*0.1*2.*3.1415927);
    for(i=2; i<12; i++)
      {
        if(i > 3) wgthmf /= 10.;
        plot[i].Scale(wgthmf,plot[i]);
      }
    wgthmf = sig/(nevhmf*0.05*3.1415927);
    for(i=16;i<22;i++) plot[i].Scale(wgthmf,plot[i]);
    plot[25].Scale(wgthmf,plot[25]);
    for(i=1; i<25; i++) plot[i].DumpToFile(i, aFile);
    plot[29].DumpToFile(29, aFile);

    delete pana;
    delete plot;
    delete pv; 
    delete help7;
    delete help8;
}


