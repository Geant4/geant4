// Test- Job simulates the experimental data given in
// H. Fesefeldt, W.Ochs and L.Stodolsky, Phys.Lett. 74B (1978) 389.
// B.Andersson and G.Gustafson, Phys.Lett. 84B (1979) 483.
// W.Ochs and T.Shimada, Z.Phys.C 4 (1980) 141.
// I.V.Ajinenko et al., Z.Phys.C 16 (1983) 291.
// I.V.Ajinenko et al., Z.Phys.C 43 (1989) 37.
// The same comparison has been done with Geant3 MC
// H. Fesefeldt, July 1998.
// This program runs Geant4 for two particles and
// three materials and produces for each a file
// "particle""Material".plot.

#include "G4HEPionPlusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEProtonInelastic.hh"
#include "G4HEPlot.hh"

int main()
{
//*******************************************************

  G4int NumberEvents = 100000;

  G4double pIncident = 250.;  
  //  G4double pIncident = 24.;
  //  G4double pIncident = 12.;
 
  G4String aFile = "Kaon250.plot";

  G4HEKaonPlusInelastic theInteraction;

  G4HEVector incidentParticle;
  incidentParticle.setDefinition("KaonPlus");

//******************************************************* 
        
  G4double A = 1.01;
  G4double Z = 1.;

  G4int maxpart = 500;
  G4HEVector* pv = new G4HEVector[maxpart];
  incidentParticle.setMomentumAndUpdate(0., 0., pIncident); 
  theInteraction.SetMaxNumberOfSecondaries(maxpart);
  theInteraction.SetVerboseLevel(2);
  G4int i = 0;

  G4HEPlot* plot = new G4HEPlot[34];
  G4HEPlot* plotE = new G4HEPlot[4]; 
  G4HEVector* pana = new G4HEVector[11];
  G4double* help7 = new G4double[maxpart];
  G4double* help8 = new G4double[maxpart];
  G4double* help9 = new G4double[maxpart]; 

  plot[ 1].Init(20, 0., 0.1);  
  plot[ 2].Init(20, 0., 0.1);
  plot[ 3].Init(20, 0., 0.1);
  plot[ 4].Init(20, 0., 0.1);
  plot[ 5].Init(20, 0., 0.1);
  plot[ 6].Init(20, 0., 0.1);
  plot[ 7].Init(20, 0., 0.1);
  plot[ 8].Init(20, 0., 0.1);
  plot[ 9].Init(20, 0., 0.1);
  plot[10].Init(20, 0., 0.1);
  plot[11].Init(20, 0., 0.1);
  plot[12].Init(20, 0., 0.1);
  plot[13].Init(20, 0., 0.1);
  plot[14].Init(20, 0., 0.1);
  plot[15].Init(20, 0., 0.1);
  plot[16].Init(20, 0., 0.1);
  plot[17].Init(20, 0., 0.1);
  plot[18].Init(20, 0., 0.1);
  plot[19].Init(20, 0., 0.1);
  plot[20].Init(20, 0., 0.1);
  plot[21].Init(20, 0., 0.1);
  plot[22].Init(40, 0.,0.25);
  plot[23].Init(40, 0.,0.25);
  plot[24].Init(40, 0.,0.25);
  plot[25].Init(40, 0.,0.25);
  plot[26].Init(40, 0., 0.05);
  plot[27].Init(40, 0., 0.05);
  plot[28].Init(40, 0., 0.05);
  plot[29].Init(40, 0., 0.05);
  plot[30].Init(40, 0., 0.05);
  plot[31].Init(40, 0., 0.05);
  plot[32].Init(40, 0., 0.05);
  plot[33].Init(40, 0., 0.05);

  plotE[1].Init(50, 0.6, 0.01);
  plotE[2].Init(50, 0.6, 0.01);
  plotE[3].Init(50, 0., 4.);


  G4int nevhmf = 0;

  pana[1] = incidentParticle;
  pana[2].setDefinition("Proton");
  pana[2].setMomentumAndUpdate(0.,0.,0.);
  pana[3].Add(pana[1], pana[2]);
  pana[4].Lor(pana[1], pana[3]);
  G4double xup = pana[1].getEnergy() + pana[1].getMomentum().z();
  G4double xlo = pana[1].getEnergy() - pana[1].getMomentum().z();
  G4double drap = log(xup/xlo)/2.;
  cout << " Beam Rapidity " << drap << G4endl;

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
       if(iteration == 1) theInteraction.SetVerboseLevel(2);
       else               theInteraction.SetVerboseLevel(0);
       theInteraction.ApplyYourself( pv, incidentParticle, A, Z);
       G4int vecLength = theInteraction.GetNumberOfSecondaries();
       G4int ngk = vecLength;
       G4int nchhmf = 0;
       G4int nhmf = ngk;
       G4double xfast = -1.; 
       G4int ifast = -1; 
       G4double xhmf = -10.;

       if(iteration < 10) cout << " iteration " << iteration << G4endl;
       for(i=0; i<ngk; i++)
         {
           if(fabs(pv[i].getCharge()) > 0.5)
             {
               nchhmf++;
               G4double phf = pv[i].Length();
               if(phf > 1.2 || pv[i].getName() != "Proton")
                 {
                   G4HEVector pTemp = pv[i];
                   if(pTemp.getCharge() > 0.5) pv[i].setDefinition("PionPlus");
                   else pv[i].setDefinition("PionMinus");
                   pv[i].setMomentum(pTemp.getMomentum());
                   pv[i].setEnergyAndUpdate(sqrt(phf*phf + pv[i].getMass()*pv[i].getMass()));
                 }
             }                         
           pana[5].Lor(pv[i], pana[3]);
           G4double cost = pana[4].CosAng(pana[5]);
           G4double sint = sqrt(1.-cost*cost);
           if(sint < 1.e-20) sint = fabs(theInteraction.normal())*0.01;
           help7[i] = (cost/sint)/pana[4].getEnergy() + 0.075*theInteraction.normal();
           help8[i] = pana[5].getEnergy()/pana[4].getEnergy();
           if(fabs(pv[i].getCharge()) > 0.5)
             {
               xhmf = pana[5].getMomentum().z()/pana[4].getMomentum().z();
               if(xhmf > xfast)
                 {
                   xfast = xhmf;
                   ifast = i;
                 }
             }
           help9[i] = pana[4].Ang(pana[5]) * 180./M_PI; 
           if(iteration < 10)
             {
               pv[i].Print(i);
               pana[5].Print(5);
               cout << " help " << help7[i] << " " << help8[i] << " " << help9[i] << G4endl;
             } 
         }
        
       plotE[3].Fill((G4double)nhmf, 1.);
       if(nhmf <= 2) continue;
       if(nchhmf <= 1) continue;
       nevhmf++;
       
       G4double etot = 0.;
       G4double ekin = 0.;
       for(i=0; i<ngk; i++)
         {
           etot += pv[i].getEnergy();
           ekin += pv[i].getKineticEnergy();
         }
       etot /= pana[3].getMomentum().z();
       ekin /= pana[3].getMomentum().z();
       plotE[1].Fill(etot,1.);
       plotE[2].Fill(ekin,1.);
      
       G4int ichhmf =G4std::min(5,(nchhmf-2)/2+1);
       for(i=0; i<ngk; i++)
	 {
           G4double rl, dl;
           G4double rlhmf = fabs(help7[i]);
           if(rlhmf <= 1.)
             {
               rl = rlhmf;
               dl = 0.1;
             }
           else
             {
               rl = 2.-1./rlhmf;
               dl = 0.1*rlhmf*rlhmf;
             }
           G4double wgte = help8[i]/dl;
           G4double wgtq = pv[i].getCharge()/dl;
           if(help7[i] > 0.)
             {
               plot[1].Fill(rl,wgte);
               plot[3].Fill(rl,wgtq);
               plot[6+ichhmf].Fill(rl,wgte);
               plot[11+ichhmf].Fill(rl,wgtq);
             }
           else
             {
               plot[2].Fill(rl,wgte);
               plot[4].Fill(rl,wgtq);
             }
         }
       pana[6].setZero();
       pana[7].setZero();
       pana[8].setZero();
       pana[9].setZero();
       G4double* echmf = new G4double[5];
       G4double* cchmf = new G4double[5];
       for(i=0; i<5; i++)
         {
           echmf[i] = 0.;
           cchmf[i] = 0.;
         }
       G4double ethmf = 0.;
       for(i=0; i<ngk; i++)
	 {
           G4double tetahf = help9[i];
           xhmf = help8[i];
           ethmf+=help8[i];
           if(tetahf > 90.) tetahf=180.-tetahf;
           if(tetahf > 30.)
             {
               pana[6].Add(pana[6],pv[i]);
               echmf[1]+=xhmf;
               cchmf[1]+=pv[i].getCharge();
             }
           if(tetahf > 45.)
             {
               pana[7].Add(pana[7],pv[i]);
               echmf[2]+=xhmf;
               cchmf[2]+=pv[i].getCharge();
             } 
           if(tetahf > 60.)
             {
               pana[8].Add(pana[8],pv[i]);
               echmf[3]+=xhmf;
               cchmf[3]+=pv[i].getCharge();
             } 
           if(tetahf > 75.)
             {
               pana[9].Add(pana[9],pv[i]);
               echmf[4]+=xhmf;
               cchmf[4]+=pv[i].getCharge();
             } 
         }
       plot[22].Fill(pana[6].getMass(),4.);
       plot[23].Fill(pana[7].getMass(),4.);
       plot[24].Fill(pana[8].getMass(),4.);
       plot[25].Fill(pana[9].getMass(),4.);
       if(iteration < 10)
         {
           cout << " Masses " << pana[6].getMass() << " " << pana[7].getMass() << " " 
		<< pana[8].getMass() << " " << pana[9].getMass() << G4endl;
           cout << echmf[1] << " " << echmf[2] << " " << echmf[3] << " " << echmf[4] << G4endl; 
           cout << cchmf[1] << " " << cchmf[2] << " " << cchmf[3] << " " << cchmf[4] << G4endl; 
         }
       for(i=1; i<5; i++)
         {
           echmf[i]/=ethmf;
           plot[25+i].Fill(echmf[i],1.);
           plot[29+i].Fill(echmf[i],cchmf[i]);
         }
       delete [] echmf;
       delete [] cchmf; 

       //       cout << " iteration " << iteration << " " << ngk << G4endl;       

       iteration++; 

     }  while (iteration < NumberEvents) ;

    cout << " number of useful events " << nevhmf << G4endl; 
   
    G4double wgthmf = 1./nevhmf;

    for(i = 1; i<=4; i++) plot[i].Scale(wgthmf, plot[i]);

    plot[5].Divide(1., 1., plot[3], plot[1]);   
    plot[6].Divide(1., 1., plot[4], plot[2]);

    for(i=1; i<=5; i++) plot[16+i].Divide(1., 1., plot[11+i], plot[6+i]); 

    wgthmf*=30.;
    for(i=22; i<=25; i++) plot[i].Scale(wgthmf,plot[i]); 

    plot[30].Divide(1., 1., plot[30], plot[26]);
    plot[31].Divide(1., 1., plot[31], plot[27]);
    plot[32].Divide(1., 1., plot[32], plot[28]);
    plot[33].Divide(1., 1., plot[33], plot[29]);

    for(i=1; i<=33; i++) plot[i].DumpToFile(i, aFile);
    for(i=1; i<=3; i++) plotE[i].DumpToFile(50+i, aFile);         

    delete pana;
    delete plot;
    delete plotE;
    delete pv; 
    delete help7;
    delete help8;
    delete help9;

}


