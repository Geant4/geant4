// Test- Job simulates the experimental data given in
// N.M.Agababyan et al., Z.Phys.C (1990)
// pi+, K+ on H, Al and Au at 250 GeV/c 
// The same comparison has been done with Geant3 MC
// H. Fesefeldt, July 1998.
// This program runs Geant4 for two particles and
// three materials and produces for each a file
// "particle""Material".plot.

#include "G4HEPionPlusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEPlot.hh"

int main()
{
//*******************************************************

  G4int NumberEvents = 100000;

  G4double A, Z;
  A = 196.97;
  Z = 79.;
  //  A = 26.98;
  //  Z = 13.;
  //  A = 1.01;
  //  Z = 1.;  
 
  G4String aFile = "PionAu.plot";
  G4String aHFile = "PionH.plot";

  G4HEPionPlusInelastic theInteraction;

  G4HEVector incidentParticle;
  incidentParticle.setDefinition("PionPlus");

//******************************************************* 
        
  G4int maxpart = 1000;
  G4HEVector* pv = new G4HEVector[maxpart];
  incidentParticle.setMomentumAndUpdate(0., 0., 250.); 
  theInteraction.SetMaxNumberOfSecondaries(maxpart);
  theInteraction.SetVerboseLevel(2);
  G4int i = 0;
  for(i=0; i<100; i++) theInteraction.testArray[i] = 0.;

  G4HEPlot* plot = new G4HEPlot[56]; 
  G4HEVector* pana = new G4HEVector[11];

  plot[ 1].Init(30, 0., 1.);
  plot[ 2].Init(50, 0., 2.);
  plot[ 3].Init(25, 0., 2.);
  plot[ 4].Init(30, 0., 0.05);
  plot[ 5].Init(22, -5.5, 0.5);
  plot[35].Init(22, -5.5, 0.5);
  plot[23].Init(22, -5.5, 0.5);
  plot[ 6].Init(22, -5.5, 0.5);
  plot[36].Init(22, -5.5, 0.5);
  plot[24].Init(22, -5.5, 0.5);
  plot[ 7].Init(22, -5.5, 0.5);
  plot[37].Init(22, -5.5, 0.5);
  plot[25].Init(22, -5.5, 0.5);
  plot[ 8].Init(100, 0., 0.05);
  plot[38].Init(100, 0., 0.05);
  plot[26].Init(100, 0., 0.05);
  plot[ 9].Init(100, 0., 0.05);
  plot[39].Init(100, 0., 0.05);
  plot[27].Init(100, 0., 0.05);
  plot[10].Init(100, 0., 0.05);
  plot[40].Init(100, 0., 0.05);
  plot[28].Init(100, 0., 0.05);
  plot[11].Init(32, -0.8, 0.05);
  plot[12].Init(32, -0.8, 0.05);
  plot[13].Init(32, -0.8, 0.05);
  plot[14].Init(20, 0., 1.);
  plot[15].Init(20, 0., 1.);
  plot[16].Init(20, 0., 1.);
  plot[17].Init(40, 0., 0.05);
  plot[47].Init(40, 0., 0.05);
  plot[29].Init(40, 0., 0.05);
  plot[18].Init(40, 0., 0.05);
  plot[48].Init(40, 0., 0.05);
  plot[30].Init(40, 0., 0.05);
  plot[31].Init( 5, 0., 1.);
  plot[32].Init( 5, 0., 1.);
  plot[33].Init( 5, 0., 1.);
  plot[19].Init(10, -5., 0.5);
  plot[20].Init(10, -5., 0.5);
  plot[21].Init(10, -5., 0.5);
  plot[22].Init(10, -5., 0.5);
  plot[51].Init(50, 0.9, 0.01);
  plot[52].Init(50, 0.6, 0.01);
  plot[53].Init(50, 0., 0.1);
  plot[54].Init(50, 0., 1.);
  plot[55].Init(50, 0., 1.);
  plot[56].Init(50, 0., 0.01);

  G4int nevhmf = 0;
  G4int nsav = 15;
  G4int nposal = 0;
  G4int npos = 0;
  G4int nneg = 0;
  G4int nneg1 = 0;
  G4int nneg2 = 0;
  G4double* aver = new G4double[51]; 
  for (i=0; i<51; i++) aver[i] = 0.;

  pana[1] = incidentParticle;
  pana[2].setMass(0.9383);
  pana[2].setKineticEnergyAndUpdate(0.);
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
         }
       if(iteration == 5) theInteraction.SetVerboseLevel(0);
       theInteraction.ApplyYourself( pv, incidentParticle, A, Z);
       G4int vecLength = theInteraction.GetNumberOfSecondaries();
       G4int ngk = vecLength;
       G4int nchhmf = 0;
       G4double ymax = -100.;
       G4double ymax2 = -100.;
       G4double xmax = -100.;
       G4int iymax = 0;
       G4int iymax2 = 0;
       G4double tmin = 100.;
       G4double* rapi = new G4double[ngk];
       G4double* xfey = new G4double[ngk];
       pana[10].setZero();
       for(i=0; i<ngk; i++)
         {
           if(fabs(pv[i].getCharge()) > 0.5)
             {
               nchhmf++;
               G4double tetahf = pv[i].Ang(pana[1]);
               G4double phf = pv[i].Length();
               G4double thmf = pana[1].getMomentum().z()*phf*tetahf*tetahf;
               if(thmf < tmin) tmin=thmf;
               if(phf > 1.2)
                 {
                   pv[i].setMass(theInteraction.PionPlus.getMass());
                   pv[i].setEnergy(sqrt(phf*phf + pv[i].getMass()*pv[i].getMass()));
                   if(pv[i].getCharge() < 0.) {pv[i].setCode(theInteraction.PionMinus.getCode());}
                   else {pv[i].setCode(theInteraction.PionPlus.getCode());}
                 }
             }                         
           pana[5].Lor(pv[i], pana[3]);
           G4double ppp = pana[5].Length();
           xup = pana[5].getEnergy() + pana[5].getMomentum().z();
           xlo = pana[5].getEnergy() - pana[5].getMomentum().z();
           xup = xup/xlo;
           rapi[i] = -10.;
           if(xup > 1.e-6) rapi[i] = log(xup)/2.;
           if(fabs(pv[i].getCharge()) > 0.5)
             {
               if(rapi[i] > ymax)
                 {
                   ymax2 = ymax;
                   iymax2 = iymax;
                   ymax = rapi[i];
                   iymax = i;
                 }
             }
           xfey[i] = pana[5].getMomentum().z()/pana[4].getMomentum().z();
           if(fabs(pv[i].getCharge()) > 0.5)
             {
               if(xfey[i] > xmax) xmax = xfey[i];
               pana[10].Add(pana[10], pv[i]);
             }           
         }               

       if(nchhmf <= 2) continue;
       
       plot[56].Fill(tmin, 1.);
       pana[10].setMomentum(pana[1].getMomentum().z() - pana[10].getMomentum().z());
       G4double pthmf = sqrt(sqr(pana[10].getMomentum().x())+sqr(pana[10].getMomentum().y()));
       plot[53].Fill(pthmf, 1.);
       plot[54].Fill(pana[10].getMomentum().z(), 1.);
       if(pthmf < 0.2 && pana[10].getMomentum().z() < 9.)
         {
           plot[55].Fill((double)nchhmf - 0.5, 1.);
         }               
       nevhmf++;
       G4double etothf = 0.;
       G4double ekinhf = 0.;
       for(i=0; i<ngk; i++)
         {
           etothf += pv[i].getEnergy();
           ekinhf += pv[i].getKineticEnergy();
         }
       etothf /= pana[1].Length();
       ekinhf /= pana[1].Length();
       plot[51].Fill(etothf, 1.);
       plot[52].Fill(ekinhf, 1.);
       
       G4int nnehmf = 0;
       for(i=0; i<ngk; i++)
         {
           pv[i].setSide(2);
           if(fabs(pv[i].getCharge()) < 0.5) continue;
           if(pv[i].getType() == "Nucleus") continue;
           G4double phmf = pv[i].Length();
           if(pv[i].getName() == "Proton")
             {
               if(phmf < 0.15 || G4UniformRand() > 3.5*phmf)
                 {
                   pv[i].setSide(2);
                   continue;
                 }
               if(phmf > 0.7+G4UniformRand()*0.3)
                 {
                   pv[i].setSide(0);
                   if(pv[i].getCharge() < -0.5) nnehmf++;
                 }
               else
                 {
                   pv[i].setSide(1);
                   pv[i].setMass(theInteraction.Proton.getMass());
                   pv[i].setEnergy(phmf*phmf+pv[i].getMass()*pv[i].getMass());
                   pv[i].setCode(theInteraction.Proton.getCode());
                 }
             }
           else
             {
               pv[i].setSide(0);
               if(pv[i].getCharge() < -0.5) nnehmf++;
             }
           G4double phmf1 = G4std::max(0., phmf*(1.+0.2*theInteraction.normal()));
           pv[i].SmulAndUpdate(pv[i], phmf1/phmf);
           pana[5].Lor(pv[i], pana[3]);
           xup = pana[5].getEnergy() + pana[5].getMomentum().z();
           xlo = pana[5].getEnergy() - pana[5].getMomentum().z();
           xup /= xlo;
           rapi[i] = -10.;
           if(xup > 1.e-6) rapi[i] = log(xup)/2.;
           xfey[i] = pana[5].getMomentum().z()/pana[4].getMomentum().z();
           if(pv[i].getName() == "Proton") plot[4].Fill(phmf, 1.);
         }
       G4double rneg = (double)nnehmf;
       for(i=0; i<ngk; i++)
	 {
           if(pv[i].getSide() == 0)
             {
               G4double ran = G4UniformRand();
               if((G4int)A == 1)
                 {
                   if(ran < 0.0085*rneg) pv[i].setSide(2);
                 }
               else
                 {
                   if(xfey[i] > 0.)
                     {
                       if(ran < 0.0085*rneg) pv[i].setSide(2);
                     }
                   else
                     {
                       if(ran < 0.2) pv[i].setSide(2);
                     }
                 }
             }
	 }     
       G4int ngrey = 0;
       G4int nshow = 0;
       G4int nshowm = 0;
       drap = ymax - 5.3;
       for (i=0; i<ngk; i++)
	 {
           pthmf = sqrt(sqr(pv[i].getMomentum().x())+sqr(pv[i].getMomentum().y()));
           if(pv[i].getSide() == 0)
             {
               nshow++;
               if(pv[i].getCharge() > 0.5)
                 {
                   nposal++;
                   npos++;
                   plot[5].Fill(rapi[i],2.);
                   plot[6].Fill(rapi[i],2.);
                   plot[8].Fill(pthmf,20.);
                   plot[9].Fill(pthmf,20.);
                 }
               if(pv[i].getCharge() < -0.5)
                 {
                   nshowm++;
                   nneg++;
                   plot[7].Fill(rapi[i],2.);
                   plot[10].Fill(pthmf,20.);
                   plot[11].Fill(xfey[i],pthmf);
                   plot[12].Fill(xfey[i],1.);
                 }
             }
           if(pv[i].getSide() == 1)
             {
               ngrey++;
               nposal++;
               plot[5].Fill(rapi[i],2.);
               plot[8].Fill(pthmf,20.);
             }
         }
       for(i=0; i<ngk; i++)
         {
           if(pv[i].getCharge() < -0.5)
             {
                pthmf = sqrt(sqr(pv[i].getMomentum().x())+sqr(pv[i].getMomentum().y()));
                G4double xhmf = (G4double)nshowm + 0.5;
                plot[14].Fill(xhmf,pthmf);
                plot[15].Fill(xhmf,1.);
                if(rapi[i] > 0.)
                  {
                    plot[17].Fill(pthmf,1.);
                    nneg1++;
                  }
                else
                  {
                    plot[18].Fill(pthmf,1.);
                    nneg2++;
                  }
             }
          }           
       plot[19].Fill(drap,2.);
       plot[20].Fill(drap,(G4double)nshowm);
       plot[21].Fill(drap,1.);

       G4double* yrap = new G4double[6];
       G4int* irap = new G4int[6];
       for (i=0; i<6; i++)
	 {
           yrap[i] = -100;
           irap[i] = 0;
         }
       for(i=0; i<ngk; i++)
         {
           if(fabs(pv[i].getCharge()) > 0.5 && pv[i].getSide() <= 1)
             {
               if(rapi[i] > yrap[1])
                 {
                    yrap[5] = yrap[4];
                    irap[5] = irap[4];
                    yrap[4] = yrap[3];
                    irap[4] = irap[3];
                    yrap[3] = yrap[2];
                    irap[3] = irap[2];
                    yrap[2] = yrap[1];
                    irap[2] = irap[1];
                    yrap[1] = rapi[i];
                    irap[1] = i;
                 }
               else if(rapi[i] > yrap[2])
                 {
                    yrap[5] = yrap[4];
                    irap[5] = irap[4];
                    yrap[4] = yrap[3];
                    irap[4] = irap[3];
                    yrap[3] = yrap[2];
                    irap[3] = irap[2];
                    yrap[2] = rapi[i];
                    irap[2] = i;
                 }
               else if(rapi[i] > yrap[3])
                 {
                    yrap[5] = yrap[4];
                    irap[5] = irap[4];
                    yrap[4] = yrap[3];
                    irap[4] = irap[3];
                    yrap[3] = rapi[i];
                    irap[3] = i;
                 }    
               else if(rapi[i] > yrap[4])
                 {
                    yrap[5] = yrap[4];
                    irap[5] = irap[4];
                    yrap[4] = rapi[i];
                    irap[4] = i;
                 }
               else if(rapi[i] > yrap[5])
                 {
                    yrap[5] = rapi[i];
                    irap[5] = i;
                 }                                
	     }
         }
       for(G4int j=1; j<6; j++)
         {
             if(irap[j] != 0)
	       {
                  i = irap[j];
                  G4double efrac = pv[i].getEnergy()/pana[1].getEnergy();
                  plot[31].Fill((G4double)j -0.5, efrac);
                  plot[32].Fill((G4double)j -0.5, 1.);
               }
         }
       plot[1].Fill((G4double)ngrey + 0.001, 1.);
       plot[2].Fill((G4double)nshow +0.001, 1.);
       plot[3].Fill((G4double)nshowm +0.001, 1.);
       aver[1]+= nshow+ngrey;
       aver[2]+= nshowm;
       aver[3]+= ngrey;
       aver[6]+= sqr(nshow+ngrey);
       aver[7]+= sqr(nshowm);
       aver[8]+= sqr(ngrey);
       ngrey = G4std::min(ngrey,19);
       ngrey = G4std::max(ngrey,0);
       aver[10+ngrey] += (G4double)nshow;
       aver[30+ngrey] += 1.;       

       delete rapi;
       delete xfey;
       delete yrap;
       delete irap;

       iteration++; 

     }  while (iteration < NumberEvents) ;

    cout << " number of useful events " << nevhmf << G4endl; 
   
  i = 0; 
    do
      {
        theInteraction.testArray[i] /= iteration;
        cout << "testArray[" << i << "] = " << theInteraction.testArray[i] << G4endl;
        i++;
      } while (i <8);
    for(i=8;i<13;i++)
      {
        theInteraction.testArray[i] /= theInteraction.testArray[13];
        cout << "testArray[" << i << "] = " << theInteraction.testArray[i] << G4endl;
      } 
    i = 13;
        cout << "testArray[" << i << "] = " << theInteraction.testArray[i] << G4endl; 
    for(i=14;i<19;i++)
      {
        theInteraction.testArray[i] /= theInteraction.testArray[19];
        cout << "testArray[" << i << "] = " << theInteraction.testArray[i] << G4endl;
      } 
    i = 19;
        cout << "testArray[" << i << "] = " << theInteraction.testArray[i] << G4endl;
    
    for(i=1; i<10; i++)
      {
        aver[i] /= nevhmf;
      }
    aver[6] -= aver[1]*aver[1];
    aver[7] -= aver[2]*aver[2];
    aver[8] -= aver[3]*aver[3];
    aver[6] = G4std::max(aver[6], 0.);
    aver[7] = G4std::max(aver[7], 0.);
    aver[8] = G4std::max(aver[8], 0.);
    aver[6] = sqrt(aver[6]);
    aver[7] = sqrt(aver[7]);
    aver[8] = sqrt(aver[8]);
    cout << " mean, dispersion " << aver[1] << " , " << aver[6] << " , " 
         << aver[2] << " , " << aver[7] << " , " << aver[3] << " , " 
         << aver[8] << G4endl;
    for(i=10; i<30; i++)
      {
         if(aver[20+i] > 0.) aver[i] /= aver[20+i];
         G4int ngrey = i-10;
         cout << " ngrey, <nshow> " << ngrey << " , " << aver[i] << " , " 
              << aver[20+i] << G4endl;
      }
    G4double wgthmf = 1./nevhmf;
    for(i = 1; i<8; i++) plot[i].Scale(wgthmf, plot[i]);
    plot[19].Scale(wgthmf, plot[19]);    
    plot[8].Scale(1./nposal, plot[8]);    
    plot[9].Scale(1./npos, plot[9]);    
    plot[10].Scale(1./nneg, plot[10]);    
    plot[17].Scale(1./nneg1, plot[17]);    
    plot[18].Scale(1./nneg2, plot[18]);    
    plot[13].Divide(1., 1., plot[11], plot[12]);   
    plot[16].Divide(1., 1., plot[14], plot[15]);   
    plot[22].Divide(1., 1., plot[20], plot[21]);   
    plot[33].Divide(1., 1., plot[31], plot[32]);   


    if(aHFile != "")
      {
        for(i=5; i<11; i++)      
          {
            cout << " plot " << i << " to plot " << 30+i << G4endl;
            plot[30+i].GetFromFile( i, aHFile);
          }
        plot[47].GetFromFile( 17, aHFile);
        cout << " plot " << 17 << " to plot " << 47 << G4endl;
        plot[48].GetFromFile( 18, aHFile);
        cout << " plot " << 18 << " to plot " << 48 << G4endl;

        plot[23].Divide(1.,1.,plot[5],plot[35]);
        plot[24].Divide(1.,1.,plot[6],plot[36]);
        plot[25].Divide(1.,1.,plot[7],plot[37]);
        plot[26].Divide(1.,1.,plot[8],plot[38]);
        plot[27].Divide(1.,1.,plot[9],plot[39]);
        plot[28].Divide(1.,1.,plot[10],plot[40]);
        plot[29].Divide(1.,1.,plot[17],plot[47]);
        plot[30].Divide(1.,1.,plot[18],plot[48]);
      }

    for (i=1; i<34; i++)
      {
        plot[i].DumpToFile( i, aFile);
      }
    for (i=51; i<57; i++)
      {
        plot[i].DumpToFile( i, aFile);
      }

    delete pana;
    delete aver;
    delete plot;
    delete pv; 
}


