// Test- Job simulates the experimental data given in
// I.V.Ajinenko et al., Z.Phys.C42, 377 (1989)
// I.V.Ajinenko et al., Z.Phys.C46, 569 (1990)
// M.Adamus et al., Z.Phys.C32, 475 (1986)
// pi+, K+, p  on H, Al and Au at 250 GeV/c 
// The same comparison has been done with Geant3 MC
// H. Fesefeldt, August 1998.
// This program runs Geant4 for three particles and
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

  G4double A, Z;
  //  A = 196.97;
  //  Z = 79.;
  //  A = 26.98;
  //  Z = 13.;
  A = 1.01;
  Z = 1.;  
 
  G4String aFile = "ProtonH.plot";
  G4String aHFile = "";

  G4HEProtonInelastic theInteraction;

  G4HEVector incidentParticle;
  incidentParticle.setDefinition("Proton");

//******************************************************* 
        
  G4int maxpart = 1000;
  G4HEVector* pv = new G4HEVector[maxpart];
  incidentParticle.setMomentumAndUpdate(0., 0., 250.); 
  theInteraction.SetMaxNumberOfSecondaries(maxpart);
  theInteraction.SetVerboseLevel(2);
  G4int i = 0;

  G4HEPlot* plot = new G4HEPlot[25]; 
  G4HEVector* pana = new G4HEVector[11];

  plot[ 1].Init(30, 0., 1.);
  plot[ 2].Init(50, 1., 2.);
  plot[ 3].Init(25, 0., 1.);
  plot[ 4].Init(30, 0., 0.05);
  plot[ 5].Init(30, 0., 1.);
  plot[ 6].Init(30, 0., 1.);
  plot[ 7].Init(30, 0., 1.);
  plot[ 8].Init(30, 0., 1.);
  plot[ 9].Init(30, 0., 1.);
  plot[10].Init(30, 0., 1.);
  plot[11].Init(30, 0., 1.);
  plot[12].Init(30, 0., 1.);
  plot[13].Init(30, 0., 1.);
  plot[14].Init(30, 0., 1.);
  plot[15].Init(30, 0., 1.);
  plot[16].Init(60, 0., 1.);
  plot[17].Init(60, 0., 1.);
  plot[18].Init(30, 0., 1.);
  plot[19].Init(30, 0., 1.);
  plot[21].Init(50, 0.6, 0.01);
  plot[22].Init(50, 0.6, 0.01);
  plot[23].Init(50, 0., 4.);

  G4int nevhmf = 0;

  pana[1] = incidentParticle;
  pana[2].setMass(0.9383);
  pana[2].setKineticEnergyAndUpdate(0.);
  pana[3].Add(pana[1], pana[2]);
  pana[4].Lor(pana[1], pana[3]);

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
       plot[1].Fill((G4double)ngk, 1.); 

       G4int nhmf = 0;
       G4int nchhmf = 0;

       G4double* rapi = new G4double[ngk];
       G4double* xfey = new G4double[ngk];

       for(i=0; i<ngk; i++)
         {
           nhmf++;
           if(fabs(pv[i].getCharge()) > 0.5)
             {
               nchhmf++;
               G4double phf = pv[i].Length(); 
               if(phf > 1.2)
                 {
                   pana[10] = pv[i];
                   if(pv[i].getCharge() < 0.) {pv[i].setDefinition("PionMinus");}
                   else {pv[i].setDefinition("PionPlus");}
                   pv[i].setMomentum(pana[10].getMomentum());
                   pv[i].setEnergy(sqrt(phf*phf + pv[i].getMass()*pv[i].getMass()));

                 }
             }                         
           pana[5].Lor(pv[i], pana[3]);
           G4double ppp = pana[5].Length();
           G4double xup = pana[5].getEnergy() + pana[5].getMomentum().z();
           G4double xlo = pana[5].getEnergy() - pana[5].getMomentum().z();
           xup = xup/xlo;
           rapi[i] = log(xup)/2.;
           xfey[i] = pana[5].getMomentum().z()/pana[4].getMomentum().z();
         }               

       if(nhmf <= 2) continue;
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
       plot[21].Fill(etothf, 1.);
       plot[22].Fill(ekinhf, 1.);
       
       G4int ngrey = 0;
       G4int nneg  = 0; 
       for(i=0; i<ngk; i++)
         {
           pv[i].setSide(2);
           if(fabs(pv[i].getCharge()) < 0.5) continue;
           if(pv[i].getType() == "Nucleus") continue;
           G4double phmf = pv[i].Length();
           if(pv[i].getName() == "Proton")
             {
               if(phmf < 0.15) continue;

               phmf *= (1.+0.2*theInteraction.normal());
               if(G4UniformRand() > 5.*phmf) continue;
               G4double ran = 0.9 + 0.5*G4UniformRand();
               if(phmf > ran) 
	         {
                   pv[i].setSide(0);
                   if(pv[i].getCharge() < -0.5) nneg++;
                 }
               else
                 {
                   plot[4].Fill(phmf, 1.);
                   if(phmf < 0.19) 
                       continue;
                   else if(phmf > 0.9)
                     {
                       pv[i].setSide(0);
                       if(pv[i].getCharge() < -0.5) nneg++;
		     }
                   else
                     { 
                       pana[10] = pv[i];
                       pv[i].setDefinition("Proton");
                       pv[i].setMomentum(pana[10].getMomentum());
                       pv[i].SmulAndUpdate(pv[i],phmf/pana[10].Length());
                       pv[i].setSide(1);
                       ngrey++;
                       continue;
                     } 
                 }
             }
           else
             {
               pv[i].setSide(0);
               if(pv[i].getCharge() < -0.5) nneg++;
             }
           G4double phmf1 = G4std::max(0., phmf*(1.+0.2*theInteraction.normal()));
           pv[i].SmulAndUpdate(pv[i], phmf1/phmf);
           pana[5].Lor(pv[i], pana[3]);
           G4double xup = pana[5].getEnergy() + pana[5].getMomentum().z();
           G4double xlo = pana[5].getEnergy() - pana[5].getMomentum().z();
           rapi[i] = log(xup/xlo)/2.;
           xfey[i] = pana[5].getMomentum().z()/pana[4].getMomentum().z();
         }

       G4double rneg = (double)nneg;
       for(i=0; i<ngk; i++)
	 {
           if(pv[i].getSide() == 0 && xfey[i] > 0.)
             {
               if(G4UniformRand() < 0.0085*rneg) pv[i].setSide(2);
             }
         }
       if(A > 1.5)
         {
           for(i=0; i<ngk; i++)
	     {
               if(pv[i].getSide() == 0 && xfey[i] > 0.)
                 {
                   if(G4UniformRand() < 0.2) pv[i].setSide(2);
                 }
             }
         } 
       G4int nall = 0;
       G4int npos = 0;
             nneg = 0;
       G4int nnegl = 0;
       G4int nnegr = 0;
       G4int nallfw = 0;
       G4int nallbw = 0;
       G4int nnegfw = 0;
       G4int nnegbw = 0;

       for (i=0; i<ngk; i++)
	 {
           G4double pthmf = sqrt(sqr(pv[i].getMomentum().x())+sqr(pv[i].getMomentum().y()));
           if(pv[i].getSide() == 0)
             {
               nall++;
               if(rapi[i] > 0.) nallfw++;
               else             nallbw++;
               if(pv[i].getCharge() < -0.5)
                 {
                   nneg++;
                   if(rapi[i] > -2.5 && rapi[i] < -1.0) nnegl++;
                   if(rapi[i] > -0.5 && rapi[i] < 1.0) nnegr++;
                   if(rapi[i] > 0.) nnegfw++;
                   else             nnegbw++;
                 }
               if(pv[i].getCharge() > 0.5) npos++;
             }
	 }
       rneg = (G4double)nneg;
       G4double rneg2 = rneg*rneg;
       G4double rgrey = (G4double)ngrey + 0.001;
       G4double rall = (G4double)nall;
       G4double rall2 = rall*rall;
       plot[1].Fill(rgrey,1.);
       plot[2].Fill(rall+0.001,1.);
       plot[3].Fill(rneg+0.001,1.);
       plot[5].Fill(rgrey,rneg);
       plot[6].Fill(rgrey,rneg2);
       plot[8].Fill(rgrey,rall);
       plot[9].Fill(rgrey,rall2);
       G4double rnegl = (G4double)nnegl;
       G4double rnegr = (G4double)nnegr;
       G4double rneglr = rnegl*rnegr;
       plot[11].Fill(rgrey,rnegl);
       plot[12].Fill(rgrey,rnegr);
       plot[13].Fill(rgrey,rneglr);
       nallbw += ngrey;
       plot[16].Fill((G4double)nallfw + 0.001,1.);
       plot[17].Fill((G4double)nallbw + 0.001,1.);
       plot[18].Fill((G4double)nnegfw + 0.001,1.);
       plot[19].Fill((G4double)nnegbw + 0.001,1.); 

       delete rapi;
       delete xfey;

       iteration++; 

     }  while (iteration < NumberEvents) ;

    cout << " number of useful events " << nevhmf << G4endl; 
   
    G4double wgthmf = 1./nevhmf;
    for(i = 1; i<=19; i++) plot[i].Scale(wgthmf, plot[i]);
    G4double wgthm1 = 1000./1.4;
    plot[4].Scale(wgthm1, plot[19]);    


    plot[5].Divide(1., 1., plot[5], plot[1]);   
    plot[6].Divide(1., 1., plot[6], plot[1]);   
    plot[7].Multiply(1., 1., plot[5], plot[5]);   
    plot[6].Add(1., -1., plot[6], plot[7]);   
    plot[6].Sqrt(1.,plot[6]);
    plot[7].Divide(1., 1., plot[6], plot[5]);
    plot[8].Divide(1., 1., plot[8], plot[1]);
    plot[9].Divide(1., 1., plot[9], plot[1]);
    plot[10].Multiply(1., 1., plot[8], plot[8]);
    plot[9].Add(1., -1., plot[9], plot[10]);
    plot[9].Sqrt(1., plot[9]);
    plot[10].Divide(1., 1., plot[9], plot[8]);
    plot[11].Divide(1., 1., plot[11], plot[1]);
    plot[12].Divide(1., 1., plot[12], plot[1]);
    plot[13].Divide(1., 1., plot[13], plot[1]);
    plot[14].Multiply(1., 1., plot[11], plot[12]);
    plot[15].Add(1., -1., plot[13], plot[14]);

        
    for (i=1; i<20; i++)
      {
        plot[i].DumpToFile( i, aFile);
      }
    for (i=21; i<24; i++)
      {
        plot[i].DumpToFile( i, aFile);
      }

    delete pana;
    delete plot;
    delete pv; 
}



