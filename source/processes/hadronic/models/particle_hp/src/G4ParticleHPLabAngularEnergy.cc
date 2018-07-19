//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 080808 Bug fix in serching mu bin and index for theBuff2b by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPLabAngularEnergy.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"

void G4ParticleHPLabAngularEnergy::Init(std::istream & aDataFile)
{
  aDataFile >> nEnergies;
  theManager.Init(aDataFile);
  theEnergies = new G4double[nEnergies];
  nCosTh = new G4int[nEnergies];
  theData = new G4ParticleHPVector * [nEnergies];
  theSecondManager = new G4InterpolationManager [nEnergies];
  for(G4int i=0; i<nEnergies; i++)
  {
    aDataFile >> theEnergies[i];
    theEnergies[i]*=eV;
    aDataFile >> nCosTh[i];
    theSecondManager[i].Init(aDataFile); 
    theData[i] = new G4ParticleHPVector[nCosTh[i]];
    G4double label;
    for(G4int ii=0; ii<nCosTh[i]; ii++)
    {
      aDataFile >> label;
      theData[i][ii].SetLabel(label);
      theData[i][ii].Init(aDataFile, eV);
    }
  }
}

G4ReactionProduct * G4ParticleHPLabAngularEnergy::Sample(G4double anEnergy, G4double massCode, G4double )
{
   G4ReactionProduct * result = new G4ReactionProduct;
   G4int Z = static_cast<G4int>(massCode/1000);
   G4int A = static_cast<G4int>(massCode-1000*Z);

   if(massCode==0)
   {
     result->SetDefinition(G4Gamma::Gamma());
   }
   else if(A==0)
   {
     result->SetDefinition(G4Electron::Electron());     
     if(Z==1) result->SetDefinition(G4Positron::Positron());
   }
   else if(A==1)
   {
     result->SetDefinition(G4Neutron::Neutron());
     if(Z==1) result->SetDefinition(G4Proton::Proton());
   }
   else if(A==2)
   {
     result->SetDefinition(G4Deuteron::Deuteron());      
   }
   else if(A==3)
   {
     result->SetDefinition(G4Triton::Triton());  
     if(Z==2) result->SetDefinition(G4He3::He3());
   }
   else if(A==4)
   {
     result->SetDefinition(G4Alpha::Alpha());
     if(Z!=2) throw G4HadronicException(__FILE__, __LINE__, "Unknown ion case 1");    
   }
   else
   {
     throw G4HadronicException(__FILE__, __LINE__, "G4ParticleHPLabAngularEnergy: Unknown ion case 2");
   }
   
   // get theta, E
   G4double cosTh, secEnergy;
   G4int i, it(0);
   // find the energy bin
   for(i=0; i<nEnergies; i++)
   {
     it = i;
     if ( anEnergy < theEnergies[i] ) break;
   }
   //080808
   //if ( it == 0 || it == nEnergies-1 ) // it marks the energy bin
   if ( it == 0 ) // it marks the energy bin
   {
     G4cout << "080808 Something unexpected is happen in G4ParticleHPLabAngularEnergy " << G4endl;
     // integrate the prob for each costh, and select theta.
     G4double * running = new G4double [nCosTh[it]];
     running[0]=0;
     for(i=0;i<nCosTh[it]; i++)
     {
       if(i!=0) running[i] = running[i-1];
       running[i]+=theData[it][i].GetIntegral(); // Does interpolated integral.
     }
     G4double random = running[nCosTh[it]-1]*G4UniformRand();
     G4int ith(0);
     for(i=0;i<nCosTh[it]; i++)
     {
       ith = i;
       if(random<running[i]) break;
     }
     //080807
     //if ( ith == 0 || ith == nCosTh[it]-1 ) //ith marks the angluar bin
     if ( ith == 0 ) //ith marks the angluar bin
     {
        cosTh = theData[it][ith].GetLabel();
        secEnergy = theData[it][ith].Sample();
        currentMeanEnergy = theData[it][ith].GetMeanX();
     }
     else
     {
       //080808
       //G4double x1 = theData[it][ith-1].GetIntegral();
       //G4double x2 = theData[it][ith].GetIntegral();
       G4double x1 = running [ ith-1 ];
       G4double x2 = running [ ith ];
       G4double x = random;
       G4double y1 = theData[it][ith-1].GetLabel();
       G4double y2 = theData[it][ith].GetLabel();
       cosTh = theInt.Interpolate(theSecondManager[it].GetInverseScheme(ith),
                                  x, x1, x2, y1, y2);
       G4ParticleHPVector theBuff1;
       theBuff1.SetInterpolationManager(theData[it][ith-1].GetInterpolationManager());
       G4ParticleHPVector theBuff2;
       theBuff2.SetInterpolationManager(theData[it][ith].GetInterpolationManager());
       x1=y1;
       x2=y2;
       G4double y, mu;
       for(i=0;i<theData[it][ith-1].GetVectorLength(); i++)
       {
         mu = theData[it][ith-1].GetX(i);
         y1 = theData[it][ith-1].GetY(i);
         y2 = theData[it][ith].GetY(mu);

         y = theInt.Interpolate(theSecondManager[it].GetScheme(ith), 
                                cosTh, x1,x2,y1,y2);
         theBuff1.SetData(i, mu, y);
       }
       for(i=0;i<theData[it][ith].GetVectorLength(); i++)
       {
         mu = theData[it][ith].GetX(i);
         y1 = theData[it][ith-1].GetY(mu);
         y2 = theData[it][ith].GetY(i);
         y = theInt.Interpolate(theSecondManager[it].GetScheme(ith), 
                                cosTh, x1,x2,y1,y2);
         theBuff2.SetData(i, mu, y);
       }
       G4ParticleHPVector theStore;
       theStore.Merge(&theBuff1, &theBuff2);
       secEnergy = theStore.Sample();
       currentMeanEnergy = theStore.GetMeanX();
     }
     delete [] running;
   }
   else // this is the small big else.
   {
     G4double x, x1, x2, y1, y2, y, tmp, E;
     // integrate the prob for each costh, and select theta.
     G4ParticleHPVector run1;
     run1.SetY(0, 0.);
     for(i=0;i<nCosTh[it-1]; i++)
     {
       if(i!=0) run1.SetY(i, run1.GetY(i-1));
       run1.SetX(i, theData[it-1][i].GetLabel());
       run1.SetY(i, run1.GetY(i)+theData[it-1][i].GetIntegral());
     }
     G4ParticleHPVector run2;
     run2.SetY(0, 0.); 
     for(i=0;i<nCosTh[it]; i++)
     {
       if(i!=0) run2.SetY(i, run2.GetY(i-1));
       run2.SetX(i, theData[it][i].GetLabel());
       run2.SetY(i, run2.GetY(i)+theData[it][i].GetIntegral());
     }
     // get the distributions for the correct neutron energy
     x = anEnergy;
     x1 = theEnergies[it-1];
     x2 = theEnergies[it];
     G4ParticleHPVector thBuff1; // to be interpolated as run1.
     thBuff1.SetInterpolationManager(theSecondManager[it-1]);
     for(i=0; i<run1.GetVectorLength(); i++)
     {
       tmp = run1.GetX(i); //theta
       y1 = run1.GetY(i); // integral
       y2 = run2.GetY(tmp);
       y = theInt.Interpolate(theManager.GetScheme(it), x, x1,x2,y1,y2);
       thBuff1.SetData(i, tmp, y);
     }
     G4ParticleHPVector thBuff2;
     thBuff2.SetInterpolationManager(theSecondManager[it]);
     for(i=0; i<run2.GetVectorLength(); i++)
     {
       tmp = run2.GetX(i); //theta
       y1 = run1.GetY(tmp); // integral
       y2 = run2.GetY(i);
       y = theInt.Lin(x, x1,x2,y1,y2);
       thBuff2.SetData(i, tmp, y);
     }
     G4ParticleHPVector theThVec;
     theThVec.Merge(&thBuff1 ,&thBuff2); // takes care of interpolation
     G4double random = (theThVec.GetY(theThVec.GetVectorLength()-1)
                        -theThVec.GetY(0))   *G4UniformRand();
     G4int ith(0);
     for(i=1;i<theThVec.GetVectorLength(); i++)
     {
       ith = i;
       if(random<theThVec.GetY(i)-theThVec.GetY(0)) break;
     }
     {
       // calculate theta
       G4double xx, xx1, xx2, yy1, yy2;
       xx =  random;
       xx1 = theThVec.GetY(ith-1)-theThVec.GetY(0); // integrals
       xx2 = theThVec.GetY(ith)-theThVec.GetY(0);
       yy1 = theThVec.GetX(ith-1); // std::cos(theta)
       yy2 = theThVec.GetX(ith);
       cosTh = theInt.Interpolate(theSecondManager[it].GetScheme(ith), 
                                  xx, xx1,xx2,yy1,yy2);
     }
     G4int i1(0), i2(0);
     // get the indixes of the vectors close to theta for low energy
     // first it-1 !!!! i.e. low in energy
     for(i=0; i<nCosTh[it-1]; i++)
     {
       i1 = i;
       if(cosTh<theData[it-1][i].GetLabel()) break;
     }
     // now get the prob at this energy for the right theta value
     x = cosTh;
     x1 = theData[it-1][i1-1].GetLabel();
     x2 = theData[it-1][i1].GetLabel();
     G4ParticleHPVector theBuff1a;
     theBuff1a.SetInterpolationManager(theData[it-1][i1-1].GetInterpolationManager());
     for(i=0;i<theData[it-1][i1-1].GetVectorLength(); i++)
     {
       E = theData[it-1][i1-1].GetX(i);
       y1 = theData[it-1][i1-1].GetY(i);
       y2 = theData[it-1][i1].GetY(E);
       y = theInt.Lin(x, x1,x2,y1,y2);
       theBuff1a.SetData(i, E, y); // wrong E, right theta.
     }
     G4ParticleHPVector theBuff2a;
     theBuff2a.SetInterpolationManager(theData[it-1][i1].GetInterpolationManager());
     for(i=0;i<theData[it-1][i1].GetVectorLength(); i++)
     {
       E = theData[it-1][i1].GetX(i);
       y1 = theData[it-1][i1-1].GetY(E);
       y2 = theData[it-1][i1].GetY(i);
       y = theInt.Lin(x, x1,x2,y1,y2);
       theBuff2a.SetData(i, E, y); // wrong E, right theta.
     }
     G4ParticleHPVector theStore1;
     theStore1.Merge(&theBuff1a, &theBuff2a); // wrong E, right theta, complete binning

     // get the indixes of the vectors close to theta for high energy
     // then it !!!! i.e. high in energy
     for(i=0; i<nCosTh[it]; i++)
     {
       i2 = i;
       if(cosTh<theData[it][i2].GetLabel()) break;
     }                  // sonderfaelle mit i1 oder i2 head on fehlen. @@@@@
     x1 = theData[it][i2-1].GetLabel();
     x2 = theData[it][i2].GetLabel();
     G4ParticleHPVector theBuff1b;
     theBuff1b.SetInterpolationManager(theData[it][i2-1].GetInterpolationManager());
     for(i=0;i<theData[it][i2-1].GetVectorLength(); i++)
     {
       E = theData[it][i2-1].GetX(i);
       y1 = theData[it][i2-1].GetY(i);
       y2 = theData[it][i2].GetY(E);
       y = theInt.Lin(x, x1,x2,y1,y2);
       theBuff1b.SetData(i, E, y); // wrong E, right theta.
     }
     G4ParticleHPVector theBuff2b;
     theBuff2b.SetInterpolationManager(theData[it][i2].GetInterpolationManager());
     //080808  i1 -> i2
     //for(i=0;i<theData[it][i1].GetVectorLength(); i++)
     for(i=0;i<theData[it][i2].GetVectorLength(); i++)
     {
       //E = theData[it][i1].GetX(i);
       //y1 = theData[it][i1-1].GetY(E);
       //y2 = theData[it][i1].GetY(i);
       E = theData[it][i2].GetX(i);
       y1 = theData[it][i2-1].GetY(E);
       y2 = theData[it][i2].GetY(i);
       y = theInt.Lin(x, x1,x2,y1,y2);
       theBuff2b.SetData(i, E, y); // wrong E, right theta.
     }
     G4ParticleHPVector theStore2;
     theStore2.Merge(&theBuff1b, &theBuff2b); // wrong E, right theta, complete binning
     // now get to the right energy.

     x = anEnergy;
     x1 = theEnergies[it-1];
     x2 = theEnergies[it];
     G4ParticleHPVector theOne1;
     theOne1.SetInterpolationManager(theStore1.GetInterpolationManager());
     for(i=0; i<theStore1.GetVectorLength(); i++)
     {
       E = theStore1.GetX(i);
       y1 = theStore1.GetY(i);
       y2 = theStore2.GetY(E);
       y = theInt.Interpolate(theManager.GetScheme(it), x, x1,x2,y1,y2);
       theOne1.SetData(i, E, y); // both correct
     }
     G4ParticleHPVector theOne2;
     theOne2.SetInterpolationManager(theStore2.GetInterpolationManager());
     for(i=0; i<theStore2.GetVectorLength(); i++)
     {
       E = theStore2.GetX(i);
       y1 = theStore1.GetY(E);
       y2 = theStore2.GetY(i);
       y = theInt.Interpolate(theManager.GetScheme(it), x, x1,x2,y1,y2);
       theOne2.SetData(i, E, y); // both correct
     }
     G4ParticleHPVector theOne;
     theOne.Merge(&theOne1, &theOne2); // both correct, complete binning

     secEnergy = theOne.Sample();
     currentMeanEnergy = theOne.GetMeanX();
   }

// now do random direction in phi, and fill the result.

   result->SetKineticEnergy(secEnergy);
   
   G4double phi = twopi*G4UniformRand();
   G4double theta = std::acos(cosTh);
   G4double sinth = std::sin(theta);
   G4double mtot = result->GetTotalMomentum(); 
   G4ThreeVector tempVector(mtot*sinth*std::cos(phi), mtot*sinth*std::sin(phi), mtot*std::cos(theta) );
   result->SetMomentum(tempVector);
   
   return result;
}
