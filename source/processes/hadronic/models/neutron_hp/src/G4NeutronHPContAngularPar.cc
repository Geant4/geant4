// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPContAngularPar.hh"
#include "G4NeutronHPLegendreStore.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4NeutronHPVector.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4NeutronHPKallbachMannSyst.hh"
#include "G4ParticleTable.hh"
 
  G4NeutronHPContAngularPar::G4NeutronHPContAngularPar()
  {
    theAngular = NULL;
    currentMeanEnergy = -2;
  }
  G4NeutronHPContAngularPar::~G4NeutronHPContAngularPar()
  {
    if(theAngular!=NULL) delete [] theAngular;
  }

  void G4NeutronHPContAngularPar::Init(ifstream & aDataFile)
  {
    aDataFile >> theEnergy >> nEnergies >> nDiscreteEnergies >> nAngularParameters;
    theEnergy *= eV;
    theAngular = new G4NeutronHPList [nEnergies];
    for(G4int i=0; i<nEnergies; i++)
    {
      G4double sEnergy;
      aDataFile >> sEnergy;
      sEnergy*=eV;
      theAngular[i].SetLabel(sEnergy);
      theAngular[i].Init(aDataFile, nAngularParameters, 1.);
    }
  }

  G4ReactionProduct * 
  G4NeutronHPContAngularPar::Sample(G4double anEnergy, G4double massCode, G4double targetMass, 
                                    G4int angularRep, G4int interpolE)
  {
    G4ReactionProduct * result = new G4ReactionProduct;
    G4int Z = massCode/1000;
    G4int A = massCode-1000*Z;
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
      if(Z!=2) G4Exception("G4NeutronHPContAngularPar: Unknown ion case 1");    
    }
    else
    {
      result->SetDefinition(G4ParticleTable::GetParticleTable()->FindIon(Z,A,0,Z));
    }
    G4int i;
    G4int it;
    G4double fsEnergy;
    G4double cosTh;
    if(angularRep==1)
    {
      G4double random = G4UniformRand();
      G4double * running = new G4double[nEnergies];
      running[0]=0;
      G4double weighted = 0;
      for(i=1; i<nEnergies; i++)
      {
        if(i!=0) 
        {
          running[i]=running[i-1];
        }
        running[i] += theInt.GetBinIntegral(theManager.GetScheme(i-1),
                             theAngular[i-1].GetLabel(), theAngular[i].GetLabel(),
                             theAngular[i-1].GetValue(0), theAngular[i].GetValue(0));
        weighted += theInt.GetWeightedBinIntegral(theManager.GetScheme(i-1),
                             theAngular[i-1].GetLabel(), theAngular[i].GetLabel(),
                             theAngular[i-1].GetValue(0), theAngular[i].GetValue(0));
      }
       // cash the mean energy in this distribution
      currentMeanEnergy = weighted/running[nEnergies-1];
      
      for(i=0; i<nEnergies; i++)
      {
        it = i;
        if(random<running[i]/running[nEnergies-1]) break;
      }
      if(it<nDiscreteEnergies||it==0) 
      {
        if(it == 0)
        {
          fsEnergy = theAngular[it].GetLabel();
          G4NeutronHPLegendreStore theStore(1);
          theStore.Init(0,fsEnergy,nAngularParameters);
          for(i=0;i<nAngularParameters;i++)
          {
            theStore.SetCoeff(0,i,theAngular[it].GetValue(i));
          }
          // use it to sample.
          cosTh = theStore.SampleMax(fsEnergy);
        }
        else
        {
          G4double e1, e2;
          e1 = theAngular[it].GetLabel();
          e2 = theAngular[it-1].GetLabel();
          fsEnergy = theInt.Interpolate(theManager.GetInverseScheme(it),
                                        random,
                                        running[it-1]/running[nEnergies-1], 
                                        running[it]/running[nEnergies-1],
                                        e1, e2);
          // fill a Legendrestore
          G4NeutronHPLegendreStore theStore(2);
          theStore.Init(0,e1,nAngularParameters);
          theStore.Init(1,e2,nAngularParameters);
          for(i=0;i<nAngularParameters;i++)
          {
            theStore.SetCoeff(0,i,theAngular[it-1].GetValue(i));
            theStore.SetCoeff(1,i,theAngular[it].GetValue(i));
          }
          // use it to sample.
          theStore.SetManager(theManager);
          cosTh = theStore.SampleMax(fsEnergy);
        }
      }
      else // continuum contribution
      {
        G4double x1 = running[it-1]/running[nEnergies-1];
        G4double x2 = running[it]/running[nEnergies-1];
        G4double y1 = theAngular[it-1].GetLabel();
        G4double y2 = theAngular[it].GetLabel();
        fsEnergy = theInt.Interpolate(theManager.GetInverseScheme(it),
                                      random,x1,x2,y1,y2);
        G4NeutronHPLegendreStore theStore(2);
        theStore.Init(0,y1,nAngularParameters);
        theStore.Init(1,y2,nAngularParameters);
        theStore.SetManager(theManager);
        for(i=0;i<nAngularParameters;i++)
        {
          theStore.SetCoeff(0,i,theAngular[it-1].GetValue(i));
          theStore.SetCoeff(1,i,theAngular[it].GetValue(i));
        }
        // use it to sample.
        cosTh = theStore.SampleMax(fsEnergy);
      }
      delete [] running;
    }
    else if(angularRep==2)
    {
      // first get the energy (already the right for this incoming energy)
      G4int i;
      G4double * running = new G4double[nEnergies];
      running[0]=0;
      G4double weighted = 0;
      for(i=1; i<nEnergies; i++)
      {
        if(i!=0) running[i]=running[i-1];
        running[i] += theInt.GetBinIntegral(theManager.GetScheme(i-1),
                             theAngular[i-1].GetLabel(), theAngular[i].GetLabel(),
                             theAngular[i-1].GetValue(0), theAngular[i].GetValue(0));
        weighted += theInt.GetWeightedBinIntegral(theManager.GetScheme(i-1),
                             theAngular[i-1].GetLabel(), theAngular[i].GetLabel(),
                             theAngular[i-1].GetValue(0), theAngular[i].GetValue(0));
      }
       // cash the mean energy in this distribution
      currentMeanEnergy = weighted/running[nEnergies-1];
      
      G4int it;
      G4double randkal = G4UniformRand();
      for(i=0; i<nEnergies; i++)
      {
        it = i;
        if(randkal<running[i]/running[nEnergies-1]) break;
      }
      
      // interpolate the secondary energy.
      G4double x, x1,x2,y1,y2;
      if(it==0) it=1;
      x = randkal*running[nEnergies-1];
      x1 = running[it-1];
      x2 = running[it];
      G4double compoundFraction;
      // interpolate energy
      y1 = theAngular[it-1].GetLabel();
      y2 = theAngular[it].GetLabel();
      fsEnergy = theInt.Interpolate(theManager.GetInverseScheme(it-1), 
                                    x, x1,x2,y1,y2);
      // for theta interpolate the compoundFractions
      G4double cLow = theAngular[it-1].GetValue(1);
      G4double cHigh = theAngular[it].GetValue(1);
      compoundFraction = theInt.Interpolate(theManager.GetScheme(it),
                                            fsEnergy, y1, y2, cLow,cHigh);
      delete [] running;
      
      // get cosTh
      G4double incidentEnergy = anEnergy;
      G4double incidentMass = G4Neutron::Neutron()->GetPDGMass();
      G4double productEnergy = fsEnergy;
      G4double productMass = result->GetMass();
      G4int targetZ = G4int(theTargetCode/1000);
      G4int targetA = G4int(theTargetCode-1000*targetZ);
      G4double targetMass = theTarget->GetMass();
      G4int residualA = targetA+1-A;
      G4int residualZ = targetZ-Z;
      G4double residualMass =  residualZ*G4Proton::Proton()->GetPDGMass();
               residualMass +=(residualA-residualZ)*G4Neutron::Neutron()->GetPDGMass();
               residualMass -= G4NucleiPropertiesTable::GetBindingEnergy(residualZ, residualA);
      G4NeutronHPKallbachMannSyst theKallbach(compoundFraction,
                                              incidentEnergy, incidentMass,
                                              productEnergy, productMass,
                                              residualMass, residualA, residualZ,
                                              targetMass, targetA, targetZ);
      cosTh = theKallbach.Sample(anEnergy);
    }
    else if(angularRep>10&&angularRep<16)
    {
      G4double random = G4UniformRand();
      G4double * running = new G4double[nEnergies];
      running[0]=0;      
      G4double weighted = 0;
      for(i=1; i<nEnergies; i++)
      {
        if(i!=0) running[i]=running[i-1];
        running[i] += theInt.GetBinIntegral(theManager.GetScheme(i-1),
                             theAngular[i-1].GetLabel(), theAngular[i].GetLabel(),
                             theAngular[i-1].GetValue(0), theAngular[i].GetValue(0));
        weighted += theInt.GetWeightedBinIntegral(theManager.GetScheme(i-1),
                             theAngular[i-1].GetLabel(), theAngular[i].GetLabel(),
                             theAngular[i-1].GetValue(0), theAngular[i].GetValue(0));
      }
       // cash the mean energy in this distribution
      currentMeanEnergy = weighted/running[nEnergies-1];
      
      for(i=0; i<nEnergies; i++)
      {
        it = i;
        if(random<running[i]/running[nEnergies-1]) break;
      }
      if(it<nDiscreteEnergies||it==0) 
      {
        if(it==0)
        {
          fsEnergy = theAngular[0].GetLabel();          
          G4NeutronHPVector theStore; 
	  G4int aCounter = 0;
          for(G4int i=1; i<nAngularParameters; i+=2) 
          {
            theStore.SetX(aCounter, theAngular[0].GetValue(i));
            theStore.SetY(aCounter, theAngular[0].GetValue(i+1));
	    aCounter++;	    
          }
          G4InterpolationManager aMan;
          aMan.Init(angularRep-10, nAngularParameters-1);
          theStore.SetInterpolationManager(aMan);
          cosTh = theStore.Sample();
        }
        else 
        {
          fsEnergy = theAngular[it].GetLabel();
          G4NeutronHPVector theStore; 
          G4InterpolationManager aMan;
          aMan.Init(angularRep-10, nAngularParameters-1);
          theStore.SetInterpolationManager(aMan); // Store interpolates f(costh)
          G4InterpolationScheme currentScheme = theManager.GetInverseScheme(it);
	  G4int aCounter = 0;
          for(G4int i=1; i<nAngularParameters; i+=2) 
          {
            theStore.SetX(aCounter, theAngular[it].GetValue(i));
            theStore.SetY(aCounter, theInt.Interpolate(currentScheme, 
                                       random,
                                       running[it-1]/running[nEnergies-1],
                                       running[it]/running[nEnergies-1],
                                       theAngular[it-1].GetValue(i+1),
                                       theAngular[it].GetValue(i+1)));
	    aCounter++;	    
          }
          cosTh = theStore.Sample();
        }
      }
      else
      {
        G4double x1 = running[it-1]/running[nEnergies-1];
        G4double x2 = running[it]/running[nEnergies-1];
        G4double y1 = theAngular[it-1].GetLabel();
        G4double y2 = theAngular[it].GetLabel();
        fsEnergy = theInt.Interpolate(theManager.GetInverseScheme(it),
                                      random,x1,x2,y1,y2);
        G4NeutronHPVector theBuff1;
        G4NeutronHPVector theBuff2;
        G4InterpolationManager aMan;
        aMan.Init(angularRep-10, nAngularParameters-1);
//        theBuff1.SetInterpolationManager(aMan); // Store interpolates f(costh)
//        theBuff2.SetInterpolationManager(aMan); // Store interpolates f(costh)
        for(i=0; i<nAngularParameters; i++) // i=1 ist wichtig!
        {
          theBuff1.SetX(i, theAngular[it-1].GetValue(i));
          theBuff1.SetY(i, theAngular[it-1].GetValue(i+1));
          theBuff2.SetX(i, theAngular[it].GetValue(i));
          theBuff2.SetY(i, theAngular[it].GetValue(i+1));
          i++;
        }
        G4NeutronHPVector theStore;
        theStore.SetInterpolationManager(aMan); // Store interpolates f(costh)        
        x1 = y1;
        x2 = y2;
        G4double x, y;
        for(i=0;i<theBuff1.GetVectorLength(); i++);
        {
          x = theBuff1.GetX(i); // costh binning identical
          y1 = theBuff1.GetY(i);
          y2 = theBuff2.GetY(i);
          y = theInt.Interpolate(theManager.GetScheme(it),
                                 fsEnergy, theAngular[it-1].GetLabel(), 
                                 theAngular[it].GetLabel(), y1, y2);
          theStore.SetX(i, x);
          theStore.SetY(i, y);
        }
        cosTh = theStore.Sample();
      }
      delete [] running;
    }
    else
    {
      G4Exception("G4NeutronHPContAngularPar::Sample: Unknown angular representation");
    }
    result->SetKineticEnergy(fsEnergy);
    G4double phi = twopi*G4UniformRand();
    G4double theta = acos(cosTh);
    G4double sinth = sin(theta);
    G4double mtot = result->GetTotalMomentum();
    G4ThreeVector tempVector(mtot*sinth*cos(phi), mtot*sinth*sin(phi), mtot*cos(theta) );
    result->SetMomentum(tempVector);
//  return the result.    
    return result;
  }
