// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPDiscreteTwoBody.hh"
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
#include "G4NeutronHPLegendreStore.hh"

  G4NeutronHPDiscreteTwoBody::G4NeutronHPDiscreteTwoBody()
  {
    theCoeff = NULL;
  }
  G4NeutronHPDiscreteTwoBody::~G4NeutronHPDiscreteTwoBody()
  {
    if(theCoeff!=NULL) delete [] theCoeff;
  }
  
G4ReactionProduct * G4NeutronHPDiscreteTwoBody::Sample(G4double anEnergy, G4double massCode, G4double mass)
{ // Interpolation still only for the most used parts; rest to be Done @@@@@
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
     if(Z!=2) G4Exception("Unknown ion case 1");    
   }
   else
   {
     G4Exception("G4NeutronHPDiscreteTwoBody: Unknown ion case 2");
   }
   
// get cosine(theta)
   G4int i, it;
   G4double cosTh;
   for(i=0; i<nEnergy; i++)
   {
     it = i;
     if(theCoeff[i].GetEnergy()>anEnergy) break;
   }
   if(it==0||it==nEnergy-1)
   {
     if(theCoeff[it].GetRepresentation()==0)
     {
       G4NeutronHPLegendreStore theStore(1);
       theStore.SetCoeff(0, theCoeff);
       theStore.SetManager(theManager);
       cosTh = theStore.SampleMax(anEnergy);
     }
     else if(theCoeff[it].GetRepresentation()==12) // means LINLIN
     {
       G4NeutronHPVector theStore; 
       G4InterpolationManager aManager;
       aManager.Init(LINLIN, theCoeff[it].GetNumberOfPoly()/2);
       theStore.SetInterpolationManager(aManager);
       for(i=0;i<theCoeff[it].GetNumberOfPoly(); i++)
       {
         theStore.SetX(i, theCoeff[it].GetCoeff(i));
         theStore.SetY(i, theCoeff[it].GetCoeff(i));
	 i++;
       }
       cosTh = theStore.Sample();
     }
     else if(theCoeff[it].GetRepresentation()==14) //this is LOGLIN
     {
       G4NeutronHPVector theStore;
       G4InterpolationManager aManager;
       aManager.Init(LOGLIN, theCoeff[it].GetNumberOfPoly()/2);
       theStore.SetInterpolationManager(aManager);
       for(i=0;i<theCoeff[it].GetNumberOfPoly(); i++)
       {
         theStore.SetX(i, theCoeff[it].GetCoeff(i));
         theStore.SetY(i, theCoeff[it].GetCoeff(i));
	 i++;
       }
       cosTh = theStore.Sample(); 
     }
     else
     {
       G4Exception("unknown representation type in Two-body scattering");
     }
   }
   else
   {
     if(theCoeff[it].GetRepresentation() == theCoeff[it-1].GetRepresentation())
     {
       if(theCoeff[it].GetRepresentation()==0)
       {
	 G4NeutronHPLegendreStore theStore(2);
	 theStore.SetCoeff(0, &(theCoeff[it-1]));
	 theStore.SetCoeff(1, &(theCoeff[it]));
         G4InterpolationManager aManager;
         aManager.Init(theManager.GetScheme(it), 2);
         theStore.SetManager(aManager);
	 cosTh = theStore.SampleMax(anEnergy);
       }
       else if(theCoeff[it].GetRepresentation()==12) // LINLIN
       {
	 G4NeutronHPVector theBuff1;
         G4InterpolationManager aManager1;
         aManager1.Init(LINLIN, theCoeff[it-1].GetNumberOfPoly()/2);
         theBuff1.SetInterpolationManager(aManager1);
	 for(i=0;i<theCoeff[it-1].GetNumberOfPoly(); i++)
	 {
           theBuff1.SetX(i, theCoeff[it-1].GetCoeff(i));
           theBuff1.SetY(i, theCoeff[it-1].GetCoeff(i));
	   i++;
	 }
	 G4NeutronHPVector theBuff2;
         G4InterpolationManager aManager2;
         aManager2.Init(LINLIN, theCoeff[it].GetNumberOfPoly()/2);
         theBuff2.SetInterpolationManager(aManager2);
	 for(i=0;i<theCoeff[it].GetNumberOfPoly(); i++)
	 {
           theBuff2.SetX(i, theCoeff[it].GetCoeff(i));
           theBuff2.SetY(i, theCoeff[it].GetCoeff(i));
	   i++;
	 }

	 G4double x1 = theCoeff[it-1].GetEnergy();
	 G4double x2 = theCoeff[it].GetEnergy();
	 G4double x = anEnergy;
	 G4double y1, y2, y, mu;

	 G4NeutronHPVector theStore1;
         theStore1.SetInterpolationManager(aManager1);
	 G4NeutronHPVector theStore2;
         theStore2.SetInterpolationManager(aManager2);
	 G4NeutronHPVector theStore;
	 
	 // for fixed mu get p1, p2 and interpolate according to x
         for(i=0; i<theBuff1.GetVectorLength(); i++)
         {
           mu = theBuff1.GetX(i);
           y1 = theBuff1.GetY(i);
           y2 = theBuff2.GetY(mu);
           y = theInt.Interpolate(theManager.GetScheme(it), x, x1,x2,y1,y2);
           theStore1.SetData(i, mu, y);
         }
         for(i=0; i<theBuff2.GetVectorLength(); i++)
         {
           mu = theBuff2.GetX(i);
           y1 = theBuff2.GetY(i);
           y2 = theBuff1.GetY(mu);
           y = theInt.Interpolate(theManager.GetScheme(it), x, x1,x2,y1,y2);
           theStore2.SetData(i, mu, y);
         }
         theStore.Merge(&theStore1, &theStore2); // merge takes care of interpolationschemes
	 cosTh = theStore.Sample();
       }
       else if(theCoeff[it].GetRepresentation()==14) 
       {
	 G4NeutronHPVector theBuff1;
         G4InterpolationManager aManager1;
         aManager1.Init(LOGLIN, theCoeff[it-1].GetNumberOfPoly()/2);
         theBuff1.SetInterpolationManager(aManager1);
	 for(i=0;i<theCoeff[it-1].GetNumberOfPoly(); i++)
	 {
           theBuff1.SetX(i, theCoeff[it-1].GetCoeff(i));
           theBuff1.SetY(i, theCoeff[it-1].GetCoeff(i));
	   i++;
	 }
	 
	 G4NeutronHPVector theBuff2;
         G4InterpolationManager aManager2;
         aManager2.Init(LOGLIN, theCoeff[it].GetNumberOfPoly()/2);
         theBuff2.SetInterpolationManager(aManager2);
	 for(i=0;i<theCoeff[it].GetNumberOfPoly(); i++)
	 {
           theBuff2.SetX(i, theCoeff[it].GetCoeff(i));
           theBuff2.SetY(i, theCoeff[it].GetCoeff(i));
	   i++;
	 }

	 G4double x1 = theCoeff[it-1].GetEnergy();
	 G4double x2 = theCoeff[it].GetEnergy();
	 G4double x = anEnergy;
	 G4double y1, y2, y, mu;

	 G4NeutronHPVector theStore1;
         theStore1.SetInterpolationManager(aManager1);
	 G4NeutronHPVector theStore2;
         theStore2.SetInterpolationManager(aManager2);
	 G4NeutronHPVector theStore;
	 
	 // for fixed mu get p1, p2 and interpolate according to x
         for(i=0; i<theBuff1.GetVectorLength(); i++)
         {
           mu = theBuff1.GetX(i);
           y1 = theBuff1.GetY(i);
           y2 = theBuff2.GetY(mu);
           y = theInt.Interpolate(theManager.GetScheme(it), x, x1,x2,y1,y2);
           theStore1.SetData(i, mu, y);
         }
         for(i=0; i<theBuff2.GetVectorLength(); i++)
         {
           mu = theBuff2.GetX(i);
           y1 = theBuff2.GetY(i);
           y2 = theBuff1.GetY(mu);
           y = theInt.Interpolate(theManager.GetScheme(it), x, x1,x2,y1,y2);
           theStore2.SetData(i, mu, y);
         }
         theStore.Merge(&theStore1, &theStore2); 
	 cosTh = theStore.Sample(); 
       }
       else
       {
         G4Exception("Two neighbouring distributions with different interpolation");
       }
     }
     else
     {
       G4Exception("unknown representation type in Two-body scattering, case 2");
     }
   }
   
// now get the energy from kinematics and Q-value.

   G4double restEnergy = anEnergy+GetQValue();
   
// assumed to be in CMS @@@@@@@@@@@@@@@@@

   G4double residualMass =   GetTarget()->GetMass() + GetNeutron()->GetMass()
                           - result->GetMass() - GetQValue();
   G4double kinE = restEnergy/(1+result->GetMass()/residualMass); // non relativistic @@
   result->SetKineticEnergy(kinE); // non relativistic @@
   G4double phi = twopi*G4UniformRand();
   G4double theta = acos(cosTh);
   G4double sinth = sin(theta);
   G4double mtot = result->GetTotalMomentum(); 
   G4ThreeVector tempVector(mtot*sinth*cos(phi), mtot*sinth*sin(phi), mtot*cos(theta) );
   result->SetMomentum(tempVector);
   
// some garbage collection
   
// return the result   
   return result;
}
