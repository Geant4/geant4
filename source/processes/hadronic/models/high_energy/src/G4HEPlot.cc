// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HEPlot.cc,v 1.3 1999-12-15 14:52:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#include "globals.hh"
#include "G4ios.hh"

//
// G4 Gheisha friend class G4GHEPlot
// last modified: H. Fesefeldt 02-July--1998

#include "G4HEPlot.hh"
  

void 
G4HEPlot::Init( G4int nbin, G4double xstart, G4double xbin ) 
   {
      Xbin = xbin;
      Nbin = nbin;
      Xstart = xstart;
      Entries = 0;
      EntriesOverflow = 0;
      EntriesUnderflow = 0;
      Weight = 0.;
      WeightOverflow = 0.;
      WeightUnderflow = 0.;
      Xvalue = new G4double[Nbin];
      Yvalue = new G4double[Nbin];    
      for(G4int i=0; i<Nbin; i++)
        {
          Xvalue[i] = Xstart + i*Xbin;
          Yvalue[i] = 0.;
        } 
      return; 
   }
 
void 
G4HEPlot::Add( G4double s1, G4double s2, 
               const G4HEPlot & p1, const G4HEPlot & p2 )
   {
     if(p1.Nbin != p2.Nbin)
       {
         G4cout << "G4HEPlot::Add: Plots must have same number of bins !" << G4endl;
         return;
       }
     if(Nbin != p1.Nbin)
       {
         G4cout << "G4HEPlot::Add: Plot must be initialized  before using it !" << G4endl;
         return;
       }  

     for(G4int i=0; i<Nbin; i++)
       { 
         Yvalue[i] = s1*p1.Yvalue[i] + s2*p2.Yvalue[i];
       }
     Entries = p1.Entries + p2.Entries;
     EntriesOverflow = p1.EntriesOverflow + p2.EntriesOverflow;
     EntriesUnderflow = p1.EntriesUnderflow + p2.EntriesUnderflow;
     Weight = s1*p1.Weight + s2*p2.Weight;
     WeightUnderflow = s1*p1.WeightUnderflow + s2*p2.WeightUnderflow;
     WeightOverflow = s1*p1.WeightOverflow + s2*p2.WeightOverflow;
     return; 
   }

void 
G4HEPlot::Multiply( G4double s1, G4double s2, 
                    const G4HEPlot & p1, const G4HEPlot & p2 )
   {
     if(p1.Nbin != p2.Nbin)
       {
         G4cout << "G4HEPlot::Multiply: Plots must have same number of bins !" << G4endl;
         return;
       }
     if(Nbin != p1.Nbin)
       {
         G4cout << "G4HEPlot::Multiply: Plot must be initialized  before using it !" << G4endl;
         return;
       }  

     for(G4int i=0; i<Nbin; i++)
       { 
         Yvalue[i] = s1*p1.Yvalue[i]* s2*p2.Yvalue[i];
       }
     Entries = p1.Entries * p2.Entries;
     EntriesOverflow = p1.EntriesOverflow * p2.EntriesOverflow;
     EntriesUnderflow = p1.EntriesUnderflow * p2.EntriesUnderflow;
     Weight = s1*p1.Weight * s2*p2.Weight;
     WeightUnderflow = s1*p1.WeightUnderflow * s2*p2.WeightUnderflow;
     WeightOverflow = s1*p1.WeightOverflow * s2*p2.WeightOverflow;
     return;
  }

void 
G4HEPlot::Divide( G4double s1, G4double s2, 
                  const G4HEPlot & p1, const G4HEPlot & p2 )
   {
     if(p1.Nbin != p2.Nbin)
       {
         G4cout << "G4HEPlot::Divide: Plots must have same number of bins !" << G4endl;
         return;
       }
     if(Nbin != p1.Nbin)
       {
         G4cout << "G4HEPlot::Divide: Plot must be defined before using it !" << G4endl;
         return;
       }  

     for(G4int i=0; i<Nbin; i++)
       { 
         if(p2.Yvalue[i] == 0.) 
           { 
             Yvalue[i] = 0.;
           }
         else
           {
             Yvalue[i] = (s1*p1.Yvalue[i]) / (s2*p2.Yvalue[i]);
           }
       } 
     if(p2.Entries > 0)
       Entries = p1.Entries / p2.Entries;
     if(p2.EntriesOverflow > 0)
       EntriesOverflow = p1.EntriesOverflow / p2.EntriesOverflow;
     if(p2.EntriesUnderflow > 0)
       EntriesUnderflow = p1.EntriesUnderflow / p2.EntriesUnderflow;
     if(p2.Weight != 0.)
       Weight = s1*p1.Weight / s2*p2.Weight;
     if(p2.WeightUnderflow != 0.)
       WeightUnderflow = s1*p1.WeightUnderflow  / s2*p2.WeightUnderflow;
     if(p2.WeightOverflow != 0.)
       WeightOverflow = s1*p1.WeightOverflow / s2*p2.WeightOverflow;
     return;
   }

void 
G4HEPlot::Scale( G4double s, const G4HEPlot & p)
   {
     if(Nbin != p.Nbin)
       {
         G4cout << "G4HEPlot::Add: Plot must be defined before using it !" << G4endl;
         return;
       }  

     for(G4int i=0; i<Nbin; i++)
       { 
         Yvalue[i] = s*p.Yvalue[i];
       }
     Entries = p.Entries ;
     EntriesOverflow = p.EntriesOverflow;
     EntriesUnderflow = p.EntriesUnderflow;
     Weight = s*p.Weight;
     WeightUnderflow = s*p.WeightUnderflow;
     WeightOverflow = s*p.WeightOverflow;
     return;
   }

void
G4HEPlot::XScale(G4double a, G4double b, const G4HEPlot & p)
   {
     G4int i;
     Entries = 0;
     EntriesOverflow = 0;
     EntriesUnderflow = 0;
     Weight = 0.;
     WeightOverflow = 0.;
     WeightUnderflow = 0.;
     Xstart = p.Xstart;
     Xbin = p.Xbin;
     Nbin = p.Nbin;
     for(i=0; i<Nbin; i++)
       {
         Xvalue[i] = p.Yvalue[i];
         Yvalue[i] = 0.;
       }
     G4double xval = Xstart - Xbin/2.;
     for(i=0; i<Nbin; i++)
       {
         xval += Xbin;
         Fill(xval*a + b, Xvalue[i]); 
       }
    for(i=0; i<Nbin; i++)
       {
         Xvalue[i] = Xstart + i*Xbin;
       } 
   }      

void 
G4HEPlot::Log( G4double s, const G4HEPlot & p)
   {
     if(Nbin != p.Nbin)
       {
         G4cout << "G4HEPlot::Log: Plot must be defined before using it !" << G4endl;
         return;
       }  

     for(G4int i=0; i<Nbin; i++)
       { 
         if(s*p.Yvalue[i] <= 0.) Yvalue[i] = 0.;
         else  Yvalue[i] = log10(s*p.Yvalue[i]);
       }
     
     Entries = p.Entries ;
     EntriesOverflow = p.EntriesOverflow;
     EntriesUnderflow = p.EntriesUnderflow;
     if(p.Weight > 0)
       Weight = log10(s*p.Weight);
     if(p.WeightUnderflow > 0)
       WeightUnderflow = log10(s*p.WeightUnderflow);
     if(p.WeightOverflow > 0)
       WeightOverflow = log10(s*p.WeightOverflow);
     return;
   }

void
G4HEPlot::Sqrt(G4double s, const G4HEPlot & p)
   {
     if(Nbin != p.Nbin)
       {
         G4cout << " G4HEPlot::Sqrt: Plot must be defined before using it !" << G4endl;
         return;
       }
     for (G4int i=0; i<Nbin; i++)
       {
         if(s*p.Yvalue[i] <= 0.) Yvalue[i] = 0.;
         else   Yvalue[i] = sqrt(s*p.Yvalue[i]);
       }
     Entries = p.Entries;
     EntriesOverflow = p.EntriesOverflow;
     EntriesUnderflow = p.EntriesUnderflow;
     if(s*p.Weight > 0.) Weight = sqrt(s*p.Weight);
     if(s*p.WeightUnderflow > 0.) WeightUnderflow = sqrt(s*p.WeightUnderflow);
     if(s*p.WeightOverflow > 0.) WeightOverflow = sqrt(s*p.WeightOverflow);
     return;
   } 

void 
G4HEPlot::Reset()
   {
     for(G4int i=0; i<Nbin; i++)
       {
         Yvalue[i] = 0.;
       }
     Entries = 0;
     EntriesOverflow = 0;
     EntriesUnderflow = 0;
     Weight = 0.;
     WeightOverflow = 0.;
     WeightUnderflow = 0.;
     return;
   }

void 
G4HEPlot::Fill(G4double x, G4double weight)
   {
     G4int i = (int)((x - Xstart)/Xbin);
     if(i < 0)
       {
         EntriesUnderflow++;
         WeightUnderflow += weight;
       }
     else if(i >= Nbin)
       {
         EntriesOverflow++;
         WeightOverflow += weight;
       }
     else
       {
         Entries++;
         Weight += weight;
         Yvalue[i] += weight;
       }
     return;
   }    

void 
G4HEPlot::Print( G4String  name, G4int iplot)
   {
     G4cout << name << iplot << " = new TH1F(" << '"' << iplot 
          << '"' << "," << '"' <<  '"' << "," << Nbin << "," 
          << Xstart << "," << Xstart+Nbin*Xbin << ");" << G4endl; 
     G4cout << "for(I=1; I<=" << Nbin <<"; I++) {" << G4endl;
     G4cout << "Y[I] = 0.;" << G4endl;
     G4cout << "}" << G4endl;
     for (G4int i=0; i<Nbin; i++)
       {
         G4cout << "Y[" << i+1 << "] = " << Yvalue[i] << ";" << G4endl;
       } 
     G4cout << "XA = " << Xstart << ";" << G4endl;
     G4cout << "STEP2 = " << Xbin/2. << ";" << G4endl;
     G4cout << "for(I=1; I<=" << Nbin << "; I++) {" << G4endl;
     G4cout << "XA = XA + STEP2;" << G4endl;
     G4cout << name << iplot << "->Fill(XA,Y[I]);" << G4endl;
     G4cout << "XA = XA + STEP2;" << G4endl;
     G4cout << "}" << G4endl;
     return;                         
   }

void
G4HEPlot::DumpToFile(G4int aPlot, G4String aName)
   {
     FILE* fz;
     fz = fopen(aName, "a");
     fprintf(fz,"%3d %3d %12.5e %12.5e \n",aPlot, Nbin, Xstart, Xbin);
     G4double xval = Xstart - Xbin/2.;
     for(G4int i=0; i<Nbin; i++)
       {
         xval += Xbin;
         fprintf(fz,"%12.5e %12.5e \n", xval, Yvalue[i]);
       }
     fclose(fz);
   }

void
G4HEPlot::GetFromFile(G4int aPlot, G4String aName)
   {
     FILE* fz;
     G4int ip, nb;
     G4double xs,xb,x,y;
     if((fz = fopen(aName, "r")) != NULL)
       {
         while(fscanf(fz,"%ld %ld %le %le", &ip, &nb, &xs, &xb)==4)
           {
             if(ip == aPlot) 
               {
                 if(Nbin != nb) Init(nb, xs, xb);
               }
             for (G4int i=0; i<nb; i++)
               {
                 fscanf(fz,"%le %le", &x, &y);
                 if(ip == aPlot) Fill(x,y);
               }
             if(ip == aPlot) break;
           }
         fclose(fz);
         return;
       }
     else
       {
         G4cout << " File " << aName << " not found " << G4endl;
       }
     return;
   }



