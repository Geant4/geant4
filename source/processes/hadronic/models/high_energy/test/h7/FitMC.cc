// Test- Job simulates the experimental data given in
// N.M.Agababyan et al., Z.Phys.C (1990)
// pi+, K+ on H, Al and Au at 250 GeV/c 
// The same comparison has been done with Geant3 MC
// H. Fesefeldt, July 1998.
// This program uses the files "particle""Material".plot
// as input and produces two files, "Fith2.plot" and
// and a file for direct input to the Root- program 
// on standard- output. 

#include "G4HEPlot.hh"
#include "globals.hh"

int main()
{
  G4int i,j;
  G4String aPIM110File = "PionMinus110.plot";
  G4String aPIP110File = "PionPlus110.plot";
  G4String aKM110File = "KaonMinus110.plot";
  G4String aKP110File = "KaonPlus110.plot";
  G4String aPIM58File = "PionMinus58.plot";
  G4String aKM58File = "KaonMinus58.plot";

  G4String aFitFile = "Fith7.plot";

  G4HEPlot* plot = new G4HEPlot[100];

  G4int iplot[] = {3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 32, 33, 34, 35,
                   56, 57, 58, 59, 60, 76, 77, 78, 79, 80};  
  plot[ 3].Init(25, 0., 0.1);
  plot[ 4].Init(25, 0., 0.1);
  plot[ 5].Init(25, 0., 0.1);
  plot[ 6].Init(25, 0., 0.1);
  plot[ 7].Init(25, 0., 0.1);
  plot[ 8].Init(25, 0., 0.1);
  plot[ 9].Init(25, 0., 0.1);
  plot[12].Init(10, 0., 0.1);
  plot[13].Init(25, 0., 0.1);
  plot[14].Init(25, 0., 0.1);
  plot[15].Init(25, 0., 0.1);
  plot[32].Init(10, 0., 0.1);
  plot[33].Init(25, 0., 0.1);
  plot[34].Init(25, 0., 0.1);
  plot[35].Init(25, 0., 0.1);
  plot[56].Init(20, 0., 0.05);
  plot[76].Init(20, 0., 0.05);
  plot[57].Init(20, 0., 0.05);
  plot[77].Init(20, 0., 0.05);
  plot[58].Init(20, 0., 0.05);
  plot[78].Init(20, 0., 0.05);
  plot[59].Init(20, 0., 0.05);
  plot[79].Init(20, 0., 0.05);
  plot[60].Init(20, 0., 0.05);
  plot[80].Init(20, 0., 0.05);

  for(i=1; i<=6; i++)
    {
      G4int ki = 1;
      if(i > 2) ki = 2;
      if(i > 4) ki=i-2;
      for(j=1; j<=20; j++) 
	{
          if(j<3) continue;
          if(j==10) continue;
          if(j==11) continue;
          if(ki == 1 && j > 15) continue;
          if(ki == 2 && (j < 12 || j > 15)) continue;
          if(ki >= 3 && j < 16) continue;
          G4int npl = (ki-1)*20+j;
          for(G4int l = 0; l<25; l++)
            {
              if(npl == iplot[l]) goto label;
            }
          cout << " wrong plot " << npl << G4endl;
          continue;
	label:
          if(i==1) plot[npl].GetFromFile(j, aPIM110File);      
          if(i==2) plot[npl].GetFromFile(j, aPIP110File);
          if(i==3) plot[npl].GetFromFile(j, aKM110File);
          if(i==4) plot[npl].GetFromFile(j, aKP110File);
          if(i==5) plot[npl].GetFromFile(j, aPIM58File);
          if(i==6) plot[npl].GetFromFile(j, aKM58File);
          if(ki <= 2) plot[npl].Scale(0.5, plot[npl]);
        }
    }
  for(i=3; i<=9; i++)
    {
      plot[i].DumpToFile(i, aFitFile);
      plot[i].Print("hist",100+i);
    }
  for(i=12; i<=15; i++)
    {
      plot[i].DumpToFile(i, aFitFile);
      plot[i].Print("hist",100+i);
    }
  for(i=32; i<=35; i++)
    {
      plot[i].DumpToFile(i, aFitFile);
      plot[i].Print("hist",100+i);
    }
  for(i=1; i<=5; i++)
    {
      plot[55+i].DumpToFile(55+i, aFitFile);
      plot[55+i].Print("hist",155+i);
      plot[75+i].DumpToFile(75+i, aFitFile);
      plot[75+i].Print("hist",175+i);
    }

  delete plot;
}


