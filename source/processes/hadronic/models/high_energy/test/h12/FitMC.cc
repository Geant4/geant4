#include "G4HEPlot.hh"
#include "globals.hh"

int main()
{
  G4int i;
  G4String aProton12File = "Proton12.plot";
  G4String aProton24File = "Proton24.plot";
  G4String aProton250File = "Proton250.plot";
  G4String aPion250File = "Pion250.plot";
  G4String aKaon250File = "Kaon250.plot";

  G4String aFitFile = "Fith12.plot";

  G4HEPlot* plot = new G4HEPlot[100]; 

  plot[ 1].Init(21,-5.5, 0.5);


  plot[ 1].GetFromFile( 1, aProton12File);
  plot[ 2].GetFromFile( 2, aProton12File);
  plot[ 3].GetFromFile( 3, aProton12File);
  plot[ 4].GetFromFile( 4, aProton12File);
  plot[ 5].GetFromFile( 5, aProton12File);
  plot[ 6].GetFromFile( 6, aProton12File);
  plot[ 7].GetFromFile(23, aProton12File);
  plot[ 9].GetFromFile(24, aProton12File);
  plot[ 8].GetFromFile(31, aProton12File);
  plot[10].GetFromFile(32, aProton12File);
  plot[21].GetFromFile( 1, aProton24File);
  plot[22].GetFromFile( 2, aProton24File);
  plot[23].GetFromFile( 3, aProton24File);
  plot[24].GetFromFile( 4, aProton24File);
  plot[25].GetFromFile( 5, aProton24File);
  plot[26].GetFromFile( 6, aProton24File);
  plot[27].GetFromFile(23, aProton24File);
  plot[29].GetFromFile(24, aProton24File);
  plot[28].GetFromFile(31, aProton24File);
  plot[30].GetFromFile(32, aProton24File);
  plot[41].GetFromFile( 1, aProton250File);
  plot[42].GetFromFile( 2, aProton250File);
  plot[43].GetFromFile( 3, aProton250File);
  plot[44].GetFromFile( 4, aProton250File);
  plot[45].GetFromFile( 5, aProton250File);
  plot[46].GetFromFile( 6, aProton250File);
  plot[71].GetFromFile(17, aPion250File);
  plot[72].GetFromFile(18, aPion250File);
  plot[73].GetFromFile(19, aPion250File);
  plot[74].GetFromFile(20, aPion250File);
  plot[75].GetFromFile(21, aPion250File);
  plot[65].GetFromFile( 5, aPion250File);
  plot[91].GetFromFile(17, aKaon250File);
  plot[92].GetFromFile(18, aKaon250File);
  plot[93].GetFromFile(19, aKaon250File);
  plot[94].GetFromFile(20, aKaon250File);
  plot[95].GetFromFile(21, aKaon250File);
  plot[85].GetFromFile( 5, aKaon250File);

  for(G4int l=1; l<=3; l++)
    {
      for(G4int i=1; i<=6; i++)
        {
          G4int npl = (l-1)*20+i; 
          plot[npl].DumpToFile(npl , aFitFile);
          plot[npl].Print("hist", 100+npl);
        } 
    }
  plot[ 7].DumpToFile( 7, aFitFile);
  plot[27].DumpToFile(27, aFitFile);
  plot[ 8].DumpToFile( 8, aFitFile);
  plot[28].DumpToFile(28, aFitFile);
  plot[ 9].DumpToFile( 9, aFitFile);
  plot[29].DumpToFile(29, aFitFile);
  plot[10].DumpToFile(10, aFitFile);
  plot[30].DumpToFile(30, aFitFile);
  plot[65].DumpToFile(65, aFitFile);
  plot[85].DumpToFile(85, aFitFile);
  plot[ 7].Print("hist", 107);
  plot[27].Print("hist", 127);
  plot[ 8].Print("hist", 108);
  plot[28].Print("hist", 128);
  plot[ 9].Print("hist", 109);
  plot[29].Print("hist", 129);
  plot[10].Print("hist", 110);
  plot[30].Print("hist", 130);
  plot[65].Print("hist", 165);
  plot[85].Print("hist", 185);

  for(l =4; l<=5; l++)
    {
      for(G4int i=11; i<=15; i++)
	{
          G4int npl=(l-1)*20+i;
          plot[npl].DumpToFile(npl, aFitFile);
          plot[npl].Print("hist", 100+npl);
        }
    }  

    delete plot;
}


