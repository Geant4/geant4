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
  G4int i;
  G4String aPionHFile = "PionH.plot";
  G4String aPionAlFile = "PionAl.plot";
  G4String aPionAuFile = "PionAu.plot";
  G4String aKaonHFile = "KaonH.plot";
  G4String aKaonAlFile = "KaonAl.plot";
  G4String aKaonAuFile = "KaonAu.plot";
  G4String aFitFile = "Fith2.plot";

  G4HEPlot* plot = new G4HEPlot[62]; 

  plot[ 1].Init(21,-5.5, 0.5);
  plot[ 2].Init(21,-5.5, 0.5);
  plot[ 3].Init(21,-5.5, 0.5);
  plot[ 4].Init(64, 0.0, 0.05);
  plot[ 5].Init(64, 0.0, 0.05);
  plot[ 6].Init(21,-5.5, 0.5);
  plot[ 7].Init(21,-5.5, 0.5);
  plot[ 8].Init(21,-5.5, 0.5);
  plot[ 9].Init(64, 0.0, 0.05);
  plot[10].Init(64, 0.0, 0.05);
  plot[31].Init(22,-5.5, 0.5);
  plot[32].Init(22,-5.5, 0.5);
  plot[33].Init(22,-5.5, 0.5);
  plot[34].Init(22,-5.5, 0.5);
  plot[35].Init(22,-5.5, 0.5);
  plot[36].Init(22,-5.5, 0.5);
  plot[37].Init(100, 0., 0.05);
  plot[38].Init(100, 0., 0.05);
  plot[39].Init(100, 0., 0.05);
  plot[40].Init(100, 0., 0.05);
  plot[41].Init(100, 0., 0.05);
  plot[42].Init(100, 0., 0.05);
  plot[43].Init(32, -0.8, 0.05);
  plot[44].Init(32, -0.8, 0.05);
  plot[45].Init(32, -0.8, 0.05);
  plot[46].Init(20, -0.5, 1.0);
  plot[47].Init(20, -0.5, 1.0);
  plot[48].Init(20, -0.5, 1.0);
  plot[49].Init(40, 0.0, 0.05);
  plot[50].Init(40, 0.0, 0.05);
  plot[51].Init(40, 0.0, 0.05);
  plot[52].Init(40, 0.0, 0.05);
  plot[53].Init(10, -5.0, 0.5);
  plot[54].Init(10, -5.0, 0.5); 
  plot[55].Init(10, -5.0, 0.5);
  plot[56].Init(10, -5.0, 0.5);
  plot[57].Init(10, -5.0, 0.5);
  plot[58].Init(10, -5.0, 0.5);
  plot[59].Init(5, 0.5, 1.0);
  plot[60].Init(5, 0.5, 1.0);
  plot[61].Init(5, 0.5, 1.0);

  plot[ 1].GetFromFile( 5, aPionAlFile);
  plot[ 1].GetFromFile( 5, aKaonAlFile);
  plot[ 2].GetFromFile( 6, aPionAlFile);
  plot[ 2].GetFromFile( 6, aKaonAlFile);
  plot[ 3].GetFromFile( 7, aPionAlFile);
  plot[ 3].GetFromFile( 7, aKaonAlFile);
  plot[ 4].GetFromFile( 9, aPionAlFile);
  plot[ 4].GetFromFile( 9, aKaonAlFile);
  plot[ 5].GetFromFile(10, aPionAlFile);
  plot[ 5].GetFromFile(10, aKaonAlFile);
  plot[ 6].GetFromFile( 5, aPionAuFile);
  plot[ 6].GetFromFile( 5, aKaonAuFile);
  plot[ 7].GetFromFile( 6, aPionAuFile);
  plot[ 7].GetFromFile( 6, aKaonAuFile);
  plot[ 8].GetFromFile( 7, aPionAuFile);
  plot[ 8].GetFromFile( 7, aKaonAuFile);
  plot[ 9].GetFromFile( 9, aPionAuFile);
  plot[ 9].GetFromFile( 9, aKaonAuFile);
  plot[10].GetFromFile(10, aPionAuFile);
  plot[10].GetFromFile(10, aKaonAuFile);

  for(i=1; i<6; i++)
    {
      plot[i].Scale(0.5, plot[i]);
    }
  for(i=6; i<11; i++)
    {
      plot[i].Scale(5., plot[i]);
    }
  
  plot[43].GetFromFile(13, aPionHFile);
  plot[43].GetFromFile(13, aKaonHFile);
  plot[46].GetFromFile(16, aPionHFile);
  plot[46].GetFromFile(16, aKaonHFile);
  plot[53].GetFromFile(19, aPionHFile);
  plot[53].GetFromFile(19, aKaonHFile);
  plot[56].GetFromFile(22, aPionHFile);
  plot[56].GetFromFile(22, aKaonHFile);
  plot[59].GetFromFile(33, aPionHFile);
  plot[59].GetFromFile(33, aKaonHFile);

  plot[31].GetFromFile(23, aPionAlFile);
  plot[31].GetFromFile(23, aKaonAlFile);  
  plot[32].GetFromFile(23, aPionAuFile);
  plot[33].GetFromFile(23, aKaonAuFile);  
  plot[33].GetFromFile(24, aPionAlFile);
  plot[33].GetFromFile(24, aKaonAlFile);  
  plot[34].GetFromFile(24, aPionAuFile);
  plot[34].GetFromFile(24, aKaonAuFile);  
  plot[35].GetFromFile(25, aPionAlFile);
  plot[35].GetFromFile(25, aKaonAlFile);  
  plot[36].GetFromFile(25, aPionAuFile);
  plot[36].GetFromFile(25, aKaonAuFile);  
  plot[37].GetFromFile(26, aPionAlFile);
  plot[37].GetFromFile(26, aKaonAlFile);  
  plot[38].GetFromFile(26, aPionAuFile);
  plot[38].GetFromFile(26, aKaonAuFile);  
  plot[39].GetFromFile(27, aPionAlFile);
  plot[39].GetFromFile(27, aKaonAlFile);  
  plot[40].GetFromFile(27, aPionAuFile);
  plot[40].GetFromFile(27, aKaonAuFile);  
  plot[41].GetFromFile(28, aPionAlFile);
  plot[41].GetFromFile(28, aKaonAlFile);  
  plot[42].GetFromFile(28, aPionAuFile);
  plot[42].GetFromFile(28, aKaonAuFile);  
  plot[44].GetFromFile(13, aPionAlFile);
  plot[44].GetFromFile(13, aKaonAlFile);  
  plot[45].GetFromFile(13, aPionAuFile);
  plot[45].GetFromFile(13, aKaonAuFile);  
  plot[47].GetFromFile(16, aPionAlFile);
  plot[47].GetFromFile(16, aKaonAlFile);  
  plot[48].GetFromFile(16, aPionAuFile);
  plot[48].GetFromFile(16, aKaonAuFile);  
  plot[49].GetFromFile(29, aPionAlFile);
  plot[49].GetFromFile(29, aKaonAlFile);  
  plot[50].GetFromFile(29, aPionAuFile);
  plot[50].GetFromFile(29, aKaonAuFile);  
  plot[51].GetFromFile(30, aPionAlFile);
  plot[51].GetFromFile(30, aKaonAlFile);  
  plot[52].GetFromFile(30, aPionAuFile);
  plot[52].GetFromFile(30, aKaonAuFile);  
  plot[54].GetFromFile(19, aPionAlFile);
  plot[54].GetFromFile(19, aKaonAlFile);  
  plot[55].GetFromFile(19, aPionAuFile);
  plot[55].GetFromFile(19, aKaonAuFile);  
  plot[57].GetFromFile(22, aPionAlFile);
  plot[57].GetFromFile(22, aKaonAlFile);  
  plot[58].GetFromFile(22, aPionAuFile);
  plot[58].GetFromFile(22, aKaonAuFile);  
  plot[60].GetFromFile(33, aPionAlFile);
  plot[60].GetFromFile(33, aKaonAlFile);  
  plot[61].GetFromFile(33, aPionAuFile);
  plot[61].GetFromFile(33, aKaonAuFile);


  for(i=31; i<=61; i++)
    {
      plot[i].Scale(0.5, plot[i]);
    }

  plot[46].XScale( 1., -0.5, plot[46]);
  plot[47].XScale( 1., -0.5, plot[47]);
  plot[48].XScale( 1., -0.5, plot[48]);
  plot[59].XScale( 1., -0.5, plot[59]);
  plot[60].XScale( 1., -0.5, plot[60]);
  plot[61].XScale( 1., -0.5, plot[61]);

  for(i=1; i<11; i++)
    {
      plot[i].DumpToFile( i, aFitFile);
      plot[i].Print("hist", 100+i);
    }
  for(i=31; i<62; i++)
    {
      plot[i].DumpToFile( i, aFitFile);
      plot[i].Print("hist", 100+i);
    }
    delete plot;
}


