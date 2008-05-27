// new file (created 20 Apr 2008), is in use in all cases of E, PlotLog.C files: 800, 1200, 1500, 1600, 3000 MeV

{
Double_t x1[3000] = {.0};
Double_t y1[3000] = {.0};
Double_t errx1[3000] = {.0};
Double_t erry1[3000] = {.0};

string file30 = "test30.txt";                
ifstream  data_file(file30.c_str(),  ios::in );
if(data_file.fail()) cout << "Error opening data file " << file30 << endl;
else                 cout << "data file " << file30 << " is opened" << endl;

TString  linetype;
Double_t angle, xcoor, dxcoor;

for (i = 0; i < npl; i++) {
  data_file >> linetype >> bin >> angle;
  if(linetype == "END") break;
  cout << bin << " " << angle << endl;
  for (j=0; j<bin; j++) {
    data_file >> xcoor >> dxcoor >> y1[j] >> erry1[j];
    x1[j]=log10(xcoor);
    errx1[j]=dxcoor/1000;
    cout  << xcoor << " " << dxcoor << " " << x1[j] << " " << y1[j] << endl;
  }
  cout << i << " " << bin << endl;
  cout << "***" << endl;
  cout << "***" << endl;
  gr[i] = new TGraphErrors(bin, x1, y1, errx1, erry1);
}
 cout << "Add Data is complete .." << endl;

}
