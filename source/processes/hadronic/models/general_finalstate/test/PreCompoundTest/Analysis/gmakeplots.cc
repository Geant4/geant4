#include "TROOT.h"
#include "TApplication.h"

#include "mkplmakeplots.h"

using namespace std;

int main(int argc, char ** argv)
{
  TApplication theApp("App", &argc, argv);

  mkplmakeplots MainWindow(gClient->GetRoot(), 712, 950);

  theApp.Run();

  return 0;
}
