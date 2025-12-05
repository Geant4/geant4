#include <map>
#include <string>
#include <vector>

void processMicrotrack(const char* filename = "microtrack.root")
{
  TFile* f = new TFile(filename, "UPDATE");
  if (!f || f->IsZombie()) {
    Error("processMicrotrack", "Cannot open file %s", filename);
    return;
  }

  // Map: histogram name -> (x axis, y axis)
  std::map<std::string, std::pair<const char*, const char*>> h1Axis = {
    {"fe", {"#varepsilon_{1} (keV)", "f(#varepsilon_{1})"}},
    {"efe", {"#varepsilon_{1} (keV)", "#varepsilon_{1} f(#varepsilon_{1})"}},
    {"e2fe", {"#varepsilon_{1} (keV)", "#varepsilon^{2}_{1} f(#varepsilon_{1})"}},
    {"fy", {"y (keV/#mum)", "f(y)"}},
    {"yfy", {"y (keV/#mum)", "y f(y)"}},
    {"y2fy", {"y (keV/#mum)", "y^{2} f(y)"}},
    {"fz", {"z_{1} (Gy)", "f(z_{1})"}},
    {"zfz", {"z_{1} (Gy)", "z_{1} f(z_{1})"}},
    {"z2fz", {"z_{1} (Gy)", "z^{2}_{1} f(z_{1})"}},
    {"Nsel", {"N_{sel}", "Counts"}},
    {"Nsite", {"N_{site}", "Counts"}},
    {"Nint", {"N_{int}", "Counts"}},
    {"KinE_in", {"T_{in} (MeV)", "Counts"}},
    {"KinE_out", {"T_{out} (MeV)", "Counts"}},
  };

  // 1) Set axis titles (if histograms already exist)
  for (const auto& kv : h1Axis) {
    if (TH1* h = (TH1*)f->Get(kv.first.c_str())) {
      h->GetXaxis()->SetTitle(kv.second.first);
      h->GetYaxis()->SetTitle(kv.second.second);
      h->Write("", TObject::kOverwrite);  // persist updates
    }
  }

  f->Close();
}

