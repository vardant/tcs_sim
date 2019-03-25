#include <TH2.h>
#include <fstream>
#include <TCanvas.h>
#include <TStyle.h>

using namespace std;

void field_map() {

  TH2F* zmap = new TH2F("zmap","g2p target field, Z component [T]",
			50, -0.5, 49.5, 50, -0.5, 49.5);
			//300, -0.5, 299.5, 300, -0.5, 299.5);
  zmap->GetXaxis()->SetTitle("Z [cm]");
  zmap->GetYaxis()->SetTitle("R [cm]");

  TH2F* rmap = (TH2F*)zmap->Clone("rmap");
  rmap->SetTitle("g2p target field, R component [T]");

  string label;
  float z, r, Bz, Br, Btot;

  ifstream ifs;
  ifs.open("g2p_hallbfield.dat");
  ifs >> label >> label >> label >> label >> label;
  for (int ir=0; ir<300; ir++)
    for (int iz=0; iz<300; iz++) {
      ifs >> z >> r >> Bz >> Br >> Btot;
      zmap->Fill(z,r,Bz);
      rmap->Fill(z,r,Br);
    }
  ifs.close();

  gStyle->SetOptStat(0);
  new TCanvas("zmap");
  zmap->Draw("colz");
  new TCanvas("rmap");
  rmap->Draw("colz");
}
