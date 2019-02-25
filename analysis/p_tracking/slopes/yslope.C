#include <TGraphErrors.h>
#include <fstream>
#include <string>
#include <iostream>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>

using namespace std;

//Fit yslope versus yslope at vertex.

void yslope(string datfile="yslope.dat") {

  vector<double> p;
  vector<double> perr;
  vector<double> pinv;
  vector<double> pinverr;
  vector<double> yslope;
  vector<double> yslerr;
  vector<double> yoffset;
  vector<double> yofferr;

  vector<double> yslope_inv;
  vector<double> yslinv_err;

  ifstream fin;
  fin.open(datfile.c_str(),ios::in);

  double P, sl, sler, of, ofer;

  while (fin >> P >> of >> ofer >> sl >> sler) {
    p.push_back(P);
    perr.push_back(0.);
    pinv.push_back(1./P);
    pinverr.push_back(0.);
    yslope.push_back(sl);
    yslerr.push_back(sler);
    yoffset.push_back(of);
    yofferr.push_back(ofer);

    yslope_inv.push_back(1./sl-1.);
    yslinv_err.push_back(sler/sl/sl);
  }

  fin.close();

  gStyle->SetOptFit();

  TGraphErrors* gslope =
    new TGraphErrors(p.size(),&p[0],&yslope[0],&perr[0],&yslerr[0]);
  gslope->SetTitle("slope versus P");
  gslope->Fit("pol0");
  gslope->SetMarkerStyle(20);
  new TCanvas("slope");
  gslope->GetXaxis()->SetTitle("P [GeV/c]");
  gslope->GetYaxis()->SetTitle("slope");
  gslope->Draw("AP");

  TGraphErrors* goffset =
    new TGraphErrors(pinv.size(),&pinv[0],&yoffset[0],&pinverr[0],&yofferr[0]);
  goffset->SetTitle("offset versus 1/P");
  goffset->Fit("pol1");
  goffset->SetMarkerStyle(20);
  new TCanvas("offset");
  goffset->GetXaxis()->SetTitle("1/P [1/(GeV/c)]");
  goffset->GetYaxis()->SetTitle("offset");
  goffset->Draw("AP");
  /*
  TGraphErrors* gslope_vs_pinv =
    new TGraphErrors(pinv.size(),&pinv[0],&yslope[0],&perr[0],&yslerr[0]);
  gslope_vs_pinv->SetTitle("slope versus Pinv");
  gslope_vs_pinv->Fit("pol1");
  gslope_vs_pinv->SetMarkerStyle(20);
  new TCanvas("slope_vs_pinv");
  gslope_vs_pinv->GetXaxis()->SetTitle("1/P [c/GeV]");
  gslope_vs_pinv->GetYaxis()->SetTitle("slope");
  gslope_vs_pinv->Draw("AP");
  */
  TGraphErrors* gslope_inv =
    new TGraphErrors(p.size(),&p[0],&yslope_inv[0],&perr[0],&yslinv_err[0]);
  gslope_inv->SetTitle("1/slope-1 versus P");
  gslope_inv->Fit("expo");
  gslope_inv->SetMarkerStyle(20);
  new TCanvas("slope_inv_vs_p");
  gslope_inv->GetXaxis()->SetTitle("P [GeV/c]");
  gslope_inv->GetYaxis()->SetTitle("1/slope - 1");
  gslope_inv->Draw("AP");

}
