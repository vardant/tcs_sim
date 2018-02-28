#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH1F.h"
#include <TStyle.h>
#include <TGaxis.h>
#include <TPaveLabel.h>
#include <iostream>
#include <sstream>   //for istringstream

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

// e- and e+ fluxes from target.

#define ENEG_IND  11
#define EPOS_IND -11

void ee_flux(string dir="ROOTfiles",
	     string phases="gamma11GeV", int maxfiles=3,
string comment="11 GeV photon beam on UVA target, e- & e+ fluxes from target") {

  const double emax=12000.;
  const int    nbin=1200;
  
  TH1D* hfeneg = new TH1D("hfeneg", "e- flux", nbin, 0., emax);
  hfeneg->GetXaxis()->SetTitle("E [MeV]");
  hfeneg->GetYaxis()->SetTitle(Form("count/s/%3.0f MeV",emax/nbin));

  TH1D *hfepos = (TH1D*)hfeneg->Clone("hfepos");
  hfepos->SetTitle("e+ flux");

  //Chain root files

  TChain ch("target");   //a chain to process Tree tree
  istringstream iss;
  iss.str(phases);
  string phase;
  while (iss >> phase) {
    cout << "phase: " << phase << endl;
    for (int ifile=0; ifile<maxfiles; ifile++) {
      string rootfile = dir + "/" + phase + "/" + Form("%d.root",ifile);
      if (FILE *file = fopen(rootfile.c_str(), "r")) {
        fclose(file);
	ch.Add(rootfile.c_str(),0);
	cout << "  Add root file: " << rootfile << endl;
      }
      else
	cout << "  File: " << rootfile << " does not exist" << endl;
    }
  }
		     
  std::vector<double> *edep_vec = 0;
  std::vector<int> *pid_vec = 0;

  TBranch *edep_br = 0;
  TBranch *pid_br = 0;

  ch.SetBranchAddress("edepcont",&edep_vec,&edep_br);
  ch.SetBranchAddress("pidcont",&pid_vec,&pid_br);

  Long64_t nentries = ch.GetEntriesFast();
  cout << "nentries = " << nentries << endl;

  for (Int_t i = 0; i < nentries; i++) {
      
    Long64_t tentry = ch.LoadTree(i);
    if (i < 10) cout << "tentry = " << tentry << "  i = " << i << endl;

    edep_br->GetEntry(tentry);
    pid_br->GetEntry(tentry);

    for (UInt_t j = 0; j < pid_vec->size(); ++j) {

      switch (pid_vec->at(j)) {
      case ENEG_IND : hfeneg->Fill(edep_vec->at(j),1.);
	break;
      case EPOS_IND : hfepos->Fill(edep_vec->at(j),1.);
	break;
      default: ;
      }

    }

  }

  // Scale histograms.

  const double CPSflux = 1.5e+12;      //photon/s
  hfeneg->Scale(CPSflux/double(nentries));
  hfepos->Scale(CPSflux/double(nentries));

  // Plot.

  //  gStyle->SetOptStat(0);

  TCanvas *ct = new TCanvas("ee_fluxes", "e-, e+ fluxes", 720, 600);

  TPaveLabel* title = new TPaveLabel(0.1,0.96,0.9,0.99,comment.c_str());
  title->SetBorderSize(0);
  title->SetFillColor(0);
  title->Draw();

  TPad* hPad = new TPad("histos","histos",0.01,0.05,0.95,0.95);
  hPad->Draw();
  hPad->cd();

  hPad->Divide(1,2);
  
  hPad->cd(1);
  gPad->SetLogy();
  hfeneg->Draw("hist");
  hPad->cd(2);
  gPad->SetLogy();
  hfepos->Draw("hist");

}
