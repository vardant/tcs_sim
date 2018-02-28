#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TStyle.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <iostream>
#include <sstream>   //for istringstream

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

void totedeps(string tree="calo", string dir="ROOTfiles",
	      string phases="gamma11GeV.tar_mag_lhe_ln2_wins",
	      int maxfiles=4, double emax=12000.,
      string comment="11 GeV photon beam on UVA target, Edep in calorimeters") {

  TH1D* hEtotPos = new TH1D("hetotpos", ("Top "+ tree).c_str(), 120, 0., emax);
  hEtotPos->GetXaxis()->SetTitle("Edep [MeV]");

  TH1D *hEtotNeg = (TH1D*)hEtotPos->Clone("hetotneg");
  hEtotNeg->SetTitle(("Bottom " + tree).c_str());

  TH1D *hEtot = (TH1D*)hEtotPos->Clone("hetot");
  hEtot->SetTitle(("Top " + tree + " Bottom " + tree).c_str());
  
  TH2D* hEtop_vs_Ebot = new TH2D("hetop_vs_ebot",
				 ("Top " + tree + " vs Bootom " + tree).c_str(),
				 120, 0., emax, 120, 0., emax);
  hEtop_vs_Ebot->GetXaxis()->SetTitle("Edep(bottom) [MeV]");
  hEtop_vs_Ebot->GetYaxis()->SetTitle("Edep(top) [MeV]");
  
  //Chain root files

  TChain ch(tree.c_str());   //a chain to process Tree tree
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

  std::vector<int> *det_vec = 0;
  //  std::vector<unsigned int> *chan_vec = 0;
  std::vector<double> *edep_vec = 0;
  //  std::vector<int> *pid_vec = 0;

  TBranch *det_br = 0;
  //  TBranch *chan_br = 0;
  TBranch *edep_br = 0;
  //  TBranch *pid_br = 0;

  ch.SetBranchAddress("detcont",&det_vec,&det_br);
  //  ch.SetBranchAddress("chancont",&chan_vec,&chan_br);
  ch.SetBranchAddress("edepcont",&edep_vec,&edep_br);
  //  ch.SetBranchAddress("pidcont",&pid_vec,&pid_br);

  Long64_t nentries = ch.GetEntriesFast();
  cout << "nentries = " << nentries << endl;

  for (Int_t i = 0; i < nentries; i++) {
      
    Long64_t tentry = ch.LoadTree(i);
    if (i < 10) cout << "tentry = " << tentry << "  i = " << i << endl;

    det_br->GetEntry(tentry);
    //    chan_br->GetEntry(tentry);
    edep_br->GetEntry(tentry);
    //    pid_br->GetEntry(tentry);

    double etotpos = 0.;
    double etotneg = 0.;

    for (UInt_t j = 0; j < det_vec->size(); ++j) {
 
      int det = det_vec->at(j);
      //      unsigned int chan = chan_vec->at(j);
      double edep = edep_vec->at(j);

      //cout << " hit " << j << ": " << det << " " << chan << " "
      //	     << edep << endl;

      if (det > 0) {
	etotpos += edep;
      }
      else {
	etotneg += edep;
      }

    }

    hEtotPos->Fill(etotpos,1.);
    hEtotNeg->Fill(etotneg,1.);
    hEtot->Fill(etotpos+etotneg,1.);
    if(etotpos>0. && etotneg>0.) hEtop_vs_Ebot->Fill(etotpos,etotneg,1.);

    //      if (det_vec->size() !=0) getchar();
  }

  // Since we passed the address of a local variable we need
  // to remove it.
  ch.ResetBranchAddresses();

  // Scale histograms.

  ////  const double CPSflux = 1.5e+12;      //photon/s
  ////  hEtotPos->Scale(CPSflux/double(nentries));
  ////  hEtotNeg->Scale(CPSflux/double(nentries));
  ////  hEtot->Scale(CPSflux/double(nentries));

   // Plot.

  gStyle->SetOptStat(0);

  TCanvas *ct = new TCanvas(("Edeps_" + tree).c_str(), "Total Edep", 900, 600);

  TPaveLabel* title = new TPaveLabel(0.1,0.96,0.9,0.99,comment.c_str());
  title->SetBorderSize(0);
  title->SetFillColor(0);
  title->Draw();

  TPad* hPad = new TPad("histos","histos",0.01,0.05,0.95,0.95);
  hPad->Draw();
  hPad->cd();

  hPad->Divide(2,2);
  
  hPad->cd(1);
  gPad->SetLogy();
  hEtotPos->Draw("hist");
  hPad->cd(2);
  gPad->SetLogy();
  hEtotNeg->Draw("hist");
  hPad->cd(3);
  gPad->SetLogy();
  hEtot->Draw("hist");
  hPad->cd(4);
  hEtop_vs_Ebot->Draw("cont");
  
  //  TCanvas *c2 = new TCanvas("EtopEbottom",
  //			    ("Etop vs Ebottom in "+tree).c_str(), 600, 600);
  //  c2->cd(1);
  //  hEtop_vs_Ebot->Draw("cont");

}
