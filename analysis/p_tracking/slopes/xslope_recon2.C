#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include <TMath.h>
#include <TDecompLU.h>
#include <iostream>
#include <fstream>
#include <TText.h>
#include <TProfile.h>
#include <TStyle.h>

#include <TVector3.h>

#include <TRandom.h>

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

#define TRACKER1_DIST 1200.   //mm, consistent with G4 coding.
#define TRACKER2_DIST 1300.
#define TRACKER3_DIST 1400.

#define TILT_ANGLE 13.835     //deg, consistent with G4 coding.
#define ROT_ANGLE  10.034

#define GEM_ACCURACY 100.E-3   //mm, conservative estimate.

struct hit {
  double x;
  double y;
  int det;
  int layer;
};

double FitStraightTrack(const vector<hit> &hitlist, const int axis,
			double &offset, double &slope);

void RotateToLab(TVector3 &v, const int ndet);

double xslope_reconstruct(const double &xslope_det, const int &quarter);

TRandom RandomGen;

double GEMSmear(const double xtrack);

//Reconstruction of track's slope in XZ plane at vertex, by use of
//track's slope from GEMs.

void xslope_recon2(string rootfile="RootFiles/p_test.root", int quarter=0,
		   string comment="test") {

  //Chain root files

  TChain chain_tracker("tracker");   //a chain to process Tree "tracker"
  TChain chain_beam("beam");
  if (FILE *file = fopen(rootfile.c_str(), "r")) {
    fclose(file);
    chain_tracker.Add(rootfile.c_str(),0);
    chain_beam.Add(rootfile.c_str(),0);
    cout << "Root file: " << rootfile << endl;
  }
  else {
    cout << "  File: " << rootfile << " does not exist" << endl;
    return;
  }

  Long64_t nentries = chain_tracker.GetEntriesFast();
  cout << "nentries = " << nentries << endl;

  // tracker branch.

  // Declaration of leaf types
  vector<double>  *xcont = 0;
  vector<double>  *ycont = 0;
  vector<double>  *edepcont = 0;
  vector<double>  *lengthcont = 0;
  vector<int>     *detcont = 0;
  vector<int>     *layercont = 0;
  vector<int>     *pidcont = 0;
  vector<int>     *pidorigcont = 0;
  vector<int>     *trackidcont = 0;
  vector<int>     *nstepcont = 0;

  // List of branches
  TBranch        *b_xcont = 0;   //!
  TBranch        *b_ycont = 0;   //!
  TBranch        *b_edepcont = 0;   //!
  TBranch        *b_lengthcont = 0;   //!
  TBranch        *b_detcont = 0;   //!
  TBranch        *b_layercont = 0;   //!
  TBranch        *b_pidcont = 0;   //!
  TBranch        *b_pidorigcont = 0;   //!
  TBranch        *b_trackidcont = 0;   //!
  TBranch        *b_nstepcont = 0;   //!

  chain_tracker.SetBranchAddress("xcont", &xcont, &b_xcont);
  chain_tracker.SetBranchAddress("ycont", &ycont, &b_ycont);
  chain_tracker.SetBranchAddress("edepcont", &edepcont, &b_edepcont);
  chain_tracker.SetBranchAddress("lengthcont", &lengthcont, &b_lengthcont);
  chain_tracker.SetBranchAddress("detcont", &detcont, &b_detcont);
  chain_tracker.SetBranchAddress("layercont", &layercont, &b_layercont);
  chain_tracker.SetBranchAddress("pidcont", &pidcont, &b_pidcont);
  chain_tracker.SetBranchAddress("pidorigcont", &pidorigcont, &b_pidorigcont);
  chain_tracker.SetBranchAddress("trackidcont", &trackidcont, &b_trackidcont);
  chain_tracker.SetBranchAddress("nstepcont", &nstepcont, &b_nstepcont);

  //Beam branch.

  Double_t        px;
  Double_t        py;
  Double_t        pz;
  TBranch        *b_px;
  TBranch        *b_py;
  TBranch        *b_pz;
  chain_beam.SetBranchAddress("px", &px, &b_px);
  chain_beam.SetBranchAddress("py", &py, &b_py);
  chain_beam.SetBranchAddress("pz", &pz, &b_pz);

  uint max_hitlist_size = 0;

  TH1F* h_dslope = new TH1F("h_dslope", "#delta(slope)", 100, -0.05, 0.05);

  TH2F* h_dslope_vs_slope;
  TH2F* h_slope_det_vs_slope;
  if (quarter==0 || quarter==3) {
    h_dslope_vs_slope = new TH2F("h_dslope_vs_slope",
	     "#delta(slope) vs slope", 100, 0., 0.40, 100, -0.05, 0.05);
    h_slope_det_vs_slope = new TH2F("h_slope_vs_slope",
	     "slope_det vs slope", 100, 0., 0.40, 100, 0., 0.40);
  }
  else {
    h_dslope_vs_slope = new TH2F("h_dslope_vs_slope",
	     "#delta(slope) vs slope", 100, -0.40, 0., 100, -0.05, 0.05);
    h_slope_det_vs_slope = new TH2F("h_slope_vs_slope",
	     "slope_det vs slope", 100, -0.40, 0., 100, -0.40, 0.);
  }

  TH2F* h_dslope_vs_p = new TH2F("h_dslope_vs_p",
	     "#delta(slope) vs p", 100, 0., 2., 100, -0.05, 0.05);

  ////  RandomGen.SetSeed(0);   //seed from comp. clock.

  for (int ientry=0; ientry<nentries; ientry++) {

    Long64_t tentry = chain_tracker.LoadTree(ientry);
    if (ientry < 10)
      cout << "tentry = " << tentry << "  ientry = " << ientry << endl;

    b_xcont->GetEntry(tentry);
    b_ycont->GetEntry(tentry);
    b_detcont->GetEntry(tentry);
    b_layercont->GetEntry(tentry);
    b_pidcont->GetEntry(tentry);
    b_pidorigcont->GetEntry(tentry);
    b_trackidcont->GetEntry(tentry);

    Long64_t bentry = chain_beam.LoadTree(ientry);
    b_px->GetEntry(bentry);
    b_py->GetEntry(bentry);
    b_pz->GetEntry(bentry);

    vector<hit> hitlist;

    //    cout << "Filling hitlist:" << endl;
    //    cout << " detcont size = " << detcont->size() << endl;
    //    cout << " pidcont size = " << pidcont->size() << endl;
    //    cout << " pidorigcont size = " << pidorigcont->size() << endl;
    //    cout << " trackidcont size = " << trackidcont->size() << endl;

    for (UInt_t j = 0; j < detcont->size(); ++j) {
      if (pidcont->at(j) == pidorigcont->at(j) &&
      	  trackidcont->at(j) == 1 && detcont->at(j) == quarter) {
	hit goodHit = {GEMSmear(xcont->at(j)), GEMSmear(ycont->at(j)),
		       detcont->at(j), layercont->at(j)};
	hitlist.push_back(goodHit);
      }
    }

    //    cout << " hitlist size = " << hitlist.size() << endl;

    bool good_hitlist = false;

    for (uint j=0; j<hitlist.size(); j++) {
      if (hitlist.at(j).layer == 0) {
	good_hitlist = true;
	break;
      }
    }

    if (good_hitlist) {

      for (uint j=0; j<hitlist.size(); j++) {
	if (!good_hitlist) break;
	for (uint k=j+1; k<hitlist.size(); k++)
	  if (hitlist.at(k).layer == hitlist.at(j).layer) {
	    good_hitlist = false;
	    break;
	  }
      }

    }

    good_hitlist = good_hitlist && hitlist.size() > 1;

    if (good_hitlist) {

      if (hitlist.size() > max_hitlist_size)
	max_hitlist_size = hitlist.size();

      double xposition = 999999.;
      double yposition = 999999.;
      for (uint j=0; j<hitlist.size(); j++) {
	if (hitlist.at(j).layer == 0) {
	  xposition = hitlist.at(j).x;
	  yposition = hitlist.at(j).y;
	  break;
	}
      }

      TVector3 vxyz(xposition, yposition, TRACKER1_DIST);

      RotateToLab(vxyz, quarter);

      double xslope_det = atan(vxyz.X()/vxyz.Z());

      double xslope_recon = xslope_reconstruct(xslope_det, quarter);

      //      cout << "slope0      = " << px/pz << endl;
      //      cout << "slope_recon = " << xslope_recon << endl;
      //      cout << "P      = " << sqrt(px*px+py*py+pz*pz)/1000. << endl;
      //      getchar();

      h_dslope->Fill(xslope_recon-atan(px/pz));
      h_dslope_vs_slope->Fill(atan(px/pz), xslope_recon-atan(px/pz));
      h_dslope_vs_p->Fill(sqrt(px*px+py*py+pz*pz)/1000.,
			  xslope_recon-atan(px/pz));
      h_slope_det_vs_slope->Fill(atan(px/pz), xslope_det);
    }

  }  //over entries

  TCanvas* cds = new TCanvas("dslope");
  h_dslope->SetTitle("");
  h_dslope->GetXaxis()->SetTitle("#delta#phi [rad]");
  h_dslope->GetYaxis()->SetTitle("count");
  h_dslope->Draw();
  TText *tds = new TText(0.25,0.80,Form("Quarter %d",quarter+1));
  tds->SetNDC();
  tds->SetTextFont(43);
  tds->SetTextAlign(22);
  tds->SetTextSize(18);
  tds->Draw();
  cds->SaveAs(Form("dslopeX_%d.gif",quarter+1));
  cds->SaveAs(Form("dslopeX_%d.pdf",quarter+1));

  TCanvas* cdss = new TCanvas("dslope_vs_slope");
  h_dslope_vs_slope->SetTitle("");
  h_dslope_vs_slope->SetStats(0);
  h_dslope_vs_slope->GetXaxis()->SetTitle("#phi [rad]");
  h_dslope_vs_slope->GetYaxis()->SetTitle("#delta#phi [rad]");
  h_dslope_vs_slope->Draw();
  TText *tdss = new TText(0.75,0.80,Form("Quarter %d",quarter+1));
  tdss->SetNDC();
  tdss->SetTextFont(43);
  tdss->SetTextAlign(22);
  tdss->SetTextSize(18);
  tdss->Draw();
  cdss->SaveAs(Form("dslopeX_vs_slopeX_%d.gif",quarter+1));
  cdss->SaveAs(Form("dslopeX_vs_slopeX_%d.pdf",quarter+1));

  TCanvas* cdsp = new TCanvas("dslope_vs_p");
  h_dslope_vs_p->SetTitle("");
  h_dslope_vs_p->SetStats(0);
  h_dslope_vs_p->GetXaxis()->SetTitle("P [GeV/c]");
  h_dslope_vs_p->GetYaxis()->SetTitle("#delta#phi [rad]");
  h_dslope_vs_p->Draw();
  tdss->Draw();
  cdsp->SaveAs(Form("dslopeX_vs_p_%d.gif",quarter+1));
  cdsp->SaveAs(Form("dslopeX_vs_p_%d.pdf",quarter+1));

  TCanvas* css = new TCanvas("slope_det_vs_slope");
  h_slope_det_vs_slope->SetTitle("");
  h_slope_det_vs_slope->SetStats(0);
  h_slope_det_vs_slope->GetXaxis()->SetTitle("#phi [rad]");
  h_slope_det_vs_slope->GetYaxis()->SetTitle("#phi(det) [rad]");
  h_slope_det_vs_slope->Draw();
  tdss->Draw();
  TH1D* h_slope_det_vs_slope_prof = h_slope_det_vs_slope->ProfileX();
  gStyle->SetOptFit();
  gStyle->SetStatY(0.55);
  double xlo = 0.03;
  double xhi = 0.30;
  switch (quarter) {
  case 0:
  case 3:
    h_slope_det_vs_slope_prof->Fit("pol1","","sames",xlo,xhi);
    break;
  case 1:
  case 2:
    h_slope_det_vs_slope_prof->Fit("pol1","","sames",-xhi,-xlo);
    break;
  default:
    cout << "*** Wrong fit limits! ***" << endl;
  }
  css->SaveAs(Form("slopeX_det_vs_slopeX_%d.gif",quarter+1));
  css->SaveAs(Form("slopeX_det_vs_slopeX_%d.pdf",quarter+1));

  cout << "Max hitlist size = " << max_hitlist_size << endl;
}

//------------------------------------------------------------------------------

double FitStraightTrack(const vector<hit> &hitlist, const int axis,
			double &offset, double &slope) {

  const double zpos[] {TRACKER1_DIST, TRACKER2_DIST, TRACKER3_DIST};

  double sum_zx = 0.;
  double sum_x = 0.;
  double sum_z = 0.;
  double sum_zz = 0.;

  uint n = hitlist.size();
  for (uint i=0; i<n; i++) {
    double z = zpos[hitlist.at(i).layer];
    double x = (axis==0 ? hitlist.at(i).x : hitlist.at(i).y);
    sum_zx += z*x;
    sum_x += x;
    sum_z += z;
    sum_zz += z*z;
  }

  slope = (n*sum_zx - sum_x*sum_z)/(n*sum_zz - sum_z*sum_z);
  offset = (sum_x - sum_z*slope)/n;

  double x2 = 0.;
  for (uint i=0; i<n; i++) {
    double z = zpos[hitlist.at(i).layer];
    double x = (axis==0 ? hitlist.at(i).x : hitlist.at(i).y);
    x2 += (x-slope*z-offset)*(x-slope*z-offset)/(GEM_ACCURACY*GEM_ACCURACY);
  };
  x2 /= 2.;

  return x2;
}

//------------------------------------------------------------------------------

void RotateToLab(TVector3 &v, const int ndet) {

  double xangle = 0.;
  double yangle = 0.;
  //  int xflip = 0;
  //  int yflip = 0;

  switch (ndet) {
    case 0:
      xangle = -TILT_ANGLE;
      yangle =  ROT_ANGLE;
      //      xflip = +1;
      //      yflip = +1;
      break;
    case 1:
      xangle = -TILT_ANGLE;
      yangle = -ROT_ANGLE;
      //      xflip = -1;
      //      yflip = +1;
      break;
    case 2:
      xangle =  TILT_ANGLE;
      yangle = -ROT_ANGLE;
      //      xflip = -1;
      //      yflip = -1;
      break;
    case 3:
      xangle =  TILT_ANGLE;
      yangle =  ROT_ANGLE;
      //      xflip = +1;
      //      yflip = -1;
      break;
  default:
    cout << "*** RotateToLab: wrong ndet = " << ndet << " ! ***" << endl;
  }

  //  cout << "RotateToLab:" << endl;
  //  cout << " ndet   = " << ndet << endl;
  //  cout << " xangle = " << xangle << endl;
  //  cout << " yangle = " << yangle << endl;
  //  cout << " xflip  = " << xflip << endl;
  //  cout << " yflip  = " << yflip << endl;

  //  cout << " vector before flip:" << endl;
  //  v.Print();

  //  v.SetX(xflip*v.X());
  //  v.SetY(yflip*v.Y());

  //  cout << " vector after flip:" << endl;
  //  v.Print();

  v.RotateX(xangle*TMath::DegToRad());

  //  cout << " vector after X-rotataion:" << endl;
  //  v.Print();

  v.RotateY(yangle*TMath::DegToRad());

  //  cout << " vector after XY-rotataion:" << endl;
  //  v.Print();
  //  getchar();
}

//------------------------------------------------------------------------------

double GEMSmear(const double xtrack) {
  return RandomGen.Gaus(xtrack, GEM_ACCURACY);
  //  return RandomGen.Gaus(xtrack, 0);
}

//------------------------------------------------------------------------------

double xslope_reconstruct(const double &xslope_det, const int &quarter) {

 double p0 = 0.;   // -0.000818, 0.0006236, -0.0004909, 0.0004326
 switch (quarter) {
 case 0:
 case 2:
   p0 = -0.00065;
   break;
 case 1:
 case 3:
   p0 = 0.00046;
   break;
 default:
   cout << "xslope_reconstruct: wrong quarter = " << quarter << " !" << endl;
 }

 double p1 = 1.;
 switch (quarter) {
 case 0:
 case 1:
   p1 = 0.963515;
   break;
 case 2:
 case 3:
   p1 = 1.06;
   break;
 default:
   cout << "xslope_reconstruct: wrong quarter = " << quarter << " !" << endl;
 }

  return (xslope_det - p0)/p1;
}
