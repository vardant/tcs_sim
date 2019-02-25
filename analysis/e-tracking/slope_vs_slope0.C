#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH1.h"
#include "TProfile.h"
#include "TStyle.h"
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <TAxis.h>

#include <TVector3.h>

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

#define TRACKER1_DIST 1200.   //mm, consistent with G4 coding.
#define TRACKER2_DIST 1300.
#define TRACKER3_DIST 1400.

#define TILT_ANGLE 13.835     //deg, consistent with G4 coding.
#define ROT_ANGLE  10.034

#define GEM_ACCURACY 100.E-3   //um

struct hit {
  double x;
  double y;
  int det;
  int layer;
};

double FitStraightTrack(const vector<hit> &hitlist, const int axis,
			double &offset, double &slope);

void RotateToLab(TVector3 &v, const int ndet);

//Slope at GEM trackers versus slope at vertex, for fixed momentum
// at vertex of 5 GeV.

void slope_vs_slope0(string rootfile="RootFiles/e-pyz1thx18thy30.root",
		     int quarter=3, string comment="test") {

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

  Double_t        py;
  Double_t        pz;
  TBranch        *b_py;
  TBranch        *b_pz;
  chain_beam.SetBranchAddress("py", &py, &b_py);
  chain_beam.SetBranchAddress("pz", &pz, &b_pz);

  uint max_hitlist_size = 0;

  string htitle = "Slope vs slope at vertex, "+rootfile;
  TH2F* h_slope_vs_slope0 = new TH2F("h_slope_vs_slope0", htitle.c_str(),
				     100, -0.5,-0.1, 100, -0.45, -0.05); //3GeV
				     //100, -0.55,-0.1, 100, -0.45,-0.05);//2GeV
                                     //100, -0.55, -0.2, 100, -0.45, 0.); //1GeV

				     //100, 0.,0.4, 100, 0., 0.5); //2GeV
				     //100, -0.1,0.4, 100, 0., 0.5); //2GeV
				     //100, -0.2,0.3, 100, 0., 0.5); //1GeV

  //				     100, -0.5,-0.1, 100, -0.45, -0.05); //3GeV
  //				     100, -0.55,-0.1, 100, -0.45, -0.05); //2GeV
  //				     100, -0.55, -0.2, 100, -0.45, 0.); //1GeV

				     //100, -0.05, 0.35, 100, 0., 0.45);//1.894GeV
				     //100, 0., 0.4, 100, 0., 0.45);//3GeV
				     //100, 0., 0.45, 100, 0., 0.45); //4.894GeV
				     //100, -0.2, 0.3, 100, 0., 0.45);//0.894GeV

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
	hit goodHit = {xcont->at(j), ycont->at(j), detcont->at(j),
		       layercont->at(j)};
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

      double xslope = 0.;
      double xoffset = 0.;
      double xChi2 = FitStraightTrack(hitlist, 0, xoffset, xslope);

      double yslope = 0.;
      double yoffset = 0.;
      double yChi2 = FitStraightTrack(hitlist, 1, yoffset, yslope);

      //      cout << "hitlist size = " << hitlist.size() << endl;
      //      cout << "xslope  = " << xslope << endl;
      //      cout << "xoffset = " << xoffset << endl;
      //      cout << "xChi2   = " << xChi2 << endl;
      //      cout << "yslope  = " << yslope << endl;
      //      cout << "yoffset = " << yoffset << endl;
      //      cout << "yChi2   = " << yChi2 << endl;
      //      getchar();

      //      TVector3 voff(xoffset, yoffset, 0.);
      //      RotateToLab(voff, ndet);

      TVector3 vdir(xslope, yslope, 1.);
      ////      TVector3 vdir(0., 0., 1.);

      RotateToLab(vdir, quarter);

      h_slope_vs_slope0->Fill(atan(py/pz), atan(vdir.Y()/vdir.Z()));

      //      if(py/pz>0.2) {
      //	cout << "slope0 = " << py/pz << endl;
      //	cout << "slope  = " << vdir.Y()/vdir.Z() << endl;
      //      }
    }

  }  //over entries

  new TCanvas("slope_vs_slope0");
  gStyle->SetOptFit();
  h_slope_vs_slope0->GetXaxis()->SetTitle("#theta_{Vertex} [rad]");
  h_slope_vs_slope0->GetYaxis()->SetTitle("#theta_{GEMs} [rad]");
  //  h_slope_vs_slope0->Draw();
  TH1D* h_prof = h_slope_vs_slope0->ProfileX();
  h_prof->GetXaxis()->SetTitle("#theta at Vertex [rad]");
  h_prof->GetYaxis()->SetTitle("#theta at GEMs [rad]");
  h_prof->Fit("pol1","","",-0.41,-0.13); //8GeV
  //  h_prof->Fit("pol1","","",-0.42,-0.13); //5GeV
  //  h_prof->Fit("pol1","","",-0.43,-0.16); //4GeV
  //  h_prof->Fit("pol1","","",-0.45,-0.18); //3GeV
  //  h_prof->Fit("pol1","","",-0.48,-0.21); //2GeV
  //  h_prof->Fit("pol1","","",-0.52,-0.31); //1GeV
  
  //  h_prof->Fit("pol1","","",0.08,0.37); //7GeV
  //  h_prof->Fit("pol1","","",0.07,0.36); //6GeV
  //  h_prof->Fit("pol1","","",0.07,0.35); //5GeV
  //  h_prof->Fit("pol1","","",0.06,0.34); //4GeV
  //  h_prof->Fit("pol1","","",0.05,0.32); //3GeV
  //  h_prof->Fit("pol1","","",0.02,0.28); //2GeV
  //  h_prof->Fit("pol1","","",-0.07,0.19); //1GeV

  //  h_prof->Fit("pol1","","",-0.41,-0.13); //8GeV
  //  h_prof->Fit("pol1","","",-0.42,-0.14); //5GeV
  //  h_prof->Fit("pol1","","",-0.44,-0.15); //4GeV
  //  h_prof->Fit("pol1","","",-0.45,-0.18); //3GeV
  //  h_prof->Fit("pol1","","",-0.48,-0.21); //2GeV
  //  h_prof->Fit("pol1","","",-0.52,-0.31); //1GeV

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
