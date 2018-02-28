#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH1.h"
#include "TH2.h"
#include <TStyle.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <iostream>
#include <sstream>   //for istringstream
#include <fstream>

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

// Histogram TCS kinematic variables of events accepted in the TCS setup.
// Events are tagged as accepted if there are significant energy depostions
// in the detectors.

#define ELEID 11
#define POSID -11
#define PROID 2212

//Thresholds on energy depositions [MeV].

#define EeleTXMIN 0.1
#define EposTXMIN 0.1
#define EproTXMIN 0.1

#define EeleTYMIN 0.1
#define EposTYMIN 0.1
#define EproTYMIN 0.1

#define EeleHXMIN 1.5
#define EposHXMIN 1.5
#define EproHXMIN 1.5

#define EeleHYMIN 1.5
#define EposHYMIN 1.5
#define EproHYMIN 1.5

/*
#define EeleTXMIN 0.12
#define EposTXMIN 0.12
#define EproTXMIN 0.18

#define EeleTYMIN 0.12
#define EposTYMIN 0.12
#define EproTYMIN 0.18

#define EeleHXMIN 1.5
#define EposHXMIN 1.5
#define EproHXMIN 2.2

#define EeleHYMIN 1.5
#define EposHYMIN 1.5
#define EproHYMIN 2.2
*/

#define EeleCALOMIN 1200
#define EposCALOMIN 1200

//Offsets.

#define TXOffset -17.5  //Horizontal offset, half of the width of detector
#define TYOffset   5.3  //Vertical offset to exclude +/- 6 deg.
#define HXOffset -50.   //Horizontal offset, half of the width of detector
#define HYOffset  15.5  //Vertical offset to exclude +/- 6 deg.
#define CXOffset -50.   //Horizontal offset, half of the width of detector
#define CYOffset  16.8  //Vertical offset to exclude +/- 6 deg.

#define TYSize 107.
#define HYSize 317.
#define CYSize 342.

//Transverse sizes of elements of detectors

#define dTX 0.1
#define dTY 0.1
#define dHX 1.
#define dHY 1.
#define dCa 2.

void GetEnergies(TChain &ch, int i,
		 TBranch* &det_br, TBranch* &edep_br, TBranch* &pid_br,
		 std::vector<int>* &det_vec,
		 std::vector<double>* &edep_vec,
		 std::vector<int>* &pid_vec,
		 double &eelepos, double &epospos, double &epropos,
		 double &eeleneg, double &eposneg, double &eproneg);

void GetEnergies(TChain &ch, int i,
		 TBranch* &det_br, TBranch* &edep_br,
		 std::vector<int>* &det_vec,
		 std::vector<double>* &edep_vec,
		 double &epos, double &eneg);

double GetCoordinate(TChain &ch, int i,
		     TBranch* &det_br, TBranch* &edep_br, TBranch* &pid_br,
		     TBranch* &chan_br,
		     std::vector<int>* &det_vec,
		     std::vector<double>* &edep_vec,
		     std::vector<int>* &pid_vec,
		     std::vector<int>* &chan_vec,
		     int pid, int side, double offset, double width);

void GetCoordinates(TChain &ch, int i,
		    TBranch* &det_br, TBranch* &edep_br, TBranch* &pid_br,
		    TBranch* &col_br, TBranch* &row_br,
		    std::vector<int>* &det_vec,
		    std::vector<double>* &edep_vec,
		    std::vector<int>* &pid_vec,
		    std::vector<int>* &col_vec, std::vector<int>* &row_vec,
		    int pid, int side, double xoffset, double yoffset,
		    double width,
		    double &x, double &y);

  void tcs_ana(string tracked_file, string primary_file, string tagged_file) {

  // TCS setup detectors' responses.

  //Chain file with tracked events

  TChain ch_tx("trackerx");   //a chain to process Tree tree
  TChain ch_ty("trackery");
  TChain ch_hx("hodox");
  TChain ch_hy("hodoy");
  TChain ch_ca("calo");
  TChain ch_kin("kin");

  if (FILE *file = fopen(tracked_file.c_str(), "r")) {
    fclose(file);
    ch_tx.Add(tracked_file.c_str(),0);
    ch_ty.Add(tracked_file.c_str(),0);
    ch_hx.Add(tracked_file.c_str(),0);
    ch_hy.Add(tracked_file.c_str(),0);
    ch_ca.Add(tracked_file.c_str(),0);
    ch_kin.Add(tracked_file.c_str(),0);
  }
  else
    cout << "  File: " << tracked_file << " does not exist" << endl;

  std::vector<int> *tx_vec = 0;
  std::vector<int> *ty_vec = 0;
  std::vector<int> *hx_vec = 0;
  std::vector<int> *hy_vec = 0;
  std::vector<int> *ca_vec = 0;

  std::vector<double> *tx_edep_vec = 0;
  std::vector<double> *ty_edep_vec = 0;
  std::vector<double> *hx_edep_vec = 0;
  std::vector<double> *hy_edep_vec = 0;
  std::vector<double> *ca_edep_vec = 0;

  std::vector<int> *tx_pid_vec = 0;
  std::vector<int> *ty_pid_vec = 0;
  std::vector<int> *hx_pid_vec = 0;
  std::vector<int> *hy_pid_vec = 0;
  std::vector<int> *ca_pid_vec = 0;

  std::vector<int> *tx_chan_vec = 0;
  std::vector<int> *ty_chan_vec = 0;
  std::vector<int> *hx_chan_vec = 0;
  std::vector<int> *hy_chan_vec = 0;
  std::vector<int> *ca_col_vec = 0;
  std::vector<int> *ca_row_vec = 0;

  TBranch *tx_br = 0;
  TBranch *ty_br = 0;
  TBranch *hx_br = 0;
  TBranch *hy_br = 0;
  TBranch *ca_br = 0;

  TBranch *tx_edep_br = 0;
  TBranch *ty_edep_br = 0;
  TBranch *hx_edep_br = 0;
  TBranch *hy_edep_br = 0;
  TBranch *ca_edep_br = 0;

  TBranch *tx_pid_br = 0;
  TBranch *ty_pid_br = 0;
  TBranch *hx_pid_br = 0;
  TBranch *hy_pid_br = 0;
  TBranch *ca_pid_br = 0;

  TBranch *tx_chan_br = 0;
  TBranch *ty_chan_br = 0;
  TBranch *hx_chan_br = 0;
  TBranch *hy_chan_br = 0;
  TBranch *ca_col_br = 0;
  TBranch *ca_row_br = 0;

  ch_tx.SetBranchAddress("detcont",&tx_vec,&tx_br);
  ch_ty.SetBranchAddress("detcont",&ty_vec,&ty_br);
  ch_hx.SetBranchAddress("detcont",&hx_vec,&hx_br);
  ch_hy.SetBranchAddress("detcont",&hy_vec,&hy_br);
  ch_ca.SetBranchAddress("detcont",&ca_vec,&ca_br);

  ch_tx.SetBranchAddress("edepcont",&tx_edep_vec,&tx_edep_br);
  ch_ty.SetBranchAddress("edepcont",&ty_edep_vec,&ty_edep_br);
  ch_hx.SetBranchAddress("edepcont",&hx_edep_vec,&hx_edep_br);
  ch_hy.SetBranchAddress("edepcont",&hy_edep_vec,&hy_edep_br);
  ch_ca.SetBranchAddress("edepcont",&ca_edep_vec,&ca_edep_br);

  ch_tx.SetBranchAddress("pidcont",&tx_pid_vec,&tx_pid_br);
  ch_ty.SetBranchAddress("pidcont",&ty_pid_vec,&ty_pid_br);
  ch_hx.SetBranchAddress("pidcont",&hx_pid_vec,&hx_pid_br);
  ch_hy.SetBranchAddress("pidcont",&hy_pid_vec,&hy_pid_br);
  ch_ca.SetBranchAddress("pidcont",&ca_pid_vec,&ca_pid_br);

  ch_tx.SetBranchAddress("chancont",&tx_chan_vec,&tx_chan_br);
  ch_ty.SetBranchAddress("chancont",&ty_chan_vec,&ty_chan_br);
  ch_hx.SetBranchAddress("chancont",&hx_chan_vec,&hx_chan_br);
  ch_hy.SetBranchAddress("chancont",&hy_chan_vec,&hy_chan_br);
  ch_ca.SetBranchAddress("colcont",&ca_col_vec,&ca_col_br);
  ch_ca.SetBranchAddress("rowcont",&ca_row_vec,&ca_row_br);

  //Open the primary root file of events from TCS geterator.
  
  // Kinematic quantities.

  Double_t psf;
  Double_t Q2;
  Double_t t;
  Double_t s;
  Double_t xi;
  Double_t tau, eta;
  Double_t phi_cm, the_cm;
  Double_t Eg;
  Double_t flux_factor;
  Double_t crs_BH;

  ch_kin.SetBranchAddress("psf", &psf);
  ch_kin.SetBranchAddress("Q2", &Q2);
  ch_kin.SetBranchAddress("t", &t);
  ch_kin.SetBranchAddress("s", &s);
  ch_kin.SetBranchAddress("xi", &xi);
  ch_kin.SetBranchAddress("tau", &tau);
  ch_kin.SetBranchAddress("eta", &eta);
  ch_kin.SetBranchAddress("the_cm", &the_cm);
  ch_kin.SetBranchAddress("phi_cm", &phi_cm);
  ch_kin.SetBranchAddress("Eg", &Eg);
  ch_kin.SetBranchAddress("flux_factor", &flux_factor);
  ch_kin.SetBranchAddress("crs_BH", &crs_BH);

  //Open the primary root file (from TCS generator).
  TFile *prim_file = new TFile(primary_file.c_str());
  TTree *primTree = (TTree*)prim_file->Get("TCS_Tree");

  TFile* tag_file = new TFile(tagged_file.c_str(),"recreate");
  TTree* taggedTree = primTree->CloneTree();

  ////  bool acceptance_flag;
  int acceptance_flag;
  TBranch *acceptance_br = taggedTree->Branch("acceptance_flag",
					      &acceptance_flag,
					      "acceptanse_flag/I");

  //Hit coordinates per particle and detector.

  double txelepos, tyelepos, hxelepos, hyelepos, cxelepos, cyelepos;
  double txeleneg, tyeleneg, hxeleneg, hyeleneg, cxeleneg, cyeleneg;

  double txpospos, typospos, hxpospos, hypospos, cxpospos, cypospos;
  double txposneg, typosneg, hxposneg, hyposneg, cxposneg, cyposneg;

  double txpropos, typropos, hxpropos, hypropos, cxpropos, cypropos;
  double txproneg, typroneg, hxproneg, hyproneg, cxproneg, cyproneg;

  //e-, + (top) det-s.
  TBranch *txelepos_b = taggedTree->Branch("tx_ele_pos",&txelepos,"txelepos/D");
  TBranch *tyelepos_b = taggedTree->Branch("ty_ele_pos",&tyelepos,"tyelepos/D");
  TBranch *hxelepos_b = taggedTree->Branch("hx_ele_pos",&hxelepos,"hxelepos/D");
  TBranch *hyelepos_b = taggedTree->Branch("hy_ele_pos",&hyelepos,"hyelepos/D");
  TBranch *cxelepos_b = taggedTree->Branch("cx_ele_pos",&cxelepos,"cxelepos/D");
  TBranch *cyelepos_b = taggedTree->Branch("cy_ele_pos",&cyelepos,"cyelepos/D");
  //e-, - (bottom) det-s.
  TBranch *txeleneg_b = taggedTree->Branch("tx_ele_neg",&txeleneg,"txeleneg/D");
  TBranch *tyeleneg_b = taggedTree->Branch("ty_ele_neg",&tyeleneg,"tyeleneg/D");
  TBranch *hxeleneg_b = taggedTree->Branch("hx_ele_neg",&hxeleneg,"hxeleneg/D");
  TBranch *hyeleneg_b = taggedTree->Branch("hy_ele_neg",&hyeleneg,"hyeleneg/D");
  TBranch *cxeleneg_b = taggedTree->Branch("cx_ele_neg",&cxeleneg,"cxeleneg/D");
  TBranch *cyeleneg_b = taggedTree->Branch("cy_ele_neg",&cyeleneg,"cyeleneg/D");
  //e+, + (top) det-s.
  TBranch *txpospos_b = taggedTree->Branch("tx_pos_pos",&txpospos,"txpospos/D");
  TBranch *typospos_b = taggedTree->Branch("ty_pos_pos",&typospos,"typospos/D");
  TBranch *hxpospos_b = taggedTree->Branch("hx_pos_pos",&hxpospos,"hxpospos/D");
  TBranch *hypospos_b = taggedTree->Branch("hy_pos_pos",&hypospos,"hypospos/D");
  TBranch *cxpospos_b = taggedTree->Branch("cx_pos_pos",&cxpospos,"cxpospos/D");
  TBranch *cypospos_b = taggedTree->Branch("cy_pos_pos",&cypospos,"cypospos/D");
  //e+, - (bottom) det-s.
  TBranch *txposneg_b = taggedTree->Branch("tx_pos_neg",&txposneg,"txposneg/D");
  TBranch *typosneg_b = taggedTree->Branch("ty_pos_neg",&typosneg,"typosneg/D");
  TBranch *hxposneg_b = taggedTree->Branch("hx_pos_neg",&hxposneg,"hxposneg/D");
  TBranch *hyposneg_b = taggedTree->Branch("hy_pos_neg",&hyposneg,"hyposneg/D");
  TBranch *cxposneg_b = taggedTree->Branch("cx_pos_neg",&cxposneg,"cxposneg/D");
  TBranch *cyposneg_b = taggedTree->Branch("cy_pos_neg",&cyposneg,"cyposneg/D");
  //proton, + (top) det-s.
  TBranch *txpropos_b = taggedTree->Branch("tx_pro_pos",&txpropos,"txpropos/D");
  TBranch *typropos_b = taggedTree->Branch("ty_pro_pos",&typropos,"typropos/D");
  TBranch *hxpropos_b = taggedTree->Branch("hx_pro_pos",&hxpropos,"hxpropos/D");
  TBranch *hypropos_b = taggedTree->Branch("hy_pro_pos",&hypropos,"hypropos/D");
  TBranch *cxpropos_b = taggedTree->Branch("cx_pro_pos",&cxpropos,"cxpropos/D");
  TBranch *cypropos_b = taggedTree->Branch("cy_pro_pos",&cypropos,"cypropos/D");
  //proton, - (bottom) det-s.
  TBranch *txproneg_b = taggedTree->Branch("tx_pro_neg",&txproneg,"txproneg/D");
  TBranch *typroneg_b = taggedTree->Branch("ty_pro_neg",&typroneg,"typroneg/D");
  TBranch *hxproneg_b = taggedTree->Branch("hx_pro_neg",&hxproneg,"hxproneg/D");
  TBranch *hyproneg_b = taggedTree->Branch("hy_pro_neg",&hyproneg,"hyproneg/D");
  TBranch *cxproneg_b = taggedTree->Branch("cx_pro_neg",&cxproneg,"cxproneg/D");
  TBranch *cyproneg_b = taggedTree->Branch("cy_pro_neg",&cyproneg,"cyproneg/D");

  // Histograms.

  TH2D* h_q2_t_det = new TH2D("h_q2_t_det", "Q'^{2} vs -t",
                              100, 0., 1., 120, 3.5, 9.5);
  TH2D* h_q2_t_vrt = (TH2D*) h_q2_t_det->Clone("h_q2_t_vrt");

  TH1D* h_t_det = new TH1D("h_t_det", "-t", 110, 0., 1.1);
  TH1D* h_q2_det = new TH1D("h_q2_det", "Q'^{2}", 120, 3.5, 9.5);
  TH1D* h_s_det = new TH1D("h_s_det", "s", 80, 8., 24.);
  TH1D* h_tau_det = new TH1D("h_tau_det", "#tau", 100, 0., 1.);
  TH1D* h_eta_det = new TH1D("h_eta_det", "#eta", 100, -0.5, 0.);

  TH2D* h_the_phi_cm_det = new TH2D("h_the_phi_cm_det",
                    "#theta_{CM} vs #phi_{CM}", 90, 0., 360., 45, 0., 180.);

  TH1D* h_phi_cm_det = new TH1D("h_phi_cm_det", "#phi_{CM}", 120, 0., 360.);
  TH1D* h_the_cm_det = new TH1D("h_the_cm_det", "#theta_{CM}", 60, 0., 180.);

  //

  // UVA HIPS source.
  // 1.5*10^12 photon flux,
  // 3 cm target length, 0.6 packing fraction,
  // 0.52*10^23 cm^-2 NH3 molecules,
  // 1.5 *10^23 cm^-2 polarized hydrogen atoms,
  // 30 days of running, 2.592*10^6 s:
  const Double_t L0 = 5.85e41;      //cm^-2

  const Double_t pbn = 1.e-36;

  std::ofstream ofs("acceptance_flag.dat");

  // Loop over entries.

  Long64_t nentries = ch_kin.GetEntriesFast();
  cout << "nentries = " << nentries << endl;

  for (Int_t i = 0; i < nentries; i++) {

    Long64_t tentry = ch_kin.LoadTree(i);
    //    det_br->GetEntry(tentry);

    ch_kin.GetEntry(i);

    //    cout << "Kimenatic quantities:" << endl;
    //    cout << "psf = " << psf << endl;
    //    cout << "Q2 = " << Q2 << endl;
    //    cout << "t = " << t << endl;
    //    cout << "s = " << s << endl;
    //    cout << "xi = " << xi << endl;
    //    cout << "tau = " << tau << endl;
    //    cout << "eta = " << eta << endl;
    //    cout << "phi_cm = " << phi_cm << endl;
    //    cout << "the_cm = " << the_cm << endl;
    //    cout << "Eg = " << Eg << endl;
    //    cout << "flux_factor = " << flux_factor << endl;
    //    cout << "crs_BH = " << crs_BH << endl;
    //    getchar();

    double weight = 0.;
    if (crs_BH<1000.) weight = flux_factor*psf*crs_BH*(pbn*L0)/nentries;
    //    if (crs_BH<1000.) weight = 1; //phase space calc.

    //Fill vertex histograms.

    h_q2_t_vrt->Fill(-t,Q2,weight);

    // Get hit coordinates in the dectors;

    txelepos=  GetCoordinate(ch_tx, i, tx_br, tx_edep_br, tx_pid_br, tx_chan_br,
			     tx_vec, tx_edep_vec, tx_pid_vec, tx_chan_vec,
			     ELEID, +1, TXOffset, dTX);
  
    tyelepos=  GetCoordinate(ch_ty, i, ty_br, ty_edep_br, ty_pid_br, ty_chan_br,
			     ty_vec, ty_edep_vec, ty_pid_vec, ty_chan_vec,
			     ELEID, +1, TYOffset, dTY);
  
    hxelepos=  GetCoordinate(ch_hx, i, hx_br, hx_edep_br, hx_pid_br, hx_chan_br,
			     hx_vec, hx_edep_vec, hx_pid_vec, hx_chan_vec,
			     ELEID, +1, HXOffset, dHX);
  
    hyelepos=  GetCoordinate(ch_hy, i, hy_br, hy_edep_br, hy_pid_br, hy_chan_br,
			     hy_vec, hy_edep_vec, hy_pid_vec, hy_chan_vec,
			     ELEID, +1, HYOffset, dHY);

    GetCoordinates(ch_ca, i, ca_br, ca_edep_br, ca_pid_br, ca_col_br, ca_row_br,
		   ca_vec, ca_edep_vec, ca_pid_vec, ca_col_vec, ca_row_vec,
		   ELEID, +1, CXOffset, CYOffset, dCa, cxelepos, cyelepos);

    txelepos_b->Fill();
    tyelepos_b->Fill();
    hxelepos_b->Fill();
    hyelepos_b->Fill();
    cxelepos_b->Fill();
    cyelepos_b->Fill();

    txeleneg=  GetCoordinate(ch_tx, i, tx_br, tx_edep_br, tx_pid_br, tx_chan_br,
			     tx_vec, tx_edep_vec, tx_pid_vec, tx_chan_vec,
			     ELEID, -1, TXOffset, dTX);
  
    tyeleneg=  GetCoordinate(ch_ty, i, ty_br, ty_edep_br, ty_pid_br, ty_chan_br,
			     ty_vec, ty_edep_vec, ty_pid_vec, ty_chan_vec,
			     ELEID, -1, -TYOffset-TYSize, dTY);
  
    hxeleneg=  GetCoordinate(ch_hx, i, hx_br, hx_edep_br, hx_pid_br, hx_chan_br,
			     hx_vec, hx_edep_vec, hx_pid_vec, hx_chan_vec,
			     ELEID, -1, HXOffset, dHX);
  
    hyeleneg=  GetCoordinate(ch_hy, i, hy_br, hy_edep_br, hy_pid_br, hy_chan_br,
			     hy_vec, hy_edep_vec, hy_pid_vec, hy_chan_vec,
			     ELEID, -1, -HYOffset-HYSize, dHY);
  
    GetCoordinates(ch_ca, i, ca_br, ca_edep_br, ca_pid_br, ca_col_br, ca_row_br,
		   ca_vec, ca_edep_vec, ca_pid_vec, ca_col_vec, ca_row_vec,
		   ELEID, -1, CXOffset, -CYOffset-CYSize, dCa,
		   cxeleneg, cyeleneg);

    txeleneg_b->Fill();
    tyeleneg_b->Fill();
    hxeleneg_b->Fill();
    hyeleneg_b->Fill();
    cxeleneg_b->Fill();
    cyeleneg_b->Fill();

    txpospos=  GetCoordinate(ch_tx, i, tx_br, tx_edep_br, tx_pid_br, tx_chan_br,
			     tx_vec, tx_edep_vec, tx_pid_vec, tx_chan_vec,
			     POSID, +1, TXOffset, dTX);
  
    typospos=  GetCoordinate(ch_ty, i, ty_br, ty_edep_br, ty_pid_br, ty_chan_br,
			     ty_vec, ty_edep_vec, ty_pid_vec, ty_chan_vec,
			     POSID, +1, TYOffset, dTY);
  
    hxpospos=  GetCoordinate(ch_hx, i, hx_br, hx_edep_br, hx_pid_br, hx_chan_br,
			     hx_vec, hx_edep_vec, hx_pid_vec, hx_chan_vec,
			     POSID, +1, HXOffset, dHX);
  
    hypospos=  GetCoordinate(ch_hy, i, hy_br, hy_edep_br, hy_pid_br, hy_chan_br,
			     hy_vec, hy_edep_vec, hy_pid_vec, hy_chan_vec,
			     POSID, +1, HYOffset, dHY);
  
    GetCoordinates(ch_ca, i, ca_br, ca_edep_br, ca_pid_br, ca_col_br, ca_row_br,
		   ca_vec, ca_edep_vec, ca_pid_vec, ca_col_vec, ca_row_vec,
		   POSID, +1, CXOffset, CYOffset, dCa, cxpospos, cypospos);

    txpospos_b->Fill();
    typospos_b->Fill();
    hxpospos_b->Fill();
    hypospos_b->Fill();
    cxpospos_b->Fill();
    cypospos_b->Fill();

    txposneg=  GetCoordinate(ch_tx, i, tx_br, tx_edep_br, tx_pid_br, tx_chan_br,
			     tx_vec, tx_edep_vec, tx_pid_vec, tx_chan_vec,
			     POSID, -1, TXOffset, dTX);
  
    typosneg=  GetCoordinate(ch_ty, i, ty_br, ty_edep_br, ty_pid_br, ty_chan_br,
			     ty_vec, ty_edep_vec, ty_pid_vec, ty_chan_vec,
			     POSID, -1, -TYOffset-TYSize, dTY);
  
    hxposneg=  GetCoordinate(ch_hx, i, hx_br, hx_edep_br, hx_pid_br, hx_chan_br,
			     hx_vec, hx_edep_vec, hx_pid_vec, hx_chan_vec,
			     POSID, -1, HXOffset, dHX);
  
    hyposneg=  GetCoordinate(ch_hy, i, hy_br, hy_edep_br, hy_pid_br, hy_chan_br,
			     hy_vec, hy_edep_vec, hy_pid_vec, hy_chan_vec,
			     POSID, -1, -HYOffset-HYSize, dHY);
  
    GetCoordinates(ch_ca, i, ca_br, ca_edep_br, ca_pid_br, ca_col_br, ca_row_br,
		   ca_vec, ca_edep_vec, ca_pid_vec, ca_col_vec, ca_row_vec,
		   POSID, -1, CXOffset, -CYOffset-CYSize, dCa,
		   cxposneg, cyposneg);

    txposneg_b->Fill();
    typosneg_b->Fill();
    hxposneg_b->Fill();
    hyposneg_b->Fill();
    cxposneg_b->Fill();
    cyposneg_b->Fill();

    txpropos=  GetCoordinate(ch_tx, i, tx_br, tx_edep_br, tx_pid_br, tx_chan_br,
			     tx_vec, tx_edep_vec, tx_pid_vec, tx_chan_vec,
			     PROID, +1, TXOffset, dTX);
  
    typropos=  GetCoordinate(ch_ty, i, ty_br, ty_edep_br, ty_pid_br, ty_chan_br,
			     ty_vec, ty_edep_vec, ty_pid_vec, ty_chan_vec,
			     PROID, +1, TYOffset, dTY);
  
    hxpropos=  GetCoordinate(ch_hx, i, hx_br, hx_edep_br, hx_pid_br, hx_chan_br,
			     hx_vec, hx_edep_vec, hx_pid_vec, hx_chan_vec,
			     PROID, +1, HXOffset, dHX);
  
    hypropos=  GetCoordinate(ch_hy, i, hy_br, hy_edep_br, hy_pid_br, hy_chan_br,
			     hy_vec, hy_edep_vec, hy_pid_vec, hy_chan_vec,
			     PROID, +1, HYOffset, dHY);
  
    GetCoordinates(ch_ca, i, ca_br, ca_edep_br, ca_pid_br, ca_col_br, ca_row_br,
		   ca_vec, ca_edep_vec, ca_pid_vec, ca_col_vec, ca_row_vec,
		   PROID, +1, CXOffset, CYOffset, dCa, cxpropos, cypropos);

    txpropos_b->Fill();
    typropos_b->Fill();
    hxpropos_b->Fill();
    hypropos_b->Fill();
    cxpropos_b->Fill();
    cypropos_b->Fill();

    txproneg=  GetCoordinate(ch_tx, i, tx_br, tx_edep_br, tx_pid_br, tx_chan_br,
			     tx_vec, tx_edep_vec, tx_pid_vec, tx_chan_vec,
			     PROID, -1, TXOffset, dTX);
  
    typroneg=  GetCoordinate(ch_ty, i, ty_br, ty_edep_br, ty_pid_br, ty_chan_br,
			     ty_vec, ty_edep_vec, ty_pid_vec, ty_chan_vec,
			     PROID, -1, -TYOffset-TYSize, dTY);
  
    hxproneg=  GetCoordinate(ch_hx, i, hx_br, hx_edep_br, hx_pid_br, hx_chan_br,
			     hx_vec, hx_edep_vec, hx_pid_vec, hx_chan_vec,
			     PROID, -1, HXOffset, dHX);
  
    hyproneg=  GetCoordinate(ch_hy, i, hy_br, hy_edep_br, hy_pid_br, hy_chan_br,
			     hy_vec, hy_edep_vec, hy_pid_vec, hy_chan_vec,
			     PROID, -1, -HYOffset-HYSize, dHY);
  
    GetCoordinates(ch_ca, i, ca_br, ca_edep_br, ca_pid_br, ca_col_br, ca_row_br,
		   ca_vec, ca_edep_vec, ca_pid_vec, ca_col_vec, ca_row_vec,
		   PROID, -1, CXOffset, -CYOffset-CYSize, dCa,
		   cxproneg, cyproneg);

    txproneg_b->Fill();
    typroneg_b->Fill();
    hxproneg_b->Fill();
    hyproneg_b->Fill();
    cxproneg_b->Fill();
    cyproneg_b->Fill();

    // Check detector signals (distinguish e+,e-,p vertices).

    double e_ele_pos_tx = 0.;
    double e_pos_pos_tx = 0.;
    double e_pro_pos_tx = 0.;
    double e_ele_neg_tx = 0.;
    double e_pos_neg_tx = 0.;
    double e_pro_neg_tx = 0.;

    GetEnergies(ch_tx, i,
		tx_br,  tx_edep_br,  tx_pid_br,
		tx_vec, tx_edep_vec, tx_pid_vec,
		e_ele_pos_tx, e_pos_pos_tx, e_pro_pos_tx,
		e_ele_neg_tx, e_pos_neg_tx, e_pro_neg_tx);

    //    cout << "Eele: " << eelepos << "  " << eeleneg << endl;
    //    cout << "Epos: " << epospos << "  " << eposneg << endl;
    //    cout << "Epro: " << epropos << "  " << eproneg << endl;
    //    getchar();

    double e_ele_pos_ty = 0.;
    double e_pos_pos_ty = 0.;
    double e_pro_pos_ty = 0.;
    double e_ele_neg_ty = 0.;
    double e_pos_neg_ty = 0.;
    double e_pro_neg_ty = 0.;

    GetEnergies(ch_ty, i,
		ty_br,  ty_edep_br,  ty_pid_br,
		ty_vec, ty_edep_vec, ty_pid_vec,
		e_ele_pos_ty, e_pos_pos_ty, e_pro_pos_ty,
		e_ele_neg_ty, e_pos_neg_ty, e_pro_neg_ty);

    double e_ele_pos_hx = 0.;
    double e_pos_pos_hx = 0.;
    double e_pro_pos_hx = 0.;
    double e_ele_neg_hx = 0.;
    double e_pos_neg_hx = 0.;
    double e_pro_neg_hx = 0.;

    GetEnergies(ch_hx, i,
		hx_br,  hx_edep_br,  hx_pid_br,
		hx_vec, hx_edep_vec, hx_pid_vec,
		e_ele_pos_hx, e_pos_pos_hx, e_pro_pos_hx,
		e_ele_neg_hx, e_pos_neg_hx, e_pro_neg_hx);

    double e_ele_pos_hy = 0.;
    double e_pos_pos_hy = 0.;
    double e_pro_pos_hy = 0.;
    double e_ele_neg_hy = 0.;
    double e_pos_neg_hy = 0.;
    double e_pro_neg_hy = 0.;

    GetEnergies(ch_hy, i,
		hy_br,  hy_edep_br,  hy_pid_br,
		hy_vec, hy_edep_vec, hy_pid_vec,
		e_ele_pos_hy, e_pos_pos_hy, e_pro_pos_hy,
		e_ele_neg_hy, e_pos_neg_hy, e_pro_neg_hy);

    double e_ele_pos_ca = 0.;
    double e_pos_pos_ca = 0.;
    double e_pro_pos_ca = 0.;
    double e_ele_neg_ca = 0.;
    double e_pos_neg_ca = 0.;
    double e_pro_neg_ca = 0.;

    GetEnergies(ch_ca, i,
		ca_br,  ca_edep_br,  ca_pid_br,
		ca_vec, ca_edep_vec, ca_pid_vec,
		e_ele_pos_ca, e_pos_pos_ca, e_pro_pos_ca,
		e_ele_neg_ca, e_pos_neg_ca, e_pro_neg_ca);

    bool good_ele_pos = e_ele_pos_tx > EeleTXMIN &&
                        e_ele_pos_ty > EeleTYMIN &&
                        e_ele_pos_hx > EeleHXMIN &&
                        e_ele_pos_hy > EeleHYMIN &&
                        e_ele_pos_ca > EeleCALOMIN;

    bool good_ele_neg = e_ele_neg_tx > EeleTXMIN &&
                        e_ele_neg_ty > EeleTYMIN &&
                        e_ele_neg_hx > EeleHXMIN &&
                        e_ele_neg_hy > EeleHYMIN &&
                        e_ele_neg_ca > EeleCALOMIN;

    bool good_pos_pos = e_pos_pos_tx > EposTXMIN &&
                        e_pos_pos_ty > EposTYMIN &&
                        e_pos_pos_hx > EposHXMIN &&
                        e_pos_pos_hy > EposHYMIN &&
                        e_pos_pos_ca > EposCALOMIN;

    bool good_pos_neg = e_pos_neg_tx > EposTXMIN &&
                        e_pos_neg_ty > EposTYMIN &&
                        e_pos_neg_hx > EposHXMIN &&
                        e_pos_neg_hy > EposHYMIN &&
                        e_pos_neg_ca > EposCALOMIN;

    bool good_pro_pos = e_pro_pos_tx > EproTXMIN &&
                        e_pro_pos_ty > EproTYMIN &&
                        e_pro_pos_hx > EproHXMIN &&
                        e_pro_pos_hy > EproHYMIN;

    bool good_pro_neg = e_pro_neg_tx > EproTXMIN &&
                        e_pro_neg_ty > EproTYMIN &&
                        e_pro_neg_hx > EproHXMIN &&
                        e_pro_neg_hy > EproHYMIN;

    if (good_ele_pos && good_ele_neg)
      cout << "*** good_ele_pos & good_ele_neg! ***" << endl;

    if (good_pos_pos && good_pos_neg)
      cout << "*** good_pos_pos & good_pos_neg! ***" << endl;

    if (good_pro_pos && good_pro_neg)
      cout << "*** good_pro_pos & good_pro_neg! ***" << endl;

    acceptance_flag = true;

    if (!good_ele_pos && !good_ele_neg) acceptance_flag = false;

    if (!good_pos_pos && !good_pos_neg) acceptance_flag = false;

    if (!good_pro_pos && !good_pro_neg) acceptance_flag = false;;

    if (!(good_ele_pos && good_pos_neg) &&
	!(good_ele_neg && good_pos_pos))
      acceptance_flag = false;

    ofs << acceptance_flag << endl;
    acceptance_br->Fill();

    if (!acceptance_flag) continue;
    
    // Fill detector histograms.

    h_q2_t_det->Fill(-t,Q2,weight);

    h_t_det->Fill(-t, weight);
    h_s_det->Fill(s, weight);
    //    h_xi_det->Fill(xi, weight);
    h_q2_det->Fill(Q2, weight);
    h_tau_det->Fill(tau, weight);
    h_eta_det->Fill(eta, weight);
    h_the_cm_det->Fill(the_cm*TMath::RadToDeg(), weight);
    h_phi_cm_det->Fill(phi_cm*TMath::RadToDeg(), weight);
    h_the_phi_cm_det->Fill(phi_cm*TMath::RadToDeg(), the_cm*TMath::RadToDeg(),
                           weight);
  }

  tag_file->cd();
  taggedTree->Write();

   // Plot.

  TCanvas *c7 = new TCanvas("kinematics", "TCS kinematics", 900, 900);
  c7->Divide(3, 3);

  // rot_ana histograms.

  TFile* rot_ana_file = new TFile(
  "/home/vardan/tcs.project/rot_ana/uva_field.right_xyz.rafo/rot_ana.root");

  c7->cd(2);
  TH1F* h_t_old = (TH1F*)rot_ana_file->Get("h_t_det");
  h_t_old->SetMaximum(h_t_det->GetMaximum()*1.1);
  h_t_old->SetLineColor(3);
  h_t_old->Draw("HIST");

  c7->cd(3);
  TH1F* h_q2_old = (TH1F*)rot_ana_file->Get("h_q2_det");
  h_q2_old->SetMaximum(h_q2_det->GetMaximum()*1.1);
  h_q2_old->SetLineColor(3);
  h_q2_old->Draw("HIST");

  c7->cd(4);
  TH1F* h_s_old = (TH1F*)rot_ana_file->Get("h_s_det");
  h_s_old->SetMaximum(h_s_det->GetMaximum()*1.1);
  h_s_old->SetLineColor(3);
  h_s_old->Draw("HIST");

  c7->cd(5);
  TH1F* h_tau_old = (TH1F*)rot_ana_file->Get("h_tau_det");
  h_tau_old->SetMaximum(h_tau_det->GetMaximum()*1.1);
  h_tau_old->SetLineColor(3);
  h_tau_old->Draw("HIST");

  c7->cd(6);
  TH1F* h_eta_old = (TH1F*)rot_ana_file->Get("h_eta_det");
  h_eta_old->SetMaximum(h_eta_det->GetMaximum()*1.1);
  h_eta_old->SetLineColor(3);
  h_eta_old->Draw("HIST");

  c7->cd(7);
  TH2F* h_the_phi_cm_old = (TH2F*)rot_ana_file->Get("h_the_phi_cm_det");
  h_the_phi_cm_old->SetMarkerColor(3);
  h_the_phi_cm_old->Draw();

  c7->cd(8);
  TH1F* h_phi_cm_old = (TH1F*)rot_ana_file->Get("h_phi_cm_det");
  h_phi_cm_old->SetMaximum(h_phi_cm_det->GetMaximum()*1.1);
  h_phi_cm_old->SetMinimum(0.);
  h_phi_cm_old->SetLineColor(3);
  h_phi_cm_old->Draw("HIST");

  c7->cd(9);
  TH1F* h_the_cm_old = (TH1F*)rot_ana_file->Get("h_the_cm_det");
  h_the_cm_old->SetMaximum(h_the_cm_det->GetMaximum()*1.1);
  h_the_cm_old->SetLineColor(3);
  h_the_cm_old->Draw("HIST");

  // G4 histograms.

  int ipad = 0;

  c7->cd(++ipad);
  h_q2_t_vrt->GetXaxis()->SetTitle("-t (GeV^{2})");
  h_q2_t_vrt->GetYaxis()->SetTitle("Q'^{2} (GeV^{2})");
  h_q2_t_vrt->SetMarkerColor(3);
  h_q2_t_vrt->Draw("HIST");
  h_q2_t_det->SetMarkerColor(2);
  h_q2_t_det->Draw("HIST same");

  c7->cd(++ipad);
  h_t_det->GetXaxis()->SetTitle("-t (GeV^{2})");
  h_t_det->SetLineColor(2);
  h_t_det->Draw("HIST same");

  c7->cd(++ipad);
  h_q2_det->GetXaxis()->SetTitle("Q'^{2} (GeV^{2})");
  h_q2_det->SetLineColor(2);
  h_q2_det->Draw("HIST same");

  c7->cd(++ipad);
  h_s_det->GetXaxis()->SetTitle("s (GeV^{2})");
  h_s_det->SetLineColor(2);
  h_s_det->Draw("HIST same");

  c7->cd(++ipad);
  h_tau_det->GetXaxis()->SetTitle("#tau");
  h_tau_det->SetLineColor(2);
  h_tau_det->Draw("HIST same");

  TLegend* legend = new TLegend(0.467,0.725,0.867,0.500);
  //  legend->SetHeader("The Legend Title","C");
  legend->AddEntry(h_tau_old,"C++/ROOT","l");
  legend->AddEntry(h_tau_det,"Geant 4","l");
  legend->Draw();

  c7->cd(++ipad);
  h_eta_det->GetXaxis()->SetTitle("#eta");
  h_eta_det->SetLineColor(2);
  h_eta_det->Draw("HIST same");

  c7->cd(++ipad);
  h_the_phi_cm_det->GetXaxis()->SetTitle("#phi_{CM} (deg)");
  h_the_phi_cm_det->GetYaxis()->SetTitle("#theta_{CM} (deg)");
  h_the_phi_cm_det->SetMarkerColor(2);
  h_the_phi_cm_det->Draw("HIST same");

  c7->cd(++ipad);
  h_phi_cm_det->GetXaxis()->SetTitle("#phi_{CM} (deg)");
  h_phi_cm_det->SetLineColor(2);
  h_phi_cm_det->Draw("HIST same");
  cout << "Integral " << h_phi_cm_det->Integral() << endl;

  c7->cd(++ipad);
  h_the_cm_det->GetXaxis()->SetTitle("#theta_{CM} (deg)");
  h_the_cm_det->SetLineColor(2);
  h_the_cm_det->Draw("HIST same");

  //  TCanvas *c11 = new TCanvas("momenta", "TCS momenta", 1500, 500);
  //  c11->Divide(3, 1);

}

//------------------------------------------------------------------------------

void GetEnergies(TChain &ch, int i,
		 TBranch* &det_br, TBranch* &edep_br, TBranch* &pid_br,
		 std::vector<int>* &det_vec,
		 std::vector<double>* &edep_vec,
		 std::vector<int>* &pid_vec,
		 double &eelepos, double &epospos, double &epropos,
		 double &eeleneg, double &eposneg, double &eproneg) {

  Long64_t tentry = ch.LoadTree(i);
  if (i < 10) cout << "tentry = " << tentry << "  i = " << i << endl;

  det_br->GetEntry(tentry);
  edep_br->GetEntry(tentry);
  pid_br->GetEntry(tentry);

  eelepos = 0.;
  epospos = 0.;
  epropos = 0.;
  eeleneg = 0.;
  eposneg = 0.;
  eproneg = 0.;

  //  cout << "GetEnergies: det_vec size = " << det_vec->size() << endl;

  for (UInt_t j = 0; j < det_vec->size(); ++j) {
 
    int det = det_vec->at(j);
    double edep = edep_vec->at(j);
    int pid = pid_vec->at(j);
    if (pid != ELEID && pid != POSID && pid != PROID)
      cout << "*** Wrong particle ID: " << pid << " ***" << endl;

    //cout<< " hit " << j << ": " << det << " " << pid << " " << edep << endl;

    if (det > 0)
      switch (pid) {
      case ELEID : eelepos += edep; break;
      case POSID : epospos += edep; break;
      case PROID : epropos += edep; break;
      }
    else
      switch (pid) {
      case ELEID : eeleneg += edep; break;
      case POSID : eposneg += edep; break;
      case PROID : eproneg += edep; break;
      }

  }

}

//------------------------------------------------------------------------------

void GetEnergies(TChain &ch, int i,
		 TBranch* &det_br, TBranch* &edep_br,
		 std::vector<int>* &det_vec,
		 std::vector<double>* &edep_vec,
		 double &epos, double &eneg) {

  Long64_t tentry = ch.LoadTree(i);
  if (i < 10) cout << "tentry = " << tentry << "  i = " << i << endl;

  det_br->GetEntry(tentry);
  edep_br->GetEntry(tentry);

  epos = 0.;
  eneg = 0.;

  //  cout << "GetEnergies: det_vec size = " << det_vec->size() << endl;

  for (UInt_t j = 0; j < det_vec->size(); ++j) {
 
    int det = det_vec->at(j);
    double edep = edep_vec->at(j);

    //cout<< " hit " << j << ": " << det << " " << pid << " " << edep << endl;

    if (det > 0)
      epos += edep;
    else
      eneg += edep;
  }

};

//------------------------------------------------------------------------------

double GetCoordinate(TChain &ch, int i,
		     TBranch* &det_br, TBranch* &edep_br, TBranch* &pid_br,
		     TBranch* &chan_br,
		     std::vector<int>* &det_vec,
		     std::vector<double>* &edep_vec,
		     std::vector<int>* &pid_vec,
		     std::vector<int>* &chan_vec,
		     int pid, int side, double offset, double width) {
  
  Long64_t tentry = ch.LoadTree(i);
  //  if (i < 10) cout << "tentry = " << tentry << "  i = " << i << endl;

  det_br->GetEntry(tentry);
  edep_br->GetEntry(tentry);
  pid_br->GetEntry(tentry);
  chan_br->GetEntry(tentry);
  
  //  cout << "GetCoordinate: det_vec size = " << det_vec->size() << endl;

  double meanX = 0.;
  double esum = 0;

  for (UInt_t j = 0; j < det_vec->size(); ++j) {

    if (det_vec->at(j)* side > 0 && pid_vec->at(j) == pid) {
      double edep = edep_vec->at(j);
      esum += edep;
      meanX += edep * (chan_vec->at(j)+0.5)*width;
    }

  }

  if (esum > 0.) {
    meanX /= esum;
    meanX += offset;
  }
  else
    meanX = 999999;

  return meanX;
}

//------------------------------------------------------------------------------

void GetCoordinates(TChain &ch, int i,
		    TBranch* &det_br, TBranch* &edep_br, TBranch* &pid_br,
		    TBranch* &col_br, TBranch* &row_br,
		    std::vector<int>* &det_vec,
		    std::vector<double>* &edep_vec,
		    std::vector<int>* &pid_vec,
		    std::vector<int>* &col_vec, std::vector<int>* &row_vec,
		    int pid, int side, double xoffset, double yoffset,
		    double width,
		    double &x, double &y) {
  
  Long64_t tentry = ch.LoadTree(i);
  //  if (i < 10) cout << "tentry = " << tentry << "  i = " << i << endl;

  det_br->GetEntry(tentry);
  edep_br->GetEntry(tentry);
  pid_br->GetEntry(tentry);
  col_br->GetEntry(tentry);
  row_br->GetEntry(tentry);
  
  //  cout << "GetCoordinate: det_vec size = " << det_vec->size() << endl;

  double meanX = 0.;
  double meanY = 0.;
  double esum = 0;

  for (UInt_t j = 0; j < det_vec->size(); ++j) {

    if (det_vec->at(j)* side > 0 && pid_vec->at(j) == pid) {
      double edep = edep_vec->at(j);
      esum += edep;
      meanX += edep * (col_vec->at(j)+0.5)*width;
      meanY += edep * (row_vec->at(j)+0.5)*width;
    }

  }

  if (esum > 0.) {
    meanX /= esum;
    meanY /= esum;
    meanX += xoffset;
    meanY += yoffset;
  }
  else {
    meanX = 999999;
    meanY = 999999;
  }

  x = meanX;
  y = meanY;
}
