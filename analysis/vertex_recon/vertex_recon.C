#define TCSTracker_cxx
#define TCSKin_cxx

#include <iostream>
#include "include/TCSTracker.h"
#include "include/TCSKin.h"

#include <TChain.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "TH1.h"
#include "TH2.h"
#include <TCanvas.h>
#include <TRandom.h>
#include <TMath.h>
#include <TStyle.h>
#include <TProfile.h>

using namespace std;

#define ELEID   11
#define POSID  -11
#define PROID 2212

#define NQUARTER 4
#define NLAYER   3

#define Mp 0.938

#define M2Pg_MIN -2.e38   //GeV^2
#define M2Pg_MAX  2.e38

TRandom RandomGen;

double edep_calo(double p);
double XSlopeReconP(double xslope_det, int quarter);
double YSlopeReconE(const double yslope, const double pyz, int particle_id);
double YSlopeReconP(const double yslope, const double pyz, int quarter);

void vertex_recon(int runlo=1000, int runhi=1010) {

  //Chain root files
  //  TChain chain_tracker("tracker");   //a chain to process Tree "tracker"
  //  TChain chain_beam("beam");
  //  TChain chain_kin("kin");
  //  int nchain = 0;
  //  for (int runno=runlo; runno<runhi; runno++) {
  //    string rootfile=Form("RootFiles.tracked/DEEPGen_%d_tracked.root",runno);
  //    if (FILE *file = fopen(rootfile.c_str(), "r")) {
  //      fclose(file);
  //      chain_tracker.Add(rootfile.c_str(),0);
  //      chain_beam.Add(rootfile.c_str(),0);
  //      chain_kin.Add(rootfile.c_str(),0);
  //      //      cout << rootfile << " chained" << endl;
  //      nchain++;
  //    }
  //    //    else
  //    //      cout << "  File: " << rootfile << " does not exist" << endl;
  //  }
  //  chain_tracker.MakeClass("TCSTracker");
  //  chain_beam.MakeClass("TCSBeam");
  //  chain_kin.MakeClass("TCSKin");

  TCSTracker tracker(runlo, runhi);
  TCSKin kin(runlo, runhi);

  Long64_t nentries = tracker.GetEntries();
  cout << "nentries = " << nentries << endl;

  const string particle_name[3] {"e-",  "e+",  "P"};
  const int particle_id[3]      {ELEID, POSID, PROID};
  for (int i=0; i<3; i++)
    cout << particle_name[i] << "  " << particle_id[i] << endl;

  ////  RandomGen.SetSeed(0);   //seed from comp. clock.

  TH1F* h_dphi[3];
  for (int i=0; i<3; i++) {
    string name = "h_dphi_"+ particle_name[i];
    string title = "#delta#phi for "+particle_name[i];
    h_dphi[i] = new TH1F(name.c_str(), title.c_str(), 100, -0.05, 0.05);
    h_dphi[i]->GetXaxis()->SetTitle("#delta#phi [rad]");
    h_dphi[i]->GetYaxis()->SetTitle("count/0.001 [rad^{-1}]");
  }

  TH1F* h_dtheta[3];
  for (int i=0; i<3; i++) {
    string name = "h_dtheta_"+ particle_name[i];
    string title = "#delta#theta for "+particle_name[i];
    double thlo = -0.05;
    double thhi =  0.05;
    if (i==2) {
      thlo = -0.25;
      thhi =  0.25;
    }
    h_dtheta[i] = new TH1F(name.c_str(), title.c_str(), 100, thlo, thhi);
    h_dtheta[i]->GetXaxis()->SetTitle("#delta#theta [rad]");
    h_dtheta[i]->GetYaxis()->SetTitle("count/0.001 [rad^{-1}]");
  }

  TH1F* h_dphi_gstar = new TH1F("h_dphi_gstar", "#delta#phi for #gamma*",
				100, -0.05, 0.05);
  TH1F* h_dtheta_gstar=new TH1F("h_dtheta_gstar","#delta#theta for #gamma*",
				100, -0.05, 0.05);
  TH1F* h_de_gstar = new TH1F("h_de_gstar", "#deltaE for #gamma*",
			      100, -0.5, 0.5);

  h_dphi_gstar->GetXaxis()->SetTitle("#delta#theta [rad]");
  h_dphi_gstar->GetYaxis()->SetTitle("count/0.001 [rad^{-1}]");
  h_dtheta_gstar->GetXaxis()->SetTitle("#delta#theta [rad]");
  h_dtheta_gstar->GetYaxis()->SetTitle("count/0.001 [rad^{-1}]");
  h_de_gstar->GetXaxis()->SetTitle("#deltaE [GeV]");
  h_de_gstar->GetYaxis()->SetTitle("count/0.001 [GeV^{-1}]");

  TH1F* h_dpx_p = new TH1F("h_dpx_p", "#deltaPx for P", 100, -0.5, 0.5);
  TH1F* h_dpy_p = new TH1F("h_dpy_p", "#deltaPy for P", 100, -0.5, 0.5);
  TH1F* h_dpz_p = new TH1F("h_dpz_p", "#deltaPz for P", 100, -0.5, 0.5);
  TH1F* h_de_p  = new TH1F("h_de_p",  "#deltaE for P", 100,  -0.5, 0.5);

  //Kinematics' histograms.

  TH1F* h_dt = new TH1F("h_dt",  "#deltat", 100,  -1, 1.);
  h_dt->GetXaxis()->SetTitle("#deltat [GeV^{2}]");

  TH2F* h_dt_vs_t = new TH2F("h_dt_vs_t", "#deltat vs -t",
			     20, 0., 2., 100, -1., 1.);
  h_dt_vs_t->GetXaxis()->SetTitle("-t_{VERTEX} [GeV^{2}]");
  h_dt_vs_t->GetYaxis()->SetTitle("#deltat [GeV^{2}]");

  TH1F* h_dq2 = new TH1F("h_dq2",  "#deltaQ'^{2}", 100,  -1., 1.);
  h_dq2->GetXaxis()->SetTitle("#deltaQ'^{2} [GeV^{2}]");

  TH2F* h_dq2_vs_q2 = new TH2F("h_dq2_vs_q2", "#deltaQ'^{2} vs Q'^{2}",
			       100, 3., 10., 100, -1., 1.);
  h_dq2_vs_q2->GetXaxis()->SetTitle("Q'^{2} [GeV^{2}]");
  h_dq2_vs_q2->GetYaxis()->SetTitle("#deltaQ'^{2}_{VETEX} [GeV^{2}]");

  TH1F* h_ds = new TH1F("h_ds",  "#deltas", 100,  -5, 5.);
  h_ds->GetXaxis()->SetTitle("#deltas [GeV^{2}]");

  TH2F* h_ds_vs_s = new TH2F("h_ds_vs_s", "#deltas vs s",
			     100, 8., 24., 100, -5., 5.);
  h_ds_vs_s->GetXaxis()->SetTitle("s_{VERTEX} [GeV^{2}]");
  h_ds_vs_s->GetYaxis()->SetTitle("#deltas [GeV^{2}]");

  TH1F* h_dxi = new TH1F("h_dxi",  "#delta#xi", 100,  -0.25, 0.25);
  h_dxi->GetXaxis()->SetTitle("#delta#xi");

  TH2F* h_dxi_vs_xi = new TH2F("h_dxi_vs_xi", "#delta#xi vs #xi",
			       100, 0., 0.6, 100, -0.25, 0.25);
  h_dxi_vs_xi->GetXaxis()->SetTitle("#xi_{VERTEX}");
  h_dxi_vs_xi->GetYaxis()->SetTitle("#delta#xi");

  TH1F* h_dtau = new TH1F("h_dtau",  "#delta#tau", 100,  -0.25, 0.25);
  h_dtau->GetXaxis()->SetTitle("#delta#tau");

  TH2F* h_dtau_vs_tau = new TH2F("h_dtau_vs_tau", "#delta#tau vs #tau",
				 100, 0.1, 0.7, 100, -0.25, 0.25);
  h_dtau_vs_tau->GetXaxis()->SetTitle("#tau_{VERTEX}");
  h_dtau_vs_tau->GetYaxis()->SetTitle("#delta#tau");

  TH1F* h_deta = new TH1F("h_deta",  "#delta#eta", 100,  -0.25, 0.25);
  h_deta->GetXaxis()->SetTitle("#delta#eta");

  TH2F* h_deta_vs_eta = new TH2F("h_deta_vs_eta", "#delta#eta vs #eta",
				 100, 0., 0.6, 100, -0.25, 0.25);
  h_deta_vs_eta->GetXaxis()->SetTitle("-#eta_{VERTEX}");
  h_deta_vs_eta->GetYaxis()->SetTitle("#delta#eta");

  TH1F* h_dthetaCM = new TH1F("h_dthetaCM","#delta#theta_{CM}", 100,-0.5,0.5);
  h_dthetaCM->GetXaxis()->SetTitle("#delta#theta_{CM} [rad]");

  TH2F* h_dthetaCM_vs_thetaCM = new TH2F("h_dthetaCM_vs_thetaCM",
  "#delta#theta_{CM} vs #theta_{CM}", 100, 0., TMath::Pi(), 100, -0.5, 0.5);
  h_dthetaCM_vs_thetaCM->GetXaxis()->SetTitle("#theta_{CM} (VERTEX) [rad]");
  h_dthetaCM_vs_thetaCM->GetYaxis()->SetTitle("#delta#theta_{CM} [rad]");

  TH1F* h_dphiCM = new TH1F("h_dphiCM","#delta#phi_{CM}", 100,-0.5,0.5);
  h_dphiCM->GetXaxis()->SetTitle("#delta#phi_{CM} [rad]");

  TH2F* h_dphiCM_vs_phiCM = new TH2F("h_dphiCM_vs_phiCM",
  "#delta#phi_{CM} vs #phi_{CM}", 100, 0., TMath::Pi(), 100, -0.5, 0.5);
  h_dphiCM_vs_phiCM->GetXaxis()->SetTitle("#phi_{CM} (VERTEX) [rad]");
  h_dphiCM_vs_phiCM->GetYaxis()->SetTitle("#delta#phi_{CM} [rad]");

  TH2F* h_phiCM_vs_phiCM = new TH2F("h_phiCM_vs_phiCM",
 "phi_{CM} vs #phi_{CM}", 100, 0., 2*3.14159, 100, -2*3.14159,2*3.14159);

  TH1F* h_eg = new TH1F("h_eg",  "E#gamma reconstructed", 100,  0., 15.);
  h_eg->GetXaxis()->SetTitle("E#gamma [GeV]");
  TH1F* h_deg = new TH1F("h_deg",  "#deltaE#gamma", 100,  -1, 1.);
  h_deg->GetXaxis()->SetTitle("#deltaE#gamma [GeV]");
  TH1F* h_m2dpg = new TH1F("h_m2dpg", "M^{2}(#deltaP#gamma)", 100,  -1, 1.);
  h_m2dpg->GetXaxis()->SetTitle("M^{2}(#deltaP#gamma) [GeV^{2}]");
  TH1F* h_m2pg = new TH1F("h_m2pg",  "M^{2}(P#gamma)", 100,  -5, 5.);
  h_m2pg->GetXaxis()->SetTitle("M^{2}(P#gamma) [GeV^{2}]");

  for (int ientry=0; ientry<nentries; ientry++) {

    TVector3 vdir[3];    //tracks' directions at GEM trackers
    int hit_quarter[3];

    int num_good_track =0;
    for (int iparticle=0; iparticle<3; iparticle++) {
      for (int quarter=0; quarter<NQUARTER; quarter++) {
	if (tracker.GoodTrack(ientry, particle_id[iparticle], iparticle+1,
			      quarter, vdir[iparticle])) {
	  hit_quarter[iparticle] = quarter;
	  num_good_track++;
	  //	  cout << "quarter = " << quarter
	  //	       << "  particle_id = " << particle_id[iparticle]
	  //	       << "  num_good_track = " << num_good_track << endl;
	  //	  getchar();
	  break;
	}
      }
    }

    //    if (num_good_track > 2) {
    //cout << "ientry = " << ientry << "   num_good_track = " << num_good_track
    //	   << endl;
    //      for (int i=0; i<3; i++)
    //	vdir[i].Print();
    //      getchar();
    //    }

    if (num_good_track == 3) {

      kin.GetEntry(ientry);
      KinVar kv = kin.GetKinVar();

      TLorentzVector Pg_vertex(0., 0., kv.Eg, kv.Eg);   //Eg in GeV
      TLorentzVector P_target(0., 0., 0., Mp);
      TLorentzVector P_vertex[3];
      TLorentzVector Pgstar_vertex;

      TLorentzVector Pg_recon;
      TLorentzVector P_recon[3];
      TLorentzVector Pgstar_recon;

      P_vertex[0].SetPx(kv.pminus[1]/1000.);
      P_vertex[0].SetPy(kv.pminus[2]/1000.);
      P_vertex[0].SetPz(kv.pminus[3]/1000.);
      P_vertex[0].SetE (kv.pminus[0]/1000.);

      P_vertex[1].SetPx(kv.pplus[1]/1000.);
      P_vertex[1].SetPy(kv.pplus[2]/1000.);
      P_vertex[1].SetPz(kv.pplus[3]/1000.);
      P_vertex[1].SetE (kv.pplus[0]/1000.);

      P_vertex[2].SetPx(kv.precoil[1]/1000.);
      P_vertex[2].SetPy(kv.precoil[2]/1000.);
      P_vertex[2].SetPz(kv.precoil[3]/1000.);
      P_vertex[2].SetE (kv.precoil[0]/1000.);

      Pgstar_vertex = P_vertex[0] + P_vertex[1];

      double phi_vertex[3];
      double theta_vertex[3];

      double phi_recon[3];
      double theta_recon[3];

      phi_vertex[0] = atan(kv.pminus[1]/kv.pminus[3]);
      phi_vertex[1] = atan(kv.pplus[1]/kv.pplus[3]);
      phi_vertex[2] = atan(kv.precoil[1]/kv.precoil[3]);

      theta_vertex[0] = atan(kv.pminus[2]/kv.pminus[3]);
      theta_vertex[1] = atan(kv.pplus[2]/kv.pplus[3]);
      theta_vertex[2] = atan(kv.precoil[2]/kv.precoil[3]);

      for (int i=0; i<2; i++)
	phi_recon[i] = atan(vdir[i].X()/vdir[i].Z());

      phi_recon[2] = XSlopeReconP(atan(vdir[2].X()/vdir[2].Z()),hit_quarter[2]);

      for (int i=0; i<2; i++) {

	double theta_det = atan(vdir[i].Y()/vdir[i].Z());
	//      double yposition_det = atan(vxyz.Y()/vxyz.Z());

	double px, py, pz;
	switch (i) {
	case 0: px = kv.pminus[1]; py = kv.pminus[2]; pz = kv.pminus[3]; break;
	case 1: px = kv.pplus[1];  py = kv.pplus[2];  pz = kv.pplus[3]; break;
	default: ;
	}

	double e_det = edep_calo(sqrt(px*px+py*py+pz*pz)/1000.);

	double pyz_det = e_det * 
	  sqrt(vdir[i].Y()*vdir[i].Y()+vdir[i].Z()*vdir[i].Z())/vdir[i].Mag();

	theta_recon[i] = YSlopeReconE(theta_det, pyz_det, particle_id[i]);

	TVector3 p_recon(tan(phi_recon[i]), tan(theta_recon[i]), 1.);
	p_recon *= e_det/p_recon.Mag();

	P_recon[i].SetE(e_det);
	P_recon[i].SetVect(p_recon);

      } //e-e+

      Pgstar_recon = P_recon[0] + P_recon[1];

      //Use co-planarity.
      TVector3 p_recon(-Pgstar_recon.Px(), -Pgstar_recon.Py(),
		       -Pgstar_recon.Px()/tan(phi_recon[2]));
                       ////      sqrt(P_vertex[2].Py()*P_vertex[2].Py()+
                       ////      P_vertex[2].Pz()*P_vertex[2].Pz())*
		       ////cos(theta_recon[2]));
		       ////	       -Pgstar_recon.Py()/tan(theta_vertex[2]));
		       ////	       -Pgstar_vertex.Px()/tan(phi_recon[2]));
		       ////	       P_vertex[2].Pz());
		       ////	       -Pgstar_recon.Px()/tan(phi_vertex[2]));

      //Iterate theta_recon, p_recon.Z of proton.

      double theta_det = atan(vdir[2].Y()/vdir[2].Z());
      double pyz = sqrt(p_recon.Y()*p_recon.Y() + p_recon.Z()*p_recon.Z());
      double theta_rec = YSlopeReconP(theta_det, pyz, hit_quarter[2]);
      p_recon.SetZ(pyz*cos(theta_rec));

      double e_recon = sqrt(Mp*Mp + p_recon.Mag2());

      P_recon[2].SetVect(p_recon);
      P_recon[2].SetE(e_recon);

      Pg_recon = Pgstar_recon + P_recon[2] - P_target;

      //      if (P_recon[2].Pz() < 0.)
      //	cout << "Pz_proton = " << P_recon[2].Pz() << endl;

      theta_recon[2] = atan(P_recon[2].Py()/P_recon[2].Pz());
      ////      theta_recon[2] = atan(Pgstar_recon.Py()/Pgstar_recon.Px()*
      ////			    tan(phi_recon[2]));

      /*
      TVector3 beta = Pgstar_recon.BoostVector();
      TLorentzVector Pcm_recon[3];
      for (int i=0; i<3; i++) {
	Pcm_recon[i] = P_recon[i];
	Pcm_recon[i].Boost(-beta);
      }

      double thetaCM = Pcm_recon[0].Angle(-Pcm_recon[2].Vect());
      //double thetaCM = TMath::Pi() - Pcm_recon[0].Angle(Pcm_recon[2].Vect());
      ////      double thetaCM = Pcm_recon[0].Angle(TVector3(0.,0.,1.));
      ////      TVector3 pcm = Pcm_recon[0].Vect();
      ////      double thetaCM = acos(pcm.Z()/
      ////sqrt(pcm.X()*pcm.X()+pcm.Y()*pcm.Y()+pcm.Z()*pcm.Z()));

      h_dthetaCM->Fill(thetaCM - kv.the_cm);
      h_dthetaCM_vs_thetaCM->Fill(kv.the_cm, thetaCM - kv.the_cm);
      ////      h_dthetaCM_vs_thetaCM->Fill(kv.the_cm, thetaCM);

      //      beta.Print();
      //      TLorentzVector pb = Pgstar_recon;
      //      pb.Boost(-beta);
      //      Pgstar_recon.Print();
      //      cout << Pgstar_recon.M() << endl;
      //      pb.Print();
      //      getchar();
      */

      //thetaCM reconstruction.

      TVector3 beta_gstar = Pgstar_recon.BoostVector();
      TLorentzVector P_gsb[3];
      for (int i=0; i<3; i++) {
	P_gsb[i] = P_recon[i];
	P_gsb[i].Boost(-beta_gstar);
      }

      double thetaCM = P_gsb[0].Angle(-P_gsb[2].Vect());

      //phiCM reconstruction.

      ////      TVector3 beta_gp = (Pg_recon + P_target).BoostVector();
      TVector3 beta_gp = (TLorentzVector(0.,0.,Pg_recon.E(),Pg_recon.E()) +
			  P_target).BoostVector();
      TLorentzVector P_gpb[3];
      for (int i=0; i<3; i++) {
	P_gpb[i] = P_recon[i];
	P_gpb[i].Boost(-beta_gp);
      }
      TLorentzVector Pt_gpb = P_target;
      Pt_gpb.Boost(-beta_gp);

      TVector3 ppp_prp = Pt_gpb.Vect().Cross(P_gpb[2].Vect());
      TVector3 pee_prp = (P_gpb[1].Vect()).Cross(P_gpb[0].Vect());
      double phiCM = acos((pee_prp*ppp_prp)/(pee_prp.Mag()*ppp_prp.Mag()));

      TVector3 px_prp = ppp_prp.Cross(P_gpb[2].Vect());
      if (pee_prp*px_prp < 0.)
	  phiCM = -phiCM +2.*TMath::Pi();

      //Fill histograms.

      bool good_event = Pg_recon.M2() > M2Pg_MIN && Pg_recon.M2() < M2Pg_MAX;

      if (good_event) {

	for (int i=0; i<3; i++) {
	  h_dphi[i]->Fill(phi_recon[i] - phi_vertex[i]);
	  h_dtheta[i]->Fill(theta_recon[i]-theta_vertex[i]);
	}

	h_dphi_gstar->Fill(Pgstar_recon.Px()/Pgstar_recon.Pz() -
			   Pgstar_vertex.Px()/Pgstar_vertex.Pz());
	h_dtheta_gstar->Fill(Pgstar_recon.Py()/Pgstar_recon.Pz() -
			     Pgstar_vertex.Py()/Pgstar_vertex.Pz());
	h_de_gstar->Fill(Pgstar_recon.E() - Pgstar_vertex.E());

	h_dpx_p->Fill(P_recon[2].Px() - P_vertex[2].Px());
	h_dpy_p->Fill(P_recon[2].Py() - P_vertex[2].Py());
	h_dpz_p->Fill(P_recon[2].Pz() - P_vertex[2].Pz());
	h_de_p->Fill (P_recon[2].E()  - P_vertex[2].E());

	double t_recon  = (P_target - P_recon[2]).M2();
	h_dt->Fill(t_recon - kv.t);
	h_dt_vs_t->Fill(-kv.t, t_recon - kv.t);

	double Q2_recon = Pgstar_recon.M2();
	h_dq2->Fill(Q2_recon - kv.Q2);
	h_dq2_vs_q2->Fill(kv.Q2, Q2_recon - kv.Q2);

	double s_recon  = (Pgstar_recon + P_recon[2]).M2();
	h_ds->Fill(s_recon - kv.s);
	h_ds_vs_s->Fill(kv.s, s_recon - kv.s);

	double xi_recon = Q2_recon/(2.*(s_recon - Mp*Mp) + t_recon - Q2_recon);
	h_dxi->Fill(xi_recon - kv.xi);
	h_dxi_vs_xi->Fill(kv.xi, xi_recon - kv.xi);

	double tau_recon = Q2_recon/(s_recon - Mp*Mp);
	h_dtau->Fill(tau_recon - kv.tau);
	h_dtau_vs_tau->Fill(kv.tau, tau_recon - kv.tau);

	double eta_recon = tau_recon/(tau_recon - 2.);
	h_deta->Fill(eta_recon - kv.eta);
	h_deta_vs_eta->Fill(-kv.eta, eta_recon - kv.eta);

	h_eg->Fill(Pg_recon.E());
	h_deg->Fill(Pg_recon.E() - kin.Eg);
	h_m2dpg->Fill((Pg_recon-Pg_vertex).M2());
	h_m2pg->Fill(Pg_recon.M2());

	h_dthetaCM->Fill(thetaCM - kv.the_cm);
	h_dthetaCM_vs_thetaCM->Fill(kv.the_cm, thetaCM - kv.the_cm);

	h_dphiCM->Fill(phiCM - kv.phi_cm);
	h_dphiCM_vs_phiCM->Fill(kv.phi_cm, phiCM - kv.phi_cm);
	h_phiCM_vs_phiCM->Fill(kv.phi_cm, phiCM);

      }   //good event

    }   //triple coin.

  }  //entries


  TCanvas* cphi = new TCanvas("cphi", "#phi residuals", 1400, 335);
  cphi->Divide(3,1);
  for (int i=0; i<3; i++) {
    cphi->cd(i+1);
    h_dphi[i]->Draw();
  }

  TCanvas* ctheta = new TCanvas("ctheta", "#theta residuals", 1400, 335);
  ctheta->Divide(3,1);
  for (int i=0; i<3; i++) {
    ctheta->cd(i+1);
    h_dtheta[i]->Draw();
  }

  TCanvas* cgstar = new TCanvas("cgstar", "gamma* residuals", 1400, 335);
  cgstar->Divide(3,1);
  cgstar->cd(1);
  h_dphi_gstar->Draw();
  cgstar->cd(2);
  h_dtheta_gstar->Draw();
  cgstar->cd(3);
  h_de_gstar->Draw();

  TCanvas* cproton = new TCanvas("cproton", "recoil proton residuals",
				 1400, 1000);
  cproton->Divide(2,2);
  cproton->cd(1);
  h_dpx_p->Draw();
  cproton->cd(2);
  h_dpy_p->Draw();
  cproton->cd(3);
  h_dpz_p->Draw();
  cproton->cd(4);
  h_de_p->Draw();

  TCanvas* cg = new TCanvas("cg", "#gamma reconstruction", 1000, 700);
  cg->Divide(2,2);
  cg->cd(1);
  h_eg->Draw();
  cg->cd(2);
  h_deg->Draw();
  cg->cd(3);
  h_m2pg->Draw();
  cg->cd(4);
  h_m2dpg->Draw();

  TCanvas* cq2 = new TCanvas("cq2", "Q'^{2} residuals", 1400, 500);
  cq2->Divide(2,1);
  cq2->cd(1);
  h_dq2->Draw();
  cq2->cd(2);
  h_dq2_vs_q2->Draw();

  TCanvas* ct = new TCanvas("ct", "-t residuals", 1400, 500);
  ct->Divide(2,1);
  ct->cd(1);
  h_dt->Draw();
  ct->cd(2);
  ////  h_dt_vs_t->Draw();
  TH1D* h_dt_prof = h_dt_vs_t->ProfileX("",1,-1,"s");   //standard deviation
  TH2D* h_dt_frame = new TH2D("frame",h_dt_vs_t->GetTitle(),
			      100, h_dt_vs_t->GetXaxis()->GetXmin(),
			      h_dt_vs_t->GetXaxis()->GetXmax(),
			      100, h_dt_vs_t->GetYaxis()->GetXmin(),
			      h_dt_vs_t->GetYaxis()->GetXmax());
  h_dt_frame->GetXaxis()->SetTitle(h_dt_vs_t->GetXaxis()->GetTitle());
  h_dt_frame->GetYaxis()->SetTitle(h_dt_vs_t->GetYaxis()->GetTitle());
  h_dt_frame->SetStats(0);
  h_dt_frame->Draw();
  //  gStyle->SetOptStat(0);
  h_dt_prof->Draw("same");

  TCanvas* cs = new TCanvas("cs", "s residuals", 1400, 500);
  cs->Divide(2,1);
  cs->cd(1);
  h_ds->Draw();
  cs->cd(2);
  h_ds_vs_s->Draw();

  TCanvas* cxi = new TCanvas("cxi", "#xi residuals", 1400, 500);
  cxi->Divide(2,1);
  cxi->cd(1);
  h_dxi->Draw();
  cxi->cd(2);
  h_dxi_vs_xi->Draw();

  TCanvas* ctau = new TCanvas("ctau", "#tau residuals", 1400, 500);
  ctau->Divide(2,1);
  ctau->cd(1);
  h_dtau->Draw();
  ctau->cd(2);
  h_dtau_vs_tau->Draw();

  TCanvas* ceta = new TCanvas("ceta", "#eta residuals", 1400, 500);
  ceta->Divide(2,1);
  ceta->cd(1);
  h_deta->Draw();
  ceta->cd(2);
  h_deta_vs_eta->Draw();

  TCanvas* cthetaCM = new TCanvas("cthetaCM", "#thetaCM residuals", 1400, 500);
  cthetaCM->Divide(2,1);
  cthetaCM->cd(1);
  h_dthetaCM->Draw();
  cthetaCM->cd(2);
  h_dthetaCM_vs_thetaCM->Draw();

  TCanvas* cphiCM = new TCanvas("cphiCM", "#phiCM residuals", 1400, 500);
  cphiCM->Divide(2,1);
  cphiCM->cd(1);
  h_dphiCM->Draw();
  cphiCM->cd(2);
  h_dphiCM_vs_phiCM->Draw();
  ////  h_phiCM_vs_phiCM->Draw();

  gStyle->SetTitleSize(0.06,"t"); 
  gStyle->SetLabelSize(0.045,"X");
  gStyle->SetLabelSize(0.045,"Y");
  gStyle->SetTitleSize(0.045,"X");
  gStyle->SetTitleSize(0.045,"Y");

  TCanvas* ckin = new TCanvas("ckin","Kin. quantities' residuals",1200,800);
  ckin->Divide(3,3);
  ckin->cd(1); h_dq2->Draw();
  ckin->cd(2); h_dt->Draw();
  ckin->cd(3); h_ds->Draw();
  ckin->cd(4); h_dxi->Draw();
  ckin->cd(5); h_deta->Draw();
  ckin->cd(6); h_dtau->Draw();
  ckin->cd(7); h_dthetaCM->Draw();
  ckin->cd(8); h_dphiCM->Draw();
  ckin->cd(9); h_deg->Draw();

}

//=============================================================================

double edep_calo(double p) {

  //HYCAL resolution, p in GeV/c.
  const double a=0.009;
  const double b=0.025;
  const double c=0.010;
  double sigma = sqrt(a*a+b*b/p+c*c/(p*p));

  return RandomGen.Gaus(p, sigma);
}

//------------------------------------------------------------------------------

double XSlopeReconP(double xslope_det, int quarter) {

  //Phi reconstruction for proton.

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

//------------------------------------------------------------------------------

double YSlopeReconE(const double yslope, const double pyz, int particle_id) {

  //average of quarters 0 and 3 for e-
  const double slope = 0.99545;
  const double offset_p0 = 0.00040304;
  const double offset_p1 = 0.18465;

  double offset = offset_p0 + offset_p1/pyz;

  double yslope_recon = 0.;
  switch (particle_id) {
  case  11: yslope_recon = (yslope-offset)/slope; break;
  case -11: yslope_recon = (yslope+offset)/slope; break;
  default: cout << "TCSTracker::YSlopeReconE: wrong particle_id = "
		<< particle_id << endl;
  }

  return yslope_recon;
}

//------------------------------------------------------------------------------

double YSlopeReconP(const double yslope, const double pyz, int quarter) {

  double sl_p0, sl_p1, of_p0, of_p1;
  switch (quarter) {
  case 0:
  case 1: sl_p0 = -1.662; sl_p1 = -2.114; of_p0 = -0.009478; of_p1 = -0.1704;
    break;
  case 2:
  case 3: sl_p0 = -2.312; sl_p1 = -2.236; of_p0 =  0.007202; of_p1 = -0.188;
    break;
  default: sl_p0 = 0.; sl_p1 = 0.; of_p0 =  0.; of_p1 = 0.;
    cout << "YSlopeReconP: wrong quarter = " << quarter << endl;
  }

  double offset = of_p0 + of_p1/pyz;
  double slope = 1./(1.+exp(sl_p0+sl_p1*pyz));

  return (yslope-offset)/slope;
}
