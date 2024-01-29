// Script to plot the fiducial regions.

#include "../includes/HCalConstants.h"

void plot_fiducial( const char* rootfile )
{
	// Load the root file and get the TTree "T" created from elastic event selection analysis scripts.
	TFile* anarootfile = new TFile(rootfile, "READ");
	TTree* T = (TTree*)anarootfile->Get("T"); 

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("e.kine.W2", 1);
	T->SetBranchStatus("sbs.hcal.pred.x", 1);
	T->SetBranchStatus("sbs.hcal.pred.y", 1);
	T->SetBranchStatus("cut.passfiducial", 1);
	T->SetBranchStatus("sbs.hcal.bestclus.x", 1);
	T->SetBranchStatus("sbs.hcal.bestclus.y", 1);
	
	double x {0.};
	double y {0.};
	bool pass_fiducial {false};
	double det_x {0.};
	double det_y {0.};
	double w2 {0.};

	T->SetBranchAddress("sbs.hcal.pred.x", &x);
	T->SetBranchAddress("sbs.hcal.pred.y", &y);
	T->SetBranchAddress("cut.passfiducial", &pass_fiducial);
	T->SetBranchAddress("sbs.hcal.bestclus.x", &det_x);
	T->SetBranchAddress("sbs.hcal.bestclus.y", &det_y);
	T->SetBranchAddress("e.kine.W2", &w2);

	// Histogram definitions //
	TH2D* h2_neutron_xy_all = new TH2D("h2_neutron_xy_all", "Neutron hyp. - ALL events; HCal Y (m); HCal X (m)", 200, -2.0, 2.0, 400, -4, 4);
	TH2D* h2_neutron_xy_fidpass = new TH2D("h2_xy_neutron_fidpass", "Neutron hyp. - events PASSED Fiducial Cut; HCal Y (m); HCal X (m)", 200, -2.0, 2.0, 400, -4, 4);
	TH2D* h2_neutron_xy_fidfail = new TH2D("h2_xy_neutron_fidfail", "Neutron hyp. - events FAILED Fiducial Cut; HCal Y (m); HCal X (m)", 200, -2.0, 2.0, 400, -4, 4);
	TH2D* h2_proton_xy_all = new TH2D("h2_proton_xy_all", "Proton hyp. - ALL events; HCal Y (m); HCal X (m)", 200, -2.0, 2.0, 400, -4, 4);
	TH2D* h2_proton_xy_fidpass = new TH2D("h2_xy_proton_fidpass", "Proton hyp. - events PASSED Fiducial Cut; HCal Y (m); HCal X (m)", 200, -2.0, 2.0, 400, -4, 4);
	TH2D* h2_proton_xy_fidfail = new TH2D("h2_xy_proton_fidfail", "Proton hyp. - events FAILED Fiducial Cut; HCal Y (m); HCal X (m)", 200, -2.0, 2.0, 400, -4, 4);
	TH2D* h2_detected_xy = new TH2D("h2_detected_xy", "Detected HCal X vs Y; HCal Y (m); HCal X (m)", 200, -2.0, 2.0, 400, -4, 4);
	const double avg_proton_deflection = -0.758239;
	const double stddev_protondis = 0.199934;


	long nevent {0};

	while ( T->GetEntry(nevent++) )
	{
		if ( w2 < 0.5 || w2 > 1.1 ) continue;
		
		h2_neutron_xy_all->Fill(y, x);
		h2_proton_xy_all->Fill(y, x + avg_proton_deflection);
		
		if ( pass_fiducial ) 
		{
			h2_neutron_xy_fidpass->Fill(y, x);
			h2_proton_xy_fidpass->Fill(y, x + avg_proton_deflection);
		}
		else
		{
			h2_neutron_xy_fidfail->Fill(y, x);
			h2_proton_xy_fidfail->Fill(y, x + avg_proton_deflection);
		}
	
		h2_detected_xy->Fill(det_y, det_x);		
	}

	TLine* l_hcalouter_xlow = new TLine( HCalConst::hcal_active_ylow, HCalConst::hcal_active_xlow, HCalConst::hcal_active_yhigh, HCalConst::hcal_active_xlow );
	TLine* l_hcalouter_xhigh = new TLine( HCalConst::hcal_active_ylow, HCalConst::hcal_active_xhigh, HCalConst::hcal_active_yhigh, HCalConst::hcal_active_xhigh );
	TLine* l_hcalouter_ylow = new TLine( HCalConst::hcal_active_ylow, HCalConst::hcal_active_xlow, HCalConst::hcal_active_ylow, HCalConst::hcal_active_xhigh );
	TLine* l_hcalouter_yhigh = new TLine( HCalConst::hcal_active_yhigh, HCalConst::hcal_active_xlow, HCalConst::hcal_active_yhigh, HCalConst::hcal_active_xhigh );


	TLine* l_hcalactive_xlow = new TLine( HCalConst::hcal_active_ylow_safe_pass2, HCalConst::hcal_active_xlow_safe_pass2, HCalConst::hcal_active_yhigh_safe_pass2, HCalConst::hcal_active_xlow_safe_pass2 );
	TLine* l_hcalactive_xhigh = new TLine( HCalConst::hcal_active_ylow_safe_pass2, HCalConst::hcal_active_xhigh_safe_pass2, HCalConst::hcal_active_yhigh_safe_pass2, HCalConst::hcal_active_xhigh_safe_pass2 );
	TLine* l_hcalactive_ylow = new TLine( HCalConst::hcal_active_ylow_safe_pass2, HCalConst::hcal_active_xlow_safe_pass2, HCalConst::hcal_active_ylow_safe_pass2, HCalConst::hcal_active_xhigh_safe_pass2 );
	TLine* l_hcalactive_yhigh = new TLine( HCalConst::hcal_active_yhigh_safe_pass2, HCalConst::hcal_active_xlow_safe_pass2, HCalConst::hcal_active_yhigh_safe_pass2, HCalConst::hcal_active_xhigh_safe_pass2 );
	

	TCanvas* C1 = new TCanvas("C1", "C1", 800, 800);
	C1->Divide(2,1);
	C1->cd(1);
	h2_neutron_xy_all->Draw("COLZ");
	l_hcalouter_xlow->Draw();
	l_hcalouter_xhigh->Draw();
	l_hcalouter_ylow->Draw();
	l_hcalouter_yhigh->Draw();
	l_hcalactive_xlow->Draw();
	l_hcalactive_xhigh->Draw();
	l_hcalactive_ylow->Draw();
	l_hcalactive_yhigh->Draw();
	C1->cd(2);
	h2_proton_xy_all->Draw("COLZ");
	l_hcalouter_xlow->Draw();
	l_hcalouter_xhigh->Draw();
	l_hcalouter_ylow->Draw();
	l_hcalouter_yhigh->Draw();
	l_hcalactive_xlow->Draw();
	l_hcalactive_xhigh->Draw();
	l_hcalactive_ylow->Draw();
	l_hcalactive_yhigh->Draw();


	TCanvas* C2 = new TCanvas("C2", "C2", 800, 800);
	C2->Divide(2,1);
	C2->cd(1);
	h2_neutron_xy_fidpass->Draw("COLZ");
	l_hcalouter_xlow->Draw();
	l_hcalouter_xhigh->Draw();
	l_hcalouter_ylow->Draw();
	l_hcalouter_yhigh->Draw();
	l_hcalactive_xlow->Draw();
	l_hcalactive_xhigh->Draw();
	l_hcalactive_ylow->Draw();
	l_hcalactive_yhigh->Draw();
	C2->cd(2);
	h2_proton_xy_fidpass->Draw("COLZ");
	l_hcalouter_xlow->Draw();
	l_hcalouter_xhigh->Draw();
	l_hcalouter_ylow->Draw();
	l_hcalouter_yhigh->Draw();
	l_hcalactive_xlow->Draw();
	l_hcalactive_xhigh->Draw();
	l_hcalactive_ylow->Draw();
	l_hcalactive_yhigh->Draw();

	TCanvas* C3 = new TCanvas("C3", "C3", 800, 800);
	C3->Divide(2,1);
	C3->cd(1);
	h2_neutron_xy_fidfail->Draw("COLZ");
	l_hcalouter_xlow->Draw();
	l_hcalouter_xhigh->Draw();
	l_hcalouter_ylow->Draw();
	l_hcalouter_yhigh->Draw();
	l_hcalactive_xlow->Draw();
	l_hcalactive_xhigh->Draw();
	l_hcalactive_ylow->Draw();
	l_hcalactive_yhigh->Draw();
	C3->cd(2);
	h2_proton_xy_fidfail->Draw("COLZ");
	l_hcalouter_xlow->Draw();
	l_hcalouter_xhigh->Draw();
	l_hcalouter_ylow->Draw();
	l_hcalouter_yhigh->Draw();
	l_hcalactive_xlow->Draw();
	l_hcalactive_xhigh->Draw();
	l_hcalactive_ylow->Draw();
	l_hcalactive_yhigh->Draw();

	TCanvas* C4 = new TCanvas("C4", "C4", 400, 800);
	h2_detected_xy->Draw("COLZ");
	l_hcalouter_xlow->Draw();
	l_hcalouter_xhigh->Draw();
	l_hcalouter_ylow->Draw();
	l_hcalouter_yhigh->Draw();
	l_hcalactive_xlow->Draw();
	l_hcalactive_xhigh->Draw();
	l_hcalactive_ylow->Draw();
	l_hcalactive_yhigh->Draw();

}