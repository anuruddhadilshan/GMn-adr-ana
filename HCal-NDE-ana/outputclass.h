#include "eventclass.h"

#ifndef RESULTSCLASS_H
#define RESULTSCLASS_H

class Output
{
private:

	Event& m_event;
	const char* m_outputfilename;
	TFile* m_fout;
	TTree* m_resultstree;

	// Define the member variables to hold data for output TTree variables.
	double m_bbtrvz{0.};
	double m_bbtrth{0.};
	double m_bbtrx{0.};
	double m_bbtrp{0.};
	double m_bbtrpx{0.};
	double m_bbtrpy{0.};
	double m_bbshe{0.};
	double m_bbpse{0.};
	double m_bbshpse{0.};
	double m_bbeoverp{0.};
	double m_sbshcalx{0.};
	double m_sbshcaly{0.};
	double m_sbshcale{0.};
	double m_adctimediffhcalsh{0.};
	double m_bbphitgt{0.};
	double m_bbpoltgt{0.};
	double m_photone{0.};
	double m_sbsneutronpx{0.};
	double m_sbsneutronpy{0.};
	double m_sbsneutronpz{0.};
	double m_sbsneutrone{0.};
	double m_sbsneutronhcalx{0.};
	double m_sbsneutronhcaly{0.};
	double m_sbsneutronhcaldx{0.};
	double m_sbsneutronhcaldy{0.};
	bool m_didpasshcalshadctimecut{false};
	bool m_didpassrcut{false};
	bool m_hcalndermthdpass{false};

	// Analysis results output histograms.
	TH1D* m_h1_bb_tr_th;
	TH2D* m_h2_bb_tr_th_vs_x;
	TH1D* m_h1_phi_tgt;
	TH1D* m_h1_pol_tgt;
	TH2D* m_h2_hcalclusX_vs_phi;
	TH1D* m_h1_photon_e;
	TH1D* m_h1_sbs_neutron_hcal_xpos;
	TH1D* m_h1_sbs_neutron_hcal_ypos;
	TH2D* m_h2_sbs_neutron_hcal_xy;
	TH2D* m_h2_sbs_neutron_hcal_dxdy;

	// Cut parameters output histograms.
	TH1D* m_h1_bb_tr_vz;
	TH1D* m_h1_adctimediff_hcalsh;
	TH1D* m_h1_hcal_e;
	TH1D* m_h1_bb_tr_p;
	TH1D* m_h1_bb_shps_e;
	TH1D* m_h1_bb_eoverp;

public:

	Output(Event& event, const char* outputfilename) : m_event{event}, m_outputfilename{outputfilename}
	{
		m_fout = new TFile(Form("%s.root",outputfilename),"RECREATE");
		m_resultstree = new TTree("T","Down beding track analysis");

		m_resultstree->Branch("bb.tr.vz", &m_bbtrvz);
		m_resultstree->Branch("bb.tr.th", &m_bbtrth);
		m_resultstree->Branch("bb.tr.x", &m_bbtrx);
		m_resultstree->Branch("bb.tr.p", &m_bbtrp);
		m_resultstree->Branch("bb.tr.px", &m_bbtrpx);
		m_resultstree->Branch("bb.tr.py", &m_bbtrpy);
		m_resultstree->Branch("bb.sh.e", &m_bbshe);
		m_resultstree->Branch("bb.ps.e", &m_bbpse);
		m_resultstree->Branch("bb.shps.e", &m_bbshpse);
		m_resultstree->Branch("bb.eoverp", &m_bbeoverp);
		m_resultstree->Branch("sbs.hcal.x", &m_sbshcalx);
		m_resultstree->Branch("sbs.hcal.y", &m_sbshcaly);
		m_resultstree->Branch("sbs.hcal.e", &m_sbshcale);
		m_resultstree->Branch("adctimediff.hcalsh", &m_adctimediffhcalsh);
		m_resultstree->Branch("bb.phitgt", &m_bbphitgt);
		m_resultstree->Branch("bb.poltgt", &m_bbpoltgt);
		m_resultstree->Branch("beam.photon.e", &m_photone);
		m_resultstree->Branch("sbs.neutronrecon.px", &m_sbsneutronpx);
		m_resultstree->Branch("sbs.neutronrecon.py", &m_sbsneutronpy);
		m_resultstree->Branch("sbs.neutronrecon.pz", &m_sbsneutronpz);
		m_resultstree->Branch("sbs.neutronrecon.e", &m_sbsneutrone);
		m_resultstree->Branch("sbs.neutronrecon.hcalx", &m_sbsneutronhcalx);
		m_resultstree->Branch("sbs.neutronrecon.hcaly", &m_sbsneutronhcaly);
		m_resultstree->Branch("sbs.neutronrecon.hcaldx", &m_sbsneutronhcaldx);
		m_resultstree->Branch("sbs.neutronrecon.hcaldy", &m_sbsneutronhcaldy);
		m_resultstree->Branch("adctimecorr.didpass", &m_didpasshcalshadctimecut);
		m_resultstree->Branch("hcalnde.rcutdidpass", &m_didpassrcut);
		m_resultstree->Branch("hcalnde.rmthddidpass", &m_hcalndermthdpass);

		//Define the output analysis results histograms.
		m_h1_bb_tr_th = new TH1D("h1_bb_tr_th0", "BB track theta (dx/dz) distribution; dx/dz", 100, -0.2, 0.8);
		m_h2_bb_tr_th_vs_x = new TH2D("h2_bb_tr_th0_vs_x", "BB track dx/dz vs x distribution; x (m); dx/dz", 160, -0.8, 0.8, 1200, -0.6, 0.6);
		m_h1_phi_tgt = new TH1D("h1_phi_tgt", "Azimuthal scattering angle (#phi) distribution; #phi degrees", 700, -35, 35);
		m_h1_pol_tgt = new TH1D("h1_pol_tgt", "Polar scattering angle (#theta) distribution; #theta degrees", 900, 0, 90);
		m_h2_hcalclusX_vs_phi = new TH2D("h2_hcalclusX_vs_phi", "HCal cluster vertical pos vs #phi distribution; #phi degrees; sbs.hcal.x (m)", 450, -10, 35, 500, -3.0, 2.0);
		m_h1_photon_e = new TH1D("h1_photon_e", "Reconstructed photon energy; #gamma Energy (GeV)", 1000, 0, 10);
		m_h1_sbs_neutron_hcal_xpos = new TH1D("h1_sbs_neutron_hcal_xpos", "Reconstructed x hit position of the neutron on HCal; X_{pos} (m)", 1000, -4, 6);
		m_h1_sbs_neutron_hcal_ypos = new TH1D("h1_sbs_neutron_hcal_ypos", "Reconstructed y hit position of the neutron on HCal; Y_{pos} (m)", 500, -3, 2);
		m_h2_sbs_neutron_hcal_xy = new TH2D("h2_sbs_neutron_hcal_xy", "Reconstructed y vs x hit position of the neutron on HCal; Y_{pos} (m); X_{pos} (m)", 500, -3, 2, 1000, -4, 6);
		m_h2_sbs_neutron_hcal_dxdy = new TH2D("h2_sbs_neutron_hcal_dxdy", "HCal dx vs dy; dy = Y_{hcal}-Y_{predicted} (m); dx = X_{hcal}-X_{predicted} (m)", 400, -2, 2, 800, -4, 4);

		//Define the output cut parameter histograms.
		m_h1_bb_tr_vz = new TH1D("h1_bb_tr_vz", "BB track vertex Z position distribution; vertex Z (m)", 600, -0.15, 0.15);
		m_h1_adctimediff_hcalsh = new TH1D("h1_adctimediff_hcalsh","HCal ADC time - SH ADC time; HCal_{ADCtime}-BBCal_{ADCtime} (ns); Entries",300,-100,200);
		m_h1_hcal_e = new TH1D("h1_hcal_e","HCal Energy Deopsited; HCal Energy (GeV); Entries",250,0,0.5);
		m_h1_bb_tr_p = new TH1D("h1_bb_tr_p", "BB track momentum distribution; Momentum (GeV/c)", 1000, 0, 10.0);
		m_h1_bb_shps_e = new TH1D("h1_bb_shps_e","Total SH and PS cluster energy sum; SH+PS Energy (GeV)' Entries",500,0,5);	
		m_h1_bb_eoverp = new TH1D("h1_bb_eoverp", "BigBite E/P distribution; E/P", 1000, -5, 5);
	}

	void copyFromEvent()
	{
		m_bbtrvz = m_event.return_BBTrVz();
		m_bbtrth = m_event.return_BBTrth();
		m_bbtrx = m_event.return_BBTrx();
		m_bbtrp = m_event.return_BBTrP();
		m_bbtrpx = m_event.return_BBTrPx();
		m_bbtrpy = m_event.return_BBtrPy();
		m_bbshe = m_event.return_BBSHe();
		m_bbpse = m_event.return_BBPSe();
		m_bbshpse = m_event.return_BBSHPSe();
		m_bbeoverp = m_event.return_EoverP();
		m_sbshcalx = m_event.return_SBSHCalx();
		m_sbshcaly = m_event.return_SBSHCaly();
		m_sbshcale = m_event.return_SBSHCale();
		m_adctimediffhcalsh = m_event.return_ADCTimeDiffHCalSH();
		m_bbphitgt = m_event.return_BBphitgt();
		m_bbpoltgt = m_event.return_BBpoltgt();
		m_photone = m_event.return_PhotonE();
		m_sbsneutronpx = m_event.return_NeutronPx();
		m_sbsneutronpy = m_event.return_NeutronPy();
		m_sbsneutronpz = m_event.return_NeutronPz();
		m_sbsneutrone = m_event.return_NeutronE();
		m_sbsneutronhcalx = m_event.return_HCalNeutronPredXpos();
		m_sbsneutronhcaly = m_event.return_HCalNeutronPredYpos(); 
		m_sbsneutronhcaldx = m_event.return_HCalNeutrondx();
		m_sbsneutronhcaldy = m_event.return_HCalNeutrondy();
		m_didpasshcalshadctimecut = m_event.return_DidPassHCalBBSH_ADCtime_Corr();
		m_didpassrcut = m_event.return_DidHitInCircle();
		m_hcalndermthdpass = m_event.return_DidPassRmthd();
	}

	void fillOutTree()
	{
		m_resultstree->Fill();
	}

	void fillHistos()
	{
		// Analysis results histos.
		m_h1_bb_tr_th->Fill(m_bbtrth);
		m_h2_bb_tr_th_vs_x->Fill(m_bbtrx, m_bbtrth);
		m_h1_phi_tgt->Fill(m_bbphitgt);
		m_h1_pol_tgt->Fill(m_bbpoltgt);
		m_h2_hcalclusX_vs_phi->Fill(m_bbphitgt, m_sbshcalx);
		m_h1_photon_e->Fill(m_photone);
		m_h1_sbs_neutron_hcal_xpos->Fill(m_sbsneutronhcalx);
		m_h1_sbs_neutron_hcal_ypos->Fill(m_sbsneutronhcaly);
		m_h2_sbs_neutron_hcal_xy->Fill(m_sbsneutronhcaly, m_sbsneutronhcalx);
		m_h2_sbs_neutron_hcal_dxdy->Fill(m_sbsneutronhcaldy, m_sbsneutronhcaldx);

		// Cut parameters histos.
		m_h1_bb_tr_vz->Fill(m_bbtrvz);
		m_h1_adctimediff_hcalsh->Fill(m_adctimediffhcalsh);
		m_h1_hcal_e->Fill(m_sbshcale);
		m_h1_bb_tr_p->Fill(m_bbtrp);
		m_h1_bb_shps_e->Fill(m_bbshpse);
		m_h1_bb_eoverp->Fill(m_bbeoverp);
	}

	void closeOutFile()
	{
		m_resultstree->Write(0, TObject::kWriteDelete, 0);
		m_h1_bb_tr_th->Write();
		m_h2_bb_tr_th_vs_x->Write();
		m_h1_phi_tgt->Write();
		m_h1_pol_tgt->Write();
		m_h2_hcalclusX_vs_phi->Write();
		m_h1_photon_e->Write();
		m_h1_sbs_neutron_hcal_xpos->Write();
		m_h1_sbs_neutron_hcal_ypos->Write();
		m_h2_sbs_neutron_hcal_xy->Write();
		m_h2_sbs_neutron_hcal_dxdy->Write();

		m_h1_bb_tr_vz->Write();
		m_h1_adctimediff_hcalsh->Write();
		m_h1_hcal_e->Write();
		m_h1_bb_tr_p->Write();
		m_h1_bb_shps_e->Write();
	}

private:

	// Define TCanvases to print output analysis histograms and make a PDF file.
	static const int m_nanacanvas{10};
	TCanvas* m_anaCan[m_nanacanvas];

	// Define TCanvases to print cut parameter output histograms.
	static const int m_ncutcanvas{6};
	TCanvas* m_cutCan[m_ncutcanvas];

public:

	void make_anapdf()
	{
		
		m_anaCan[0] = new TCanvas("Track dx/dz distribution");;
		m_h1_bb_tr_th->Draw();

		m_anaCan[1] = new TCanvas("dx/dz vs track x pos distribution");
		m_h2_bb_tr_th_vs_x->Draw("COLZ");

		m_anaCan[2] = new TCanvas("Azimuthal scattering (#phi) angle distribution");
		m_h1_phi_tgt->Draw();

		m_anaCan[3] = new TCanvas("Pola scattering (#theta) angle distribution");
		m_h1_pol_tgt->Draw();

		m_anaCan[4] = new TCanvas("HCal cluster vertical pos vs #phi distribution");
		m_h2_hcalclusX_vs_phi->Draw("COLZ");

		m_anaCan[5] = new TCanvas("Reconstructed photon energy distribution");
		m_h1_photon_e->Draw();

		m_anaCan[6] = new TCanvas("Reconstructed neutron X position on HCal");
		m_h1_sbs_neutron_hcal_xpos->Draw();

		m_anaCan[7] = new TCanvas("Reconstructed neutron Y position on HCal");
		m_h1_sbs_neutron_hcal_ypos->Draw();

		m_anaCan[8] = new TCanvas("Reconstructed neutron Y vs X position on HCal");
		m_h2_sbs_neutron_hcal_xy->Draw("COLZ");

		m_anaCan[9] = new TCanvas("HCal dx vs dy");
		m_h2_sbs_neutron_hcal_dxdy->Draw("COLZ");

		TString pdffilename = Form("%s_anahistos.pdf", m_outputfilename);
		TString openfilename = pdffilename+"(";
		TString closefilename = pdffilename+")";

		double lmargin=0.15;
	  	double rmargin=0.15;
	    double bmargin=0.15;
	    double tmargin=0.09;

	    for (int icanvas = 0; icanvas < m_nanacanvas; icanvas++)
		{
			if(icanvas == 0) m_anaCan[icanvas]->Print(openfilename);
			else	if (icanvas == m_nanacanvas-1) m_anaCan[icanvas]->Print(closefilename);
			else m_anaCan[icanvas]->Print(pdffilename);
		}
	}

	void make_cutpdf()
	{
		
		m_cutCan[0] = new TCanvas("BB track vertex distribution");;
		m_h1_bb_tr_vz->Draw();

		m_cutCan[1] = new TCanvas("HCal and BBCal ADC time difference");
		m_h1_adctimediff_hcalsh->Draw();

		m_cutCan[2] = new TCanvas("HCal cluster energy distribution");
		m_h1_hcal_e->Draw();

		m_cutCan[3] = new TCanvas("BB track momentum distribution");
		m_h1_bb_tr_p->Draw();

		m_cutCan[4] = new TCanvas("SH cluster energy + PS cluster energy distribution");
		m_h1_bb_shps_e->Draw();

		m_cutCan[5] = new TCanvas("BigBite E/P distribution");
		m_h1_bb_eoverp->Draw();

		TString pdffilename = Form("%s_cuthistos.pdf", m_outputfilename);
		TString openfilename = pdffilename+"(";
		TString closefilename = pdffilename+")";

		double lmargin=0.15;
	  	double rmargin=0.15;
	    double bmargin=0.15;
	    double tmargin=0.09;

	    for (int icanvas = 0; icanvas < m_ncutcanvas; icanvas++)
		{
			if(icanvas == 0) m_cutCan[icanvas]->Print(openfilename);
			else	if (icanvas == m_ncutcanvas-1) m_cutCan[icanvas]->Print(closefilename);
			else m_cutCan[icanvas]->Print(pdffilename);
		}
	}
};

#endif