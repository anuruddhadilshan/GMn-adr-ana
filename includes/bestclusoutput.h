#include "elasticevent.h"
#include "besthcalclus.h"

#ifndef BESTCLUSOUTPUT_H
#define BESTCLUSOUTPUT_H

class BestClusOutput
{
private:

	ElasticEvent& m_event;
	BestHCalClus& m_bestHCalClus;
	const char* m_outputfilename;
	TFile* m_fout;
	TTree* m_resultstree;

	// Define the member variables to hold data for output TTree variables.
	double m_Q2 {0.};
	double m_W2 {0.};
	double m_predhcal_xpos {0.};
	double m_predhcal_ypos {0.};
	double m_bbtrvz {0.};
	// double m_bbtrth {0.};
	// double m_bbtrx {0.};
	double m_bbtrp{0.};
	// double m_bbtrpx{0.};
	// double m_bbtrpy{0.};
	double m_bbshe{0.};
	double m_bbpse{0.};
	double m_bbeoverp{0.};

	double m_sbshcalx{0.};
	double m_sbshcaly{0.};
	double m_sbshcale{0.};
	double m_sbshcalatimeblk{0.};
	// double m_bbphitgt{0.};
	// double m_bbpoltgt{0.};
	
	const static int m_MAXHCALCLUS {10};
	double* m_sbshcalcluse; // Ptr to "sbs.hcal.clus.e" array.
	double* m_sbshcalclusatime; // Ptr to "sbs.hcal.clus.atime" array.
	double* m_sbshcalclusx; // Ptr to "sbs.hcal.clus.x" array.
	double* m_sbshcalclusy; // Ptr to "sbs.hcal.clus.y" array.
	double* m_sbshcal_heclus_e; // Ptr to sorted sbs.hcal.clus.e array
	double* m_sbshcal_heclus_atime; // Ptr to sorted sbs.hcal.clus.atime array
	double* m_sbshcal_heclus_x; // Ptr to sorted sbs.hcal.clus.x array
	double* m_sbshcal_heclus_y; // Ptr to sorted sbs.hcal.clus.y array
	double m_shadctime {0.};
	int m_sbshcal_bestclus_indx {0}; // Index of the best cluster from Atime+HighestEnergy method.
	double m_sbshcal_bestclus_e {0.};
	double m_sbshcal_bestclus_atime {0.};
	double m_sbshcal_bestclus_x {0.};
	double m_sbshcal_bestclus_y {0.};

	// Analysis results output histograms.
	TH1D* m_h1_Q2;
	TH1D* m_h1_W2;
	TH1D* m_h1_predhcal_xpos;
	TH1D* m_h1_predhcal_ypos;
	TH2D* m_h2_predhcal_xy;
	TH1D* m_h1_heblkclus_dx;
	TH1D* m_h1_heclus_dx;
	TH1D* m_h1_heclus_1_dx;
	TH1D* m_h1_heclus_2_dx;
	TH1D* m_h1_bestclus_dx;
	TH2D* m_h2_bestclus_dxdy;
	TH1D* m_h1_bestclus_1_dx;
	TH1D* m_h1_bestclus_2_dx;
	TH1D* m_h1_bestclus_3_dx;
	TH1D* m_h1_bestclus_4_dx;
	//TH2D* m_h2_sbs_neutron_hcal_dxdy;

	// Diagnostics output histograms.
	TH1D* m_h1_bb_tr_vz;
	TH1D* m_h1_bb_tr_p;
	TH1D* m_h1_bb_ps_e;
	TH1D* m_h1_bb_eoverp;
	TH1D* m_h1_hcalheclu_sh_atimediff;

public:

	BestClusOutput(ElasticEvent& event, BestHCalClus& bestHCalClus, const char* outputfilename) : m_event{event}, m_bestHCalClus{bestHCalClus}, m_outputfilename{outputfilename},
	m_sbshcalcluse{m_event.return_HCalClusEArryPtr()}, m_sbshcalclusatime{m_event.return_HCalClusAtimeArryPtr()}, m_sbshcalclusx{m_event.return_HCalClusXArryPtr()}, m_sbshcalclusy{m_event.return_HCalClusYArryPtr()},
	m_sbshcal_heclus_e{m_event.return_HCalHEClusEArryPtr()}, m_sbshcal_heclus_atime{m_event.return_HCalHEClusAtimeArryPtr()}, m_sbshcal_heclus_x{m_event.return_HCalHEClusXArryPtr()}, m_sbshcal_heclus_y{m_event.return_HCalHEClusYArryPtr()}
	{
		m_fout = new TFile(Form("%s.root",outputfilename),"RECREATE");
		m_resultstree = new TTree("T","Best HCal cluster analysis");

		m_resultstree->Branch("ekine.Q2", &m_Q2);
		m_resultstree->Branch("ekine.W2", &m_W2);
		m_resultstree->Branch("bb.tr.vz", &m_bbtrvz);
		// m_resultstree->Branch("bb.tr.th", &m_bbtrth);
		// m_resultstree->Branch("bb.tr.x", &m_bbtrx);
		m_resultstree->Branch("bb.tr.p", &m_bbtrp);
		// m_resultstree->Branch("bb.tr.px", &m_bbtrpx);
		// m_resultstree->Branch("bb.tr.py", &m_bbtrpy);
		m_resultstree->Branch("bb.sh.e", &m_bbshe);
		m_resultstree->Branch("bb.ps.e", &m_bbpse);
		m_resultstree->Branch("bb.eoverp", &m_bbeoverp);
		m_resultstree->Branch("bb.sh.adctime", &m_shadctime);
		m_resultstree->Branch("sbs.hcal.pred.x", &m_predhcal_xpos);
		m_resultstree->Branch("sbs.hcal.pred.y", &m_predhcal_ypos);
		m_resultstree->Branch("sbs.hcal.x", &m_sbshcalx);
		m_resultstree->Branch("sbs.hcal.y", &m_sbshcaly);
		m_resultstree->Branch("sbs.hcal.e", &m_sbshcale);
		m_resultstree->Branch("sbs.hcal.atimeblk", &m_sbshcalatimeblk);
		// m_resultstree->Branch("bb.phitgt", &m_bbphitgt);
		// m_resultstree->Branch("bb.poltgt", &m_bbpoltgt);
		// m_resultstree->Branch();
		m_resultstree->Branch("sbs.hcal.clus.e", m_sbshcalcluse, "m_sbshcalcluse[10]/D");
		m_resultstree->Branch("sbs.hcal.clus.atime", m_sbshcalclusatime, "m_sbshcalclusatime[10]/D");
		m_resultstree->Branch("sbs.hcal.clus.x", m_sbshcalclusx, "m_sbshcalclusx[10]/D");
		m_resultstree->Branch("sbs.hcal.clus.y", m_sbshcalclusy, "m_sbshcalclusy[10]/D");
		m_resultstree->Branch("sbs.hcal.heclus.e", m_sbshcal_heclus_e, "m_sbshcal_heclus_e[10]/D");
		m_resultstree->Branch("sbs.hcal.heclus.atime", m_sbshcal_heclus_atime, "m_sbshcal_heclus_atime[10]/D");
		m_resultstree->Branch("sbs.hcal.heclus.x", m_sbshcal_heclus_x, "m_sbshcal_heclus_x[10]/D");
		m_resultstree->Branch("sbs.hcal.heclus.y", m_sbshcal_heclus_y, "m_sbshcal_heclus_y[10]/D");
		m_resultstree->Branch("sbs.hcal.bestclus.indx", &m_sbshcal_bestclus_indx);
		m_resultstree->Branch("sbs.hcal.bestclus.e", &m_sbshcal_bestclus_e);
		m_resultstree->Branch("sbs.hcal.bestclus.atime", &m_sbshcal_bestclus_atime);
		m_resultstree->Branch("sbs.hcal.bestclus.x", &m_sbshcal_bestclus_x);
		m_resultstree->Branch("sbs.hcal.bestclus.y", &m_sbshcal_bestclus_y);
		
		// Analysis results output histograms.
		m_h1_Q2 = new TH1D("h1_Q2", "Q^{2}/(GeV^{2}/c^{2})", 300,0, 15);
		m_h1_W2 = new TH1D("h1_W2", "W^{2}/(GeV^{2}/c^{2})", 200, 0, 5);
		m_h1_predhcal_xpos = new TH1D("h1_predhcal_xpos", "Neutron hypothesis reconstructed x hit position on HCal; X_{pos} (m)", 1000, -4, 6);
		m_h1_predhcal_ypos = new TH1D("h1_predhcal_ypos", "Neutron hypothesis reconstructed y hit position on HCal; Y_{pos} (m)", 500, -3, 2);
		m_h2_predhcal_xy = new TH2D("h2_predhcal_xy", "Neutron hypothesis reconstructed x vs y hit position on HCal; Y_{pos} (m); X_{pos} (m)", 500, -3, 2, 1000, -4, 6);
		m_h1_heblkclus_dx = new TH1D("h1_heblkclus_dx", "Highest Energy Block HCal Cluster dx; dx (m)", 800, -4, 4);
		m_h1_heclus_dx = new TH1D("h1_heclus_dx", "Highest Energy HCal Cluster dx; dx (m)", 800, -4, 4);
		m_h1_heclus_1_dx = new TH1D("h1_heclus_1_dx", "2nd Highest Energy HCal Cluster dx; dx (m)", 800, -4, 4);
		m_h1_heclus_2_dx = new TH1D("h1_heclus_2_dx", "3rd Highest Energy HCal Cluster dx; dx (m)", 800, -4, 4);
		m_h1_bestclus_dx = new TH1D("h1_bestclus_dx", "Best HCal Cluster from Timing+HE algorithm; dx (m)", 800, -4, 4);
		m_h1_bestclus_1_dx = new TH1D("h1_bestclus_1_dx", "Best HCal Cluster -1; dx (m)", 800, -4, 4);
		m_h1_bestclus_2_dx = new TH1D("h1_bestclus_2_dx", "Best HCal Cluster -2; dx (m)", 800, -4, 4);
		m_h1_bestclus_3_dx = new TH1D("h1_bestclus_3_dx", "Best HCal Cluster -3; dx (m)", 800, -4, 4);
		m_h1_bestclus_4_dx = new TH1D("h1_bestclus_4_dx", "Best HCal Cluster -4; dx (m)", 800, -4, 4);
		m_h2_bestclus_dxdy = new TH2D("h2_bestclus_dxdy", "HCal Best Clus dx vs dy; dy = Y_{hcal}-Y_{predicted} (m); dx = X_{hcal}-X_{predicted} (m)", 400, -2, 2, 800, -4, 4);

		// Diagnostic variables output histograms.
		m_h1_bb_tr_vz = new TH1D("h1_bb_tr_vz", "BB track vertex Z position distribution; vertex Z (m)", 600, -0.15, 0.15);
		// m_h1_hcal_e = new TH1D("h1_hcal_e","HCal Energy Deopsited; HCal Energy (GeV); Entries",250,0,0.5);
		m_h1_bb_tr_p = new TH1D("h1_bb_tr_p", "BB track momentum distribution; Momentum (GeV/c)", 1000, 0, 10.0);
		m_h1_bb_ps_e = new TH1D("h1_bb_ps_e", "BB pre-shower energy", 300,0, 3);
		m_h1_bb_eoverp = new TH1D("h1_bb_eoverp", "BigBite E/P distribution; E/P", 300, -1, 2);
		//m_h1_adctimediff_hcalsh = new TH1D("h1_adctimediff_hcalsh","HCal ADC time - SH ADC time; HCal_{ADCtime}-BBCal_{ADCtime} (ns); Entries",300,-100,200);
		m_h1_hcalheclu_sh_atimediff = new TH1D("h1_hcalheclus_sh_atimediff","HCal highest energy cluster time - SH time; HCal_{ADCtime}-BBCal_{ADCtime} (ns); Entries",3000,-100,200);
	}

	void copyFromEvent()
	{
		m_Q2 = m_event.return_Q2();
		m_W2 = m_event.return_W2();
		m_predhcal_xpos = m_event.return_nHypthsPredx();
		m_predhcal_ypos = m_event.return_nHypthsPredy();
		m_bbtrvz = m_event.return_BBTrVz();
		// m_bbtrth = m_event.return_BBTrth();
		// m_bbtrx = m_event.return_BBTrx();
		m_bbtrp = m_event.return_BBTrP();
		// m_bbtrpx = m_event.return_BBTrPx();
		// m_bbtrpy = m_event.return_BBtrPy();
		m_bbshe = m_event.return_BBSHe();
		m_bbpse = m_event.return_BBPSe();
		m_bbeoverp = m_event.return_EoverP();
		m_shadctime = m_event.return_SHADCTime();
		m_sbshcalx = m_event.return_SBSHCalx();
		m_sbshcaly = m_event.return_SBSHCaly();
		m_sbshcalatimeblk = m_event.return_SBSHCalAtimeBlk();
		m_sbshcale = m_event.return_SBSHCale();
		m_sbshcal_bestclus_indx = m_bestHCalClus.return_BestHCalClusIndx();
		m_sbshcal_bestclus_e = m_bestHCalClus.return_BestHCalClusE();
		m_sbshcal_bestclus_atime = m_bestHCalClus.return_BestHCalClusAtime();
		m_sbshcal_bestclus_x = m_bestHCalClus.return_BestHCalClusX();
		m_sbshcal_bestclus_y = m_bestHCalClus.return_BestHCalClusY();

		// m_adctimediffhcalsh = m_event.return_ADCTimeDiffHCalSH();
		// m_bbphitgt = m_event.return_BBphitgt();
		// m_bbpoltgt = m_event.return_BBpoltgt();		
		// m_clusepoint = m_event.return_HCalClusE();	
	}

	void fillOutTree()
	{
		m_resultstree->Fill();
	}

	void fillHistos()
	{
		// Analysis results histos.
		m_h1_Q2->Fill(m_Q2);
		m_h1_W2->Fill(m_W2);
		m_h1_predhcal_xpos->Fill(m_predhcal_xpos);
		m_h1_predhcal_ypos->Fill(m_predhcal_ypos);
		m_h2_predhcal_xy->Fill(m_predhcal_ypos, m_predhcal_xpos);

		m_h1_heblkclus_dx->Fill(m_sbshcalclusx[0] - m_predhcal_xpos);	
		m_h1_heclus_dx->Fill(m_sbshcal_heclus_x[0] - m_predhcal_xpos);
		m_h1_heclus_1_dx->Fill(m_sbshcal_heclus_x[1] - m_predhcal_xpos);
		m_h1_heclus_2_dx->Fill(m_sbshcal_heclus_x[2] - m_predhcal_xpos);
		m_h1_bestclus_dx->Fill(m_sbshcal_bestclus_x - m_predhcal_xpos);
		m_h2_bestclus_dxdy->Fill(m_sbshcal_bestclus_y - m_predhcal_ypos, m_sbshcal_bestclus_x - m_predhcal_xpos); 
		fillNextBestClusHistos();

		// Cut parameters histos.
		m_h1_bb_tr_vz->Fill(m_bbtrvz);
		// m_h1_adctimediff_hcalsh->Fill(m_adctimediffhcalsh);
		//m_h1_hcal_e->Fill(m_sbshcale);
		m_h1_bb_tr_p->Fill(m_bbtrp);
		m_h1_bb_ps_e->Fill(m_bbpse);
		m_h1_bb_eoverp->Fill(m_bbeoverp);
		m_h1_hcalheclu_sh_atimediff->Fill(m_sbshcal_heclus_atime[0] - m_shadctime);
	}

	void fillNextBestClusHistos()
	{
		switch ( m_sbshcal_bestclus_indx )
		{
			case 0:
			{
				m_h1_bestclus_1_dx->Fill(m_sbshcal_heclus_x[1] - m_predhcal_xpos);
				m_h1_bestclus_2_dx->Fill(m_sbshcal_heclus_x[2] - m_predhcal_xpos);
				m_h1_bestclus_3_dx->Fill(m_sbshcal_heclus_x[3] - m_predhcal_xpos);
				m_h1_bestclus_4_dx->Fill(m_sbshcal_heclus_x[4] - m_predhcal_xpos);
				return;
			}

			case 1:
			{
				m_h1_bestclus_1_dx->Fill(m_sbshcal_heclus_x[0] - m_predhcal_xpos);
				m_h1_bestclus_2_dx->Fill(m_sbshcal_heclus_x[2] - m_predhcal_xpos);
				m_h1_bestclus_3_dx->Fill(m_sbshcal_heclus_x[3] - m_predhcal_xpos);
				m_h1_bestclus_4_dx->Fill(m_sbshcal_heclus_x[4] - m_predhcal_xpos);
				return;
			}

			case 2:
			{
				m_h1_bestclus_1_dx->Fill(m_sbshcal_heclus_x[0] - m_predhcal_xpos);
				m_h1_bestclus_2_dx->Fill(m_sbshcal_heclus_x[1] - m_predhcal_xpos);
				m_h1_bestclus_3_dx->Fill(m_sbshcal_heclus_x[3] - m_predhcal_xpos);
				m_h1_bestclus_4_dx->Fill(m_sbshcal_heclus_x[4] - m_predhcal_xpos);
				return;
			}

			case 3:
			{
				m_h1_bestclus_1_dx->Fill(m_sbshcal_heclus_x[0] - m_predhcal_xpos);
				m_h1_bestclus_2_dx->Fill(m_sbshcal_heclus_x[1] - m_predhcal_xpos);
				m_h1_bestclus_3_dx->Fill(m_sbshcal_heclus_x[2] - m_predhcal_xpos);
				m_h1_bestclus_4_dx->Fill(m_sbshcal_heclus_x[4] - m_predhcal_xpos);
				return;
			}

			case 4:
			{
				m_h1_bestclus_1_dx->Fill(m_sbshcal_heclus_x[0] - m_predhcal_xpos);
				m_h1_bestclus_2_dx->Fill(m_sbshcal_heclus_x[1] - m_predhcal_xpos);
				m_h1_bestclus_3_dx->Fill(m_sbshcal_heclus_x[2] - m_predhcal_xpos);
				m_h1_bestclus_4_dx->Fill(m_sbshcal_heclus_x[3] - m_predhcal_xpos);
				return;
			}

			default:
			{
				m_h1_bestclus_1_dx->Fill(m_sbshcal_heclus_x[1] - m_predhcal_xpos);
				m_h1_bestclus_2_dx->Fill(m_sbshcal_heclus_x[2] - m_predhcal_xpos);
				m_h1_bestclus_3_dx->Fill(m_sbshcal_heclus_x[3] - m_predhcal_xpos);
				m_h1_bestclus_4_dx->Fill(m_sbshcal_heclus_x[4] - m_predhcal_xpos);
				return;	
			}
		}
	}

	void makePlots()
	{
		TCanvas* C1 = new TCanvas();
		m_h2_bestclus_dxdy->Draw("colz");

		TCanvas* C2 = new TCanvas();
		m_h1_heblkclus_dx->SetLineColor(kBlue);
		m_h1_heclus_dx->SetLineColor(kGreen);
		m_h1_bestclus_dx->SetLineColor(kRed);
		m_h1_bestclus_dx->Draw();
		m_h1_heclus_dx->Draw("same");
		m_h1_heblkclus_dx->Draw("same");

		TCanvas* C3 = new TCanvas();
		m_h1_bestclus_dx->Draw();
		m_h1_bestclus_1_dx->SetLineColor(kGreen);
		m_h1_bestclus_1_dx->Draw("same");
		m_h1_bestclus_2_dx->SetLineColor(kBlue);
		m_h1_bestclus_2_dx->Draw("same");
		m_h1_bestclus_3_dx->SetLineColor(kMagenta);
		m_h1_bestclus_3_dx->Draw("same");
		m_h1_bestclus_4_dx->SetLineColor(kOrange);
		m_h1_bestclus_4_dx->Draw("same");

		TCanvas* C4 = new TCanvas();
		m_h1_heclus_dx->SetLineColor(kRed);
		m_h1_heclus_1_dx->SetLineColor(kGreen);
		m_h1_heclus_2_dx->SetLineColor(kBlue);
		m_h1_heclus_dx->Draw();
		m_h1_heclus_1_dx->Draw("same");
 	}

	void closeOutFile()
	{
		m_resultstree->Write(0, TObject::kWriteDelete, 0);
		m_h1_Q2->Write();
		m_h1_W2->Write(); 
		m_h1_predhcal_xpos->Write();
		m_h1_predhcal_ypos->Write();
		m_h2_predhcal_xy->Write();

		m_h2_bestclus_dxdy->Write();
		m_h1_heblkclus_dx->Write();
		m_h1_heclus_dx->Write();
		m_h1_heclus_1_dx->Write();
		m_h1_heclus_2_dx->Write();
		m_h1_bestclus_dx->Write();
		m_h1_bestclus_1_dx->Write();
		m_h1_bestclus_2_dx->Write();
		m_h1_bestclus_3_dx->Write();
		m_h1_bestclus_4_dx->Write();		

		m_h1_bb_tr_vz->Write();
		//m_h1_adctimediff_hcalsh->Write();
		//m_h1_hcal_e->Write();
		m_h1_bb_tr_p->Write();
		m_h1_bb_ps_e->Write();
		m_h1_bb_eoverp->Write();
		m_h1_hcalheclu_sh_atimediff->Write();

		delete[] m_sbshcalcluse;
		delete[] m_sbshcalclusatime;
		delete[] m_sbshcalclusx;
		delete[] m_sbshcalclusy;
		delete[] m_sbshcal_heclus_e;
		delete[] m_sbshcal_heclus_atime;
		delete[] m_sbshcal_heclus_x;
		delete[] m_sbshcal_heclus_y;
	}



// private:

// 	// Define TCanvases to print output analysis histograms and make a PDF file.
// 	static const int m_nanacanvas{10};
// 	TCanvas* m_anaCan[m_nanacanvas];

// 	// Define TCanvases to print cut parameter output histograms.
// 	static const int m_ncutcanvas{6};
// 	TCanvas* m_cutCan[m_ncutcanvas];

// public:

// 	void make_anapdf()
// 	{
		
// 		m_anaCan[0] = new TCanvas("Track dx/dz distribution");;
// 		m_h1_bb_tr_th->Draw();

// 		m_anaCan[1] = new TCanvas("dx/dz vs track x pos distribution");
// 		m_h2_bb_tr_th_vs_x->Draw("COLZ");

// 		m_anaCan[2] = new TCanvas("Azimuthal scattering (#phi) angle distribution");
// 		m_h1_phi_tgt->Draw();

// 		m_anaCan[3] = new TCanvas("Pola scattering (#theta) angle distribution");
// 		m_h1_pol_tgt->Draw();

// 		m_anaCan[4] = new TCanvas("HCal cluster vertical pos vs #phi distribution");
// 		m_h2_hcalclusX_vs_phi->Draw("COLZ");

// 		m_anaCan[5] = new TCanvas("Reconstructed photon energy distribution");
// 		m_h1_photon_e->Draw();

// 		m_anaCan[6] = new TCanvas("Reconstructed neutron X position on HCal");
// 		m_h1_sbs_neutron_hcal_xpos->Draw();

// 		m_anaCan[7] = new TCanvas("Reconstructed neutron Y position on HCal");
// 		m_h1_sbs_neutron_hcal_ypos->Draw();

// 		m_anaCan[8] = new TCanvas("Reconstructed neutron Y vs X position on HCal");
// 		m_h2_sbs_neutron_hcal_xy->Draw("COLZ");

// 		m_anaCan[9] = new TCanvas("HCal dx vs dy");
// 		m_h2_sbs_neutron_hcal_dxdy->Draw("COLZ");

// 		TString pdffilename = Form("%s_anahistos.pdf", m_outputfilename);
// 		TString openfilename = pdffilename+"(";
// 		TString closefilename = pdffilename+")";

// 		double lmargin=0.15;
// 	  	double rmargin=0.15;
// 	    double bmargin=0.15;
// 	    double tmargin=0.09;

// 	    for (int icanvas = 0; icanvas < m_nanacanvas; icanvas++)
// 		{
// 			if(icanvas == 0) m_anaCan[icanvas]->Print(openfilename);
// 			else	if (icanvas == m_nanacanvas-1) m_anaCan[icanvas]->Print(closefilename);
// 			else m_anaCan[icanvas]->Print(pdffilename);
// 		}
// 	}

	// void make_cutpdf()
	// {
		
	// 	m_cutCan[0] = new TCanvas("BB track vertex distribution");;
	// 	m_h1_bb_tr_vz->Draw();

	// 	m_cutCan[1] = new TCanvas("HCal and BBCal ADC time difference");
	// 	m_h1_adctimediff_hcalsh->Draw();

	// 	m_cutCan[2] = new TCanvas("HCal cluster energy distribution");
	// 	m_h1_hcal_e->Draw();

	// 	m_cutCan[3] = new TCanvas("BB track momentum distribution");
	// 	m_h1_bb_tr_p->Draw();

	// 	m_cutCan[4] = new TCanvas("SH cluster energy + PS cluster energy distribution");
	// 	m_h1_bb_shps_e->Draw();

	// 	m_cutCan[5] = new TCanvas("BigBite E/P distribution");
	// 	m_h1_bb_eoverp->Draw();

	// 	TString pdffilename = Form("%s_cuthistos.pdf", m_outputfilename);
	// 	TString openfilename = pdffilename+"(";
	// 	TString closefilename = pdffilename+")";

	// 	double lmargin=0.15;
	//   	double rmargin=0.15;
	//     double bmargin=0.15;
	//     double tmargin=0.09;

	//     for (int icanvas = 0; icanvas < m_ncutcanvas; icanvas++)
	// 	{
	// 		if(icanvas == 0) m_cutCan[icanvas]->Print(openfilename);
	// 		else	if (icanvas == m_ncutcanvas-1) m_cutCan[icanvas]->Print(closefilename);
	// 		else m_cutCan[icanvas]->Print(pdffilename);
	// 	}
	// }
};

#endif