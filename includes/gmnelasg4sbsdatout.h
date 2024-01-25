#include "elasticeventg4sbs.h"
#include "besthcalclusforg4sbs.h"

#ifndef GMNELASG4SBSDATOUT_H
#define GMNELASG4SBSDATOUT_H

class OutputG4SBS
{
private:

	ElasticEventG4SBS& m_event;
	BestHCalClusForG4SBS& m_bestHCalClus;
	const TString m_outputdirpath;
	const char* m_outputfilename;
	const bool m_is_simc;
	TFile* m_fout;
	TTree* m_resultstree;
	TTree* m_simcinfotree;

	// Define the member variables to hold data for output TTree variables.
	double m_Q2 {0.};
	double m_W2 {0.};
	double m_mc_weight {0.}; // Event weight for the G4SBS generator.
	double m_simc_weight {0.}; // Cross section & spectral function weight (ub/MeV/sr^2) for SIMC generator.
	double m_simc_finalweight {0.}; // Final event weight for the SIMC generator.
	double m_predhcal_xpos {0.};
	double m_predhcal_ypos {0.};
	bool m_pass_fiducialcut {false};
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
	double m_bbphitgt{0.};
	double m_bbpoltgt{0.};
	
	const static int m_MAXHCALCLUS {20};
	double* m_sbshcalcluse; // Ptr to "sbs.hcal.clus.e" array.
	double* m_sbshcalclusatime; // Ptr to "sbs.hcal.clus.atime" array.
	double* m_sbshcalclusx; // Ptr to "sbs.hcal.clus.x" array.
	double* m_sbshcalclusy; // Ptr to "sbs.hcal.clus.y" array.
	double* m_sbshcal_heclus_e; // Ptr to sorted sbs.hcal.clus.e array
	double* m_sbshcal_heclus_atime; // Ptr to sorted sbs.hcal.clus.atime array
	double* m_sbshcal_heclus_x; // Ptr to sorted sbs.hcal.clus.x array
	double* m_sbshcal_heclus_y; // Ptr to sorted sbs.hcal.clus.y array
	double m_sbshcal_heclus_dx;
	double m_sbshcal_heclus_dy;
	double m_shadctime {0.};
	int m_sbshcal_bestclus_indx {0}; // Index of the best cluster from Atime+HighestEnergy method.
	double m_sbshcal_bestclus_e {0.};
	double m_sbshcal_bestclus_atime {0.};
	double m_sbshcal_bestclus_x {0.};
	double m_sbshcal_bestclus_y {0.};
	double m_sbshcal_bestclus_dx {0.};
	double m_sbshcal_bestclus_dy {0.};

	// Analysis results output histograms.
	TH1D* m_h1_Q2;
	TH1D* m_h1_W2;
	TH1D* m_h1_predhcal_xpos;
	TH1D* m_h1_predhcal_ypos;
	TH2D* m_h2_predhcal_xy;
	TH1D* m_h1_heclus_dx;
	TH1D* m_h1_heclus_dy;
	TH2D* m_h2_heclus_dxdy;
	TH1D* m_h1_bestclus_dx;
	TH1D* m_h1_bestclus_dy;
	TH2D* m_h2_bestclus_dxdy;	

	// Diagnostics output histograms.
	TH1D* m_h1_bb_tr_vz;
	TH1D* m_h1_bb_tr_p;
	TH1D* m_h1_bb_ps_e;
	TH1D* m_h1_bb_eoverp;
	TH1D* m_h1_hcalheclus_sh_atimediff;
	TH1D* m_h1_hcalbestclus_e;
	TH1D* m_h1_hcalbestclus_sh_atimediff; 

	// Variables for the SIMC simulation run info.
	double m_simc_ntried {0.};
	double m_simc_genvol {0.};
	double m_simc_luminosity {0.};
	double m_simc_charge {0.};

public:

	OutputG4SBS(ElasticEventG4SBS& event, BestHCalClusForG4SBS& bestHCalClus, TString path_outputdir, const char* outputfilename, const bool is_simc = false)
	: m_event{event}, m_bestHCalClus{bestHCalClus}, m_outputdirpath{path_outputdir}, m_outputfilename{outputfilename}, m_is_simc{is_simc},
	m_sbshcalcluse{m_event.return_HCalClusEArryPtr()}, m_sbshcalclusatime{m_event.return_HCalClusAtimeArryPtr()}, m_sbshcalclusx{m_event.return_HCalClusXArryPtr()}, m_sbshcalclusy{m_event.return_HCalClusYArryPtr()},
	m_sbshcal_heclus_e{m_event.return_HCalHEClusEArryPtr()}, m_sbshcal_heclus_atime{m_event.return_HCalHEClusAtimeArryPtr()}, m_sbshcal_heclus_x{m_event.return_HCalHEClusXArryPtr()}, m_sbshcal_heclus_y{m_event.return_HCalHEClusYArryPtr()}
	{
		m_fout = new TFile(Form("%s/%s.root", m_outputdirpath.Data(), m_outputfilename),"RECREATE");
		m_resultstree = new TTree("T","Event Tree");

		m_resultstree->Branch("e.kine.Q2", &m_Q2);
		m_resultstree->Branch("e.kine.W2", &m_W2);
		if (!is_simc) m_resultstree->Branch("MC.mc_weight", &m_mc_weight);
		if (is_simc) m_resultstree->Branch("MC.simc_weight", &m_simc_weight);
		if (is_simc) m_resultstree->Branch("MC.simc_finalweight", &m_simc_finalweight);
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
		m_resultstree->Branch("cut.passfiducial", &m_pass_fiducialcut);
		m_resultstree->Branch("sbs.hcal.x", &m_sbshcalx);
		m_resultstree->Branch("sbs.hcal.y", &m_sbshcaly);
		m_resultstree->Branch("sbs.hcal.e", &m_sbshcale);
		m_resultstree->Branch("sbs.hcal.atimeblk", &m_sbshcalatimeblk);
		// m_resultstree->Branch("bb.phitgt", &m_bbphitgt);
		// m_resultstree->Branch("bb.poltgt", &m_bbpoltgt);
		m_resultstree->Branch("sbs.hcal.clus.e", m_sbshcalcluse, "m_sbshcalcluse[10]/D");
		m_resultstree->Branch("sbs.hcal.clus.atime", m_sbshcalclusatime, "m_sbshcalclusatime[10]/D");
		m_resultstree->Branch("sbs.hcal.clus.x", m_sbshcalclusx, "m_sbshcalclusx[10]/D");
		m_resultstree->Branch("sbs.hcal.clus.y", m_sbshcalclusy, "m_sbshcalclusy[10]/D");
		m_resultstree->Branch("sbs.hcal.heclus.e", m_sbshcal_heclus_e, "m_sbshcal_heclus_e[10]/D");
		m_resultstree->Branch("sbs.hcal.heclus.atime", m_sbshcal_heclus_atime, "m_sbshcal_heclus_atime[10]/D");
		m_resultstree->Branch("sbs.hcal.heclus.x", m_sbshcal_heclus_x, "m_sbshcal_heclus_x[10]/D");
		m_resultstree->Branch("sbs.hcal.heclus.y", m_sbshcal_heclus_y, "m_sbshcal_heclus_y[10]/D");
		m_resultstree->Branch("sbs.hcal.heclus.dx", &m_sbshcal_heclus_dx);
		m_resultstree->Branch("sbs.hcal.heclus.dy", &m_sbshcal_heclus_dy);
		m_resultstree->Branch("sbs.hcal.bestclus.indx", &m_sbshcal_bestclus_indx);
		m_resultstree->Branch("sbs.hcal.bestclus.e", &m_sbshcal_bestclus_e);
		m_resultstree->Branch("sbs.hcal.bestclus.atime", &m_sbshcal_bestclus_atime);
		m_resultstree->Branch("sbs.hcal.bestclus.x", &m_sbshcal_bestclus_x);
		m_resultstree->Branch("sbs.hcal.bestclus.y", &m_sbshcal_bestclus_y);
		m_resultstree->Branch("sbs.hcal.bestclus.dx", &m_sbshcal_bestclus_dx);
		m_resultstree->Branch("sbs.hcal.bestclus.dy", &m_sbshcal_bestclus_dy);

		// Analysis results output histograms.
		m_h1_Q2 = new TH1D("h1_Q2", "Q^{2}/(GeV^{2}/c^{2})", 300,0, 15);
		m_h1_W2 = new TH1D("h1_W2", "W^{2}/(GeV^{2}/c^{2})", 200, 0, 5);
		m_h1_predhcal_xpos = new TH1D("h1_predhcal_xpos", "Neutron hypothesis reconstructed x hit position on HCal; X_{pos} (m)", 1000, -4, 6);
		m_h1_predhcal_ypos = new TH1D("h1_predhcal_ypos", "Neutron hypothesis reconstructed y hit position on HCal; Y_{pos} (m)", 500, -3, 2);
		m_h2_predhcal_xy = new TH2D("h2_predhcal_xy", "Neutron hypothesis reconstructed x vs y hit position on HCal; Y_{pos} (m); X_{pos} (m)", 500, -3, 2, 1000, -4, 6);
		m_h1_heclus_dx = new TH1D("h1_heclus_dx", "dx - From highest energy HCal cluster; dx (m)", 400, -4, 4);
		m_h1_heclus_dy = new TH1D("h1_heclus_dy", "dy - From highest energy HCal cluster; dy (m)", 200, -2, 2);
		m_h2_heclus_dxdy = new TH2D("h2_heclus_dxdy", "dx vs dy - From highest energy HCal cluster; dy = Y_{hcal}-Y_{predicted} (m); dx = X_{hcal}-X_{predicted} (m)", 400, -2, 2, 800, -4, 4);
		m_h1_bestclus_dx = new TH1D("h1_bestclus_dx", "dx - From timing+HEclus algorithm; dx (m)", 400, -4, 4);
		m_h1_bestclus_dy = new TH1D("h1_bestclus_dy", "dy - From timing+HEclus algorithm; dy (m)", 200, -2, 2);
		m_h2_bestclus_dxdy = new TH2D("h2_bestclus_dxdy", "dx vs dy - From timing+HEclus algorithm; dy = Y_{hcal}-Y_{predicted} (m); dx = X_{hcal}-X_{predicted} (m)", 400, -2, 2, 800, -4, 4);

		// Diagnostic variables output histograms.
		m_h1_bb_tr_vz = new TH1D("h1_bb_tr_vz", "BB track vertex Z position distribution; vertex Z (m)", 600, -0.15, 0.15);
		m_h1_bb_tr_p = new TH1D("h1_bb_tr_p", "BB track momentum distribution; Momentum (GeV/c)", 1000, 0, 10.0);
		m_h1_bb_ps_e = new TH1D("h1_bb_ps_e", "BB pre-shower energy", 300,0, 3);
		m_h1_bb_eoverp = new TH1D("h1_bb_eoverp", "BigBite E/P distribution; E/P", 300, -1, 2);
		m_h1_hcalheclus_sh_atimediff = new TH1D("h1_hcalheclus_sh_atimediff","HCal highest energy cluster time - SH time; HCal_{ADCtime}-BBSH_{ADCtime} (ns); Entries",3000,-100,200);
		m_h1_hcalbestclus_e = new TH1D("h_hcalbestclus_e", "HCal time+HEclus energy; Cluster Energy (GeV)", 200, 0, 2.0);
		m_h1_hcalbestclus_sh_atimediff = new TH1D("h1_hcalbestclus_sh_atimediff","HCal time+HEclus time - SH time; HCal_{ADCtime}-BBSH_{ADCtime} (ns); Entries",3000,-100,200);   

		// Initiate the TTree to store simulation run specific info.
		if (is_simc) m_simcinfotree = new TTree("S", "SIMC Simulation Info Tree");
		
		if (is_simc) m_simcinfotree->Branch("MC.simc.ntried", &m_simc_ntried);
		if (is_simc) m_simcinfotree->Branch("MC.simc.genvol", &m_simc_genvol);
		if (is_simc) m_simcinfotree->Branch("MC.simc.luminosity", &m_simc_luminosity);
		if (is_simc) m_simcinfotree->Branch("MC.simc.charge", &m_simc_charge);
	}

	void copyFromEvent()
	{
		m_Q2 = m_event.return_Q2();
		m_W2 = m_event.return_W2();
		if (!m_is_simc) m_mc_weight = m_event.return_mcWeight();
		if (m_is_simc) m_simc_weight = m_event.return_MCSimcWeight();
		if (m_is_simc) m_simc_finalweight = m_event.return_SimcFinalWeight();
		m_predhcal_xpos = m_event.return_nHypthsPredx();
		m_predhcal_ypos = m_event.return_nHypthsPredy();
		m_pass_fiducialcut = m_event.return_passFiducialCut();
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
		m_sbshcal_heclus_dx = m_sbshcal_heclus_x[0] - m_predhcal_xpos;
		m_sbshcal_heclus_dy = m_sbshcal_heclus_y[0] - m_predhcal_ypos;
		m_sbshcal_bestclus_dx = m_sbshcal_bestclus_x - m_predhcal_xpos;
		m_sbshcal_bestclus_dy = m_sbshcal_bestclus_y - m_predhcal_ypos;
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

		if (m_is_simc) 
		{
			m_h1_heclus_dx->Fill(m_sbshcal_heclus_dx, m_simc_finalweight);
			m_h1_heclus_dy->Fill(m_sbshcal_heclus_dy, m_simc_finalweight);
			m_h2_heclus_dxdy->Fill(m_sbshcal_heclus_dy, m_sbshcal_heclus_dx, m_simc_finalweight);
			m_h1_bestclus_dx->Fill(m_sbshcal_bestclus_dx, m_simc_finalweight);
			m_h1_bestclus_dy->Fill(m_sbshcal_bestclus_dy, m_simc_finalweight);
			m_h2_bestclus_dxdy->Fill(m_sbshcal_bestclus_dy, m_sbshcal_bestclus_dx, m_simc_finalweight); 
		}
		else
		{
			m_h1_heclus_dx->Fill(m_sbshcal_heclus_dx, m_mc_weight);
			m_h1_heclus_dy->Fill(m_sbshcal_heclus_dy, m_mc_weight);
			m_h2_heclus_dxdy->Fill(m_sbshcal_heclus_dy, m_sbshcal_heclus_dx, m_mc_weight);
			m_h1_bestclus_dx->Fill(m_sbshcal_bestclus_dx, m_mc_weight);
			m_h1_bestclus_dy->Fill(m_sbshcal_bestclus_dy, m_mc_weight);
			m_h2_bestclus_dxdy->Fill(m_sbshcal_bestclus_dy, m_sbshcal_bestclus_dx, m_mc_weight);
		}
		 		
		// Diagnostic histos.
		m_h1_bb_tr_vz->Fill(m_bbtrvz);
		m_h1_bb_tr_p->Fill(m_bbtrp);
		m_h1_bb_ps_e->Fill(m_bbpse);
		m_h1_bb_eoverp->Fill(m_bbeoverp);
		m_h1_hcalheclus_sh_atimediff->Fill(m_sbshcal_heclus_atime[0] - m_shadctime);
		m_h1_hcalbestclus_e->Fill(m_sbshcal_bestclus_e);
		m_h1_hcalbestclus_sh_atimediff->Fill(m_sbshcal_bestclus_atime - m_shadctime);
	}

	void closeOutFile()
	{
		m_resultstree->Write(0, TObject::kWriteDelete, 0);
		m_h1_Q2->Write();
		m_h1_W2->Write(); 
		m_h1_predhcal_xpos->Write();
		m_h1_predhcal_ypos->Write();
		m_h2_predhcal_xy->Write();
		m_h1_heclus_dx->Write();
		m_h1_heclus_dy->Write();
		m_h2_heclus_dxdy->Write();
		m_h1_bestclus_dx->Write();
		m_h1_bestclus_dy->Write();
		m_h2_bestclus_dxdy->Write();		
		
		m_h1_bb_tr_vz->Write();
		m_h1_bb_tr_p->Write();
		m_h1_bb_ps_e->Write();
		m_h1_bb_eoverp->Write();
		m_h1_hcalheclus_sh_atimediff->Write();
		m_h1_hcalbestclus_e->Write();
		m_h1_hcalbestclus_sh_atimediff->Write(); 

		if (m_is_simc) m_simcinfotree->Write();

		delete[] m_sbshcalcluse;
		delete[] m_sbshcalclusatime;
		delete[] m_sbshcalclusx;
		delete[] m_sbshcalclusy;
		delete[] m_sbshcal_heclus_e;
		delete[] m_sbshcal_heclus_atime;
		delete[] m_sbshcal_heclus_x;
		delete[] m_sbshcal_heclus_y;

		m_fout->Close();

		delete m_fout;
	}

	void copySIMCrunInfo() // *SIMC only function* Copy common run information about the SIMC event generation.
	{
		m_simc_ntried = m_event.return_SimcNtried();
		m_simc_genvol = m_event.return_SimcGenvol();
		m_simc_luminosity = m_event.return_SimcLuminosity();
		m_simc_charge = m_event.return_SimcCharge();		
	}

	void fillSIMCrunInfoTree() // *SIMC only function*
	{
		m_simcinfotree->Fill();
	}
	
};

#endif