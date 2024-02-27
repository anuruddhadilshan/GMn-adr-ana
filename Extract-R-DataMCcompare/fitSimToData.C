// This script fits SIMC+G4SBS generated data to real experiment data of GMn, and then extract *R*, the differential cross-sction ratio of n(e,e') and p(e,e').

// The order of operations:
// // 1) Read real data ROOT file and make the dx plot with reasonable elastic cuts. Apply charge normalization (with live-time correction?).
// // 2) Read d(e,e'n) and d(e,e'p) simulation ROOT files and make the simulation dx plots, with the same elastic cuts as above. Apply charge normalization.
// // 3) Overlay the data and simulation histograms and inspect --> This is before any background unfolding/fitting is done for the data histograms.
// // 4) Now apply any background estimation techniques.

#include "TMinuit.h"
#include "../includes/read_fitSimToData_config.h"

constexpr int n_backg_pol = 3; // Order of the background polynomial to be used to fit real data.
double backg_pol_par [n_backg_pol+1] {0.}; // Parameter array for background polynomial from the fit results.
constexpr int dx_nbins = 150; // Number of bins used for ALL the dx histograms.
constexpr double dx_low = -3.0; // Lower limit used for ALL the dx histograms.
constexpr double dx_high = 2.5; // Higher limit use for ALL the dx histograms.

double experhist_lowlimit {0.}; // First and last dx values of the experiment data dx histogram with non-zero bins.
double experhist_highlimit {0.}; 

void fill_ExperimentDatHistos( const TString experimentdat_rootfile, TH1D* h1_experiment_dx, const double w2_min, const double w2_max, const double dy_min, const double dy_max );
void fill_SimDatHistos( const TString simdat_rootfile, TH1D* h1_simdat_dx, const double w2_min, const double w2_max, const double dy_min, const double dy_max, const double dx_offset );
void fill_G4SBSinelHistos( const TString g4sbs_inel_rootfile, TH1D* h1_g4sbsinel_dx, const double w2_min, const double w2_max, const double dy_min, const double dy_max , double p_count, double backg_to_p_ratio);
void fit_ExperimentDatHistos( TH1D* h1_exp_dx, double pol_par [n_backg_pol+1], double& experhist_lowlimit, double& experhist_highlimit, double& backg_to_p_ratio ); //Function to fit experiment dx histogram and extract the parameters for the background histo.
void make_PolBackgHistoForSim( TH1D* h1_backg_dx, double p_count, double backg_to_p_ratio, double pol_par [n_backg_pol+1], const double low_limit, const double high_limit );
void make_PolBackgHistoForSim1( TH1D* h1_backg_dx, double par0, double par1, double par2, double par3 );
void Chi2_PolBackgAddToSim(int& npar, double* deriv, double& result, double* par, int flag);
void Chi2_PolBackgAddToSim1(int& npar, double* deriv, double& result, double* par, int flag);
void Chi2_G4SBSbackg(int& npar, double* deriv, double& result, double* par, int flag);
double return_Chi2();
void substract_PolBackgHistoFromData( TH1D* h1_exp_dx, TH1D* h1_backg_dx, double pol_par [n_backg_pol+1], const double low_limit, const double high_limit );
void plotOverlayDataAndSimHistos( TH1D* datHist, TH1D* simDeen, TH1D* simDeep, TH1D* simBackg, TH1D* SimFullHist, const char* dataTitle, const char* simTitle, const char* canvasTitle );
void makeResidualPlot( TH1D* h1_exp_dx, TH1D* h1_sim_dx );

TH1D* h1_experiment_dx;
TH1D* h1_simc_original_deen_dx;
TH1D* h1_simc_original_deep_dx;
TH1D* h1_simc_original_deeN_dx;
TH1D* h1_simc_scaled_deen_dx;
TH1D* h1_simc_scaled_deep_dx;
TH1D* h1_simc_scaled_deeN_dx;
TH1D* h1_polbackg_dx;
TH1D* h1_polbackg_scaled_dx;
TH1D* h1_g4sbsinelbackg_dx;
TH1D* h1_g4sbsinelbackg_scaled_dx;
TH1D* h1_dx_temp;



void fitSimToData( const char* configfilename, const int nbackg_tech = 0 ) 
// nbackg_tech = 0 ==> Fit real dat. dx backg with polN. Add a scaled background to the SIMC simulated data dx, distributed according to the polN from fitting real data. Minimization algorithm will ONLY scale the background preserving its shape.
// nbackg_tech = 1 ==> Fit real dat. dx backg with polN. Add a scaled background to the SIMC simulated data dx, distributed according to the polN from fitting real data. Let all the parameters in the polN change in the minimzation algorithm. 
// nbackg_tech = 2 ==> Fit real dat. dx backg with polN. Substract that background from the real data dx itself, using the polN fit results.
// nbackg_tech = 3 ==> 
{
	ConfigfileFitSimToData configfile{configfilename};

	const TString experimentdat_rootfile = configfile.return_ExprDatFile();
	const TString simc_deen_rootfile = configfile.return_SIMCdeenFile();
	const TString simc_deep_rootfile = configfile.return_SIMCdeepFile();
	const TString g4sbs_inel_rootfile = configfile.return_G4SBSinelFile();
	const double cut_W2_min = configfile.return_W2CutMin();
	const double cut_W2_max = configfile.return_W2CutMax();
	const double cut_dy_min = configfile.return_dyCutMin();
	const double cut_dy_max = configfile.return_dyCutMax();
	const double simoffset_dx_deen = configfile.return_simDeenDxOffset();
	const double simoffset_dx_deep = configfile.return_simDeepDxOffset();

	h1_experiment_dx = new TH1D("h1_experiment_dx", "dx - Experiment Data; dx (m)", dx_nbins, dx_low, dx_high);
	h1_simc_original_deen_dx = new TH1D("h1_simc_original_deen_dx", "dx - SIMC deen; dx (m)", dx_nbins, dx_low, dx_high);
	h1_simc_original_deep_dx = new TH1D("h1_simc_original_deep_dx", "dx - SIMC deep; dx (m)", dx_nbins, dx_low, dx_high);
	h1_simc_original_deeN_dx = new TH1D("h1_simc_original_deeN_dx", "dx - SIMC deeN; dx (m)", dx_nbins, dx_low, dx_high); // Sum of deen and deep histograms.
	h1_polbackg_dx = new TH1D("h1_polbackg_dx", "dx - pol background; dx (m)", dx_nbins, dx_low, dx_high);
	h1_g4sbsinelbackg_dx = new TH1D("h1_g4sbsinel_dx", "dx - inelastic backg. from G4SBS; dx (m)", dx_nbins, dx_low, dx_high);
	h1_dx_temp = new TH1D("h1_dx_temp", "dx - temp; dx (m)", dx_nbins, dx_low, dx_high);

	fill_ExperimentDatHistos(experimentdat_rootfile, h1_experiment_dx, cut_W2_min, cut_W2_max, cut_dy_min, cut_dy_max);
	fill_SimDatHistos(simc_deen_rootfile, h1_simc_original_deen_dx, cut_W2_min, cut_W2_max, cut_dy_min, cut_dy_max, simoffset_dx_deen);
	fill_SimDatHistos(simc_deep_rootfile, h1_simc_original_deep_dx, cut_W2_min, cut_W2_max, cut_dy_min, cut_dy_max, simoffset_dx_deep);

	double fn_result {0.};
	double fp_result {0.};
	double fb_result {0.};
	double fn_error {0.};
	double fp_error {0.};
	double fb_error {0.};

	if ( nbackg_tech == 0 ) //Fit real dat. dx backg with polN. Add a scaled background to the SIMC simulated data dx, distributed according to the polN from fitting real data.
	{
		double backg_to_p_ratio {0.}; // Background to proton ratio for scalling of background when adding to the simulation.
				
		// Start making dx hitograms and adding background to the simulated data histograms //
		fit_ExperimentDatHistos(h1_experiment_dx, backg_pol_par, experhist_lowlimit, experhist_highlimit, backg_to_p_ratio);

		make_PolBackgHistoForSim(h1_polbackg_dx, h1_simc_original_deep_dx->Integral(), backg_to_p_ratio, backg_pol_par, experhist_lowlimit, experhist_highlimit);

		double tot_int {0.};
		tot_int = h1_simc_original_deen_dx->Integral() + h1_simc_original_deep_dx->Integral() + h1_polbackg_dx->Integral();

		h1_simc_original_deen_dx->Scale( 1.0 / tot_int ); // Normalize to 1.0
		h1_simc_original_deep_dx->Scale( 1.0 / tot_int ); // Normalize to 1.0
		h1_polbackg_dx->Scale( 1.0 / tot_int );

		h1_dx_temp->Add( h1_simc_original_deen_dx, h1_simc_original_deep_dx );
		h1_simc_original_deeN_dx->Add( h1_dx_temp, h1_polbackg_dx );

		h1_experiment_dx->Scale( 1.0 / h1_experiment_dx->Integral() );
		// Done making the histograms and adding background //

		h1_simc_scaled_deen_dx = (TH1D*)h1_simc_original_deen_dx->Clone("h1_simc_scaled_deen_dx");
		h1_simc_scaled_deep_dx = (TH1D*)h1_simc_original_deep_dx->Clone("h1_simc_scaled_deep_dx");
		h1_polbackg_scaled_dx = (TH1D*)h1_polbackg_dx->Clone("h1_polbackg_scaled_dx");
		h1_simc_scaled_deeN_dx = (TH1D*)h1_simc_original_deeN_dx->Clone("h1_simc_scaled_deeN_dx");

		plotOverlayDataAndSimHistos( h1_experiment_dx, h1_simc_original_deen_dx, h1_simc_original_deep_dx, h1_polbackg_dx, h1_simc_original_deeN_dx, "Data", "Simulation pre-fitting", "Data and Simulation comparison before fitting" );		

		// Start the simulation fit to data process //
		// // Creating an instance of TMinuit
		TMinuit minuit(3); 

		// // Set the user-defined function to be minimized.
		minuit.SetFCN(Chi2_PolBackgAddToSim);

		// // Set starting values for parameters
		double start_fn = 1.0;
		double start_fp = 1.0;
		double start_fb = 1.0;

		int ierflg = 0;

		minuit.mnparm(0, "fn", start_fn, 0.001, 0.5, 1.5, ierflg);
		minuit.mnparm(1, "fp", start_fp, 0.001, 0.5, 1.5, ierflg);
		minuit.mnparm(2, "fb", start_fb, 0.001, 0.5, 1.5, ierflg);

		// // Set minimization strategy and other options if needed.
		double arglist[10];

		arglist[0] = 2; // Use MIGRAD algorithm.

		minuit.mnexcm("SET STR", arglist, 1, ierflg);

		// // Minimize the function
		arglist[0] = 5000; // Maximum number of function calls.
		arglist[1] = 1.0;  // Tolerence.

		minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
		// End the simulation fit tot data process //

		// // Get results
		minuit.GetParameter(0, fn_result, fn_error);
		minuit.GetParameter(1, fp_result, fp_error);
		minuit.GetParameter(2, fb_result, fb_error);
		//double min_chi2 = minuit.fAmin();
		// // Apply the optimal scale
		h1_simc_scaled_deen_dx->Scale(fn_result);
		h1_simc_scaled_deep_dx->Scale(fp_result);
		h1_polbackg_scaled_dx->Scale(fb_result);
	
		h1_dx_temp->Add( h1_simc_scaled_deen_dx, h1_simc_scaled_deep_dx );
		h1_simc_scaled_deeN_dx->Add(h1_dx_temp, h1_polbackg_scaled_dx);

		plotOverlayDataAndSimHistos( h1_experiment_dx, h1_simc_scaled_deen_dx, h1_simc_scaled_deep_dx, h1_polbackg_scaled_dx, h1_simc_scaled_deeN_dx, "Data", "Simulation", "Data and Simulation comparison after fitting" );		
	}
	else if ( nbackg_tech == 1 ) // Fit real dat. dx backg with polN. Add a scaled background to the SIMC simulated data dx, distributed according to the polN from fitting real data. Let all the parameters in the polN change in the minimzation algorithm. 
	{
		double backg_to_p_ratio {0.}; // Background to proton ratio for scalling of background when adding to the simulation.
		// double backg_pol_par [n_backg_pol+1] {0.}; // Parameter array for background polynomial from the fit results.
		
		// Start making dx hitograms and adding background to the simulated data histograms //
		fit_ExperimentDatHistos(h1_experiment_dx, backg_pol_par, experhist_lowlimit, experhist_highlimit, backg_to_p_ratio);

		make_PolBackgHistoForSim1(h1_polbackg_dx, backg_pol_par[0], backg_pol_par[1], backg_pol_par[2], backg_pol_par[3]);

		double tot_int {0.};
		tot_int = h1_simc_original_deen_dx->Integral() + h1_simc_original_deep_dx->Integral() + h1_polbackg_dx->Integral();

		h1_simc_original_deen_dx->Scale( 1.0 / tot_int ); // Normalize to 1.0
		h1_simc_original_deep_dx->Scale( 1.0 / tot_int ); // Normalize to 1.0
		h1_polbackg_dx->Scale( 1.0 / tot_int );

		h1_dx_temp->Add( h1_simc_original_deen_dx, h1_simc_original_deep_dx );
		h1_simc_original_deeN_dx->Add( h1_dx_temp, h1_polbackg_dx );

		h1_experiment_dx->Scale( 1.0 / h1_experiment_dx->Integral() );
		// Done making the histograms and adding background //

		plotOverlayDataAndSimHistos( h1_experiment_dx, h1_simc_original_deen_dx, h1_simc_original_deep_dx, h1_polbackg_dx, h1_simc_original_deeN_dx, "Data", "Simulation pre-fitting", "Data and Simulation comparison before fitting" );	

		h1_simc_scaled_deen_dx = (TH1D*)h1_simc_original_deen_dx->Clone("h1_simc_scaled_deen_dx");
		h1_simc_scaled_deep_dx = (TH1D*)h1_simc_original_deep_dx->Clone("h1_simc_scaled_deep_dx");
		h1_polbackg_scaled_dx = (TH1D*)h1_polbackg_dx->Clone("h1_polbackg_scaled_dx");
		h1_simc_scaled_deeN_dx = (TH1D*)h1_simc_original_deeN_dx->Clone("h1_simc_scaled_deeN_dx");

		// for (const auto& val : backg_pol_par)
		// {
		// 	std::cout << val << '\n';
		// }
		// Start the simulation fit to data process //
		// // Creating an instance of TMinuit
		TMinuit minuit(6); 

		// // Set the user-defined function to be minimized.
		minuit.SetFCN(Chi2_PolBackgAddToSim1);

		// // Set starting values for parameters
		double start_fn = 1.0;
		double start_fp = 1.0;
		double start_fb = 1.0;
		double start_pol_p0 = backg_pol_par[0];
		double start_pol_p1 = backg_pol_par[1];
		double start_pol_p2 = backg_pol_par[2];
		double start_pol_p3 = backg_pol_par[3];

		int ierflg = 0;

		minuit.mnparm(0, "fn", start_fn, 0.001, 0.5, 1.5, ierflg);
		minuit.mnparm(1, "fp", start_fp, 0.001, 0.5, 1.5, ierflg);
		minuit.mnparm(2, "pol_p0", start_pol_p0, 1, 0.5, 1000.0, ierflg);
		minuit.mnparm(3, "pol_p1", start_pol_p1, 1, -100, 100, ierflg);
		minuit.mnparm(4, "pol_p2", start_pol_p2, 1, -100, 100, ierflg);
		minuit.mnparm(5, "pol_p3", start_pol_p3, 1, -100, 100, ierflg);

		// // Set minimization strategy and other options if needed.
		double arglist[10];

		arglist[0] = 2; // Use MIGRAD algorithm.

		minuit.mnexcm("SET STR", arglist, 1, ierflg);

		// // Minimize the function
		arglist[0] = 5000; // Maximum number of function calls.
		arglist[1] = 1.0;  // Tolerence.

		minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
		// End the simulation fit tot data process //

		// // Get results
		double pol_p0_result {0.};
		double pol_p0_error {0.};
		double pol_p1_result {0.};
		double pol_p1_error {0.};
		double pol_p2_result {0.};
		double pol_p2_error {0.};
		double pol_p3_result {0.};
		double pol_p3_error {0.};
		minuit.GetParameter(0, fn_result, fn_error);
		minuit.GetParameter(1, fp_result, fp_error);
		minuit.GetParameter(2, pol_p0_result, pol_p0_error);
		minuit.GetParameter(3, pol_p1_result, pol_p1_error);
		minuit.GetParameter(4, pol_p2_result, pol_p2_error);
		minuit.GetParameter(5, pol_p3_result, pol_p3_error);

		make_PolBackgHistoForSim1(h1_polbackg_scaled_dx, pol_p0_result, pol_p1_result, pol_p2_result, pol_p3_result);
	}
	else if ( nbackg_tech == 2 )
	{
		double backg_to_p_ratio {0.}; // Background to proton ratio for scalling of background when adding to the simulation.
				
		// Start making dx hitograms and adding background to the simulated data histograms //
		fit_ExperimentDatHistos(h1_experiment_dx, backg_pol_par, experhist_lowlimit, experhist_highlimit, backg_to_p_ratio);

		fill_G4SBSinelHistos(g4sbs_inel_rootfile, h1_g4sbsinelbackg_dx, cut_W2_min, cut_W2_max, cut_dy_min, cut_dy_max, h1_simc_original_deep_dx->Integral(), backg_to_p_ratio);

		double tot_int {0.};
		tot_int = h1_simc_original_deen_dx->Integral() + h1_simc_original_deep_dx->Integral() + h1_g4sbsinelbackg_dx->Integral();

		h1_simc_original_deen_dx->Scale( 1.0 / tot_int ); // Normalize to 1.0
		h1_simc_original_deep_dx->Scale( 1.0 / tot_int ); // Normalize to 1.0
		h1_g4sbsinelbackg_dx->Scale( 1.0 / tot_int );

		h1_dx_temp->Add( h1_simc_original_deen_dx, h1_simc_original_deep_dx );
		h1_simc_original_deeN_dx->Add( h1_dx_temp, h1_g4sbsinelbackg_dx );

		h1_experiment_dx->Scale( 1.0 / h1_experiment_dx->Integral() );
		// Done making the histograms and adding background //

		h1_simc_scaled_deen_dx = (TH1D*)h1_simc_original_deen_dx->Clone("h1_simc_scaled_deen_dx");
		h1_simc_scaled_deep_dx = (TH1D*)h1_simc_original_deep_dx->Clone("h1_simc_scaled_deep_dx");
		h1_g4sbsinelbackg_scaled_dx = (TH1D*)h1_g4sbsinelbackg_dx->Clone("h1_polbackg_scaled_dx");
		h1_simc_scaled_deeN_dx = (TH1D*)h1_simc_original_deeN_dx->Clone("h1_simc_scaled_deeN_dx");

		plotOverlayDataAndSimHistos( h1_experiment_dx, h1_simc_original_deen_dx, h1_simc_original_deep_dx, h1_g4sbsinelbackg_dx, h1_simc_original_deeN_dx, "Data", "Simulation pre-fitting", "Data and Simulation comparison before fitting" );		

		// Start the simulation fit to data process //
		// // Creating an instance of TMinuit
		TMinuit minuit(3); 

		// // Set the user-defined function to be minimized.
		minuit.SetFCN(Chi2_G4SBSbackg);

		// // Set starting values for parameters
		double start_fn = 1.0;
		double start_fp = 1.0;
		double start_fb = 1.0;

		int ierflg = 0;

		minuit.mnparm(0, "fn", start_fn, 0.001, 0.5, 1.5, ierflg);
		minuit.mnparm(1, "fp", start_fp, 0.001, 0.5, 1.5, ierflg);
		minuit.mnparm(2, "fb", start_fb, 0.001, 0.5, 1.5, ierflg);

		// // Set minimization strategy and other options if needed.
		double arglist[10];

		arglist[0] = 2; // Use MIGRAD algorithm.

		minuit.mnexcm("SET STR", arglist, 1, ierflg);

		// // Minimize the function
		arglist[0] = 5000; // Maximum number of function calls.
		arglist[1] = 1.0;  // Tolerence.

		minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
		// End the simulation fit tot data process //

		// // Get results
		minuit.GetParameter(0, fn_result, fn_error);
		minuit.GetParameter(1, fp_result, fp_error);
		minuit.GetParameter(2, fb_result, fb_error);
		//double min_chi2 = minuit.fAmin();
		// // Apply the optimal scale
		h1_simc_scaled_deen_dx->Scale(fn_result);
		h1_simc_scaled_deep_dx->Scale(fp_result);
		h1_g4sbsinelbackg_scaled_dx->Scale(fb_result);
	
		h1_dx_temp->Add( h1_simc_scaled_deen_dx, h1_simc_scaled_deep_dx );
		h1_simc_scaled_deeN_dx->Add(h1_dx_temp, h1_g4sbsinelbackg_scaled_dx);

		plotOverlayDataAndSimHistos( h1_experiment_dx, h1_simc_scaled_deen_dx, h1_simc_scaled_deep_dx, h1_g4sbsinelbackg_scaled_dx, h1_simc_scaled_deeN_dx, "Data", "Simulation", "Data and Simulation comparison after fitting" );		

	}

	makeResidualPlot( h1_experiment_dx, h1_simc_scaled_deeN_dx );

	// Print results
	std::cout << "Optimal fn: " << fn_result << " +/- " << fn_error << '\n';
	std::cout << "Optimal fp: " << fp_result << " +/- " << fp_error << '\n';
	std::cout << "Optimal fb: " << fb_result << " +/- " << fb_error << '\n';
	double fn_div_fp = fn_result / fp_result;
	std::cout << "fn/fp: " << fn_div_fp << " +/- " << fn_div_fp * sqrt( pow((fn_error/fn_result),2) + pow((fp_error/fp_result),2) ) << '\n';
}




void fill_ExperimentDatHistos(const TString experimentdat_rootfile, TH1D* h1_experiment_dx, const double w2_min, const double w2_max, const double dy_min, const double dy_max)
{	
	TFile* rootfile = new TFile(experimentdat_rootfile.Data(),"READ");
	TTree* T = (TTree*)rootfile->Get("T");

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("e.kine.W2", 1);
	T->SetBranchStatus("cut.passfiducial", 1);
	T->SetBranchStatus("cut.passcointime", 1);
	T->SetBranchStatus("sbs.hcal.bestclus.dx", 1);
	T->SetBranchStatus("sbs.hcal.bestclus.dy", 1);

	double w2 {0.};
	bool pass_fiducial {false};
	bool pass_cointime {false};
	double dx {0.};
	double dy {0.};

	T->SetBranchAddress("e.kine.W2", &w2);
	T->SetBranchAddress("cut.passfiducial", &pass_fiducial);
	T->SetBranchAddress("cut.passcointime", &pass_cointime);
	T->SetBranchAddress("sbs.hcal.bestclus.dx", &dx);
	T->SetBranchAddress("sbs.hcal.bestclus.dy", &dy);

	long nevent {0};

	while( T->GetEntry(nevent++) )
	{
		if ( w2 < w2_min || w2 > w2_max || dy < dy_min || dy > dy_max || !pass_fiducial || !pass_cointime ) continue;

		h1_experiment_dx->Fill(dx);
	}

	//h1_experiment_dx->Scale( 1.0 / h1_experiment_dx->Integral() );
}

void fill_SimDatHistos(const TString simdat_rootfile, TH1D* h1_simdat_dx, const double w2_min, const double w2_max, const double dy_min, const double dy_max, const double dx_offset)
{
	TFile* rootfile = new TFile(simdat_rootfile.Data(),"READ");
	TTree* T = (TTree*)rootfile->Get("T");

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("e.kine.W2", 1);
	T->SetBranchStatus("cut.passfiducial", 1);
	T->SetBranchStatus("sbs.hcal.bestclus.dx", 1);
	T->SetBranchStatus("sbs.hcal.bestclus.dy", 1);
	T->SetBranchStatus("MC.simc_finalweight", 1);
	T->SetBranchStatus("sbs.hcal.pred.x", 1);
	T->SetBranchStatus("sbs.hcal.pred.y", 1);
	T->SetBranchStatus("sbs.hcal.clus.x", 1);
	T->SetBranchStatus("sbs.hcal.clus.y", 1);

	double w2 {0.};
	bool pass_fiducial {false};
	double dx {0.};
	double dy {0.};
	double simc_finalweight {0.};
	double pred_x {0.};
	double pred_y {0.};
	double clus_x [10] {0.};
	double clus_y [10] {0.};

	T->SetBranchAddress("e.kine.W2", &w2);
	T->SetBranchAddress("cut.passfiducial", &pass_fiducial);
	T->SetBranchAddress("sbs.hcal.bestclus.dx", &dx);
	T->SetBranchAddress("sbs.hcal.bestclus.dy", &dy);
	T->SetBranchAddress("MC.simc_finalweight", &simc_finalweight);
	T->SetBranchAddress("sbs.hcal.pred.x", &pred_x);
	T->SetBranchAddress("sbs.hcal.pred.y", &pred_y);
	T->SetBranchAddress("sbs.hcal.clus.x", clus_x);
	T->SetBranchAddress("sbs.hcal.clus.y", clus_y);

	long nevent {0};

	while( T->GetEntry(nevent++) )
	{
		if ( w2 < w2_min || w2> w2_max || dy < dy_min || dy > dy_max || !pass_fiducial ) continue;

		double dx_clus = clus_x[0] - pred_x;

		h1_simdat_dx->Fill(dx_clus + dx_offset, simc_finalweight);
	}

	TTree* S = (TTree*)rootfile->Get("S");
	S->SetBranchStatus("MC.simc.charge", 1);
	double simc_totcharge {0.};
	S->SetBranchAddress("MC.simc.charge", &simc_totcharge);
	S->GetEntry(0);

	//h1_simdat_dx->Scale(1/simc_totcharge);
}

void fill_G4SBSinelHistos( const TString g4sbs_inel_rootfile, TH1D* h1_g4sbsinel_dx, const double w2_min, const double w2_max, const double dy_min, const double dy_max, double p_count, double backg_to_p_ratio )
{
	TFile* rootfile = new TFile(g4sbs_inel_rootfile.Data(),"READ");
	TTree* T = (TTree*)rootfile->Get("T");

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("e.kine.W2", 1);
	T->SetBranchStatus("cut.passfiducial", 1);
	T->SetBranchStatus("sbs.hcal.bestclus.dx", 1);
	T->SetBranchStatus("sbs.hcal.bestclus.dy", 1);
	T->SetBranchStatus("MC.mc_weight", 1);
	T->SetBranchStatus("sbs.hcal.pred.x", 1);
	T->SetBranchStatus("sbs.hcal.pred.y", 1);
	T->SetBranchStatus("sbs.hcal.clus.x", 1);
	T->SetBranchStatus("sbs.hcal.clus.y", 1);

	double w2 {0.};
	bool pass_fiducial {false};
	double dx {0.};
	double dy {0.};
	double mc_weight {0.};
	double pred_x {0.};
	double pred_y {0.};
	double clus_x [10] {0.};
	double clus_y [10] {0.};

	T->SetBranchAddress("e.kine.W2", &w2);
	T->SetBranchAddress("cut.passfiducial", &pass_fiducial);
	T->SetBranchAddress("sbs.hcal.bestclus.dx", &dx);
	T->SetBranchAddress("sbs.hcal.bestclus.dy", &dy);
	T->SetBranchAddress("MC.mc_weight", &mc_weight);
	T->SetBranchAddress("sbs.hcal.pred.x", &pred_x);
	T->SetBranchAddress("sbs.hcal.pred.y", &pred_y);
	T->SetBranchAddress("sbs.hcal.clus.x", clus_x);
	T->SetBranchAddress("sbs.hcal.clus.y", clus_y);

	long nevent {0};

	while( T->GetEntry(nevent++) )
	{
		//if ( w2 < w2_min || w2> w2_max || dy < dy_min || dy > dy_max || !pass_fiducial ) continue;

		double dx_clus = clus_x[0] - pred_x;

		h1_g4sbsinel_dx->Fill(dx_clus, mc_weight);
	}

	double scaleFactor_toMathchData =  ( p_count * backg_to_p_ratio ) / ( h1_g4sbsinel_dx->Integral() );

	h1_g4sbsinel_dx->Scale( scaleFactor_toMathchData );
}

void fit_ExperimentDatHistos( TH1D* h1_exp_dx, double pol_par [n_backg_pol+1], double& experhist_lowlimit, double& experhist_highlimit, double& backg_to_p_ratio )
{
	TH1D* h1_exp_clone_dx = (TH1D*)h1_exp_dx->Clone("h1_exp_clone_dx");

	TCanvas* C1 = new TCanvas("C1", "dx Fitting", 1200, 900);
	C1->Divide(2,1);
	C1->cd(1);
	h1_exp_clone_dx->Draw("HIST");
	C1->Update();

	// 1) Let's do a first order Gaus fitting for n and p peaks.
	std::cout << "Enter an approximate demarcation value between the proton and the neutron peak: ";
	double demarcation_val {0.};
	std::cin >> demarcation_val;
	std::cout << '\n';

	double maxbincontent_dxhist = h1_exp_clone_dx->GetMaximum();
	TLine* demarcation_line = new TLine(demarcation_val,0.0,demarcation_val,maxbincontent_dxhist);
	demarcation_line->SetLineStyle(2);
	demarcation_line->Draw("SAME");

	// 1.1) Get the maximum bin content's x pos values of the proton and neutron side of the demarcation line.
	int demarcation_bin = h1_exp_clone_dx->GetXaxis()->FindBin(demarcation_val); // Bin number of the bin that the demarcation line goes through.
	int min_bin = h1_exp_clone_dx->GetXaxis()->GetFirst(); // First bin.
	int max_bin = h1_exp_clone_dx->GetXaxis()->GetLast(); // Last bin.
	double ppeak_maxbincontent = h1_exp_clone_dx->GetBinContent(min_bin);
	double ppeak_maxbin = min_bin;

	for ( int i = min_bin+1; i < demarcation_bin; i++ )
	{
		if ( h1_exp_clone_dx->GetBinContent(i) > ppeak_maxbincontent )
		{
			ppeak_maxbincontent = h1_exp_clone_dx->GetBinContent(i);
			ppeak_maxbin = i;
		}
	}

	double npeak_maxbincontent = h1_exp_clone_dx->GetBinContent(demarcation_bin);
	double npeak_maxbin = demarcation_bin;

	for ( int i = demarcation_bin; i <= max_bin; i++ )
	{
		if ( h1_exp_clone_dx->GetBinContent(i) > npeak_maxbincontent )
		{
			npeak_maxbincontent = h1_exp_clone_dx->GetBinContent(i);
			npeak_maxbin = i;
		}
	}

	// 1.2) Now generate the ranges for the first order fitting. 
	// p peak
	double ppeak_1_lowlimit = h1_exp_clone_dx->GetXaxis()->GetBinCenter(ppeak_maxbin) - 0.4;
	double ppeak_1_highlimit = demarcation_val;

	// n peak
	double npeak_1_lowlimit = demarcation_val;
	double npeak_1_highlimit = h1_exp_clone_dx->GetXaxis()->GetBinCenter(npeak_maxbin) + 0.4;

	// 1.3) Now do a first order gaus fit to the proton and neutron peak.
	TF1* pfit_1 = new TF1("pfit_1", "gaus", ppeak_1_lowlimit, ppeak_1_highlimit);
	h1_exp_clone_dx->Fit(pfit_1, "RQN");
	pfit_1->Draw("SAME");

	TF1* nfit_1 = new TF1("nfit_1", "gaus", npeak_1_lowlimit, npeak_1_highlimit);
	h1_experiment_dx->Fit(nfit_1, "RQN");
	nfit_1->Draw("SAME");
	C1->Update();

	// 2) Performing second order Gasu fits for n and p Gaus fits.
	// 2.1) Now generate more tighter ranges for a second order fitting.
	// p peak
	double pfit_1_par [3] {0.};
	pfit_1->GetParameters(&pfit_1_par[0]);
	double ppeak_2_lowlimit = pfit_1_par[1] - 1.0*pfit_1_par[2];
	double ppeak_2_highlimit = pfit_1_par[1] + 0.8*pfit_1_par[2];

	// n peak
	double nfit_1_par [3] {0.};
	nfit_1->GetParameters(&nfit_1_par[0]);
	double npeak_2_lowlimit = nfit_1_par[1] - 0.5*nfit_1_par[2];
	double npeak_2_highlimit = nfit_1_par[1] + 1.0*nfit_1_par[2];

	// 2.2) Define two second order fit functions for the proton and the neutron peaks and initilize their paramters using the paremters of the first order fit functions.
	C1->cd(2);
	h1_exp_clone_dx->Draw("HIST");

	TF1* pfit_2 = new TF1("pfit_2", "gaus", ppeak_2_lowlimit, ppeak_2_highlimit);
	pfit_2->SetParameters(pfit_1_par);
	h1_exp_clone_dx->Fit(pfit_2, "RQN");
	pfit_2->Draw("SAME");

	TF1* nfit_2 = new TF1("nfit_2", "gaus", npeak_2_lowlimit, npeak_2_highlimit);
	nfit_2->SetParameters(nfit_1_par);
	h1_exp_clone_dx->Fit(nfit_2, "RQN");
	nfit_2->Draw("SAME");

	// 3.0) Now perform a first order polN fit for the background.
	// // First let's find the first and the last non-zero bin of the experiment data dx histogram, so we don't go beyond that range.
	int firstNonZeroBin = -1;

	for ( int i = 1; i <= h1_exp_clone_dx->GetNbinsX(); i++ )
	{
		if ( h1_exp_clone_dx->GetBinContent(i) != 0.0 )
		{
			firstNonZeroBin = i;
			break;
		}
	}

	int lastNonZeroBin = -1;

	for ( int i = h1_exp_clone_dx->GetNbinsX(); i >= 1; i-- )
	{
		if ( h1_exp_clone_dx->GetBinContent(i) != 0.0 )
		{
			lastNonZeroBin = i;
			break;
		}
	}

	experhist_lowlimit = h1_exp_clone_dx->GetXaxis()->GetBinCenter(firstNonZeroBin);
	experhist_highlimit = h1_exp_clone_dx->GetXaxis()->GetBinCenter(lastNonZeroBin);

	TF1* backgfit_1 = new TF1("backgfit_1", Form("pol%i", n_backg_pol), experhist_lowlimit, experhist_highlimit); //, h1_dx_withelastcut->GetXaxis()->GetBinCenter(min_bin), h1_dx_withelastcut->GetXaxis()->GetBinCenter(max_bin));
	//backgfit_1->SetRange(h1_dx_withelastcut->GetXaxis()->GetBinCenter(min_bin), pfit_1_par[1] - 2.5*pfit_1_par[2]);
	//backgfit_1->SetRange(nfit_1_par[1] + 2.5*nfit_1_par[2], h1_dx_withelastcut->GetXaxis()->GetBinCenter(max_bin), "R+");
	h1_exp_clone_dx->Fit(backgfit_1, "RQN");
	backgfit_1->Draw("SAME");
	C1->Update();

	// 4.0) Define a total fit function and with two Gasians for the n and p and a polN for the background. 
	TF1* totalfit = new TF1("totalfit", Form("gaus(0)+gaus(3)+pol%i(6)", n_backg_pol), experhist_lowlimit, experhist_highlimit);
	
	// 4.1) Get the parameters from the two second order Gaus fits for n and p and the first order polN fit for the background.
	double totalfit_1_par [6+n_backg_pol+1] {0.};
	nfit_2->GetParameters(&totalfit_1_par[0]);
	pfit_2->GetParameters(&totalfit_1_par[3]);
	backgfit_1->GetParameters(&totalfit_1_par[6]);

	totalfit->SetParameters(totalfit_1_par);

	// 4.2) Apply the total fit.
	TCanvas* C2 = new TCanvas("C2", "Final Fits", 1200, 900);
	C2->Divide(2,1);
	C2->cd(1);
	h1_exp_clone_dx->Draw("HIST");
	h1_exp_clone_dx->Fit(totalfit, "RQN");
	totalfit->Draw("SAME");

	TLegend* legend_totalfitresultsplot = new TLegend(0.65,0.75,0.9,0.9);
	legend_totalfitresultsplot->AddEntry(h1_exp_clone_dx, "dx histogram", "l");
	legend_totalfitresultsplot->AddEntry(totalfit, "Total Fit", "l");
	legend_totalfitresultsplot->Draw("SAME");
	C2->Update();

	// 4.3) Create n and p Gaus and backg pol4 TF1 functions using the fit result parameters and draw them.
	double totalfit_result_par [6+n_backg_pol+1] {0.};
	totalfit->GetParameters(&totalfit_result_par[0]);

	// n function
	TF1* n_func = new TF1("n_func", "gaus", experhist_lowlimit, experhist_highlimit);
	n_func->SetParameters(&totalfit_result_par[0]);
	n_func->SetLineColor(kViolet);

	// p function
	TF1* p_func = new TF1("p_func", "gaus", experhist_lowlimit, experhist_highlimit);
	p_func->SetParameters(&totalfit_result_par[3]);
	p_func->SetLineColor(kOrange+7);

	// Background function
	TF1* backg_func = new TF1("backg_func", Form("pol%i", n_backg_pol), h1_exp_clone_dx->GetXaxis()->GetBinCenter(firstNonZeroBin), h1_exp_clone_dx->GetXaxis()->GetBinCenter(lastNonZeroBin));
	backg_func->SetParameters(&totalfit_result_par[6]);
	backg_func->SetLineColor(kGreen);

	C2->cd(2);
	h1_exp_clone_dx->Draw("HIST");
	n_func->Draw("SAME");
	p_func->Draw("SAME");
	backg_func->Draw("SAME");

	TLegend* legend_individfitfuncplot = new TLegend(0.65,0.75,0.9,0.9);
	legend_individfitfuncplot->AddEntry(h1_exp_clone_dx, "dx histogram", "l");
	legend_individfitfuncplot->AddEntry(n_func, "Neutron gaus", "l");
	legend_individfitfuncplot->AddEntry(p_func, "Proton gaus", "l");
	legend_individfitfuncplot->AddEntry(backg_func, Form("Background pol%i", n_backg_pol), "l");
	legend_individfitfuncplot->Draw("SAME");

	for ( int i = 0; i < n_backg_pol+1; i++ )
	{
		pol_par[i] = totalfit_result_par[i+6];
	}

	backg_to_p_ratio = backg_func->Integral( experhist_lowlimit, experhist_highlimit ) / p_func->Integral( experhist_lowlimit, experhist_highlimit);
}

void make_PolBackgHistoForSim( TH1D* h1_backg_dx, double p_count, double backg_to_p_ratio, double pol_par [n_backg_pol+1], const double low_limit, const double high_limit )
{
	int entries = (int)(p_count*backg_to_p_ratio);
	std::cout << "Backg events: " << entries << '\n';

	TF1* backg_pol = new TF1(Form("backg_pol%i", n_backg_pol), Form("pol%i", n_backg_pol), low_limit, high_limit);

	backg_pol->SetParameters(&pol_par[0]);

	for ( int i = 0; i < entries; i++ )
	{
		double randVal = backg_pol->GetRandom();

		h1_backg_dx->Fill(randVal);
	}
}

double return_Chi2() // Return the Chi2 between the experiment and scaled sumulation dx histogram.
{
	double chi2 {0.};

	double min_bin = h1_experiment_dx->GetXaxis()->FindBin(experhist_lowlimit);
	double max_bin = h1_experiment_dx->GetXaxis()->FindBin(experhist_highlimit);
	
	for ( int i = min_bin; i <= max_bin; i++ )
	{
		double real_dat = h1_experiment_dx->GetBinContent(i);
		double sim_dat = h1_simc_scaled_deeN_dx->GetBinContent(i);
		double residual = real_dat - sim_dat;

		if ( sim_dat != 0 )
		{
			chi2 += pow(residual,2)/sim_dat;
		}
	}

	return chi2;
}

// User-defined function to be minimized
void Chi2_PolBackgAddToSim(int& npar, double* deriv, double& result, double* par, int flag)
{
	h1_simc_scaled_deen_dx->Scale(par[0]);
	h1_simc_scaled_deep_dx->Scale(par[1]);
	h1_polbackg_scaled_dx->Scale(par[2]);

	h1_dx_temp->Add( h1_simc_scaled_deen_dx, h1_simc_scaled_deep_dx );
	h1_simc_scaled_deeN_dx->Add( h1_dx_temp, h1_polbackg_scaled_dx );

	// Calculate the chi-square value as the quantity to be minimized
	//result = h1_simc_scaled_deeN_dx->Chi2Test(h1_experiment_dx, "WW CHI2/NDF");

	result = return_Chi2();

	// Reset the scalling for further iterations.
	h1_simc_scaled_deen_dx->Scale(1.0/par[0]);
	h1_simc_scaled_deep_dx->Scale(1.0/par[1]);	
	h1_polbackg_scaled_dx->Scale(1.0/par[2]);

	std::cout << "Iteration Chi-square value = " << result << '\n';
}

void make_PolBackgHistoForSim1(TH1D* h1_backg_dx, double par0, double par1, double par2, double par3)
{
	TF1* backg_pol = new TF1("backg_pol3", "pol3", experhist_lowlimit, experhist_highlimit);

	backg_pol->SetParameter(0, par0);
	backg_pol->SetParameter(1, par1);
	backg_pol->SetParameter(2, par2);
	backg_pol->SetParameter(3, par3);

	double min_bin = h1_experiment_dx->GetXaxis()->FindBin(experhist_lowlimit);
	double max_bin = h1_experiment_dx->GetXaxis()->FindBin(experhist_highlimit);

	for ( int i = min_bin; i <= max_bin; i++ )
	{
		h1_backg_dx->SetBinContent( i, backg_pol->Eval(h1_experiment_dx->GetXaxis()->GetBinCenter(i)) );
	}
}

void Chi2_G4SBSbackg(int& npar, double* deriv, double& result, double* par, int flag)
{
	h1_simc_scaled_deen_dx->Scale(par[0]);
	h1_simc_scaled_deep_dx->Scale(par[1]);
	h1_g4sbsinelbackg_scaled_dx->Scale(par[2]);

	h1_dx_temp->Add( h1_simc_scaled_deen_dx, h1_simc_scaled_deep_dx );
	h1_simc_scaled_deeN_dx->Add( h1_dx_temp, h1_g4sbsinelbackg_scaled_dx );

	// Calculate the chi-square value as the quantity to be minimized
	//result = h1_simc_scaled_deeN_dx->Chi2Test(h1_experiment_dx, "WW CHI2/NDF");

	result = return_Chi2();

	// Reset the scalling for further iterations.
	h1_simc_scaled_deen_dx->Scale(1.0/par[0]);
	h1_simc_scaled_deep_dx->Scale(1.0/par[1]);	
	h1_g4sbsinelbackg_scaled_dx->Scale(1.0/par[2]);

	std::cout << "Iteration Chi-square value = " << result << '\n';
}

void Chi2_PolBackgAddToSim1(int& npar, double* deriv, double& result, double* par, int flag)
{
	h1_simc_scaled_deen_dx->Scale(par[0]);
	h1_simc_scaled_deep_dx->Scale(par[1]);
	
	make_PolBackgHistoForSim1( h1_polbackg_scaled_dx, par[2], par[3], par[4], par[5] );	

	h1_dx_temp->Add( h1_simc_scaled_deen_dx, h1_simc_scaled_deep_dx );
	h1_simc_scaled_deeN_dx->Add( h1_dx_temp, h1_polbackg_scaled_dx );

	// Calculate the chi-square value as the quantity to be minimized
	//result = h1_simc_scaled_deeN_dx->Chi2Test(h1_experiment_dx, "WW CHI2/NDF");

	result = return_Chi2();

	// Reset the scalling for further iterations.
	h1_simc_scaled_deen_dx->Scale(1.0/par[0]);
	h1_simc_scaled_deep_dx->Scale(1.0/par[1]);	
	//h1_polbackg_scaled_dx->Scale(1.0/par[2]);

	std::cout << "Iteration Chi-square value = " << result << '\n';
}


void substract_PolBackgHistoFromData( TH1D* h1_exp_dx, TH1D* h1_backg_dx, double pol_par [n_backg_pol+1], const double low_limit, const double high_limit )
{
	// Background function
	TF1* backg_func = new TF1("backg_func", Form("pol%i", n_backg_pol), low_limit, high_limit);

	backg_func->SetParameters(&pol_par[0]);

	double min_bin = h1_backg_dx->GetXaxis()->FindBin(low_limit);
	double max_bin = h1_backg_dx->GetXaxis()->FindBin(high_limit);

	for ( int i = min_bin; i <= max_bin; i++ )
	{
		double x = h1_backg_dx->GetBinCenter(i);

		double y = backg_func->Eval(x);

		double bacgk_substracted_dx = h1_exp_dx->GetBinContent(i) - y;

		if ( bacgk_substracted_dx >= 0.0 )
		{
			h1_backg_dx->SetBinContent(i, y);
			h1_exp_dx->SetBinContent(i, bacgk_substracted_dx);
		}
		else
		{
			h1_backg_dx->SetBinContent(i, 0);	
			h1_exp_dx->SetBinContent(i, 0);		
		}
	}
}

void plotOverlayDataAndSimHistos(TH1D* dataHist, TH1D* simDeen, TH1D* simDeep, TH1D* simBackg, TH1D* simFullHist, const char* dataTitle, const char* simTitle, const char* canvasTitle)
{
	TCanvas* canvas = new TCanvas(canvasTitle, canvasTitle, 800, 600);

	// // Data histogram. 
    dataHist->SetLineColor(kBlue+4);
    dataHist->SetMarkerColor(kBlue+4);
    dataHist->SetMarkerStyle(0); // Circle marker
    dataHist->SetTitle("");
    dataHist->SetStats(0);
    dataHist->Draw("HIST");
    // //

    // // Simulation Deen histogram.
    simDeen->SetFillColor(kCyan);
    simDeen->SetFillStyle(3003);
    simDeen->SetLineColorAlpha(kCyan+1, 0.5);
    simDeen->SetTitle("");
    simDeen->SetStats(0);
    simDeen->Draw("HIST+SAME");
    // //

    // // Simulation Deep histogram.
    simDeep->SetFillColor(kViolet);
    simDeep->SetFillStyle(3003);
    simDeep->SetLineColorAlpha(kViolet+1, 0.5);
    simDeep->SetTitle("");
    simDeep->SetStats(0);
    simDeep->Draw("HIST+SAME");
    // //

    // // Simulation background histogram.
    simBackg->SetFillColor(kGreen);
    simBackg->SetFillStyle(3003);
    simBackg->SetLineColorAlpha(kGreen+1, 0.5);
    simBackg->SetTitle("");
    simBackg->SetStats(0);
    simBackg->Draw("HIST+SAME");
    // //

    // // Simulation full histogram.
    // Draw the simulation histogram with a dashed red line, filled area, or different marker
    // simHist->SetLineColor(kRed);
    // SimFullHist->SetMarkerColor(kRed);
    // SimFullHist->SetMarkerStyle(21); // Square marker
    //simFullHist->SetFillColor(kRed); // Set fill color
    //simFullHist->SetFillStyle(3004); // Set fill pattern
    simFullHist->SetLineColor(kRed); // Set line color with transparency
    simFullHist->SetTitle("");
    simFullHist->SetStats(0);
    simFullHist->Draw("HIST+SAME");
    // //

    // Create a legend to label the histograms
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the position as needed
    legend->AddEntry(dataHist, dataTitle, "l");
    legend->AddEntry(simDeen, "Simulated deen events", "l");
    legend->AddEntry(simDeep, "Simulated deep events", "l");
   	legend->AddEntry(simBackg, "Simulation background", "l"); 
    legend->AddEntry(simFullHist, simTitle, "l");
    legend->Draw();

    // Update the canvas
    canvas->Draw();
}

void makeResidualPlot( TH1D* h1_exp_dx, TH1D* h1_sim_dx )
{
	// TGraph* graph = new TGraph();
	// graph->GetXaxis()->SetTitle("dx (m)");
	// graph->GetYaxis()->SetTitle("Experiment - Simulation");
	// graph->GetXaxis()->SetRangeUser(dx_low, dx_high);

	TH1D* h1_residuals_dx = new TH1D("h1_residuals_dx", "Residuals; dx (m); Experiment - Simulation", dx_nbins, dx_low, dx_high);

	double min_bin = h1_exp_dx->GetXaxis()->FindBin(experhist_lowlimit);
	double max_bin = h1_exp_dx->GetXaxis()->FindBin(experhist_highlimit);

	for ( int i = min_bin; i <= max_bin; i++ )
	{
		double residual = h1_exp_dx->GetBinContent(i) - h1_sim_dx->GetBinContent(i);

		//graph->SetPoint( i - min_bin, h1_exp_dx->GetXaxis()->GetBinCenter(i), residual );
		h1_residuals_dx->SetBinContent(i, residual);
	}

	TCanvas* canvas = new TCanvas("Residuals", "Residuals", 800, 900);

	canvas->Divide(1, 2);

	canvas->cd(1);
	h1_exp_dx->Draw("HIST");
	h1_sim_dx->Draw("HIST+SAME");
	TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the position as needed
	legend->AddEntry(h1_exp_dx, "Data", "l");
	legend->AddEntry(h1_sim_dx, "Simulation", "l");
	legend->Draw();

	canvas->cd(2);
	h1_residuals_dx->SetStats(0);
	h1_residuals_dx->SetMarkerStyle(21);
	h1_residuals_dx->Draw("");

	canvas->Draw();
}