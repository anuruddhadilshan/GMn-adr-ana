// Macro to do Gaus fits for p and n peaks of LD2 data dx plots and extract their mean positions. 
// Written for the specific task of calibrating the G4SBS sbsfield scale factor.



void subrange_ld2dxplots_gausfits( const char* rootfile, int dat_type = 1 ) // Real data = 0, Simulation data = 1.
{

	// Load the root file and get the TTree object "T"
	TFile* anarootfile = new TFile(rootfile, "READ");
	TTree* T = (TTree*)anarootfile->Get("T");

	// Apply W2 and dy cuts to clean-up background from the dx plot.
	const double cut_w2_low {-1};
	const double cut_w2_high {1.5};
	const double cut_dy_low {-0.5};
	const double cut_dy_high {0.5};

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("e.kine.W2", 1);
	T->SetBranchStatus("sbs.hcal.bestclus.dx", 1);
	T->SetBranchStatus("sbs.hcal.bestclus.dy", 1);
	if ( dat_type ) T->SetBranchStatus("MC.mc_weight", 1);

	double w2 {0.};
	double dx {0.};
	double dy {0.};
	double mc_weight {0.};

	T->SetBranchAddress("e.kine.W2", &w2);
	T->SetBranchAddress("sbs.hcal.bestclus.dx", &dx);
	T->SetBranchAddress("sbs.hcal.bestclus.dy", &dy);
	if ( dat_type ) T->SetBranchAddress("MC.mc_weight", &mc_weight);

	TH1D* h1_W2 = new TH1D("h1_W2", "W^{2}/(GeV^{2}/c^{2})", 175, -2, 5);
	TH1D* h1_dy = new TH1D("h1_dy", "dy; dy (m)", 75, -1.5, 1.5);
	const int nbins_dx {250};
	TH1D* h1_dx = new TH1D("h1_dx", "dx; dx (m)", 125, -2.5, 2.5);
		
	long nevent {0}; // Keep track of the event number in the loop.

	while ( T->GetEntry(nevent) )
	{
		nevent++;

		//////////////////////////////////////////////////////////////
		if ( w2 < cut_w2_low || w2 > cut_w2_high ) continue; // W2 cut
		if ( dy < cut_dy_low || dy > cut_dy_high ) continue; // dy cut
		//////////////////////////////////////////////////////////////
		
		h1_W2->Fill(w2);

		if ( dat_type ) 
		{
			h1_dx->Fill(dx, mc_weight);
			h1_dy->Fill(dy, mc_weight);
		}
		else 
		{
			h1_dx->Fill(dx);
			h1_dy->Fill(dy);
		}		
		
	}

	TCanvas* C1 = new TCanvas("C1", "C1", 800, 800);
	C1->Divide(2,1);
	C1->cd(1);
	h1_W2->Draw();
	C1->cd(2);
	h1_dy->Draw("HIST");
	TF1* dy_fit = new TF1("dy_fit", "gaus", -0.2, 0.2);
	h1_dy->Fit(dy_fit, "R");
	dy_fit->Draw("SAME");
	C1->Update();

	TCanvas* C2 = new TCanvas("C2", "C2", 800, 800);
	h1_dx->Draw("HIST");
	C2->Update();

	std::cout << "Enter an approximate demarcation value between the proton and the neutron peak: ";
	double demarcation_val {0.};
	std::cin >> demarcation_val;
	std::cout << '\n';

	double max_bincontent = h1_dx->GetMaximum();
	TLine* demarcation_line = new TLine(demarcation_val,0.0,demarcation_val,max_bincontent);
	demarcation_line->SetLineStyle(2);
	demarcation_line->Draw("SAME");

	//// 1) Get the maximum bin content's x pos values of the proton and neutron side of the demarcation line.
	int demarcation_bin = h1_dx->GetXaxis()->FindBin(demarcation_val); // Bin number of the bin that the demarcation line goes through.
	int min_bin = h1_dx->GetXaxis()->GetFirst(); // First non-zero bin.
	int max_bin = h1_dx->GetXaxis()->GetLast(); // Last non-zero bin.

	double ppeak_maxbincontent = h1_dx->GetBinContent(min_bin);
	double ppeak_maxbin = min_bin;

	for ( int i = min_bin+1; i < demarcation_bin; i++ )
	{
		if ( h1_dx->GetBinContent(i) > ppeak_maxbincontent )
		{
			ppeak_maxbincontent = h1_dx->GetBinContent(i);
			ppeak_maxbin = i;
		}
	}

	double npeak_maxbincontent = h1_dx->GetBinContent(demarcation_bin);
	double npeak_maxbin = demarcation_bin;

	for ( int i = demarcation_bin; i <= max_bin; i++ )
	{
		if ( h1_dx->GetBinContent(i) > npeak_maxbincontent )
		{
			npeak_maxbincontent = h1_dx->GetBinContent(i);
			npeak_maxbin = i;
		}
	} 

	std::cout << ppeak_maxbin << "   " << ppeak_maxbincontent << '\n';
	std::cout << npeak_maxbin << "   " << npeak_maxbincontent << '\n';

	//// 2) Now generate the ranges for the first order fitting. //
	// p peak
	double firstorder_ppeak_lowlimit = h1_dx->GetXaxis()->GetBinCenter(ppeak_maxbin) - 0.4;
	double firstorder_ppeak_highlimit = h1_dx->GetXaxis()->GetBinCenter(ppeak_maxbin) + 0.4;

	// n peak
	double firstorder_npeak_lowlimit = h1_dx->GetXaxis()->GetBinCenter(npeak_maxbin) - 0.4;
	double firstorder_npeak_highlimit = h1_dx->GetXaxis()->GetBinCenter(npeak_maxbin) + 0.4;

	//// 3) Now do a first order gaus fit to the proton and neutron peak.
	TF1* firstorder_pfit = new TF1("firstorder_pfit", "gaus", firstorder_ppeak_lowlimit, firstorder_ppeak_highlimit);
	h1_dx->Fit(firstorder_pfit, "R");
	firstorder_pfit->Draw("SAME");

	TF1* firstorder_nfit = new TF1("firstorder_nfit", "gaus", firstorder_npeak_lowlimit, firstorder_npeak_highlimit);
	h1_dx->Fit(firstorder_nfit, "R+");
	firstorder_nfit->Draw("SAME");
	C2->Update();

	//// 4) Now generate more tighter ranges for a second order fitting.
	// p peak
	double firstorder_pfit_par [3];
	firstorder_pfit->GetParameters(&firstorder_pfit_par[0]);
	double secondorder_ppeak_lowlimit = firstorder_pfit_par[1] - 1.0*firstorder_pfit_par[2];
	double secondorder_ppeak_highlimit = firstorder_pfit_par[1] + 0.8*firstorder_pfit_par[2];

	// n peak
	double firstorder_nfit_par [3];
	firstorder_nfit->GetParameters(&firstorder_nfit_par[0]);
	double secondorder_npeak_lowlimit = firstorder_nfit_par[1] - 0.8*firstorder_nfit_par[2];
	double secondorder_npeak_highlimit = firstorder_nfit_par[1] + 1.0*firstorder_nfit_par[2];

	
	//// 5) Define a two second order fit functions for the proton and the neutron peaks and initilize their paramters using the paremters of the first order fit functions.
	
	TCanvas* C3 = new TCanvas("C3", "C3", 800, 800);
	h1_dx->Draw("HIST");

	TF1* secondorder_pfit = new TF1("secondorder_pfit", "gaus", secondorder_ppeak_lowlimit, secondorder_ppeak_highlimit);
	secondorder_pfit->SetParameters(firstorder_pfit_par);
	h1_dx->Fit(secondorder_pfit, "R+");
	secondorder_pfit->Draw("SAME");

	TF1* secondorder_nfit = new TF1("secondorder_nfit", "gaus", secondorder_npeak_lowlimit, secondorder_npeak_highlimit);
	secondorder_nfit->SetParameters(firstorder_nfit_par);
	h1_dx->Fit(secondorder_nfit, "R+");
	secondorder_nfit->Draw("SAME");

	C3->Update();

	// Printout to the screen with the results.
	std::cout << "\nProton peak dx position: " << secondorder_pfit->GetParameter(1) << " +/- " << secondorder_pfit->GetParError(1) << " meters\n";
	std::cout << "\nNeutron peak dx position: " << secondorder_nfit->GetParameter(1) << " +/- " << secondorder_nfit->GetParError(1) << " meters\n";
	std::cout << "\ndx proton peak std.dev: " << secondorder_pfit->GetParameter(2) << " +/- " << secondorder_pfit->GetParError(2) << " meters\n";
	std::cout << "\ndy peak std.dev: " << dy_fit->GetParameter(2) << " +/- " << dy_fit->GetParError(2) << " meters\n";

}