// This script applies final elastic cuts on W2 and dy. 
// Makes a dx plot. Fit n, p, and backg with two Gausians and a 4th order polynomial.
// Extract the n and p counts from the final fit integrals and calculates R.
// This R will be of course very approximate. The elastic yields from the step above can be used to calculate the statistical uncertainty.


void fit_dx_getR( const char* rootfile, const double low_w2 = -0.593086, const double high_w2 = 2.35375, const double low_dy = -0.562647, const double high_dy = 0.568397 )
{
	// Load the root file and get the TTree "T" created from elastic event selection analysis scripts.
	TFile* anarootfile = new TFile(rootfile, "READ");
	TTree* T = (TTree*)anarootfile->Get("T"); 

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("e.kine.W2", 1);
	T->SetBranchStatus("sbs.hcal.bestclus.dx", 1);
	T->SetBranchStatus("sbs.hcal.bestclus.dy", 1);

	double w2 {0.};
	double dx {0.};
	double dy {0.};

	T->SetBranchAddress("e.kine.W2", &w2);
	T->SetBranchAddress("sbs.hcal.bestclus.dx", &dx);
	T->SetBranchAddress("sbs.hcal.bestclus.dy", &dy);

	// Histogram definitions //
	TH1D* h1_w2_noelascut = new TH1D("h1_w2_noelascut", "W2 - no elastic cuts; W^{2} (GeV^{2})", 150, -2, 5);
	TH1D* h1_dy_noelastcut = new TH1D("h1_dy_noelastcut", "dy - no elastic cuts; dy (m)", 150, -3.5, 3.5);
	TH2D* h2_dxw2_noelastcut = new TH2D("h2_dxw2_noelastcut", "dx vs W^{2} - no elastic cuts; W^{2} (GeV^{2}); dx (m)", 175, -2, 5, 150, -3.5, 2.5); // Shows the elastic/inelastic distributions throughout the dx plots.
	TH1D* h1_dx_noelastcut = new TH1D("h1_dx_noelastcut", "dx - no elastic cuts; dx (m)", 150, -3.5, 2.5); // dx plot before removing inelastics outside of the respective elastic peaks in W2 and dy plots.

	TH1D* h1_w2_withelascut = new TH1D("h1_w2_withelascut", "W2 - with elastic cuts; W^{2} (GeV^{2})", 150, -2, 5);
	TH1D* h1_dy_withelastcut = new TH1D("h1_dy_withelastcut", "dy - with elastic cuts; dy (m)", 150, -3.5, 3.5);
	TH2D* h2_dxw2_withelastcut = new TH2D("h2_dxw2_withelastcut", "dx vs W^{2} - with elastic cuts; W^{2} (GeV^{2}); dx (m)", 175, -2, 5, 150, -3.5, 2.5); // Shows the elastic/inelastic distributions throughout the dx plots.
	TH1D* h1_dx_withelastcut = new TH1D("h1_dx_withelastcut", "dx - with final elastic cuts; dx (m)", 150, -3.5, 2.5); // dx plot after removing inelastics outside of the respective elastic peaks in W2 and dy plots.

	
	long nevent {0};

	while ( T->GetEntry(nevent++) )
	{
		h1_w2_noelascut->Fill(w2);
		h1_dy_noelastcut->Fill(dy);
		h2_dxw2_noelastcut->Fill(w2,dx);
		h1_dx_noelastcut->Fill(dx);

		if ( w2 < low_w2 || w2 > high_w2 ) continue; // W2 cut.
		if ( dy < low_dy || dy > high_dy ) continue; // dy cut.

		h1_w2_withelascut->Fill(w2);
		h1_dy_withelastcut->Fill(dy);
		h2_dxw2_withelastcut->Fill(w2,dx);
		h1_dx_withelastcut->Fill(dx);
	}


	TCanvas* C1 = new TCanvas("C1", "Show Elastic Cuts on W^{2} and dy", 1200, 900);
	C1->Divide(1,2);
	C1->cd(1);
	TPad* subPad1 = (TPad*)gPad;
	subPad1->SetFillColor(kGreen-10);
	subPad1->SetFillStyle(3003);
	subPad1->Divide(2,1);
	subPad1->cd(1);
	h1_w2_noelascut->Draw();
	subPad1->cd(2);
	h1_w2_withelascut->Draw();
	C1->cd(2);
	TPad* subPad2 = (TPad*)gPad;
	subPad2->SetFillColor(kBlue-10);
	subPad2->SetFillStyle(3003);
	subPad2->Divide(2,1);
	subPad2->cd(1);
	h1_dy_noelastcut->Draw();
	subPad2->cd(2);
	h1_dy_withelastcut->Draw();
	C1->Update();

	TCanvas* C2 = new TCanvas("C2", "Compare Inelastic Background", 1200, 900);
	C2->Divide(1,2);
	C2->cd(1);
	TPad* subPad3 = (TPad*)gPad;
	subPad3->SetFillColor(kYellow-10);
	subPad3->SetFillStyle(3003);
	subPad3->Divide(2,1);
	subPad3->cd(1);
	h2_dxw2_noelastcut->Draw("COLZ");
	subPad3->cd(2);
	h2_dxw2_withelastcut->Draw("COLZ");
	C2->cd(2);
	TPad* subPad4 = (TPad*)gPad;
	subPad4->SetFillColor(kRed-10);
	subPad4->SetFillStyle(3003);
	subPad4->Divide(2,1);
	subPad4->cd(1);
	h1_dx_noelastcut->Draw();
	subPad4->cd(2);
	h1_dx_withelastcut->Draw();
	C2->Update();

	//// Begin fitting analysis to extract R ////

	TCanvas* C3 = new TCanvas("C3", "dx Fitting", 1200, 900);
	C3->Divide(2,1);
	C3->cd(1);
	h1_dx_withelastcut->Draw();
	C3->Update();

	// 1) Let's do a first order Gaus fitting for n and p peaks.
	std::cout << "Enter an approximate demarcation value between the proton and the neutron peak: ";
	double demarcation_val {0.};
	std::cin >> demarcation_val;
	std::cout << '\n';

	double maxbincontent_dxhist = h1_dx_withelastcut->GetMaximum();
	TLine* demarcation_line = new TLine(demarcation_val,0.0,demarcation_val,maxbincontent_dxhist);
	demarcation_line->SetLineStyle(2);
	demarcation_line->Draw("SAME");

	// 1.1) Get the maximum bin content's x pos values of the proton and neutron side of the demarcation line.
	int demarcation_bin = h1_dx_withelastcut->GetXaxis()->FindBin(demarcation_val); // Bin number of the bin that the demarcation line goes through.
	int min_bin = h1_dx_withelastcut->GetXaxis()->GetFirst(); // First non-zero bin.
	int max_bin = h1_dx_withelastcut->GetXaxis()->GetLast(); // Last non-zero bin.
	double ppeak_maxbincontent = h1_dx_withelastcut->GetBinContent(min_bin);
	double ppeak_maxbin = min_bin;

	for ( int i = min_bin+1; i < demarcation_bin; i++ )
	{
		if ( h1_dx_withelastcut->GetBinContent(i) > ppeak_maxbincontent )
		{
			ppeak_maxbincontent = h1_dx_withelastcut->GetBinContent(i);
			ppeak_maxbin = i;
		}
	}

	double npeak_maxbincontent = h1_dx_withelastcut->GetBinContent(demarcation_bin);
	double npeak_maxbin = demarcation_bin;

	for ( int i = demarcation_bin; i <= max_bin; i++ )
	{
		if ( h1_dx_withelastcut->GetBinContent(i) > npeak_maxbincontent )
		{
			npeak_maxbincontent = h1_dx_withelastcut->GetBinContent(i);
			npeak_maxbin = i;
		}
	} 

	// 1.2) Now generate the ranges for the first order fitting. 
	// p peak
	double ppeak_1_lowlimit = h1_dx_withelastcut->GetXaxis()->GetBinCenter(ppeak_maxbin) - 0.5;
	double ppeak_1_highlimit = h1_dx_withelastcut->GetXaxis()->GetBinCenter(ppeak_maxbin) + 0.25;

	// n peak
	double npeak_1_lowlimit = h1_dx_withelastcut->GetXaxis()->GetBinCenter(npeak_maxbin) - 0.25;
	double npeak_1_highlimit = h1_dx_withelastcut->GetXaxis()->GetBinCenter(npeak_maxbin) + 0.5;

	// 1.3) Now do a first order gaus fit to the proton and neutron peak.
	TF1* pfit_1 = new TF1("pfit_1", "gaus", ppeak_1_lowlimit, ppeak_1_highlimit);
	h1_dx_withelastcut->Fit(pfit_1, "RQN");
	pfit_1->Draw("SAME");

	TF1* nfit_1 = new TF1("nfit_1", "gaus", npeak_1_lowlimit, npeak_1_highlimit);
	h1_dx_withelastcut->Fit(nfit_1, "RQN");
	nfit_1->Draw("SAME");
	C3->Update();


	// 2) Performing second order Gasu fits for n and p Gaus fits.
	// 2.1) Now generate more tighter ranges for a second order fitting.
	// p peak
	double pfit_1_par [3] {0.};
	pfit_1->GetParameters(&pfit_1_par[0]);
	double ppeak_2_lowlimit = pfit_1_par[1] - 0.8*pfit_1_par[2];
	double ppeak_2_highlimit = pfit_1_par[1] + 0.75*pfit_1_par[2];

	// n peak
	double nfit_1_par [3] {0.};
	nfit_1->GetParameters(&nfit_1_par[0]);
	double npeak_2_lowlimit = nfit_1_par[1] - 0.45*nfit_1_par[2];
	double npeak_2_highlimit = nfit_1_par[1] + 0.8*nfit_1_par[2];

	// 2.2) Define two second order fit functions for the proton and the neutron peaks and initilize their paramters using the paremters of the first order fit functions.
	C3->cd(2);
	h1_dx_withelastcut->Draw();

	TF1* pfit_2 = new TF1("pfit_2", "gaus", ppeak_2_lowlimit, ppeak_2_highlimit);
	pfit_2->SetParameters(pfit_1_par);
	h1_dx_withelastcut->Fit(pfit_2, "RQN");
	pfit_2->Draw("SAME");

	TF1* nfit_2 = new TF1("nfit_2", "gaus", npeak_2_lowlimit, npeak_2_highlimit);
	nfit_2->SetParameters(nfit_1_par);
	h1_dx_withelastcut->Fit(nfit_2, "RQN");
	nfit_2->Draw("SAME");

	// 3.0) Now perform a first order pol4 fit for the background.
	TF1* backgfit_1 = new TF1("backgfit_1", "pol4"); //, h1_dx_withelastcut->GetXaxis()->GetBinCenter(min_bin), h1_dx_withelastcut->GetXaxis()->GetBinCenter(max_bin));
	backgfit_1->SetRange(h1_dx_withelastcut->GetXaxis()->GetBinCenter(min_bin), pfit_1_par[1] - 2.5*pfit_1_par[2]);
	backgfit_1->SetRange(nfit_1_par[1] + 2.5*nfit_1_par[2], h1_dx_withelastcut->GetXaxis()->GetBinCenter(max_bin), "R+");
	h1_dx_withelastcut->Fit(backgfit_1, "RQN");
	backgfit_1->Draw("SAME");
	C3->Update();

	// 4.0) Define a total fit function and with two Gasians for the n and p and a pol4 for the background. 
	TF1* totalfit = new TF1("totalfit", Form("gaus(0)+gaus(3)+pol4(6)"), h1_dx_withelastcut->GetXaxis()->GetBinCenter(min_bin), h1_dx_withelastcut->GetXaxis()->GetBinCenter(max_bin));
	
	// 4.1) Get the parameters from the two second order Gaus fits for n and p and the first order pol4 fit for the background.
	double totalfit_1_par [11] {0.};
	nfit_2->GetParameters(&totalfit_1_par[0]);
	pfit_2->GetParameters(&totalfit_1_par[3]);
	backgfit_1->GetParameters(&totalfit_1_par[6]);

	totalfit->SetParameters(totalfit_1_par);

	// 4.2) Apply the total fit.
	TCanvas* C4 = new TCanvas("C4", "Final Fits", 1200, 900);
	C4->Divide(2,1);
	C4->cd(1);
	h1_dx_withelastcut->Draw();
	h1_dx_withelastcut->Fit(totalfit, "RQN");
	totalfit->Draw("SAME");

	TLegend* legend_totalfitresultsplot = new TLegend(0.65,0.75,0.9,0.9);
	legend_totalfitresultsplot->AddEntry(h1_dx_withelastcut, "dx histogram", "l");
	legend_totalfitresultsplot->AddEntry(totalfit, "Total Fit", "l");
	legend_totalfitresultsplot->Draw("SAME");
	C4->Update();

	// 4.2) Create n and p Gaus and backg pol4 TF1 functions using the fit result parameters and draw them.
	double totalfit_result_par [11] {0.};
	totalfit->GetParameters(&totalfit_result_par[0]);

	// n function
	TF1* n_func = new TF1("n_func", "gaus", h1_dx_withelastcut->GetXaxis()->GetBinCenter(min_bin), h1_dx_withelastcut->GetXaxis()->GetBinCenter(max_bin));
	n_func->SetParameters(&totalfit_result_par[0]);
	n_func->SetLineColor(kGreen);

	// p function
	TF1* p_func = new TF1("p_func", "gaus", h1_dx_withelastcut->GetXaxis()->GetBinCenter(min_bin), h1_dx_withelastcut->GetXaxis()->GetBinCenter(max_bin));
	p_func->SetParameters(&totalfit_result_par[3]);
	p_func->SetLineColor(kOrange+7);

	// Background function
	TF1* backg_func = new TF1("backg_func", "pol4", h1_dx_withelastcut->GetXaxis()->GetBinCenter(min_bin), h1_dx_withelastcut->GetXaxis()->GetBinCenter(max_bin));
	backg_func->SetParameters(&totalfit_result_par[6]);
	backg_func->SetLineColor(kViolet);

	C4->cd(2);
	h1_dx_withelastcut->Draw();
	n_func->Draw("SAME");
	p_func->Draw("SAME");
	backg_func->Draw("SAME");

	TLegend* legend_individfitfuncplot = new TLegend(0.65,0.75,0.9,0.9);
	legend_individfitfuncplot->AddEntry(h1_dx_withelastcut, "dx histogram", "l");
	legend_individfitfuncplot->AddEntry(n_func, "Neutron gaus", "l");
	legend_individfitfuncplot->AddEntry(p_func, "Proton gaus", "l");
	legend_individfitfuncplot->AddEntry(backg_func, "Background pol4", "l");
	legend_individfitfuncplot->Draw("SAME");

	TCanvas* C5 = new TCanvas("C5", "Histograms from fit functions", 1200, 900);
	h1_dx_withelastcut->Draw();
	TH1* hist_from_totfitfunc = (TH1*)totalfit->GetHistogram();
	hist_from_totfitfunc->Draw("SAME");
	TH1* hist_from_pfunc = (TH1*)p_func->GetHistogram();
	hist_from_pfunc->Draw("SAME");
	TH1* hist_from_nfunc = (TH1*)n_func->GetHistogram();
	hist_from_nfunc->Draw("SAME");
	TLegend* legend_histfromfitfuncplot = new TLegend(0.65,0.75,0.9,0.9);
	legend_histfromfitfuncplot->AddEntry(h1_dx_withelastcut, "dx histogram", "l");
	legend_histfromfitfuncplot->AddEntry(hist_from_totfitfunc, "Hist from tot-fit func", "l");
	legend_histfromfitfuncplot->AddEntry(hist_from_nfunc, "Hist from n func", "l");
	legend_histfromfitfuncplot->AddEntry(hist_from_pfunc, "Hist from p func", "l");
	legend_histfromfitfuncplot->Draw("SAME");

	std::cout << "\n*** Elastic yield and R calculations from the histograms derived from the fit results ***\n";
	double nNeutrons_fromhist = hist_from_nfunc->Integral();
	std::cout << "Elastic neutron yield = " << nNeutrons_fromhist << '\n';
	double nProtons_fromhist = hist_from_pfunc->Integral();
	std::cout << "Elastic proton  yield = " << nProtons_fromhist << '\n';
	std::cout << "Total elastic yield = " << nNeutrons_fromhist + nProtons_fromhist << '\n';
	double nOverp_yieldratio_fromhists = nNeutrons_fromhist / nProtons_fromhist;
	std::cout << " * The n/p yield ratio  R = " << nOverp_yieldratio_fromhists << '\n';
	
	std::cout << "\n*** R calculations directly from the fit results ***\n";
	double area_nPeak = n_func->Integral(h1_dx_withelastcut->GetXaxis()->GetBinCenter(min_bin), h1_dx_withelastcut->GetXaxis()->GetBinCenter(max_bin));
	double area_pPeak = p_func->Integral(h1_dx_withelastcut->GetXaxis()->GetBinCenter(min_bin), h1_dx_withelastcut->GetXaxis()->GetBinCenter(max_bin));
	double nOverp_yieldratio_fromfitfunc = area_nPeak / area_pPeak;
	std::cout << "Neutron peak Gaussian's integral = " << area_nPeak << '\n';
	std::cout << "Proton  peak Gaussian's integral = " << area_pPeak << '\n';
	std::cout << "Reduced chi2 (chi2/NDF) of the total fit = " <<  totalfit->GetChisquare() / totalfit->GetNDF() << '\n';
	std::cout << " * The n/p yield ratio  R = " << nOverp_yieldratio_fromfitfunc << '\n';	

}