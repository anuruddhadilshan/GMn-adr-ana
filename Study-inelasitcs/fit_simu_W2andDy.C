// This script fits 'W2' and 'dy' distributions of simulated data and obtain lover and upper bounds of elastic peak distributions in W2 and dy plots.
// The idea is to more accurately select most of the elastic events available when analyzing real data and reject unnecessary inelastic background distributions that live outside of the range of the true elastic distributions.

void fit_simu_W2andDy( const char* rootfile )
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
	TH1D* h1_w2 = new TH1D("h1_w2", "W2 (GeV^{2}; W^{2}; W^{2} (GeV^{2})", 150, -2, 5);
	TH1D* h1_dy = new TH1D("h1_dy", "dy; dy (m)", 150, -3.5, 3.5);
	
	long nevent {0};

	while ( T->GetEntry(nevent) )
	{
		nevent++;

		h1_w2->Fill(w2);
		h1_dy->Fill(dy);
	}

	gStyle->SetOptFit();

	// 1) Perform a first order Gaus fit for W2 and dy.
	TF1* firstorder_W2fit = new TF1("firstorder_W2fit", "gaus", -2, 5);
	h1_w2->Fit(firstorder_W2fit, "R");
	double firstorder_W2fit_par [3] {0.};
	firstorder_W2fit->GetParameters(&firstorder_W2fit_par[0]);

	TF1* firstorder_Dyfit = new TF1("firstorder_Dyfit", "gaus", -3.5, 3.5);
	h1_dy->Fit(firstorder_Dyfit, "R");
	double firstorder_Dyfit_par [3] {0.};
	firstorder_Dyfit->GetParameters(&firstorder_Dyfit_par[0]);


	// 2) Derive the sub-ranges to do a second order fitting.
	double lowlimit_w2 = firstorder_W2fit_par[1] - firstorder_W2fit_par[2]*1.0;
	double highlimit_w2 = firstorder_W2fit_par[1] + firstorder_W2fit_par[2]*1.0;
	double lowlimit_dy = firstorder_Dyfit_par[1] - firstorder_Dyfit_par[2]*1.0;
	double highlimit_dy = firstorder_Dyfit_par[1] + firstorder_Dyfit_par[2]*1.0;


	// 3) Perform a second order Gaus fit for W2 and dy after initilizing the fit parameters using the first order fit.
	TF1* secondorder_W2fit = new TF1("secondorder_W2fit", "gaus", lowlimit_w2, highlimit_w2 );
	secondorder_W2fit->SetParameters(firstorder_W2fit_par);
	h1_w2->Fit(secondorder_W2fit, "R");
	double secondorder_W2fit_par [3] {0.};
	secondorder_W2fit->GetParameters(&secondorder_W2fit_par[0]);

	TF1* secondorder_Dyfit = new TF1("secondorder_Dyfit", "gaus", lowlimit_dy, highlimit_dy );
	secondorder_Dyfit->SetParameters(firstorder_Dyfit_par);
	h1_dy->Fit(secondorder_Dyfit, "R");
	double secondorder_Dyfit_par [3] {0.};
	secondorder_Dyfit->GetParameters(&secondorder_Dyfit_par[0]);


	// 4) Get W2 and dy cut boundaries.
	double low_elastic_w2 = secondorder_W2fit_par[1] - secondorder_W2fit_par[2]*3.5;
	double high_elastic_w2 = secondorder_W2fit_par[1] + secondorder_W2fit_par[2]*3.5;
	double low_elastic_dy = secondorder_Dyfit_par[1] - secondorder_Dyfit_par[2]*4.5;
	double high_elastic_dy = secondorder_Dyfit_par[1] + secondorder_Dyfit_par[2]*4.5;

	TCanvas* C1 = new TCanvas("C1", "Generate elastic-cut boundaries from simulated data", 800, 600);
	C1->Divide(2,1);

	C1->cd(1);
	h1_w2->Draw();
	secondorder_W2fit->Draw("SAME");
	TLine* l_w2low = new TLine(low_elastic_w2, 0, low_elastic_w2, h1_w2->GetMaximum()*1.05);
	l_w2low->SetLineWidth(2);
	l_w2low->SetLineColor(kRed);
	l_w2low->SetLineStyle(kDashed);
	l_w2low->Draw("SAME");
	TLine* l_w2high = new TLine(high_elastic_w2, 0, high_elastic_w2, h1_w2->GetMaximum()*1.05);
	l_w2high->SetLineWidth(2);
	l_w2high->SetLineColor(kRed);
	l_w2high->SetLineStyle(kDashed);
	l_w2high->Draw("SAME");

	C1->cd(2);
	h1_dy->Draw();
	secondorder_Dyfit->Draw("SAME");
	TLine* l_dylow = new TLine(low_elastic_dy, 0, low_elastic_dy, h1_dy->GetMaximum()*1.05);
	l_dylow->SetLineWidth(2);
	l_dylow->SetLineColor(kRed);
	l_dylow->SetLineStyle(kDashed);
	l_dylow->Draw("SAME");
	TLine* l_dyhigh = new TLine(high_elastic_dy, 0, high_elastic_dy, h1_dy->GetMaximum()*1.05);
	l_dyhigh->SetLineWidth(2);
	l_dyhigh->SetLineColor(kRed);
	l_dyhigh->SetLineStyle(kDashed);
	l_dyhigh->Draw("SAME");

	std::cout << "\n*** W2 Elastic Cut region: " << low_elastic_w2 << " < W2 < " << high_elastic_w2 << '\n';
	std::cout << "\n*** dy Elastic Cut region: " << low_elastic_dy << " < dy < " << high_elastic_dy << '\n';

}