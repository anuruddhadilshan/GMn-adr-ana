// Script to show-case the inelastic contamination in the real experiment data.

void show_inelastic_backg( const char* rootfile )
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
	TH1D* h1_w2 = new TH1D("h1_w2", "W2 (GeV^{2})", 175, -2, 5);
	TH1D* h1_dy = new TH1D("h1_dy", "dy (m)", 175, -3.5, 3.5);
	TH2D* h2_dxw2 = new TH2D("h2_dxw2", "dx vs W^{2}; W^{2} (GeV^{2}); dx (m)", 175, -2, 5, 150, -3.5, 2.5); // Shows the elastic/inelastic distributions throughout the dx plots.
	TH1D* h1_dx = new TH1D("h1_dx", "dx (m)", 150, -3.5, 2.5); // dx plot before removing inelastics outside of the respective elastic peaks in W2 and dy plots.

	long nevent {0};

	while ( T->GetEntry(nevent) )
	{
		nevent++;

		h1_w2->Fill(w2);
		h1_dy->Fill(dy);
		h2_dxw2->Fill(w2, dx);
		h1_dx->Fill(dx);
	}

	TCanvas* C1 = new TCanvas("C1", "Show Inelastic Background", 1200, 900);
	C1->Divide(2,2);

	C1->cd(1);
	h1_w2->Draw();
	TLine* l_W2_h1w2 = new TLine( 0.88, 0, 0.88, h1_w2->GetMaximum()*1.05 );
	l_W2_h1w2->SetLineWidth(2);
	l_W2_h1w2->SetLineColor(kRed);
	l_W2_h1w2->SetLineStyle(kDashed);
	l_W2_h1w2->Draw("SAME");
	

	C1->cd(2);
	h1_dy->Draw();
	TLine* l_dy_h1dy = new TLine( 0.0, 0, 0.0, h1_dy->GetMaximum()*1.05 );
	l_dy_h1dy->SetLineWidth(2);
	l_dy_h1dy->SetLineColor(kMagenta);
	l_dy_h1dy->SetLineStyle(kDashed);
	l_dy_h1dy->Draw("SAME");

	C1->cd(3);
	h2_dxw2->Draw("COLZ");
	TLine* l_W2_h2dxw2 = new TLine( 0.88, -3.5, 0.88, 2.5 );
	l_W2_h2dxw2->SetLineWidth(2);
	l_W2_h2dxw2->SetLineColor(kRed);
	l_W2_h2dxw2->SetLineStyle(kDashed);
	l_W2_h2dxw2->Draw("SAME");

	C1->cd(4);
	h1_dx->Draw();

}