// Script to compare the photon energy spectra with BB cuts and BB+HCal cuts, to study the ratio near the Brehmsstrahlung end-point.
#include <iostream>

constexpr int nplots = 4;
bool didpass_neutron_spotcut( const double rcut, const double dx, const double dy );
bool didpass_hcalsh_corrcut( const double corrtime, const double corrtime_lowcut, const double corrtime_highcut );
void make_ratio_graph(const TH1D* h1_photonE_denom, const TH1D* h1_photonE_numer, TCanvas* C);


void photonEspectra_compare( const char* inrootfile, const double rcut = 0.5, const double corrtime_lowcut = 43, const double corrtime_highcut = 60 )
{
	TFile* infile = new TFile( inrootfile, "READ" );
	TTree* T = (TTree*)infile->Get("T");
	
	T->SetBranchStatus( "*", 0 );
	T->SetBranchStatus( "beam.photon.e", 1 );
	T->SetBranchStatus( "adctimediff.hcalsh", 1 );
	T->SetBranchStatus( "sbs.neutronrecon.hcaldx", 1 );
	T->SetBranchStatus( "sbs.neutronrecon.hcaldy", 1 );

	double recon_photonE{0.};
	double adctimediff_hcalsh{0.};
	double dx{0.};
	double dy{0.};

	T->SetBranchAddress( "beam.photon.e", &recon_photonE );
	T->SetBranchAddress( "adctimediff.hcalsh", &adctimediff_hcalsh );
	T->SetBranchAddress( "sbs.neutronrecon.hcaldx", &dx );
	T->SetBranchAddress( "sbs.neutronrecon.hcaldy", &dy );
	
	//Plots
	TH2D* h2_dxdy_bigbitecuts = new TH2D("h2_dxdy_bigbitecuts", "HCal dx vs dy; dy = Y_{hcal}-Y_{predicted} (m); dx = X_{hcal}-X_{predicted} (m)", 400, -2, 2, 800, -4, 4);
	TH2D* h2_dxdy_bigbiteplushcalcuts = new TH2D("h2_dxdy_bigbiteplushcalcuts", "HCal dx vs dy; dy = Y_{hcal}-Y_{predicted} (m); dx = X_{hcal}-X_{predicted} (m)", 400, -2, 2, 800, -4, 4);
	const int nbins_photonespectra = 50;
	TH1D* h1_photonE_denom = new TH1D("h1_photonE_denom", "Recon. photon energy; E_{#gamma} (GeV)", nbins_photonespectra, 2.2, 4.2); //Histogram with photon E after BB cuts.
	TH1D* h1_photonE_numer = new TH1D("h1_photonE_numer", "Recon. photon energy; E_{#gamma} (GeV)", nbins_photonespectra, 2.2, 4.2); //Histogram with photon E after BB plus HCal cuts.

	long nevent = 0;
	
	while( T->GetEntry(nevent++) )
	{
		h2_dxdy_bigbitecuts->Fill(dy,dx);
		h1_photonE_denom->Fill(recon_photonE);

		bool pass_rcut = didpass_neutron_spotcut(rcut, dx, dy);
		bool pass_hcalsh_corrcut = didpass_hcalsh_corrcut(adctimediff_hcalsh, corrtime_lowcut, corrtime_highcut);
		bool both_cuts = pass_rcut&&pass_hcalsh_corrcut;
	
		if ( !both_cuts ) continue; //Neutron spot cut + HCal and BB SH ADC time corr cut.
					
		h2_dxdy_bigbiteplushcalcuts->Fill(dy,dx);
		h1_photonE_numer->Fill(recon_photonE);
	}

	TCanvas* C[nplots];

	C[0] = new TCanvas("c1", "dx vs dy with only BB cuts", 400, 800);
	h2_dxdy_bigbitecuts->Draw("COLZ");

	C[1] = new TCanvas("c2", "dx vs dy with BB+HCal cuts", 400, 800);
	h2_dxdy_bigbiteplushcalcuts->Draw("COLZ");

	C[2] = new TCanvas("c3", "Photon energies superimposed", 500, 800);
	h1_photonE_denom->SetLineColor(4);
	h1_photonE_denom->Draw();
	h1_photonE_numer->SetLineColor(8);
	h1_photonE_numer->SetFillColor(kGreen);
	h1_photonE_numer->SetFillStyle(3144);
	h1_photonE_numer->Draw("SAME");
	auto* legend = new TLegend(0.2, 0.8, 0.4, 1);
	legend->AddEntry(h1_photonE_denom, "denominator", "l");
	legend->AddEntry(h1_photonE_numer, "numerator", "l");
	legend->Draw();

	C[3] = new TCanvas("c4", "Ratios", 600, 500);
	make_ratio_graph(h1_photonE_denom, h1_photonE_numer, C[3]);
}

bool didpass_neutron_spotcut( const double rcut, const double dx, const double dy )
{
	double r = std::sqrt( std::pow(dx,2) + std::pow(dy,2) );

	if ( r <= rcut )
	{
		return true;
	}
		
	return false;
}

bool didpass_hcalsh_corrcut( const double corrtime, const double corrtime_lowcut, const double corrtime_highcut )
{
	if ( corrtime_lowcut < corrtime && corrtime < corrtime_highcut )
	{
		return true;
	}

	return false;
}

void make_ratio_graph(const TH1D* h1_photonE_denom, const TH1D* h1_photonE_numer, TCanvas* C)
{
	const int nbins = h1_photonE_denom->GetNbinsX();

	double ratios[nbins];
	double bincenter[nbins];

	for (int i = 0; i < nbins; i++)
	{
		if( h1_photonE_numer->GetBinContent(i+1) != 0 && h1_photonE_denom->GetBinContent(i+1) != 0 )
		{
			ratios[i] = h1_photonE_numer->GetBinContent(i+1) / h1_photonE_denom->GetBinContent(i+1);	
		}
		else
		{
			std::cout << h1_photonE_numer->GetBinContent(i+1) << " / " << h1_photonE_denom->GetBinContent(i+1) << '\n';
			ratios[i] = 0;
		}
		
		bincenter[i] = h1_photonE_denom->GetBinCenter(i+1);

		std::cout << i << '\n';
	}

	TGraph* g1_photonespectra_ratios = new TGraph(nbins, bincenter, ratios);
	C->cd();
	g1_photonespectra_ratios->SetTitle("Bin content ratios");
	g1_photonespectra_ratios->GetXaxis()->SetTitle("E_{#gamma} GeV");
	g1_photonespectra_ratios->GetYaxis()->SetTitle("Ratio");
	g1_photonespectra_ratios->Draw("AL*");
}

