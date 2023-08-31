#include <iostream>
#include <string>
#include "../includes/beam_variables.h"
#include "../includes/HCalConstants.h"
#include "read_dwnbndana_config.h"
#include "eventclass.h"
#include "outputclass.h"

void print_analysis_percentage(double cuurentevent_ana_percentage, int& previousevent_ana_percentage_int);
// void make_pdf(TCanvas* c[nCanvas], TString output_PDF_filename);

void ana_dwnbnd(const char* configfilename, const char* outputfilename = "dwnbendanaplots")
{
	Configfile configfile{};
	configfile.readin_dwnbndana_configfile(configfilename);
	TChain* C = configfile.return_TChain();

	int num_SBSkine = configfile.return_SBSKineNum();
	// Copying cut thresholds to local variables //
	double cut_MinGEMHitsOnTrack = configfile.return_MinGEMHitsOnTrack();
	double cut_VzCutUpStream = configfile.return_VzCutUpStream();
	double cut_VzCutDwnStream = configfile.return_VzCutDwnStream();
	double cut_HCalE = configfile.return_HCalECut();
	double cut_CoinCutLow = configfile.return_CoinCutLow();
	double cut_CoinCutHigh = configfile.return_CoinCutHigh();
	double cut_RCut = configfile.return_RCut();
		
	// Define an object from the event class.
	Event event{ C, num_SBSkine, cut_MinGEMHitsOnTrack, cut_VzCutUpStream, cut_VzCutDwnStream, cut_HCalE, cut_CoinCutLow, cut_CoinCutHigh, cut_RCut };

	// Define an object fomr the Output class.
	Output output{ event, outputfilename };
	
	//Variables to keep track of printout to the terminal.
	int previousevent_ana_percentage_int{0};
	long nevents {C->GetEntries()};
	std::cout << "Number of events in the TChain: " << nevents << '\n';

	std::cout << "\n--- Beginning down-bending track analysis ---\n";

	long nevent{0};

	while (event.getEntry(nevent++)) 
	{
		double ana_percentage{(nevent/(double)nevents)*100}; //Percentage of events analyzed in the Event List.
		print_analysis_percentage(ana_percentage, previousevent_ana_percentage_int);
		
		event.calc_BBTrackAngles(); // Calculates the polar and azimuthal scattering angle of the BigBite track.		
		event.calc_PhotonE(); // Calculates the photon energy from the reaction: photon + p --> pi+ + n , using the measured BigBite track momentum.
		event.calc_NeutronKin(); // Calculates the kinematics of the neutron for the reaction: photon + p --> pi+ + n, using bb.tr momentum information and reconstructed photon energy.
		event.calc_NeutronHCalIntersect(); // Calculates the "predicted hit position" of the neutron using the calculated neutron kinematics.
		event.calc_NeutronHCaldxdy(); // Calculates the "detected-predicted" hit positions on HCal.

		// Applying Cuts... //
		if ( !event.passBigBiteCuts() ) continue;	

		event.eval_HCalNDE_Rmthd();	
		
		output.copyFromEvent();
		output.fillOutTree();
		output.fillHistos();
	}

	event.print_HCalNDE_Rmthd();

	output.closeOutFile();		
	output.make_anapdf();	
	output.make_cutpdf();	
}

void print_analysis_percentage(double cuurentevent_ana_percentage, int& previousevent_ana_percentage_int)
{
	int cuurentevent_ana_percentage_int{(int)cuurentevent_ana_percentage};
	if (cuurentevent_ana_percentage_int%10==0 && cuurentevent_ana_percentage_int>previousevent_ana_percentage_int)
	{
		std::cout << cuurentevent_ana_percentage_int <<"%\n";
	}

	previousevent_ana_percentage_int = cuurentevent_ana_percentage_int;
}