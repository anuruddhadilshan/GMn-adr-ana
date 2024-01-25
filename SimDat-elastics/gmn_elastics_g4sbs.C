// This script is made to make elstic cuts and find elastic events from g4sbs simulation data.

#include <iostream>
#include <string>
#include "../includes/read_gmnelastics_g4sbs_config.h"
#include "../includes/elasticeventg4sbs.h"
#include "../includes/besthcalclusforg4sbs.h"
#include "../includes/gmnelasg4sbsdatout.h"


void print_analysis_percentage(double cuurentevent_ana_percentage, int& previousevent_ana_percentage_int);


void gmn_elastics_g4sbs( const char* configfilename, const char* outputfilename = "gmn_elastics_simdat" )
{

	ConfigfileG4SBS configfile{configfilename};

	TChain* const C = configfile.return_TChain();
	const int num_SBSkine = configfile.return_SBSKineNum();
	const TString target = configfile.return_Target();
	const double mc_beamcurrmicroamp = configfile.return_MCBeamCurrMicroAmp();
	const int mc_neventsgen = configfile.return_MCNeventsGen();
	const double num_SBSFieldScale = configfile.return_SBSFieldScale();
	// Copying cut thresholds to local variables //
	const double cut_MinPreShE = configfile.return_MinPreShECut();
	const double cut_MinGEMHitsOnTrack = configfile.return_MinGEMHitsOnTrack();
	const double cut_MinBBTrP = configfile.return_MinBBTrPCut();
	const double cut_MaxTrackChi2divNdof = configfile.return_MaxTrackChi2divNdof();
	const double cut_VzCutUpStream = configfile.return_VzCutUpStream();
	const double cut_VzCutDwnStream = configfile.return_VzCutDwnStream();
	const double cut_CoinCutMean = configfile.return_CoinCutMean(); 
	const double cut_CoinCutSigma = configfile.return_CoinCutSigma();	
	const double cut_CoinCutSigmaFactor = configfile.return_CoinCutSigmaFactor();
	const TString path_OutputDir = configfile.return_OutputDir();

	// Initialize an object from the "ElasticEvent" class.
	ElasticEventG4SBS event{ C, num_SBSkine, target, mc_beamcurrmicroamp, mc_neventsgen, num_SBSFieldScale, cut_MinPreShE, cut_MinGEMHitsOnTrack, cut_MinBBTrP, cut_MaxTrackChi2divNdof, cut_VzCutUpStream, cut_VzCutDwnStream };

	// Initialize an object from the "BestHCalClus" class to find the "best" HCal cluster.
	BestHCalClusForG4SBS bestHCalClus{ event, cut_CoinCutMean, cut_CoinCutSigma, cut_CoinCutSigmaFactor };

	// Initialize an object from the "Output" class for analysis outputs.
	OutputG4SBS output{ event, bestHCalClus, path_OutputDir, outputfilename };

	//Variables to keep track of printout to the terminal.
	int previousevent_ana_percentage_int{0};
	long nevents {C->GetEntries()};

	std::cout << "Number of events in the TChain: " << nevents << '\n';

	std::cout << "\n--- Beginning elastic-event-selection analysis of simulated data ---\n";

	long nevent {0};

	while ( event.getEntry(nevent++) )
	{

		double ana_percentage{(nevent/(double)nevents)*100}; //Percentage of events analyzed in the Event List.
		print_analysis_percentage(ana_percentage, previousevent_ana_percentage_int);

		if ( !event.passGoodElectronCuts() ) continue;

		event.calcBBTrackAngles();
		event.calcQ2andW2();
		event.calcNeutronHypthsHCalIntersect();
		event.checkFiducialCut();
		//if ( !event.passFiducialCut() ) continue;

		event.calcMCweight();

		bestHCalClus.findBestHCalClus();

		output.copyFromEvent();
		output.fillOutTree();
		output.fillHistos();
	}

	//output.makePlots();
	output.closeOutFile();

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