// Script to select elastic events SIMC+G4SBS simulated data. 

#include <iostream>
#include <string>
#include <cctype>
#include "../includes/read_gmnelastics_simc_config.h"
#include "../includes/elasticeventg4sbs.h" // ElasticEventG4SBS class is used with functionality added to handle SIMC 
#include "../includes/besthcalclusforg4sbs.h"
#include "../includes/gmnelasg4sbsdatout.h"


void print_analysis_percentage(double cuurentevent_ana_percentage, int& previousevent_ana_percentage_int);
int getJobNum(TString rootfile);
double readJLabHPC_CVSfile(const char* csvfile_constchar, int jobnum, int colnum );
int getSIMCweightPar(const int jobnum, const char* csvfile_constchar, double& nTried, double& genvol, double& luminosity, double& charge);
void makeSingleROOTfile(TString outputdirpath, const char* outfilenamepreint);


void gmn_elastics_simc( const char* configfilename, TString outputfilenamepreint = "gmn_elastics_simcdat" )
{

	ConfigfileSIMC configfile{configfilename};

	const std::vector<TString>& rootfile_vector = configfile.return_ROOTFileVector();
	const TString csvfile = configfile.return_CSVFile();
	const char* csvfile_constchar = csvfile.Data();
	const TString simc_process = configfile.return_SIMCprocess();
	outputfilenamepreint.Append('_');
	const TString outputfilenamepreint_withprocess = outputfilenamepreint.Append(simc_process);
	const int num_SBSkine = configfile.return_SBSKineNum();
	const TString target = configfile.return_Target();
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


	const int num_rootfiles = rootfile_vector.size();
	std::cout << "\n#" << num_rootfiles << " ROOT files were loaded for analysis.\n";

	std::cout << "\n--- Beginning elastic-event-selection analysis of simulated data from the *SIMC* generator ---\n";
	
	int i_rootfile {0};

	for ( const auto& filename : rootfile_vector ) // Loop through the replayed ROOT files and analyze them.
	{
		
		std::cout << "\n#" << ++i_rootfile << " ROOT file, begin analysis... " << "(" << filename << ")\n";
		std::cout << '\n';

		
		//// Read the corresponding .csv file and extract the parameters, luminosity, genvol, and Ntried for SIMC 'event weight' calculations. ///
		double nTried {0.};
		double genvol {0.};
		double luminosity {0.};
		double charge {0.};
						
		int jobnum = getJobNum(filename);

		if ( getSIMCweightPar(jobnum, csvfile_constchar, nTried, genvol, luminosity, charge) != 0 ) continue; //Skip the ROOT file if return val is not 0. The reason will be given within the function body.
		//// Complete reading the .csv file ////


		TChain* const C = new TChain("T");
		C->Add(filename);

		// Initialize an object from the "ElasticEventG4SBS" class. Use the constructor intended for the SIMC analysis.
		ElasticEventG4SBS event{ C, num_SBSkine, target, nTried, genvol, luminosity, charge, num_SBSFieldScale, cut_MinPreShE, cut_MinGEMHitsOnTrack, cut_MinBBTrP, cut_MaxTrackChi2divNdof, cut_VzCutUpStream, cut_VzCutDwnStream };

		// Initialize an object from the "BestHCalClus" class to find the "best" HCal cluster.
		BestHCalClusForG4SBS bestHCalClus{ event, cut_CoinCutMean, cut_CoinCutSigma, cut_CoinCutSigmaFactor };

		// Initialize an object from the "OutputG4SBS" class for analysis outputs. The last bool input variable must be set to "true" to indicate analysis with SIMC generated simulation data.
		const char* currentrun_outfilenamepreint = Form("%s_job_%i", outputfilenamepreint_withprocess.Data(), jobnum);
		OutputG4SBS output{ event, bestHCalClus, path_OutputDir, currentrun_outfilenamepreint, true };

		//Variables to keep track of printout to the terminal.
		int previousevent_ana_percentage_int{0};
		long nevents {C->GetEntries()};


		std::cout << "\nNumber of events in the TChain: " << nevents << '\n';

		long nevent {0};

		while ( event.getEntry(nevent++) )
		{
			double ana_percentage{(nevent/(double)nevents)*100}; //Percentage of events analyzed in the Event List.
			print_analysis_percentage(ana_percentage, previousevent_ana_percentage_int);

			if ( !event.passGoodElectronCuts() ) continue; // Pre-shower E and all the track quality cuts applied using E arm data. 

			// *** Calculating the kinematics of the scattered electron and the hadron *** //
			event.calcBBTrackAngles();
			event.calcQ2andW2();
			event.calcNeutronHypthsHCalIntersect(); 
			event.checkFiducialCut();
			//if ( !event.passFiducialCut() ) continue;

			event.calcSIMCweight();

			bestHCalClus.findBestHCalClus();

			output.copyFromEvent();
			output.fillOutTree();
			output.fillHistos();

		}

		output.copySIMCrunInfo();
		output.fillSIMCrunInfoTree();

		output.closeOutFile();

	}

	makeSingleROOTfile( path_OutputDir, outputfilenamepreint_withprocess );

	std::cout << "\n--- Analysis Completed! ---\n";
	std::cout << '\n';	
	
}

////

void print_analysis_percentage(double cuurentevent_ana_percentage, int& previousevent_ana_percentage_int)
{
	int cuurentevent_ana_percentage_int{(int)cuurentevent_ana_percentage};
	if (cuurentevent_ana_percentage_int%10==0 && cuurentevent_ana_percentage_int>previousevent_ana_percentage_int)
	{
		std::cout << cuurentevent_ana_percentage_int <<"%\n";
	}

	previousevent_ana_percentage_int = cuurentevent_ana_percentage_int;

}

int getJobNum(const TString rootfile)
{
	// Find the position of the text sequence
    TString searchText = "_job_";
    int index = rootfile.Index(searchText);

    if (index != kNPOS) 
    {
        // Extract the substring after searchText
        TString numberString = rootfile(index + searchText.Length(), rootfile.Length());

        // Find the position of the first non-numeric character
        const char* numberStringconstchar = numberString.Data(); // Convert TString to const char*
        std::string stdNumberString(numberStringconstchar); // Convert const char* to std::string

        int endPos = 0;

        for (Int_t i = 0; i < stdNumberString.length(); ++i) 
        {
            // Check if the character is not a digit
            if (!isdigit(stdNumberString[i]))
            {
                endPos = i;
                break;
            }
        }

        numberString = numberString(0, endPos); // Extract only the numeric part

        int number = numberString.Atoi();

        return number;
    } 
    else 
    {
        return -1; 
    }
}

double readJLabHPC_CVSfile(const char* csvfile_constchar, const int jobnum, const int colnum )
{
	std::ifstream csvfile(csvfile_constchar);

	if (!csvfile.is_open()) 
	{
        std::cerr << "*Error opening file: " << csvfile_constchar << '\n';
        return -2;
    }

    TString currentline;

    while ( currentline.ReadLine( csvfile ) )
    {
    	TObjArray* tokens = currentline.Tokenize(",");

    	TString firstentry = ( (TObjString*)(*tokens)[0] )->GetString(); // Retrieve the first entry.

    	if ( !firstentry.IsDigit() ) 
    	{	
    		delete tokens;
    		continue; // Reject and skip if non-digit. Ex. the first line.
    	}

    	int firstentry_val = firstentry.Atoi();

    	if ( firstentry_val != jobnum ) 
    	{
    		delete tokens;
    		continue; // Reject and skip if the first entry is not the job number we want.
    	}

    	TString requiredentry = ( (TObjString*)(*tokens)[colnum] )->GetString(); // Retrieve the parameter string we want.

    	csvfile.close();

    	delete tokens; // Clean up the tokens array

    	double requiredentry_val = requiredentry.Atof();

    	if ( requiredentry_val != 0 || (requiredentry_val == 0.0 && requiredentry == "0") ) 
    	{
    		return requiredentry_val;
    	}
    	else
    	{
    		return -1; // Conversion to double was unsuccessful for some reason. 
    	}

    }

    csvfile.close();

    std::cerr << "*Error finding the line with the job number " << jobnum << " in the CSV file: " << csvfile_constchar << '\n';

    return -3;
}

int getSIMCweightPar(const int jobnum, const char* csvfile_constchar, double& nTried, double& genvol, double& luminosity, double& charge)
{
		if ( jobnum == -1 )
		{
			std::cout << "*Error: A job number could not be found for the current ROOT file from the ROOT file name. Skipping the ROOT file.\n";
			return -1; 
		}

		// Extract the above quantities from the .csv file for that particular job number. 
		// readJLabHPC_CVSfile(): Function to retrieve the necessary parameters from the JLab-HPC generated CSV file.
		nTried = readJLabHPC_CVSfile(csvfile_constchar, jobnum, 2);
		if ( nTried == -2 || nTried == -3 ) 
		{
			std::cout << "Skipping the ROOT file.\n";
			return -1;
		}

		if ( nTried == -1 )
		{
			std::cout << "*Error reading the CSV file: 'nTried' could not be read. Skipping the ROOT file.\n";
			return -1; 
		}
			
		genvol = readJLabHPC_CVSfile(csvfile_constchar, jobnum, 3);
		if ( genvol == -1 )
		{
			std::cout << "*Error reading the CSV file: 'genvol' could not be read. Skipping the ROOT file.\n";
			return -1; 
		}

		luminosity = readJLabHPC_CVSfile(csvfile_constchar, jobnum, 4);
		if ( luminosity == -1 )
		{
			std::cout << "*Error reading the CSV file: 'luminosity' could not be read. Skipping the ROOT file.\n";
			return -1; 
		}

		charge = readJLabHPC_CVSfile(csvfile_constchar, jobnum, 6);
		if ( charge == -1 )
		{
			std::cout << "*Error reading the CSV file: 'charge' could not be read. Skipping the ROOT file.\n";
			return -1; 
		}	

		return 0;
}

void makeSingleROOTfile( TString outputdirpath, const char* outfilenamepreint ) // Makes a single ROOT file from all jobs analyzed.
{
	TChain* CT = new TChain("T");
	TChain* CS = new TChain("S");

	const char* rootfilenamestring = Form("%s/%s_job_*.root", outputdirpath.Data(), outfilenamepreint);

	CT->Add(rootfilenamestring);
	CS->Add(rootfilenamestring);

	// Charge is the only parameter that will be useful to have in a final ROOT file where all the events from multiple jobs are included?
	// Thus, only that will be preserved.

	CS->SetBranchStatus("*", 0);
	CS->SetBranchStatus("MC.simc.charge", 1);

	double charge_simc {0.};

	CS->SetBranchAddress("MC.simc.charge", &charge_simc);

	int jobnum {0};
	double total_charge {0.};
	
	while ( CS->GetEntry(jobnum++) )
	{
		total_charge += charge_simc;
	}

	TFile* finaloutrootfile = new TFile(Form("%s/%s_all.root", outputdirpath.Data(), outfilenamepreint), "RECREATE");

	TTree* T;
	T = CT->CloneTree();

	TTree* S = new TTree("S", "SIMC Simulation Info Tree");

	S->Branch("MC.simc.charge", &total_charge);
	S->Fill();	

	S->Write();
	T->Write(0, TObject::kWriteDelete, 0);

	//finaloutrootfile->Write();
	finaloutrootfile->Close();
	delete finaloutrootfile;
}