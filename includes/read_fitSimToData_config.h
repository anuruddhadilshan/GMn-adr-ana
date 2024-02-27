#include <iostream>
#include <string>

#ifndef READ_FITSIMTODATA_CONFIG_H
#define READ_FITSIMTODATA_CONFIG_H

class ConfigfileFitSimToData
{
private:
	//Defining the variables read in by the configuration file.
	// TChain* m_C = new TChain("T");
	// int m_SBSKineNum {0}; 
	// TString m_Target {""};
	// double m_MCBeamCurrMicroAmp {0.};
	// int m_MCNeventsGen {0};
	// double m_SBSFieldScale {0.};
	// double m_MinGEMHitsOnTrack {0.}; 
	// double m_MaxTrackChi2divNdof {0.};
	// double m_MinBBTrP {0.};
	// double m_VzCutUpStream {0.};
	// double m_VzCutDwnStream {0.};
	// double m_MinPreShE {0.};
	// double m_CoinCutMean {0.};
	// double m_CoinCutSigma {0.};
	// double m_CoinCutSigmaFactor {0.};	
	TString m_ExperimentDatROOTFile {""};
	TString m_SIMC_deen_ROOTfile {""};
	TString m_SIMC_deep_ROOTfile {""};
	TString m_G4SBS_inel_ROOTfile {""};
	TString m_OutputDir {""};
	double m_W2Cut_min {0.};
	double m_W2Cut_max {0.};
	double m_dyCut_min {0.};
	double m_dyCut_max {0.};	
	double m_simc_deen_dx_offset {0.};
	double m_simc_deep_dx_offset {0.};
	
public:

	ConfigfileFitSimToData() = default;	

	ConfigfileFitSimToData( const char* configfilename )
	{
		readin_configfile(configfilename);
	}

	// Function that reads in the configuration file and copy the information into the member varibales.
	void readin_configfile( const char* configfilename )
	{
		ifstream configfile(configfilename);
	  	TString currentline;

	  	std::cout <<'\n'<<"--- Reading configuration file: " << configfilename << " --- \n";

	  	//Loop to read-in cut thresholds.
	  	while( currentline.ReadLine( configfile ) )
	  	{
	  		if ( !currentline.BeginsWith("#") )
	  		{
	  			TObjArray *tokens = currentline.Tokenize(" "); 
		    	int ntokens = tokens->GetEntries();
		    	if( ntokens > 1 )
		    	{
		    		TString skey = ( (TObjString*)(*tokens)[0] )->GetString();

		    		if ( skey == "ExperimentDat_ROOTfile" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_ExperimentDatROOTFile = sval;
		    			std::cout << "Experiment data ROOT file: " << m_ExperimentDatROOTFile << '\n';
		    		}

		    		if ( skey == "SIMC_deen_ROOTfile" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_SIMC_deen_ROOTfile = sval;
		    			std::cout << "SIMC deen ROOT file: " << m_SIMC_deen_ROOTfile << '\n';
		    		}

		    		if ( skey == "SIMC_deep_ROOTfile" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_SIMC_deep_ROOTfile = sval;
		    			std::cout << "SIMC deep ROOT file: " << m_SIMC_deep_ROOTfile << '\n';
		    		}

		    		if ( skey == "G4SBS_inel_ROOTfile" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_G4SBS_inel_ROOTfile = sval;
		    			std::cout << "G4SBS inelastics ROOT file: " << m_G4SBS_inel_ROOTfile << '\n';
		    		}

		    		if ( skey == "W2_cut_min" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_W2Cut_min = sval.Atof();
		    			std::cout << "W2 cut minimum: " << m_W2Cut_min << " GeV^2\n";
		    		}

		    		if ( skey == "W2_cut_max" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_W2Cut_max = sval.Atof();
		    			std::cout << "W2 cut maximum: " << m_W2Cut_max << " GeV^2\n";
		    		}

		    		if ( skey == "dy_cut_min" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_dyCut_min = sval.Atof();
		    			std::cout << "dy cut minimum: " << m_dyCut_min << " m\n";
		    		}

		    		if ( skey == "dy_cut_max" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_dyCut_max = sval.Atof();
		    			std::cout << "dy cut maximum: " << m_dyCut_max << " m\n";
		    		}

		    		if ( skey == "sim_deen_dx_offset" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_simc_deen_dx_offset = sval.Atof();
		    			std::cout << "dx offset applied for SIMC deen events: " << m_simc_deen_dx_offset << " m\n";
		    		}

		    		if ( skey == "sim_deep_dx_offset" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_simc_deep_dx_offset = sval.Atof();
		    			std::cout << "dx offset applied for SIMC deep events: " << m_simc_deep_dx_offset << " m\n";
		    		}


		   //  		if ( skey == "SBSKineNum" )
		   //  		{
		   //  			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		   //  			m_SBSKineNum = sval.Atoi();
		   //  			std::cout << "SBS kinematic setting number: " << m_SBSKineNum << '\n';
		   //  		}

					// if ( skey == "Target" )
		   //  		{
		   //  			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		   //  			m_Target = sval;
		   //  			std::cout << "Target: " << m_Target << '\n';
		   //  		}

		   //  		if ( skey == "MCBeamCurrMicroAmp" )
		   //  		{
		   //  			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		   //  			m_MCBeamCurrMicroAmp = sval.Atof();
		   //  			std::cout << "Beam current used for the G4SBS simulation in microampheres: " << m_MCBeamCurrMicroAmp << '\n';
		   //  		}

		   //  		if ( skey == "MCNeventsGen" )
		   //  		{
		   //  			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		   //  			m_MCNeventsGen = sval.Atoi();
		   //  			std::cout << "Number of events generated in the G4SBS simulation: " << m_MCNeventsGen << '\n';
		   //  		}


		   //  		if ( skey == "SBSFieldScale" )
		   //  		{
		   //  			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		   //  			m_SBSFieldScale = sval.Atoi();
		   //  			std::cout << "SBS field scale: " << m_SBSFieldScale << '\n';
		   //  		}

		   //  		if ( skey == "MinGEMHitsOnTrack" )
		   //  		{
		   //  			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		   //  			m_MinGEMHitsOnTrack = sval.Atof();
		   //  			std::cout << "Min. # of GEM hits on track: " << m_MinGEMHitsOnTrack << '\n';
		   //  		}

		   //  		if( skey == "MinBBTrackP" )
		   //  		{
		   //  			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		   //  			m_MinBBTrP = sval.Atof();
		   //  			std::cout << "Minimum BigBite track momentum cut (GeV/c): " << m_MinBBTrP << '\n';
		   //  		}

		   //  		if( skey == "VzCutUpStream" )
		   //  		{
		   //  			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		   //  			m_VzCutUpStream = sval.Atof();
		   //  			std::cout << "BigBite track Vertex Z position upstream cut (m): " << m_VzCutUpStream << '\n';
		   //  		}

		   //  		if( skey == "VzCutDwnStream" )
		   //  		{
		   //  			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		   //  			m_VzCutDwnStream = sval.Atof();
		   //  			std::cout << "BigBite track Vertex Z position downstream cut (m): " << m_VzCutDwnStream << '\n';
		   //  		}	

		   //  		if( skey == "MinPreShE" )
		   //  		{
		   //  			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		   //  			m_MinPreShE = sval.Atof();
		   //  			std::cout << "Minimum BigBite Pre-Shower energy cut (GeV): " << m_MinPreShE << '\n';
		   //  		}

		   //  		if( skey == "MaxTrackChi2divNdof" )
		   //  		{
		   //  			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		   //  			m_MaxTrackChi2divNdof = sval.Atof();
		   //  			std::cout << "Maximum BigBite Track chi2/ndof: " << m_MaxTrackChi2divNdof << '\n';
		   //  		}

		   //  		if( skey == "CoinCutMean" )
		   //  		{
		   //  			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		   //  			m_CoinCutMean = sval.Atof();
		   //  			std::cout << "Mean of the gaus fit to the peak, sbs.hcal.heclus.atime[0]-bb.sh.atime (ns): " << m_CoinCutMean << '\n';
		   //  		}

		   //  		if( skey == "CoinCutSigma" )
		   //  		{
		   //  			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		   //  			m_CoinCutSigma = sval.Atof();
		   //  			std::cout << "Sigma of the gaus fit to the peak, sbs.hcal.heclus.atime[0]-bb.sh.atime (ns): " << m_CoinCutSigma << '\n';
		   //  		}

		   //  		if( skey == "CoinCutSigmaFactor" )
		   //  		{
		   //  			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		   //  			m_CoinCutSigmaFactor = sval.Atof();
		   //  			std::cout << "Sigma factor to be used in the HCal and SH, ADC conicidence time cut: " << m_CoinCutSigmaFactor << '\n';
		   //  		}

		    		if ( skey == "OutputDir" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_OutputDir = sval;
		    			std::cout << "Output directory path: " << m_OutputDir << '\n';
		    		}

		    	}

	  		}
	  	}

	  	std::cout <<"--- Finished reading configuration file --- \n";
	}

	// Return functions.

	TString return_ExprDatFile()
	{
		return m_ExperimentDatROOTFile;
	}

	TString return_SIMCdeenFile()
	{
		return m_SIMC_deen_ROOTfile;
	}

	TString return_SIMCdeepFile()
	{
		return m_SIMC_deep_ROOTfile;
	}

	TString return_G4SBSinelFile()
	{
		return m_G4SBS_inel_ROOTfile;
	}

	double return_W2CutMin()
	{
		return m_W2Cut_min;
	}

	double return_W2CutMax()
	{
		return m_W2Cut_max;
	}

	double return_dyCutMin()
	{
		return m_dyCut_min;
	}

	double return_dyCutMax()
	{
		return m_dyCut_max;
	}

	double return_simDeenDxOffset()
	{
		return m_simc_deen_dx_offset;
	}

	double return_simDeepDxOffset()
	{
		return m_simc_deep_dx_offset;
	}

	// TChain* return_TChain()
	// {
	// 	return m_C;
	// }

	// int return_SBSKineNum()
	// {
	// 	return m_SBSKineNum;
	// }

	// TString return_Target()
	// {
	// 	return m_Target;
	// } 

	// double return_MCBeamCurrMicroAmp()
	// {
	// 	return m_MCBeamCurrMicroAmp;
	// }

	// int return_MCNeventsGen()
	// {
	// 	return m_MCNeventsGen;
	// }

	// double return_SBSFieldScale()
	// {
	// 	return m_SBSFieldScale;
	// }

	// double return_MinGEMHitsOnTrack()
	// {
	// 	return m_MinGEMHitsOnTrack;
	// }

	// double return_VzCutDwnStream()
	// {
	// 	return m_VzCutDwnStream;
	// }

	// double return_VzCutUpStream()
	// {
	// 	return m_VzCutUpStream;
	// }

	// double return_MinBBTrPCut()
	// {
	// 	return m_MinBBTrP;
	// }	

	// double return_MinPreShECut()
	// {
	// 	return m_MinPreShE;
	// }

	// double return_MaxTrackChi2divNdof()
	// {
	// 	return m_MaxTrackChi2divNdof;
	// }

	// double return_CoinCutMean()
	// {
	// 	return m_CoinCutMean;
	// }

	// double return_CoinCutSigma()
	// {
	// 	return m_CoinCutSigma;
	// }

	// double return_CoinCutSigmaFactor()
	// {
	// 	return m_CoinCutSigmaFactor;
	// }

	TString return_OutputDir()
	{
		return m_OutputDir;
	}

};

#endif