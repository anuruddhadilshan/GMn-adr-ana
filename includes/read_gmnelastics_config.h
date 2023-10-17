#include <iostream>
#include <string>

#ifndef READ_GMNELASTICS_CONFIG_H
#define READ_GMNELASTICS_CONFIG_H

class Configfile 
{
private:
	//Defining the variables read in by the configuration file.
	TChain* m_C = new TChain("T");
	int m_ReplayPassNum {0};
	int m_SBSKineNum {0}; 
	TString m_Target {""};
	double m_SBSFieldScale {0.};
	double m_MinGEMHitsOnTrack {0.}; 
	double m_MaxTrackChi2divNdof {0.};
	double m_MinBBTrP {0.};
	double m_VzCutUpStream {0.};
	double m_VzCutDwnStream {0.};
	double m_MinPreShE {0.};
	double m_CoinCutMean {0.};
	double m_CoinCutSigma {0.};
	double m_CoinCutSigmaFactor {0.};	
	TString m_OutputDir {""};
	
	
public:

	Configfile() = default;	

	Configfile( const char* configfilename )
	{
		readin_dwnbndana_configfile(configfilename);
	}

	// Function that reads in the configuration file and copy the information into the member varibales.
	void readin_dwnbndana_configfile( const char* configfilename )
	{
		ifstream configfile(configfilename);
	  	TString currentline;

	  	std::cout <<'\n'<<"--- Reading configuration file: " << configfilename << " --- \n";

	  	//Loop to read-in input ROOT files.
	  	while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") )
	    {
	    	if( !currentline.BeginsWith("#") )
	    	{
	      		m_C->Add(currentline);
      			std::cout << "Loaded root file: " << currentline << '\n';
	    	}    
	  	}

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

		    		if ( skey == "ReplayPassNum" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_ReplayPassNum = sval.Atoi();
		    			std::cout << "Replay pass number: " << m_ReplayPassNum << '\n';
		    		}

		    		if ( skey == "SBSKineNum" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_SBSKineNum = sval.Atoi();
		    			std::cout << "SBS kinematic setting number: " << m_SBSKineNum << '\n';
		    		}

					if ( skey == "Target" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_Target = sval;
		    			std::cout << "Target: " << m_Target << '\n';
		    		}

		    		if ( skey == "SBSFieldScale" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_SBSFieldScale = sval.Atoi();
		    			std::cout << "SBS field scale: " << m_SBSFieldScale << '\n';
		    		}

		    		if ( skey == "MinGEMHitsOnTrack" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_MinGEMHitsOnTrack = sval.Atof();
		    			std::cout << "Min. # of GEM hits on track: " << m_MinGEMHitsOnTrack << '\n';
		    		}

		    		// if( skey == "HCalECut" )
		    		// {
		    		// 	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    		// 	m_HCalECut = sval.Atof();
		    		// 	std::cout << "HCal energy cut (GeV): " << m_HCalECut << '\n';
		    		// }

		    		// if( skey == "CoinCutLow" )
		    		// {
		    		// 	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    		// 	m_CoinCutLow = sval.Atof();
		    		// 	std::cout << "Low limit of 'HCal time - BBCal time' cut (ns): " << m_CoinCutLow << '\n';
		    		// }

		    		// if( skey == "CoinCutHigh" )
		    		// {
		    		// 	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    		// 	m_CoinCutHigh = sval.Atof();
		    		// 	std::cout << "High limit of 'HCal time - BBCal time' cut (ns): " << m_CoinCutHigh << '\n';
		    		// }

		    		if( skey == "MinBBTrackP" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_MinBBTrP = sval.Atof();
		    			std::cout << "Minimum BigBite track momentum cut (GeV/c): " << m_MinBBTrP << '\n';
		    		}

		    		if( skey == "VzCutUpStream" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_VzCutUpStream = sval.Atof();
		    			std::cout << "BigBite track Vertex Z position upstream cut (m): " << m_VzCutUpStream << '\n';
		    		}

		    		if( skey == "VzCutDwnStream" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_VzCutDwnStream = sval.Atof();
		    			std::cout << "BigBite track Vertex Z position downstream cut (m): " << m_VzCutDwnStream << '\n';
		    		}	

		    		if( skey == "MinPreShE" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_MinPreShE = sval.Atof();
		    			std::cout << "Minimum BigBite Pre-Shower energy cut (GeV): " << m_MinPreShE << '\n';
		    		}

		    		if( skey == "MaxTrackChi2divNdof" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_MaxTrackChi2divNdof = sval.Atof();
		    			std::cout << "Maximum BigBite Track chi2/ndof: " << m_MaxTrackChi2divNdof << '\n';
		    		}

		    		if( skey == "CoinCutMean" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_CoinCutMean = sval.Atof();
		    			std::cout << "Mean of the gaus fit to the peak, sbs.hcal.heclus.atime[0]-bb.sh.atime (ns): " << m_CoinCutMean << '\n';
		    		}

		    		if( skey == "CoinCutSigma" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_CoinCutSigma = sval.Atof();
		    			std::cout << "Sigma of the gaus fit to the peak, sbs.hcal.heclus.atime[0]-bb.sh.atime (ns): " << m_CoinCutSigma << '\n';
		    		}

		    		if( skey == "CoinCutSigmaFactor" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_CoinCutSigmaFactor = sval.Atof();
		    			std::cout << "Sigma factor to be used in the HCal and SH, ADC conicidence time cut: " << m_CoinCutSigmaFactor << '\n';
		    		}

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

	TChain* return_TChain()
	{
		return m_C;
	}

	int return_ReplayPassNum()
	{
		return m_ReplayPassNum;
	}

	int return_SBSKineNum()
	{
		return m_SBSKineNum;
	}

	TString return_Target()
	{
		return m_Target;
	} 

	double return_SBSFieldScale()
	{
		return m_SBSFieldScale;
	}

	double return_MinGEMHitsOnTrack()
	{
		return m_MinGEMHitsOnTrack;
	}

	double return_VzCutDwnStream()
	{
		return m_VzCutDwnStream;
	}

	double return_VzCutUpStream()
	{
		return m_VzCutUpStream;
	}

	// double return_HCalECut()
	// {
	// 	return m_HCalECut;
	// }

	// double return_CoinCutLow()
	// {
	// 	return m_CoinCutLow;
	// }

	// double return_CoinCutHigh()
	// {
	// 	return m_CoinCutHigh;
	// }

	double return_MinBBTrPCut()
	{
		return m_MinBBTrP;
	}	

	double return_MinPreShECut()
	{
		return m_MinPreShE;
	}

	double return_MaxTrackChi2divNdof()
	{
		return m_MaxTrackChi2divNdof;
	}

	double return_CoinCutMean()
	{
		return m_CoinCutMean;
	}

	double return_CoinCutSigma()
	{
		return m_CoinCutSigma;
	}

	double return_CoinCutSigmaFactor()
	{
		return m_CoinCutSigmaFactor;
	}

	TString return_OutputDir()
	{
		return m_OutputDir;
	}

};

#endif