#include <iostream>
#include <string>

#ifndef READ_DWNBNDANA_CONFIG_H
#define READ_DWNBNDANA_CONFIG_H

class Configfile 
{
private:
	//Defining the variables read in by the configuration file.
	TChain* m_C = new TChain("T");
	int m_SBSKine{9}; // Default to SBS 9.
	double m_MinGEMHitsOnTrack{0.}; 
	double m_VzCutUpStream{0.};
	double m_VzCutDwnStream{0.};
	double m_CoinCutMean {0.};
	double m_CoinCutSigma {0.};
	double m_CoinCutSigmaFactor {0.};
	double m_NeutronConeSigmaX {0.};
	double m_NeutronConeSigmaY {0.};
	double m_NeutronConeSigmaFactor {0.};
	double m_RCut{0.};

public:

	Configfile() = default;	

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

		    		if ( skey == "SBSKine" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_SBSKine = sval.Atoi();
		    			std::cout << "SBS kinematic setting number: " << m_SBSKine << '\n';
		    		}

		    		if ( skey == "MinGEMHitsOnTrack" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_MinGEMHitsOnTrack = sval.Atof();
		    			std::cout << "Min. # of GEM hits on track: " << m_MinGEMHitsOnTrack << '\n';
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

					if( skey == "NeutronConeSigmaX" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_NeutronConeSigmaX = sval.Atof();
		    			std::cout << "Neutron cone cut sigma x: " << m_NeutronConeSigmaX << '\n';
		    		}

		    		if( skey == "NeutronConeSigmaY" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_NeutronConeSigmaY = sval.Atof();
		    			std::cout << "Neutron cone cut sigma y: " << m_NeutronConeSigmaY << '\n';
		    		}

		    		if( skey == "NeutronConeSigmaFactor" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_NeutronConeSigmaFactor = sval.Atof();
		    			std::cout << "Neutron cone cut sigma factor to be used: " << m_NeutronConeSigmaFactor << '\n';
		    		}

		    		if( skey == "RCut" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_RCut = sval.Atof();
		    			std::cout << "Radius of the HCal search circle for hits (m): " << m_RCut << '\n';
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

	int return_SBSKineNum()
	{
		return m_SBSKine;
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

	double return_NeutronConeCutSigmaX()
	{
		return m_NeutronConeSigmaX;
	}

	double return_NeutronConeCutSigmaY()
	{
		return m_NeutronConeSigmaY;
	}

	double return_NeutronConeCutSigmaFactor()
	{
		return m_NeutronConeSigmaFactor;
	}

	double return_RCut()
	{
		return m_RCut;
	}

};

#endif