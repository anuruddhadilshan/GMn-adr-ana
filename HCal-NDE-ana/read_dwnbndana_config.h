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
	double m_HCalECut{0.};
	double m_CoinCutLow{0.};
	double m_CoinCutHigh{0.};
	double m_HCalAdcTimeCutLow{0.};
	double m_HCalAdcTimeCutHigh{0.};
	// double m_BBTrPCut{0.};
	double m_BBTrVzCut{0.};
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

		    		if( skey == "HCalECut" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_HCalECut = sval.Atof();
		    			std::cout << "HCal energy cut (GeV): " << m_HCalECut << '\n';
		    		}

		    		if( skey == "CoinCutLow" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_CoinCutLow = sval.Atof();
		    			std::cout << "Low limit of 'HCal time - BBCal time' cut (ns): " << m_CoinCutLow << '\n';
		    		}

		    		if( skey == "CoinCutHigh" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_CoinCutHigh = sval.Atof();
		    			std::cout << "High limit of 'HCal time - BBCal time' cut (ns): " << m_CoinCutHigh << '\n';
		    		}

		    		// if( skey == "BBTrPCut" )
		    		// {
		    		// 	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    		// 	m_BBTrPCut = sval.Atof();
		    		// 	std::cout << "BigBite track momentum cut (GeV/c): " << m_BBTrPCut << '\n';
		    		// }

		    		if( skey == "BBTrVzCut" )
		    		{
		    			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		    			m_BBTrVzCut = sval.Atof();
		    			std::cout << "BigBite track Vertex Z position cut (m): " << m_BBTrVzCut << '\n';
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

	double return_BBTrVzCut()
	{
		return m_BBTrVzCut;
	}

	double return_HCalECut()
	{
		return m_HCalECut;
	}

	double return_CoinCutLow()
	{
		return m_CoinCutLow;
	}

	double return_CoinCutHigh()
	{
		return m_CoinCutHigh;
	}

	// double return_BBTrPCut()
	// {
	// 	return m_BBTrPCut;
	// }

	double return_RCut()
	{
		return m_RCut;
	}

};

#endif