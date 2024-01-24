#include "elasticevent.h"

#ifndef BESTHCALCLUS_H
#define BESTHCALCLUS_H


class BestHCalClus
{
private:

	ElasticEvent& m_event; // Class for elastic event selection.
	
	const double m_coincut_mean;
	const double m_coincut_sigma;
	const double m_coincut_sigmafactor;
	const double m_coincut_low;
	const double m_coincut_high;

	double m_shadctime {0.};

	const static int m_MAXHCALCLUS {10};		
	double* m_sbshcal_heclus_e; 
	double* m_sbshcal_heclus_atime; 
	double* m_sbshcal_heclus_x; 
	double* m_sbshcal_heclus_y; 

	int m_besthcalclus_indx {0};
	bool m_passCoinCut {true}; // Will only be false in the case where even th 5th highest energy cluster fails the timing cut.


	void swap(double& a, double&b) // Function to swap two double values.
	{
		double temp = a;
		a = b;
		b = temp;
	}

	int findMaxByMagnitude(int start) // Function to find the index of the maximum element by magnitude
	{
		int maxIndex = start;

		for ( int i = start + 1; i < m_MAXHCALCLUS; i++ )
		{
			if ( m_sbshcal_heclus_e[i] > m_sbshcal_heclus_e[maxIndex] )
			{
				maxIndex = i;
			}
		}

		return maxIndex;
	}

	void makeHighestEClusArrys() // Create arrays  m_sbshcal_heclus_e, m_sbshcal_heclus_atime, m_sbshcal_heclus_x, and m_sbshcal_heclus_y, which are sorted by the descending order of cluster energy.
	{
		for ( int i = 0; i < m_MAXHCALCLUS -1; i++ )
		{ 
			int maxIndex = findMaxByMagnitude(i); // Gets the largest element's indx of the array starting from the i'th to the last elemnt of the array, sbs.hcal.clus.e

			if ( maxIndex != i ) 
			{
				swap( m_sbshcal_heclus_e[i], m_sbshcal_heclus_e[maxIndex] );
				swap( m_sbshcal_heclus_atime[i], m_sbshcal_heclus_atime[maxIndex] );
				swap( m_sbshcal_heclus_x[i], m_sbshcal_heclus_x[maxIndex] );
				swap( m_sbshcal_heclus_y[i], m_sbshcal_heclus_y[maxIndex] );
			}			
		}
	}

	void getSHtime() // Function to get the shower ADC time.
	{
		m_shadctime = m_event.return_SHADCTime();		
	}

public:

	BestHCalClus( ElasticEvent& event, const double coincut_mean, const double coincut_sigma, const double coincut_sigmafactor ) : m_event{event}, m_sbshcal_heclus_e{m_event.return_HCalHEClusEArryPtr()}, m_sbshcal_heclus_atime{m_event.return_HCalHEClusAtimeArryPtr()}, m_sbshcal_heclus_x{m_event.return_HCalHEClusXArryPtr()}, m_sbshcal_heclus_y{m_event.return_HCalHEClusYArryPtr()},
	m_coincut_mean{coincut_mean}, m_coincut_sigma{coincut_sigma}, m_coincut_sigmafactor{coincut_sigmafactor}, m_coincut_low{coincut_mean - coincut_sigmafactor*coincut_sigma}, m_coincut_high{coincut_mean + coincut_sigmafactor*coincut_sigma}	
	{
	}
	

	void findBestHCalClus() // Function to return Best HCal cluster's index, using the Atime+HighestEnergy method.
	{
		
		makeHighestEClusArrys();

		getSHtime();

		const double cointime_heclus0 = m_sbshcal_heclus_atime[0] - m_shadctime;
		const double cointime_heclus1 = m_sbshcal_heclus_atime[1] - m_shadctime;
		const double cointime_heclus2 = m_sbshcal_heclus_atime[2] - m_shadctime;
		const double cointime_heclus3 = m_sbshcal_heclus_atime[3] - m_shadctime;
		const double cointime_heclus4 = m_sbshcal_heclus_atime[4] - m_shadctime;

		m_passCoinCut = true;

		if ( m_coincut_low < cointime_heclus0 && cointime_heclus0 < m_coincut_high ) m_besthcalclus_indx = 0;
		else if ( m_coincut_low < cointime_heclus1 && cointime_heclus1 < m_coincut_high ) m_besthcalclus_indx = 1;
		else if ( m_coincut_low < cointime_heclus2 && cointime_heclus2 < m_coincut_high ) m_besthcalclus_indx = 2;
		else if ( m_coincut_low < cointime_heclus3 && cointime_heclus3 < m_coincut_high ) m_besthcalclus_indx = 3;
		else if ( m_coincut_low < cointime_heclus4 && cointime_heclus4 < m_coincut_high ) m_besthcalclus_indx = 4;
		else
		{
			m_besthcalclus_indx = 0; // Default to the 0th element. Not attemting past the fifth higest energy cluster.
			m_passCoinCut = false;
		}
	}

	int return_BestHCalClusIndx() // Index of the "best cluster" in the array sorted by descending order of cluster energy.
	{
		return m_besthcalclus_indx;
	}

	double return_BestHCalClusE()
	{
		return m_sbshcal_heclus_e[m_besthcalclus_indx];
	}  

	double return_BestHCalClusAtime()
	{
		return m_sbshcal_heclus_atime[m_besthcalclus_indx];
	}

	double return_BestHCalClusX()
	{
		return m_sbshcal_heclus_x[m_besthcalclus_indx];
	}

	double return_BestHCalClusY()
	{
		return m_sbshcal_heclus_y[m_besthcalclus_indx];
	}

	bool return_DidPassCoinCut()
	{
		return m_passCoinCut;
	}

};

#endif