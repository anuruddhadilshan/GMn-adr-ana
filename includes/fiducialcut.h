#ifndef FIDUCIALCUT_H
#define FIDUCIALCUT_H

#include "HCalConstants.h"

namespace FiducialCutConsts
{
	const std::vector<std::vector<double>> protonDeflection_simdat = { //{kinematicsetting(sbs #), target(0=LD2, 1=LH2), sbs fieldscale, avg proton deflection, std.dev x, std.dev y.
		{4, 0, 0, 0, 0, 0},
		{4, 1, 0, 0, 0, 0},
		{4, 1, 30.0, 0, 0, 0},
		{4, 0, 30.0, 0, 0, 0},
		{4, 0, 50.0, 0, 0, 0},
		//{7, 0, 85.0, 0, 0, 0},
		{7, 0, 85.0, -0.694537, 0.09877, 0.127543},
		{8, 0, 0, 0, 0, 0},
		{8, 1, 0, 0, 0, 0},
		{8, 0, 50.0, 0, 0, 0},
		{8, 0, 70.0, 0, 0, 0},
		{8, 0, 100.0, 0, 0, 0},
		{11, 0, 0, 0, 0, 0},
		{11, 0, 100.0, 0, 0, 0},
		{14, 0, 0, 0, 0, 0},
		{14, 0, 70.0, 0, 0, 0}
	};

	//// Average proton deflections for different replay passes. ////
	const std::vector<std::vector<double>> protonDeflection_pass0and1 = { //{kinematicsetting(sbs #), target(0=LD2, 1=LH2), sbs fieldscale, avg proton deflection, std.dev x, std.dev y.
		{4, 0, 0, 0, 0, 0},
		{4, 1, 0, 0, 0, 0},
		{4, 1, 30.0, 0, 0, 0},
		{4, 0, 30.0, -0.66863820, 0, 0},
		{4, 0, 50.0, -1.1193552, 0, 0},
		{7, 0, 85.0, -0.690, 0.101592, 0.168226},
		{8, 0, 0, 0, 0, 0},
		{8, 1, 0, 0, 0, 0},
		{8, 0, 50.0, -0.661306, 0, 0},
		{8, 0, 70.0, -0.920279, 0, 0},
		{8, 0, 100.0, -1.29957, 0, 0},
		{9, 0, 70.0, -0.811794, 0.180786, 0.198251},
		{11, 0, 0, 0, 0, 0},
		{11, 0, 100.0, -0.7144345, 0, 0},
		{14, 0, 0, 0, 0, 0},
		{14, 0, 70.0, -0.8045086, 0, 0}
	};

	const std::vector<std::vector<double>> protonDeflection_pass2 = { //{kinematicsetting(sbs #), target(0=LD2, 1=LH2), sbs fieldscale, avg proton deflection, std.dev x, std.dev y.
		{4, 0, 0, 0, 0, 0},
		{4, 1, 0, 0, 0, 0},
		{4, 1, 30.0, 0, 0, 0},
		{4, 0, 30.0, -0.758239, 0.199934, 0.0},
		//{4, 0, 30.0, 0.0, 0.0, 0.0},
		{4, 0, 50.0, -1.2097, 0.214125, 0.0},
		{7, 0, 85.0, -0.798306, 0.144905, 0.0},
		{8, 0, 0, 0, 0, 0},
		{8, 1, 0, 0, 0, 0},
		{8, 0, 50.0, -0.661306, 0, 0},
		{8, 0, 70.0, -0.920279, 0, 0},
		{8, 0, 100.0, -1.29957, 0, 0},
		{9, 0, 70.0, -0.811794, 0.180786, 0.198251},
		{11, 0, 0, 0, 0, 0},
		{11, 0, 100.0, -0.7144345, 0, 0},
		{14, 0, 0, 0, 0, 0},
		{14, 0, 70.0, -0.8045086, 0, 0}
	};
	
}


class SBSprotonDeflection // Class to get the average and std.dev of prton deflection distribution by the SBS magnet.
{

private:

	const int m_replaypassnum;
	const int m_sbskinenum;
	const TString m_target;
	const double m_sbsfieldscale;

	double m_targetnum {0}; // 0 = LD2, 1 = LH2. By default we will assume it is LD2 here.
	// double m_LD2_mean_protondeflection {0.};
	// double m_LD2_stddev_protondeflection {0.};
	// double m_LH2_mean_protondeflection {0.};
	// double m_LH2_stddev_protondeflection {0.};

public:

	SBSprotonDeflection( const int replayPassNum, const int sbsKineNum, const TString target, const double sbsFieldScale ) : m_replaypassnum{replayPassNum}, m_sbskinenum{sbsKineNum}, m_target{target}, m_sbsfieldscale{sbsFieldScale}
	{
		if ( target == "LH2" ) m_targetnum = 1;
	}

	SBSprotonDeflection( const int sbsKineNum, const TString target, const double sbsFieldScale ) : m_replaypassnum{-1}, m_sbskinenum{sbsKineNum}, m_target{target}, m_sbsfieldscale{sbsFieldScale} //Constructor to analyze simulated data.
	{
		if ( target == "LH2" ) m_targetnum = 1;
	}

	double return_targetnum()
	{
		return m_targetnum;
	}

	double get_proton_deflection()
	{
		if ( m_replaypassnum == -1 ) // For simulated data
		{
			for ( const auto& row : FiducialCutConsts::protonDeflection_simdat )
			{
				if ( row[0] == m_sbskinenum && row[1] == m_targetnum && row[2] == m_sbsfieldscale )
				{
					std::cout << "\nAverage proton deflection by the SBS magnet (m): " << row[3] << '\n';
					std::cout << "### IMPORTANT: The above deflection will be used to estimate the proton deflection for the fiducial cut ###\n";	
					return row[3];
				}
	     	}
		}
		else if ( m_replaypassnum == 0 || m_replaypassnum == 1 )
		{	
			for ( const auto& row : FiducialCutConsts::protonDeflection_pass0and1 )
			{
				if ( row[0] == m_sbskinenum && row[1] == m_targetnum && row[2] == m_sbsfieldscale )
				{
					std::cout << "\nAverage proton deflection by the SBS magnet (m): " << row[3] << '\n';
					std::cout << "### IMPORTANT: The above deflection will be used to estimate the proton deflection for the fiducial cut ###\n";	
					return row[3];
				}
	     	}
    	 }
    	else if ( m_replaypassnum == 2 )
		{	
			for ( const auto& row : FiducialCutConsts::protonDeflection_pass2 )
			{
				if ( row[0] == m_sbskinenum && row[1] == m_targetnum && row[2] == m_sbsfieldscale )
				{
					std::cout << "\nAverage proton deflection by the SBS magnet (m): " << row[3] << '\n';
					std::cout << "### IMPORTANT: The above deflection will be used to estimate the proton deflection for the fiducial cut ###\n";	
					return row[3];
				}
	     	}
    	 }

     	std::cout << "\n### ERROR: Either the mass replay pass number, kinematic setting number, or the percentage SBS magnet field scale entered are invalid ###\n";
		return 0;
	}

	double get_proton_deflection_stddev_x()
	{
		
		if ( m_replaypassnum == -1 )
		{
			for ( const auto& row : FiducialCutConsts::protonDeflection_simdat )
			{
				if ( row[0] == m_sbskinenum && row[1] == m_targetnum && row[2] == m_sbsfieldscale )
				{
					std::cout << "Standard deviation of the proton deflection distribution in x(dispersive) direction (m): " << row[4] << '\n';
					return row[4];
				}
	     	}
      	}
		else if ( m_replaypassnum == 0 || m_replaypassnum == 1 )
		{
			for ( const auto& row : FiducialCutConsts::protonDeflection_pass0and1 )
			{
				if ( row[0] == m_sbskinenum && row[1] == m_targetnum && row[2] == m_sbsfieldscale )
				{
					std::cout << "Standard deviation of the proton deflection distribution in x(dispersive) direction (m): " << row[4] << '\n';
					return row[4];
				}
	     	}
      	}
      	else if ( m_replaypassnum == 2 )
		{
			for ( const auto& row : FiducialCutConsts::protonDeflection_pass2 )
			{
				if ( row[0] == m_sbskinenum && row[1] == m_targetnum && row[2] == m_sbsfieldscale )
				{
					std::cout << "Standard deviation of the proton deflection distribution in x(dispersive) direction (m): " << row[4] << '\n';
					return row[4];
				}
	     	}
      	}

     	std::cout << "\n### ERROR: Either the mass replay pass number, kinematic setting number, or the percentage SBS magnet field scale entered are invalid ###\n";
		return 0;
	}

	double get_proton_deflection_stddev_y()
	{
		if ( m_replaypassnum == -1 )
		{
			for ( const auto& row : FiducialCutConsts::protonDeflection_simdat )
			{ 
				if ( row[0] == m_sbskinenum && row[1] == m_targetnum && row[2] == m_sbsfieldscale )
				{
					std::cout << "Standard deviation of the proton deflection distribution in y(non-dispersive) direction (m): " << row[5] << '\n';
					return row[5];
				}
	     	}
      	}
		else if ( m_replaypassnum == 0 || m_replaypassnum == 1 )
		{
			for ( const auto& row : FiducialCutConsts::protonDeflection_pass0and1 )
			{ 
				if ( row[0] == m_sbskinenum && row[1] == m_targetnum && row[2] == m_sbsfieldscale )
				{
					std::cout << "Standard deviation of the proton deflection distribution in y(non-dispersive) direction (m): " << row[5] << '\n';
					return row[5];
				}
	     	}
      	}
      	else if ( m_replaypassnum == 2 )
		{
			for ( const auto& row : FiducialCutConsts::protonDeflection_pass2 )
			{ 
				if ( row[0] == m_sbskinenum && row[1] == m_targetnum && row[2] == m_sbsfieldscale )
				{
					std::cout << "Standard deviation of the proton deflection distribution in y(non-dispersive) direction (m): " << row[5] << '\n';
					return row[5];
				}
	     	}
      	}

     	std::cout << "\n### ERROR: Either the mass replay pass number, kinematic setting number, or the percentage SBS magnet field scale entered are invalid ###\n";
		return 0;
	}

};


class FiducialCut
{

private: 

	SBSprotonDeflection m_protonDeflection;
	const double m_avg_proton_deflection;
	const double m_stddev_x;
	const double m_stddev_y;

	// Constant factor to multiply std.dev for safetly margins.
	const double m_stddev_factor_x;
	const double m_stddev_factor_y;

	// HCal boundaries with "double safety margins": Excluded one HCal block and a constat factor of std.dev of the dx and dy distributions. 
	double m_hcal_xlow_ssafe {0.};
	double m_hcal_xhigh_ssafe {0.};
	double m_hcal_ylow_ssafe {0.};
	double m_hcal_yhigh_ssafe {0.};
	
	const int m_targetintnum;

public:

	FiducialCut( const int replayPassNum, const int sbsKineNum, const TString target, const double sbsFieldScale, const double fidcutSafeFactorX = 0.0, const double fidcutSafeFactorY = 0.0 ) 
	: m_protonDeflection{replayPassNum,sbsKineNum,target,sbsFieldScale}, m_avg_proton_deflection{m_protonDeflection.get_proton_deflection()},
	m_stddev_x{m_protonDeflection.get_proton_deflection_stddev_x()}, m_stddev_y{m_protonDeflection.get_proton_deflection_stddev_y()}, 
	m_stddev_factor_x{fidcutSafeFactorX}, m_stddev_factor_y{fidcutSafeFactorY}, m_targetintnum{static_cast<int>(m_protonDeflection.return_targetnum())}
	{
		if ( replayPassNum == 0 || replayPassNum == 1 ) // 
		{
			m_hcal_xlow_ssafe = HCalConst::hcal_active_xlow_safe_pass1 + m_stddev_factor_x*m_stddev_x;
			m_hcal_xhigh_ssafe = HCalConst::hcal_active_xhigh_safe_pass1 - m_stddev_factor_x*m_stddev_x;
			m_hcal_ylow_ssafe = HCalConst::hcal_active_ylow_safe_pass1 + m_stddev_factor_y*m_stddev_y;
			m_hcal_yhigh_ssafe = HCalConst::hcal_active_yhigh_safe_pass1 - m_stddev_factor_y*m_stddev_y;
		}
		else if ( replayPassNum == 2 || replayPassNum == -1 ) // 
		{
			m_hcal_xlow_ssafe = HCalConst::hcal_active_xlow_safe_pass2 + m_stddev_factor_x*m_stddev_x;
			m_hcal_xhigh_ssafe = HCalConst::hcal_active_xhigh_safe_pass2 - m_stddev_factor_x*m_stddev_x;
			m_hcal_ylow_ssafe = HCalConst::hcal_active_ylow_safe_pass2 + m_stddev_factor_y*m_stddev_y;
			m_hcal_yhigh_ssafe = HCalConst::hcal_active_yhigh_safe_pass2 - m_stddev_factor_y*m_stddev_y;
		}
	} 

	int return_targetintnum()
	{
		return m_targetintnum;
	}

	bool pass_HitHCalCut( const double x_neutronhypths, const double y_neutronhypths ) // For LH2 real/simulated data analysis only. See whether the protons are hitting the active area of HCal.
	{
		double x_protonPosOnHCal = x_neutronhypths + m_avg_proton_deflection;
		double y_protonPosOnHCal = y_neutronhypths;

		if ( x_protonPosOnHCal < m_hcal_xlow_ssafe || x_protonPosOnHCal > m_hcal_xhigh_ssafe || y_protonPosOnHCal < m_hcal_ylow_ssafe || y_protonPosOnHCal > m_hcal_yhigh_ssafe ) 
		{
			return false;
		}	

		return true;

	}

	bool pass_FiducialCut( const double x_neutronhypths, const double y_neutronhypths ) // For LD2 real/simulated data. Matches the neutron and proton acceptances.
	{
		// First we see whether the neutron falls within the HCal. If not we return a "false" value as we do not need to proceed if the neutron already is outside the region we consider.
		if ( x_neutronhypths < m_hcal_xlow_ssafe || x_neutronhypths > m_hcal_xhigh_ssafe || y_neutronhypths < m_hcal_ylow_ssafe || y_neutronhypths > m_hcal_yhigh_ssafe )
		{
			return false;
		}

		// proton coordinates: Just add the expected average deflection expected to the neutron's x pos value.
		double x_protonPosOnHCal = x_neutronhypths + m_avg_proton_deflection;
		double y_protonPosOnHCal = y_neutronhypths;

		if ( x_protonPosOnHCal < m_hcal_xlow_ssafe || x_protonPosOnHCal > m_hcal_xhigh_ssafe || y_protonPosOnHCal < m_hcal_ylow_ssafe || y_protonPosOnHCal > m_hcal_yhigh_ssafe )
		{
			return false;
		}

		// We are here only if both the neutron and proton are within the HCal region considered. So we return a "true" value only at this stage, if it comes to this point.
		return true;
	}

};

#endif