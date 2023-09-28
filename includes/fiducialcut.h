#ifndef FIDUCIALCUT_H
#define FIDUCIALCUT_H

#include "HCalConstants.h"

namespace FiducialCutConsts
{
	//// Average proton deflections for different replay passes. ////
	const std::vector<std::vector<double>> protonDeflection_pass0and1 = { //{kinematicsetting(sbs #), target(0=LD2, 1=LH2), sbs fieldscale, avg proton deflection, std.dev x, std.dev y.
		{4, 0, 0, 0, 0, 0},
		{4, 1, 0, 0, 0, 0},
		{4, 1, 30.0, 0, 0, 0},
		{4, 0, 30.0, -0.66863820, 0, 0},
		{4, 0, 50.0, -1.1193552, 0, 0},
		{7, 0, 85.0, -0.693155, 0, 0},
		{8, 0, 0, 0, 0, 0},
		{8, 1, 0, 0, 0, 0},
		{8, 0, 50.0, -0.661306, 0, 0},
		{8, 0, 70.0, -0.920279, 0, 0},
		{8, 0, 100.0, -1.29957, 0, 0},
		{11, 0, 0, 0, 0, 0},
		{11, 0, 100.0, -0.7144345, 0, 0},
		{14, 0, 0, 0, 0, 0},
		{14, 0, 70.0, -0.8045086, 0, 0}
	};

	//// HCal boundaries with safety margins in HCal coordinate system. ////
	// Safety margin = Exclude 1.5 rows/columns from the edges - Andrew suggestion.
	const double m_hcal_active_xlow_safe =  -0.75 - 10.5*HCalConst::hcalblk_h;
	const double m_hcal_active_xhigh_safe = -0.75 + 10.5*HCalConst::hcalblk_h;
	const double m_hcal_active_ylow_safe = -4.5*HCalConst::hcalblk_w;
	const double m_hcal_active_yhigh_safe = 4.5*HCalConst::hcalblk_w;
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

	double return_targetnum()
	{
		return m_targetnum;
	}

	double get_proton_deflection()
	{

		if ( m_replaypassnum == 0 || m_replaypassnum == 1 )
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

     	std::cout << "\n### ERROR: Either the mass replay pass number, kinematic setting number, or the percentage SBS magnet field scale entered are invalid ###\n";
		return 0;
	}

	double get_proton_deflection_stddev_x()
	{
		if ( m_replaypassnum == 0 || m_replaypassnum == 1 )
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

     	std::cout << "\n### ERROR: Either the mass replay pass number, kinematic setting number, or the percentage SBS magnet field scale entered are invalid ###\n";
		return 0;
	}

	double get_proton_deflection_stddev_y()
	{
		if ( m_replaypassnum == 0 || m_replaypassnum == 1 )
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
	const double m_stddev_factor;

	// HCal boundaries with "double safety margins": Excluded one HCal block and a constat factor of std.dev of the dx and dy distributions. 
	const double m_hcal_xlow_ssafe;
	const double m_hcal_xhigh_ssafe;
	const double m_hcal_ylow_ssafe;
	const double m_hcal_yhigh_ssafe;
	
	const int m_targetintnum;

public:

	FiducialCut( const int replayPassNum, const int sbsKineNum, const TString target, const double sbsFieldScale ) 
	: m_protonDeflection{replayPassNum,sbsKineNum,target,sbsFieldScale}, m_avg_proton_deflection{m_protonDeflection.get_proton_deflection()}, m_stddev_x{m_protonDeflection.get_proton_deflection_stddev_x()}, m_stddev_y{m_protonDeflection.get_proton_deflection_stddev_y()}, m_stddev_factor{0.0},
	m_hcal_xlow_ssafe{FiducialCutConsts::m_hcal_active_xlow_safe + m_stddev_factor*m_stddev_x}, m_hcal_xhigh_ssafe{FiducialCutConsts::m_hcal_active_xhigh_safe - m_stddev_factor*m_stddev_x},
	m_hcal_ylow_ssafe{FiducialCutConsts::m_hcal_active_ylow_safe + m_stddev_factor*m_stddev_y}, m_hcal_yhigh_ssafe{FiducialCutConsts::m_hcal_active_yhigh_safe - m_stddev_factor*m_stddev_y},
	m_targetintnum{static_cast<int>(m_protonDeflection.return_targetnum())}
	{
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

// Find Average Proton Deflection from the SBS magnet for a given kinematic setting and SBS magnet field scale.
// Found by seeing how much the LH2 data gets shifted by the SBS magnet for the given kinematic setting and SBS magnet field scale.
// std::vector<std::vector<double>> avg_proton_deflection_onHCal = { //{kinematicsetting(sbs #), sbs fieldscale, avg proton deflection from the analysis of the LH2 data}
// 	{4, 0, 0},
// 	{4, 30.0, -0.66863820},
// 	{4, 50.0, -1.1193552},
// 	{7, 85.0, -0.693155},
// 	{8, 0, 0},
// 	{8, 50.0, -0.661306},
// 	{8, 70.0, -0.920279},
// 	{8, 100.0, -1.29957},
// 	{11, 0, 0},
// 	{11, 100.0, -0.7144345},
// 	{14, 0, 0},
// 	{14, 70.0, -0.8045086}
// };

// double getavg_proton_deflection(int kine_num, double sbsfieldscale)
// {
// 	for (int i=0; i<avg_proton_deflection_onHCal.size(); ++i)
// 	{
// 		if( avg_proton_deflection_onHCal[i][0] == kine_num && avg_proton_deflection_onHCal[i][1] == sbsfieldscale )
// 		{
// 			double average_proton_deflection = avg_proton_deflection_onHCal[i][2];
// 			std::cout << "\nAverage proton deflection by the SBS magnet, calculated using LH2 data from SBS" << kine_num << " kinematic setting and the SBS field scale " << sbsfieldscale << "% = " << average_proton_deflection << " m\n";
// 			std::cout << "### IMPORTANT: The above deflection will be used to estimate the proton deflection for the fiducial cut ###\n";
// 			return average_proton_deflection;
// 		}
// 	}

// 	std::cout << "\n### ERROR: Either the kinematic setting number or the percentage SBS magnet field scale entered are invalid ####\n";
// 	return 0;
// }
////

//// Fiducial Cut ////

// Define the boundaries of HCal dimensions with respect to the HCal coordinate system.
// const double hcal_active_xlow = HCalConst::hcal_topXpos;
// const double hcal_active_xhigh = HCalConst::hcal_botXpos;
// const double hcal_active_ylow = HCalConst::hcal_rightYpos;
// const double hcal_active_yhigh = HCalConst::hcal_leftYpos;

// // Define the boundaries of HCal active area with some "safety margins" included.
// // Safety margin = Exclude the two outer columns and two outer blocks.
// const double hcal_active_xlow_safe =  hcal_active_xlow + HCalConst::hcalblk_h;
// const double hcal_active_xhigh_safe = hcal_active_xhigh - HCalConst::hcalblk_h;
// const double hcal_active_ylow_safe = hcal_active_ylow + HCalConst::hcalblk_w;
// const double hcal_active_yhigh_safe = hcal_active_yhigh - HCalConst::hcalblk_w;

//  bool fiducial_cut(double avg_proton_deflection, double xexpected_hcal, double yexpected_hcal)
// {
// 	// neutron coordinates: The kinematic calculation to find "xexpected_hcal" and "yexpected_hcal" does not take into account the SBS magnet deflection.
// 	// So essentially what we get is the neutron position from that calculation.
// 	double neutron_xpos_hcal{xexpected_hcal};
// 	double neutron_ypos_hcal{yexpected_hcal};
// 	// First we see whether neutron falls within the HCal. If not we return a "false" value as we do not need to proceed if the neutron already is outside the region we consider.
// 	if( neutron_xpos_hcal>=hcal_active_xhigh_safe || neutron_xpos_hcal<=hcal_active_xlow_safe || neutron_ypos_hcal>=hcal_active_yhigh_safe || neutron_ypos_hcal<= hcal_active_ylow_safe ) return false;

// 	// proton coordinates: Just add the expected average deflection expected to the neutron's x pos value.
// 	double proton_xpos_hcal{xexpected_hcal+avg_proton_deflection};
// 	double proton_ypos_hcal{yexpected_hcal};  
// 	// We are here only if we have established already above that the neutron is within the HCal region we consider.
// 	// Now see whether the proton will fall within the considered region in HCal.
// 	if( proton_xpos_hcal>=hcal_active_xhigh_safe || proton_xpos_hcal<=hcal_active_xlow_safe || proton_ypos_hcal>=hcal_active_yhigh_safe || proton_ypos_hcal<= hcal_active_ylow_safe ) return false;

// 	// We are here only if both the neutron and proton are within the HCal region considered. So we return a "true" value only at this stage.
// 	return true;
// }

////

// void draw_fiducialcut(TH2D* h2_dxdy,const char* h2_dxdy_filename)
// {
// 	TCanvas* C = new TCanvas();

// 	TLine* hcal_active_horizontal_low = new TLine(hcal_active_ylow,hcal_active_xlow,hcal_active_yhigh,hcal_active_xlow);
// 	TLine* hcal_active_horizontal_high = new TLine(hcal_active_ylow,hcal_active_xhigh,hcal_active_yhigh,hcal_active_xhigh);
// 	TLine* hcal_active_vertical_left = new TLine(hcal_active_yhigh,hcal_active_xlow,hcal_active_yhigh,hcal_active_xhigh);
// 	TLine* hcal_active_vertical_right = new TLine(hcal_active_ylow,hcal_active_xlow,hcal_active_ylow,hcal_active_xhigh);
// 	TLine* hcal_active_horizontalsafemarg_low = new TLine(hcal_active_ylow_safe,hcal_active_xlow_safe,hcal_active_yhigh_safe,hcal_active_xlow_safe);
// 	TLine* hcal_active_horizontalsafemarg_high = new TLine(hcal_active_ylow_safe,hcal_active_xhigh_safe,hcal_active_yhigh_safe,hcal_active_xhigh_safe);
// 	TLine* hcal_active_verticalsafemarg_left = new TLine(hcal_active_yhigh_safe,hcal_active_xlow_safe,hcal_active_yhigh_safe,hcal_active_xhigh_safe);
// 	TLine* hcal_active_verticalsafemarg_right = new TLine(hcal_active_ylow_safe,hcal_active_xlow_safe,hcal_active_ylow_safe,hcal_active_xhigh_safe);

// 	hcal_active_horizontal_low->SetLineColorAlpha(kGreen,0);
// 	hcal_active_horizontal_low->SetLineWidth(2);
// 	hcal_active_horizontal_high->SetLineColorAlpha(kGreen,0);
// 	hcal_active_horizontal_high->SetLineWidth(2);
// 	hcal_active_vertical_left->SetLineColorAlpha(kGreen,0);
// 	hcal_active_vertical_left->SetLineWidth(2);
// 	hcal_active_vertical_right->SetLineColorAlpha(kGreen,0);
// 	hcal_active_vertical_right->SetLineWidth(2);

// 	hcal_active_horizontalsafemarg_low->SetLineColorAlpha(kRed,0);
// 	hcal_active_horizontalsafemarg_low->SetLineWidth(2);
// 	hcal_active_horizontalsafemarg_low->SetLineStyle(9);
// 	hcal_active_horizontalsafemarg_high->SetLineColorAlpha(kRed,0);
// 	hcal_active_horizontalsafemarg_high->SetLineWidth(2);
// 	hcal_active_horizontalsafemarg_high->SetLineStyle(9);
// 	hcal_active_verticalsafemarg_left->SetLineColorAlpha(kRed,0);
// 	hcal_active_verticalsafemarg_left->SetLineWidth(2);
// 	hcal_active_verticalsafemarg_left->SetLineStyle(9);
// 	hcal_active_verticalsafemarg_right->SetLineColorAlpha(kRed,0);
// 	hcal_active_verticalsafemarg_right->SetLineWidth(2);
// 	hcal_active_verticalsafemarg_right->SetLineStyle(9);

// 	h2_dxdy->Draw("COLZ");
// 	hcal_active_horizontal_low->Draw();
// 	hcal_active_horizontal_high->Draw();
// 	hcal_active_vertical_left->Draw();
// 	hcal_active_vertical_right->Draw();
// 	hcal_active_horizontalsafemarg_low->Draw();
// 	hcal_active_horizontalsafemarg_high->Draw();
// 	hcal_active_verticalsafemarg_right->Draw();
// 	hcal_active_verticalsafemarg_left->Draw();

// 	//C->SaveAs(h2_dxdy_filename);
// }

#endif