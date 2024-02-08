#ifndef HCALCONSTANTS_H
#define HCALCONSTANTS_H

namespace HCalConst 
{
	const double hcalblk_h = 0.15494; // Height of all HCAL blocks in m from MC database
	const double hcalblk_w = 0.15875; // Width of all HCAL blocks in m from MC database
	const int Nhcal_rows = 24;
	const int Nhcal_columns = 12;
	const double hcal_width = hcalblk_w*Nhcal_columns; //1.905; // m
	const double hcal_height = hcalblk_h*Nhcal_rows; // 3.71856 m
	const double hcal_topXpos = -2.355005; // Distance to the top of the HCal w.r.t HCal origin.
	const double hcal_botXpos = 1.454995; // Distance to the bottom of the HCal w.r.t HCal origin.
	const double hcal_leftYpos = 0.92964; // Distance to the left side of HCal (when looking from up-stream to down-stream direction) w.r.t HCal origin.
	const double hcal_rightYpos = -0.92964; // Distance to the right side of HCal (when looking from up-stream to down-stream direction) w.r.t HCal origin.
	const double hcal_height_abovebeamline = -0.2897; // Vertical distance (X) of the HCal origin above beamline.

	//// HCal boundaries with safety margins in HCal coordinate system, for *Pass0 and 1* ////
	// Safety margin = Exclude 1.5 rows/columns from the edges - Andrew suggestion.
	const double hcal_active_xlow_safe_pass1 =  hcal_topXpos + 1.5*hcalblk_h;
	const double hcal_active_xhigh_safe_pass1 = hcal_botXpos - 1.5*hcalblk_h;
	const double hcal_active_ylow_safe_pass1 = -4.5*hcalblk_w;
	const double hcal_active_yhigh_safe_pass1 = 4.5*hcalblk_w;	
	
	//// Actual HCal boundaries for Pass2.
	const double hcal_active_xlow = -0.75 - 12*hcalblk_h;
	const double hcal_active_xhigh = -0.75 + 12*hcalblk_h;
	const double hcal_active_ylow = -6*hcalblk_w;
	const double hcal_active_yhigh = 6*hcalblk_w;

	//// HCal boundaries with safety margins in HCal coordinate system, for *Pass2 and simulation* ////
	// Safety margin = Exclude 1.5 rows/columns from the edges - Andrew suggestion.
	const double hcal_active_xlow_safe_pass2 =  -0.75 - 10.5*hcalblk_h;
	const double hcal_active_xhigh_safe_pass2 = -0.75 + 10.5*hcalblk_h;
	const double hcal_active_ylow_safe_pass2 = -4.5*hcalblk_w;
	const double hcal_active_yhigh_safe_pass2 = 4.5*hcalblk_w;	
}


#endif