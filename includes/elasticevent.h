#include "beam_variables.h"
#include "constants.h"
#include "calc_HCalintersect.h"
#include "HCalConstants.h"
#include "fiducialcut.h"
#include <Math/Vector4D.h>
#include "TVector3.h"
#include "TMath.h"

#ifndef ELASTICEVENT_H
#define ELASTICEVENT_H

class ElasticEvent
{
private:
	
	// Read-in from the main script.
	TChain* const m_C; // TChain that holds the "T" TTree(s) of input replayed ROOT file(s).
	const int m_ReplayPassNum;
	const int m_SBSKineNum; 
	const TString m_target;
	const double m_SBSFieldScale;
	const double m_cut_MinPreShE;
	const double m_cut_MinGEMHitsOnTrack;
	const double m_cut_MinBBTrP;
	const double m_cut_MaxTrackChi2divNdof;
	const double m_cut_VzCutUpStream;
	const double m_cut_VzCutDwnStream;
	const double m_cut_CoinCutMean;
	const double m_cut_CoinCutSigam;
	const double m_cut_CoinCutSigamFactor;
	const double m_cut_CoinCutLowLimit;
	const double m_cut_CoinCutHighLimit;
		
	// Variables needed to copy information from "T" tree.
	double m_bbtrn{0};
	const static int m_MAXNTRACKS {100};
	const static int m_MAXHCALCLUS {20}; // #### MAY NEED TO INCREASE THIS IN PASS-2 ANALYSSIS! ###
	double m_bbtrvx[m_MAXNTRACKS];
	double m_bbtrvy[m_MAXNTRACKS];
	double m_bbtrvz[m_MAXNTRACKS];
	double m_bbtrth[m_MAXNTRACKS];
	double m_bbtrx[m_MAXNTRACKS];
	double m_bbtrp[m_MAXNTRACKS];
	double m_bbtrpx[m_MAXNTRACKS];
	double m_bbtrpy[m_MAXNTRACKS];
	double m_bbtrpz[m_MAXNTRACKS];
	double m_bbtrrth[m_MAXNTRACKS];
	double m_bbtrchi2[m_MAXNTRACKS];
	double m_bbtrndof[m_MAXNTRACKS];
	double m_bbgemtrnhits[m_MAXNTRACKS];
	double m_sbshcale {0.};
	double m_sbshcalx {0.};
	double m_sbshcaly {0.};
	double m_sbshcalatimeblk {0.};
	double* m_sbshcalcluse = new double[m_MAXHCALCLUS];
	double* m_sbshcal_heclus_e = new double[m_MAXHCALCLUS];
	double* m_sbshcalclusatime = new double[m_MAXHCALCLUS];
	double* m_sbshcal_heclus_atime = new double[m_MAXHCALCLUS];
	double* m_sbshcalclusx = new double[m_MAXHCALCLUS];
	double* m_sbshcal_heclus_x = new double[m_MAXHCALCLUS];
	double* m_sbshcalclusy = new double[m_MAXHCALCLUS];
	double* m_sbshcal_heclus_y = new double[m_MAXHCALCLUS];
	double m_bbshe{0.};
	double m_bbpse{0.};
	double m_shadctime{0.};
	double m_adctimediff_hcalsh{0.};
	
	// Kinematic parameters.
	GMnKinInfo m_kininfo;
	const double m_Ebeam;
	const double m_HCaldist;
	const double m_HCalangle;

	ROOT::Math::PxPyPzEVector m_Pbeam; // Four momentum of the beam electrons.
	const double m_targetmass {0.5*(Constants::n_mass + Constants::p_mass)};
	ROOT::Math::PxPyPzEVector m_Ptarg; // Four momentum of the target.
	ROOT::Math::PxPyPzEVector m_q; //Four momentum trasnferred to the scattered nucleon.

	HCalVectors m_hcalvect; // Calss to calculate intersection point of the hadron on HCal under nueutron hypothesis.

	FiducialCut m_fiducialcut; // Class to apply HCal avtive area and fiducial cuts.
	const int m_targetintnum;
		
public:

	ElasticEvent(TChain* const C, const int num_ReplayPassNum, const int num_SBSKine, const TString target, const double num_SBSFieldScale, const double cut_MinPreShE, const double cut_MinGEMHitsOnTrack, const double cut_MinBBTrP, const double cut_MaxTrackChi2divNdof, const double cut_VzCutUpStream, const double cut_VzCutDwnStream, const double cut_CoinCutMean, const double cut_CoinCutSigam, const double cut_CoinCutSigmaFactor) 
	: m_C{C}, m_ReplayPassNum{num_ReplayPassNum}, m_SBSKineNum{num_SBSKine}, m_target{target}, m_SBSFieldScale{num_SBSFieldScale}, m_cut_MinPreShE{cut_MinPreShE}, m_cut_MinGEMHitsOnTrack{cut_MinGEMHitsOnTrack}, m_cut_MinBBTrP{cut_MinBBTrP}, m_cut_MaxTrackChi2divNdof{cut_MaxTrackChi2divNdof}, m_cut_VzCutUpStream{cut_VzCutUpStream}, m_cut_VzCutDwnStream{cut_VzCutDwnStream},
	m_cut_CoinCutMean{cut_CoinCutMean}, m_cut_CoinCutSigam{cut_CoinCutSigam}, m_cut_CoinCutSigamFactor{cut_CoinCutSigmaFactor}, m_cut_CoinCutLowLimit{cut_CoinCutMean-(cut_CoinCutSigam*cut_CoinCutSigmaFactor)}, m_cut_CoinCutHighLimit{cut_CoinCutMean+(cut_CoinCutSigam*cut_CoinCutSigmaFactor)},
	m_kininfo{num_SBSKine}, m_Ebeam{m_kininfo.return_BeamEnergy()}, m_HCaldist{m_kininfo.return_HCalDis()}, m_HCalangle{m_kininfo.return_HCalTheta()},
	m_fiducialcut{m_ReplayPassNum,m_SBSKineNum,m_target,m_SBSFieldScale}, m_targetintnum{m_fiducialcut.return_targetintnum()}
	{
		m_C->SetBranchStatus("*", 0);
		m_C->SetBranchStatus("bb.tr.n", 1);
		m_C->SetBranchStatus("bb.tr.vx", 1);
		m_C->SetBranchStatus("bb.tr.vy", 1);
		m_C->SetBranchStatus("bb.tr.vz", 1);
		m_C->SetBranchStatus("bb.tr.th", 1);
		m_C->SetBranchStatus("bb.tr.x", 1);
		m_C->SetBranchStatus("bb.tr.p", 1);
		m_C->SetBranchStatus("bb.tr.px", 1);
		m_C->SetBranchStatus("bb.tr.py", 1);
		m_C->SetBranchStatus("bb.tr.pz", 1);
		m_C->SetBranchStatus("bb.tr.r_th", 1);
		m_C->SetBranchStatus("bb.tr.chi2", 1);
		m_C->SetBranchStatus("bb.tr.ndof", 1);
		m_C->SetBranchStatus("bb.gem.track.nhits", 1);
		m_C->SetBranchStatus("sbs.hcal.x", 1);
		m_C->SetBranchStatus("sbs.hcal.y", 1);
		m_C->SetBranchStatus("sbs.hcal.e", 1);
		m_C->SetBranchStatus("sbs.hcal.atimeblk", 1);
		m_C->SetBranchStatus("sbs.hcal.clus.e", 1);
		m_C->SetBranchStatus("sbs.hcal.clus.atime", 1);
		// m_C->SetBranchStatus("Ndata.sbs.hcal.clus.atime", 1);
		m_C->SetBranchStatus("sbs.hcal.clus.x", 1);
		// m_C->SetBranchStatus("Ndata.sbs.hcal.clus.x", 1);
		m_C->SetBranchStatus("sbs.hcal.clus.y", 1);
		// m_C->SetBranchStatus("Ndata.sbs.hcal.clus.y", 1);
		m_C->SetBranchStatus("bb.sh.e", 1);
		m_C->SetBranchStatus("bb.ps.e", 1);
		m_C->SetBranchStatus("bb.sh.atimeblk", 1);

		m_C->SetBranchAddress("bb.tr.n", &m_bbtrn);
		m_C->SetBranchAddress("bb.tr.vx", m_bbtrvx);
		m_C->SetBranchAddress("bb.tr.vy", m_bbtrvy);
		m_C->SetBranchAddress("bb.tr.vz", m_bbtrvz);
		m_C->SetBranchAddress("bb.tr.th", m_bbtrth);
		m_C->SetBranchAddress("bb.tr.x", m_bbtrx);
		m_C->SetBranchAddress("bb.tr.p", m_bbtrp);
		m_C->SetBranchAddress("bb.tr.px", m_bbtrpx); 
		m_C->SetBranchAddress("bb.tr.py", m_bbtrpy); 
		m_C->SetBranchAddress("bb.tr.pz", m_bbtrpz);
		m_C->SetBranchAddress("bb.tr.r_th", m_bbtrrth);
		m_C->SetBranchAddress("bb.tr.chi2", m_bbtrchi2);
		m_C->SetBranchAddress("bb.tr.ndof", m_bbtrndof);
		m_C->SetBranchAddress("bb.gem.track.nhits", m_bbgemtrnhits);
		m_C->SetBranchAddress("sbs.hcal.x", &m_sbshcalx);
		m_C->SetBranchAddress("sbs.hcal.y", &m_sbshcaly);
		m_C->SetBranchAddress("sbs.hcal.e", &m_sbshcale);
		m_C->SetBranchAddress("sbs.hcal.atimeblk", &m_sbshcalatimeblk);
		m_C->SetBranchAddress("sbs.hcal.clus.e", m_sbshcalcluse);
		m_C->SetBranchAddress("sbs.hcal.clus.atime", m_sbshcalclusatime);
		m_C->SetBranchAddress("sbs.hcal.clus.x", m_sbshcalclusx);
		m_C->SetBranchAddress("sbs.hcal.clus.y", m_sbshcalclusy);
		m_C->SetBranchAddress("bb.sh.e", &m_bbshe);
		m_C->SetBranchAddress("bb.ps.e", &m_bbpse);
		m_C->SetBranchAddress("bb.sh.atimeblk", &m_shadctime);

		m_Pbeam.SetPxPyPzE(0,0,m_Ebeam,m_Ebeam);
		m_Ptarg.SetPxPyPzE(0,0,0,m_targetmass);

		if ( m_ReplayPassNum == 2 ) m_hcalvect.make_HCal_vectors_NoVerticalOffset(m_HCaldist, m_HCalangle);
		else if ( m_ReplayPassNum == 0 || m_ReplayPassNum == 1 ) m_hcalvect.make_HCal_vectors_pass0and1(m_HCaldist, m_HCalangle);
	}

	int getEntry(int n) //Copies the enries of the "T" to the above defined member variables.
	{
		int haveData = m_C->GetEntry(n); // Return 0 if there is no data, i.e. End of the TTree.

		for ( int i = 0; i < m_MAXHCALCLUS; i++ )
		{
			m_sbshcal_heclus_e[i] = m_sbshcalcluse[i];
			m_sbshcal_heclus_atime[i] = m_sbshcalclusatime[i];
			m_sbshcal_heclus_x[i] = m_sbshcalclusx[i];
			m_sbshcal_heclus_y[i] = m_sbshcalclusy[i];
		}

		return haveData;
	}

	
	// --- BigBite cuts --- //

	bool passTrackChi2Cut() // Called by the "passBigBiteCuts() function below"
	{
		m_bbtrchi2divndof = m_bbtrchi2[0] / m_bbtrndof[0];

		if ( m_bbtrchi2divndof > m_cut_MaxTrackChi2divNdof ) return false;

		return true;
	}

	bool passGoodElectronCuts()
	{	
		// ***Track quality cuts*** //
		if ( m_bbpse < m_cut_MinPreShE ) return false; // Pre-shower energy cut for pion rejection.
		if ( m_bbgemtrnhits[0] < m_cut_MinGEMHitsOnTrack ) return false; // Min. number of gem hits on track.
		if ( m_bbtrp[0] < m_cut_MinBBTrP ) return false; // Min. BigBite track momentum as per the kinematic setting and the Q^2 of interest.
		if ( !passTrackChi2Cut() ) return false; // Use only tracks with good chi2.
		if ( m_bbtrvz[0] < m_cut_VzCutUpStream || m_bbtrvz[0] > m_cut_VzCutDwnStream  ) return false; // BigBite track vertex Z cut.
		//// 

		return true;
	}

private:

	double m_neutronhypthspred_x{0.}; // Calculated x hit position of the hadron under neutron hypothesis in HCal local detector coordinates.
	double m_nuetronhypthspred_y{0.}; // Calculated y hit position of the hadron under neutron hypothesis in HCal local detector coordinates.
	bool m_passfiducialcut {false};

public:

	void calcNeutronHypthsHCalIntersect() // Calculates the point of intersection of the hadron, under the hyphothesis that it is neutron in HCal local coordinates.
	{
		m_hcalvect.calc_expected_xyonHCal( m_q, m_bbtrvz, m_bbtrvx, m_bbtrvy);
		m_neutronhypthspred_x = m_hcalvect.return_xexpected();
		m_nuetronhypthspred_y = m_hcalvect.return_yexpected();
	}

	bool passFiducialCut() // "Fiducial/p&n acceptance matching" cut for LD2 data analysis.
	{		
		switch ( m_targetintnum )
		{
			case 0: // 0 = LD2, appy the fiducial cut.
				return m_fiducialcut.pass_FiducialCut( m_neutronhypthspred_x, m_nuetronhypthspred_y );

			case 1: // 1 = LH2, apply the HCal active area cut.
				return m_fiducialcut.pass_HitHCalCut( m_neutronhypthspred_x, m_nuetronhypthspred_y );

			default: // If m_targetintnum is neither 0 nor 1.
				return true; // Let the Fiducial cut pass so the event can proceed in analysis.
		}
	}

	void checkFiducialCut() // Check whther the fiducial cut passed or failed and update the variable.
	{
		m_passfiducialcut = passFiducialCut();
	}

		
private: 	

	// Variables to hold the secondary calculation results.
	double m_bbtrchi2divndof {0.}; // Track chi2/ndof.
	double m_bbshpse{0.}; // bb.sh.e + bb.ps.e = Total energy deposited in the sh and ps.
	double m_eoverp{0.};
	double m_bbphitgt{0.}; // Azimuthal angle of pi+ detected by the BigBite spectrometer.
	double m_bbpoltgt{0.}; // Polar angle of pi+ detected by the BigBite spectrometer.
	double m_Q2 {0.}; // Four-momentum trasfer squared.
	double m_W2 {0.}; // Invariant mass squared.
	int m_sbshcal_bestclus_indx {0}; // Best cluster ID.
	ROOT::Math::XYZVector m_q3vector; // 3 dimentional spatial vecotr of q.
	ROOT::Math::XYZVector m_p3vector; // Vector from the Hall origin.

	double m_sbsneutronhcaldx{0.}; // X coordinate difference between the measured and prediceted neutron position values. 
	double m_sbsneutronhcaldy{0.}; // Y coordinate difference between the measured and prediceted neutron position values.	

	double m_hcalshcointime{0.}; // Coin. time/time diff. between HCal best cluster and SH cluster.
	bool m_passcointcut {false}; // Fass/Fail coint cut. 

	double m_thetapq_rad {0.};
	double m_thetapq_deg {0.};

	double m_p_perp {0.}; // Missing momentum perpendicular component approximation.
	
public: 
 
	void calcBBTrackAngles() //Function to calculate bb track polar and azimuthal scattering angles.
	{
		m_bbphitgt = atan2(m_bbtrpy[0], m_bbtrpx[0]); //Azimuthal scattering angle of pi+ in radians.
		m_bbpoltgt = acos(m_bbtrpz[0]/m_bbtrp[0]); //Polar scattering angle of pi+ in radians.
	}

	void calcQ2andW2()
	{
		// Calculating q 
		ROOT::Math::PxPyPzEVector kprime; //Four vector of the scattered electron.
		kprime.SetPxPyPzE(m_bbtrpx[0],m_bbtrpy[0],m_bbtrpz[0],m_bbtrp[0]);
		m_q = m_Pbeam - kprime; 
		m_Q2 = -m_q.M2();
		m_W2 = (m_Ptarg+m_q).M2(); // Calculates the invariant mass squared (W^2) of the virtual photon - nucleon system.		
	}

	void getBestHCalClusIndx( int bestclus_indx )
	{
		m_sbshcal_bestclus_indx = bestclus_indx;
	}

	void calcCoinTime()
	{
		m_hcalshcointime = m_sbshcal_heclus_atime[m_sbshcal_bestclus_indx] - m_shadctime;
	}

	void checkCoinCut()
	{
		calcCoinTime(); // First calculate coin. time.

		if ( m_hcalshcointime > m_cut_CoinCutLowLimit && m_hcalshcointime < m_cut_CoinCutHighLimit )
		{
			m_passcointcut = true;
		}
		else
		{
			m_passcointcut = false;
		}
	}

	void make_p3Vector()
	{
		ROOT::Math::XYZVector hcal_origintocluspos( m_sbshcalclusy[m_sbshcal_bestclus_indx], -m_sbshcalclusx[m_sbshcal_bestclus_indx], 0 );
		m_p3vector = m_hcalvect.return_VertextoHCalOrigin() + hcal_origintocluspos;
	}

	void make_q3Vector()
	{
		m_q3vector = m_q.Vect();
	}
	
	void calcThetapq()
	{
		make_p3Vector();
		make_q3Vector();

		double cos_thetapq = (m_p3vector.Dot(m_q3vector)) / (m_p3vector.r()*m_q3vector.r());

		m_thetapq_rad = acos(cos_thetapq);
		m_thetapq_deg = m_thetapq_rad*TMath::RadToDeg();
	}

	void calcPperp()
	{
		m_p_perp = m_q3vector.r() * m_thetapq_rad;
	}


	// --- Return variables to the output TTree --- //

	double return_BBTrGEMNhits()
	{
		return m_bbgemtrnhits[0];
	}

	double return_BBTrChi2divNdof()
	{
		return m_bbtrchi2divndof;
	}

	double return_BBTrn()
	{
		return m_bbtrn;
	}

	double return_BBTrVz() 
	{
		return m_bbtrvz[0];
	}
	
	double return_BBTrth()
	{
		return m_bbtrth[0];
	}

	double return_BBTrx() 
	{
		return m_bbtrx[0];
	}
	
	double return_BBTrP()
	{
		return m_bbtrp[0];
	}
	
	double return_BBTrPx()
	{
		return m_bbtrpx[0];
	}

	double return_BBtrPy()
	{
		return m_bbtrpy[0];
	}
	
	double return_BBSHe()
	{
		return m_bbshe;
	}
	double return_BBPSe()
	{
		return m_bbpse;
	}

	double return_BBSHPSe()
	{
		m_bbshpse = m_bbshe + m_bbpse;
		return m_bbshpse;
	}

	double return_EoverP() 
	{
		m_bbshpse = m_bbshe + m_bbpse;
		m_eoverp = m_bbshpse / m_bbtrp[0];
		return m_eoverp;
	}

	double return_BBphitgt()
	{
		double bbphitgt_deg = m_bbphitgt*TMath::RadToDeg();
		return bbphitgt_deg;
	}	 

	double return_BBpoltgt()
	{
		double bbpoltgt_deg = m_bbpoltgt*TMath::RadToDeg();
		return bbpoltgt_deg;
	}

	double return_Q2()
	{
		return m_Q2;
	}

	double return_W2()
	{
		return m_W2;
	}

	double return_thetapq_rad()
	{
		return m_thetapq_rad;
	}

	double return_thetapq_deg()
	{
		return m_thetapq_deg;
	}

	double return_pperp()
	{
		return m_p_perp;
	}

	double return_nHypthsPredx()
	{
		return m_neutronhypthspred_x;
	}

	double return_nHypthsPredy()
	{
		return m_nuetronhypthspred_y;
	}

	bool return_passFiducialCut()
	{
		return m_passfiducialcut;
	}

	double return_CoinTime()
	{
		return m_hcalshcointime;
	}

	bool return_passCoinCut()
	{
		return m_passcointcut;
	}

	double return_SBSHCalx()
	{
		return m_sbshcalx;
	}

	double return_SBSHCaly()
	{
		return m_sbshcaly;
	}

	double return_SBSHCale()
	{
		return m_sbshcale;
	}

	double return_SBSHCalAtimeBlk()
	{
		return m_sbshcalatimeblk;
	}

	double return_ADCTimeDiffHCalSH()
	{
		return m_adctimediff_hcalsh;
	}

	double return_SHADCTime()
	{
		return m_shadctime;
	}

	double* return_HCalClusEArryPtr()
	{
		return m_sbshcalcluse;
	}

	double* return_HCalClusAtimeArryPtr()
	{
		return m_sbshcalclusatime;
	}

	double* return_HCalClusXArryPtr()
	{
		return m_sbshcalclusx;
	}

	double* return_HCalClusYArryPtr()
	{
		return m_sbshcalclusy;
	}

	double* return_HCalHEClusEArryPtr()
	{
		return m_sbshcal_heclus_e;
	}

	double* return_HCalHEClusAtimeArryPtr()
	{
		return m_sbshcal_heclus_atime;
	}

	double* return_HCalHEClusXArryPtr()
	{
		return m_sbshcal_heclus_x;
	}

	double* return_HCalHEClusYArryPtr()
	{
		return m_sbshcal_heclus_y;
	}

	// ---- **** ---- //

};

#endif