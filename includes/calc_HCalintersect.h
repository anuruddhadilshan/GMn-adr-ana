#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "HCalConstants.h"

#ifndef CALC_HCALINTERSECT_H
#define CALC_HCALINTERSECT_H

// A single class that does all the HCal vector calculations.

class HCalVectors
{

	ROOT::Math::XYZVector m_hcal_origin; // Vector pointing to the origin of HCal starting from the origin of the Hall coordinate system.
	ROOT::Math::XYZVector m_hcalXaxis;				
	ROOT::Math::XYZVector m_hcalYaxis;
	ROOT::Math::XYZVector m_hcalZaxis;
	
public:
	void make_HCal_vectors_pass0and1(double hcaldist, double hcaltheta)
	{
		m_hcalXaxis.SetXYZ(0,-1,0);
		m_hcalZaxis.SetXYZ(-sin(hcaltheta),0,cos(hcaltheta));
		m_hcalYaxis = m_hcalZaxis.Cross(m_hcalXaxis).Unit();

		//3D vector defining the position of the HCal origin w.r.t the Hall origin. Hall coordinate system = x-beam left, y-verically up, z-beam downstream.
		m_hcal_origin = hcaldist*m_hcalZaxis + HCalConst::hcal_height_abovebeamline*m_hcalXaxis;
	}

	void make_HCal_vectors_NoVerticalOffset(double hcaldist, double hcaltheta)
	{
		m_hcalXaxis.SetXYZ(0,-1,0);
		m_hcalZaxis.SetXYZ(-sin(hcaltheta),0,cos(hcaltheta));
		m_hcalYaxis = m_hcalZaxis.Cross(m_hcalXaxis).Unit();

		//3D vector defining the position of the HCal origin w.r.t the Hall origin. Hall coordinate system = x-beam left, y-verically up, z-beam downstream.
		m_hcal_origin = hcaldist*m_hcalZaxis;
	}

	void make_HCal_vectors_forsim(double hcaldist, double hcaltheta)
	{
		m_hcalXaxis.SetXYZ(0,-1,0);
		m_hcalZaxis.SetXYZ(-sin(hcaltheta),0,cos(hcaltheta));
		m_hcalYaxis = m_hcalZaxis.Cross(m_hcalXaxis).Unit();

		//3D vector defining the position of the HCal origin w.r.t the Hall origin. Hall coordinate system = x-beam left, y-verically up, z-beam downstream.
		m_hcal_origin = hcaldist*m_hcalZaxis;
	}

private: 
	ROOT::Math::XYZVector m_vertex;
	ROOT::Math::XYZVector m_vertextoHCalorigin;
	double m_hcalZdistancefromvertex {0.};
	ROOT::Math::XYZVector m_qdirection;
	ROOT::Math::XYZVector m_neutronvecdir;
	double m_hadronvectormagnitude {0.};
	double m_neutronvectormagnitude {0.};
	ROOT::Math::XYZVector m_hcalintersect;
	double m_xexpected_hcal {0.};
	double m_yexpected_hcal {0.};

public:
	void calc_expected_xyonHCal(ROOT::Math::PxPyPzEVector& q, double vz[100], double vx[100], double vy[100])
	{
		m_vertex.SetXYZ(vx[0],vy[0],vz[0]); // Vertex vector (the vector going from the origin of the Hall coordinate system(~middle of the target) to the (e,e'p/n) interaction point)
		m_vertextoHCalorigin = m_hcal_origin - m_vertex;
		m_hcalZdistancefromvertex = m_vertextoHCalorigin.Dot(m_hcalZaxis);
		m_qdirection = q.Vect().Unit(); //Unit vector defining the direction of the spacial component of the q 4-vector..
		m_hadronvectormagnitude = m_hcalZdistancefromvertex/m_qdirection.Dot(m_hcalZaxis);
		m_hcalintersect = m_hadronvectormagnitude*m_qdirection + m_vertex;

		m_xexpected_hcal = (m_hcalintersect-m_hcal_origin).Dot(m_hcalXaxis);
		m_yexpected_hcal = (m_hcalintersect-m_hcal_origin).Dot(m_hcalYaxis);
	}

	void calc_expected_NeutronxyonHCal(ROOT::Math::PxPyPzEVector& pNeutron, double vz[100])
	{
		m_vertex.SetXYZ(0,0,vz[0]); // Vertex vector (the vector going from the origin of the Hall coordinate system(~middle of the target) to the (e,e'p/n) interaction point)
		m_vertextoHCalorigin = m_hcal_origin - m_vertex;
		m_hcalZdistancefromvertex = m_vertextoHCalorigin.Dot(m_hcalZaxis);
		m_neutronvecdir = pNeutron.Vect().Unit(); //Unit vector defining the direction of the spacial component of the pNeutron 4-vector.
		m_neutronvectormagnitude = m_hcalZdistancefromvertex/m_neutronvecdir.Dot(m_hcalZaxis);
		m_hcalintersect = m_neutronvectormagnitude*m_neutronvecdir + m_vertex;

		m_xexpected_hcal = (m_hcalintersect-m_hcal_origin).Dot(m_hcalXaxis);
		m_yexpected_hcal = (m_hcalintersect-m_hcal_origin).Dot(m_hcalYaxis);
	}

	double return_xexpected()
	{
		return m_xexpected_hcal;
	}

	double return_yexpected()
	{
		return m_yexpected_hcal;
	}

	ROOT::Math::XYZVector return_q3directionvect()
	{
		return m_qdirection;
	}

	ROOT::Math::XYZVector return_HCalZaxis()
	{
		return m_hcalZaxis;
	}

	ROOT::Math::XYZVector return_VertextoHCalOrigin()
	{
		return m_vertextoHCalorigin;
	}

	double return_hcalZdistancefromvertex()
	{
		return m_hcalZdistancefromvertex;
	}

};

#endif