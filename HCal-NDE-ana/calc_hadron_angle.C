// Calculates the hadron angle for a given Q2 and beam energy for both e(p,p')e' and gamma(p,n)pi+ reactions.

constexpr double piplus_mass = 0.13957039; //GeV
constexpr double pizero_mas = 0.1349768;   //GeV
constexpr double pimin_mass = 0.13957039;  //GeV
constexpr double n_mass = 0.93956542052;   //GeV
constexpr double p_mass = 0.93827208816;   //GeV

void calc_hadron_angle(double Q2, double E_beam, double theta_deg)
{
	double theta = theta_deg*TMath::DegToRad();
	double E_gamma = E_beam;

	// Define the four momentum of the target proton.
	TLorentzVector Pproton(0,0,0,p_mass);


	//// e(p,p')e' ////
	std::cout << '\n' << "*** Electron-Proton elastic scattering kinematics ***" << '\n';
	// Define the four momentum of the incident electron (beam)
	TLorentzVector Pelectron(0,0,E_beam,E_beam);
	// Scattered electron energy under the ultrarelativistic approxiamtion.
	double E_prime = Q2 / (2*E_beam*(1-cos(theta)));
	// Define the four momentum of the scatterred electron.
	TLorentzVector Pprimeelectron(E_prime*sin(theta),0,E_prime*cos(theta),E_prime);
	// Therefore, the four momentum of the scattered proton.
	TLorentzVector Pprimeproton;
	Pprimeproton = Pelectron + Pproton - Pprimeelectron;

	// Direction vector of the scattered proton 3 vector.
	TVector3 pprimeproton_vec_dirn = Pprimeproton.Vect().Unit();
	// Define the Z axis unit vector.
	TVector3 z_hat(0,0,1);
 
	// Let the angle between the Z axis and the pprimeproton_vec_dirn be "proton_angle"
	double cos_proton_angle = pprimeproton_vec_dirn.Dot(z_hat);
	double proton_angle = acos(cos_proton_angle)*TMath::RadToDeg();
	double proton_3momentum_mag = Pprimeproton.Vect().Mag();
	std::cout << "The scattered proton angle: " << std::fixed << std::setprecision(1) << proton_angle << " degrees\n";
	std::cout << "The scattered proton momentum: " << std::fixed << std::setprecision(2) << proton_3momentum_mag << " GeV/c\n";


	//// gamma(p,n)pi+ ////
	std::cout << '\n' << "*** photon + proton -> neutron + pi+ ****" << '\n';
	// Define the four momentum of the incident photon 
	TLorentzVector Pgamma(0,0,E_gamma,E_gamma);
	// Scattered pi+ energy under ultrarelativistic approximation.
	double E_piplus = ( std::pow(piplus_mass,2) + Q2 ) / (2*E_gamma*(1-cos(theta)));
	// Define the four momentum of the scattered pi+
	TLorentzVector Ppiplus(E_piplus*sin(theta),0,E_piplus*cos(theta),E_piplus);
	// Thereofore, the four momentum of the scattered neutron.
	TLorentzVector Pneutron;
	Pneutron = Pgamma + Pproton - Ppiplus;

	// Direction of the scattered neutron 3 vector.
	TVector3 pneutron_vec_dirn = Pneutron.Vect().Unit();
	// Let the angle between the Z axis and the pneutron_vec_dirn be "piplus_angle"
	double cos_neutron_angle = pneutron_vec_dirn.Dot(z_hat);
	double neutron_angle = acos(cos_neutron_angle)*TMath::RadToDeg();
	double neutron_3momentum_mag = Pneutron.Vect().Mag();
	std::cout << "The scattered neutron angle: " << std::fixed << std::setprecision(1) << neutron_angle << " degrees\n";
	std::cout << "The scattered neutron momentum: " << std::fixed << std::setprecision(2) << neutron_3momentum_mag << " GeV/c\n";

}