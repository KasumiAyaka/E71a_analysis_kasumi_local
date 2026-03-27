// 2025/11/05
// kasumi

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}

struct Thickness_of_Materials {
	double water = 2300;//um
	double iron = 500;//um
	double emulsion = 70 * 2;//um
	double psbase = 210;//um
	double ALpack = 109;//um
};
void Measere_Range_Momentum(std::vector<Momentum_recon::Event_information>& momch, std::string output);

int main(int argc, char** argv) {
	if (argc < 3) {
		fprintf(stderr, "usage: in.momch out.txt [out.momch]");
		exit(1);
	}
	std::string file_in_momch_ori = argv[1];//input momch
	std::string file_out_txt = argv[2];
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch_ori);
	Measere_Range_Momentum(momch, file_out_txt);

	if (argc == 4) {
		std::string file_out_momch = argv[3];
		Momentum_recon::Write_Event_information_extension(file_out_momch, momch);
	}

}
void Measere_Range_Momentum(std::vector<Momentum_recon::Event_information>& momch, std::string output) {
	Momentum_recon::Mom_chain chain;
	int cnt;
	int vpl, upl, dpl;
	std::ofstream ofs(output);

	int fe, wa, em, ps;
	double d = 20000.;// um

	double water = 0.2300;//cm
	double iron = 0.0500;//cm
	double emulsion = 0.0070 * 2;//cm
	double psbase = 0.0210;//cm
	double tan;

	double rng, R, rfe, rwa, rem, rps;
	for (auto& ev : momch) {
		for (auto& c : ev.chains) {
			fe = 0; wa = 0; em = 0; ps = 0;

			vpl = ev.vertex_pl;
			upl = c.base.rbegin()->pl;
			dpl = c.base.begin()->pl;
			ofs << std::right << std::fixed
				<< std::setw(5) << std::setprecision(0) << ev.groupid << " "
				<< std::setw(3) << std::setprecision(0) << ev.vertex_pl << " "
				<< std::setw(5) << std::setprecision(0) << c.chainid << " "
				<< std::setw(3) << std::setprecision(0) << c.base.size() << " "
				<< std::setw(3) << std::setprecision(0) << c.base.rbegin()->pl - c.base.begin()->pl + 1 << " "
				<< std::setw(3) << std::setprecision(0) << c.particle_flg << " ";

			if (dpl - vpl > 0) {
				//bwd
				ofs << std::right << std::fixed
					<< std::setw(3) << std::setprecision(0) << c.base.begin()->pl << " "
					<< std::setw(8) << std::setprecision(4) << c.base.begin()->ax << " "
					<< std::setw(8) << std::setprecision(4) << c.base.begin()->ay << " "
					<< std::setw(5) << std::setprecision(0) << c.base.begin()->m[0].ph + c.base.begin()->m[1].ph << " ";

				tan = sqrt(c.base.begin()->ax * c.base.begin()->ax + c.base.begin()->ay * c.base.begin()->ay);
				if (c.particle_flg != 13) {
					c.stop_flg = 2;
					if (c.base.rbegin()->x < d || c.base.rbegin()->x > 250000 - d || c.base.rbegin()->y < d || c.base.rbegin()->y > 250000 - d) {
						//side out
						c.stop_flg = 0;
					}
					if (c.base.rbegin()->pl > 130) {
						c.stop_flg = 0;
					}
				}
			}
			else
			{
				//fwd
				tan = sqrt(c.base.rbegin()->ax * c.base.rbegin()->ax + c.base.rbegin()->ay * c.base.rbegin()->ay);
				ofs << std::right << std::fixed
					<< std::setw(3) << std::setprecision(0) << c.base.rbegin()->pl << " "
					<< std::setw(8) << std::setprecision(4) << c.base.rbegin()->ax << " "
					<< std::setw(8) << std::setprecision(4) << c.base.rbegin()->ay << " "
					<< std::setw(5) << std::setprecision(0) << c.base.rbegin()->m[0].ph + c.base.rbegin()->m[1].ph << " ";

				if (c.particle_flg != 13) {
					c.stop_flg = 2;
					if (c.base.begin()->x < d || c.base.begin()->x > 250000 - d || c.base.begin()->y < d || c.base.begin()->y > 250000 - d) {
						//side out
						c.stop_flg = 0;
					}
					if (c.base.begin()->pl < 5) {
						c.stop_flg = 0;
					}
				}
			}

			if (dpl < 16) {
				// water ecc
				if (upl % 2 == 0) {//Šï”–‡
					fe = (upl - 15 + 1) / 2;
					wa = fe - 1;
					em = upl - 15 + 1;
					ps = upl - 15 + 1;

					fe = fe + 15 - dpl + 1;
					em = em + 15 - dpl;
					ps = ps + 15 - dpl;
				}
				else {
					fe = (upl - 15) / 2;
					wa = fe - 1;
					em = upl - 15;
					ps = upl - 15;

					fe = fe + 15 - dpl + 1;
					em = em + 15 - dpl;
					ps = ps + 15 - dpl;

				}
			}
			else {
				if (upl - dpl + 1 % 2 == 1) {// npl%2==1
					fe = (upl - dpl) / 2;
					wa = fe;
					em = upl - dpl + 1;
					ps = upl - dpl + 1;
				}
				else {
					if (upl % 2 == 0) {
						fe = (upl - dpl) / 2;
						wa = fe + 1;
						em = upl - dpl + 1;
						ps = upl - dpl + 1;
					}
					else {
						fe = (upl - dpl + 1) / 2;
						wa = fe - 1;
						em = upl - dpl + 1;
						ps = upl - dpl + 1;
					}
				}

			}

			//for (auto& b : c.base) {
			//	if (b.pl < 16) {
			//		fe++;
			//	}
			//	else if(b.pl)
			//	if(b)
			//	if(b.pl%2<)
			//}
			//proton
			rfe = std::pow(fe * tan / 0.0003, 1 / 1.75);
			rwa = std::pow(wa * tan / 0.0022, 1 / 1.75);
			rem = std::pow(em * tan / 0.001, 1 / 1.75);
			rps = std::pow(ps * tan / 0.0015, 1 / 1.75);
			rng = rfe + rwa + rem + rwa;
			c.ecc_range_mom[1] = rng;

			//pion
			rfe = std::pow(fe * tan / 0.0004, 1 / 1.75);
			rwa = std::pow(wa * tan / 0.001, 1 / 1.75);
			rem = std::pow(em * tan / 0.0007, 1 / 1.75);
			rps = std::pow(ps * tan / 0.0012, 1 / 1.75);
			rng = rfe + rwa + rem + rwa;
			c.ecc_range_mom[0] = rng;


			if (c.particle_flg == 2212) {
				rfe = std::pow(fe * tan / 0.0003, 1 / 1.75);
				rwa = std::pow(wa * tan / 0.0022, 1 / 1.75);
				rem = std::pow(em * tan / 0.001, 1 / 1.75);
				rps = std::pow(ps * tan / 0.0015, 1 / 1.75);
				rng = rfe + rwa + rem + rwa;
				c.ecc_range_mom[1] = rng;
				// judge ecc stop
				if (c.base.rbegin()->x < d || c.base.rbegin()->x>250000 - d || c.base.rbegin()->y < d || c.base.rbegin()->y>250000 - d) {
					c.ecc_range_mom[1] = 0;
				}
				if (c.base.begin()->x < d || c.base.begin()->x>250000 - d || c.base.begin()->y < d || c.base.begin()->y>250000 - d) {
					c.ecc_range_mom[1] = 0;
				}
				ofs << std::right << std::fixed
					<< std::setw(8) << std::setprecision(1) << c.ecc_mcs_mom[1] << " "
					<< std::setw(8) << std::setprecision(1) << rng << " "
					//<< std::setw(8) << std::setprecision(1) << c.bm_range_mom << " "
					<< std::setw(8) << std::setprecision(1) << c.ecc_mcs_mom_error[1][0] << " "
					<< std::setw(8) << std::setprecision(1) << c.ecc_mcs_mom_error[1][1] << " "
					<< std::endl;

			}
			else {
				rfe = std::pow(fe * tan / 0.0004, 1 / 1.75);
				rwa = std::pow(wa * tan / 0.001, 1 / 1.75);
				rem = std::pow(em * tan / 0.0007, 1 / 1.75);
				rps = std::pow(ps * tan / 0.0012, 1 / 1.75);
				rng = rfe + rwa + rem + rwa;
				c.ecc_range_mom[0] = rng;
				// judge ecc stop
				if (c.base.rbegin()->x < d || c.base.rbegin()->x > 250000 - d || c.base.rbegin()->y < d || c.base.rbegin()->y > 250000 - d) {
					c.ecc_range_mom[0] = 0;
				}
				if (c.base.begin()->x < d || c.base.begin()->x > 250000 - d || c.base.begin()->y < d || c.base.begin()->y > 250000 - d) {
					c.ecc_range_mom[0] = 0;
				}
				ofs << std::right << std::fixed
					<< std::setw(8) << std::setprecision(1) << c.ecc_mcs_mom[0] << " "
					<< std::setw(8) << std::setprecision(1) << rng << " "
					//<< std::setw(8) << std::setprecision(1) << c.bm_range_mom << " "
					<< std::setw(8) << std::setprecision(1) << c.ecc_mcs_mom_error[0][0] << " "
					<< std::setw(8) << std::setprecision(1) << c.ecc_mcs_mom_error[0][1] << " "
					<< std::endl;

			}
		}

	}

}
