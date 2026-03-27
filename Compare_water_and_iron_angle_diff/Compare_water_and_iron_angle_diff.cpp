// 2024/05/09
// kasumi
// based on "Check_upstream_base.cpp

// 22025/03/08
// kasumi
// à¯êîÇ3Ç∆9ÇçÇêèÍçáÇí«â¡

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#include <iomanip>

class output_format {
public:
	int groupid, chainid, pl, dpl;
	double sigma[2], chi2[4], md, dz;
};

void Calc_divide_angle_lateral(Momentum_recon::Mom_chain& c, int pl, double& dal, double& rms, double& dal_w, double& rms_w);
void Calc_divide_angle_radial(Momentum_recon::Mom_chain& c, int pl, double& dar, double& rms, double& dar_w, double& rms_w);
double Calc_md(Momentum_recon::Mom_chain& c, int pl, double& dz);
void output_file(std::string filename, std::vector<Momentum_recon::Event_information>& mom);

int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "usage1:file-in-momch file-out-txt\n");
		exit(1);
	}

	std::string file_in_momch = argv[1];
	std::string file_out_txt = argv[2];


	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch);
	output_file(file_out_txt, momch);

}
void output_file(std::string filename, std::vector<Momentum_recon::Event_information>& mom) {
	std::ofstream ofs(filename);

	double mean[2], mean_err[2], sd[2], dal, rms_l, dar, rms_r, chi2[4], error, sigma, vph, hitnum;
	double dal_w, rms_l_w, dar_w, rms_r_w;
	double mcs, rng;
	int count[2];
	int pl = 0;
	for (auto& ev : mom) {
		pl = ev.vertex_pl;

		for (auto& c : ev.chains) {
			Calc_divide_angle_lateral(c, pl, dal, rms_l, dal_w, rms_l_w);
			Calc_divide_angle_radial(c, pl, dar, rms_r, dar_w, rms_r_w);

			//if (!isfinite(dal) || !isfinite(dar) || !isfinite(rms_l) || !isfinite(rms_r))continue;
			chi2[0] = pow(dal_w / rms_l_w, 2);
			chi2[1] = pow(dar_w / rms_r_w, 2);
			chi2[2] = pow(dal / rms_l, 2);
			chi2[3] = pow(dar / rms_r, 2);

			mcs = c.ecc_mcs_mom[0];
			rng = c.ecc_range_mom[0];

			if (c.particle_flg == 2212) {
				mcs = c.ecc_mcs_mom[1];
				rng = c.ecc_range_mom[1];
			}
			ofs << std::fixed << std::right
				<< std::setw(5) << std::setprecision(0) << ev.groupid << " "
				<< std::setw(3) << std::setprecision(0) << c.chainid << " "
				<< std::setw(4) << std::setprecision(0) << c.particle_flg << " "
				<< std::setw(3) << std::setprecision(0) << c.base.size() << " "
				<< std::setw(3) << std::setprecision(0) << c.base.rbegin()->pl - c.base.begin()->pl + 1 << " "
				<< std::setw(3) << std::setprecision(0) << mcs << " "
				<< std::setw(3) << std::setprecision(0) << rng << " "
				<< std::setw(3) << std::setprecision(0) << pl << " "
				<< std::setw(8) << std::setprecision(5) << dal_w << " "
				<< std::setw(8) << std::setprecision(5) << rms_l_w << " "
				<< std::setw(10) << std::setprecision(4) << chi2[0] << " "
				<< std::setw(8) << std::setprecision(5) << dal << " "
				<< std::setw(8) << std::setprecision(5) << rms_l << " "
				<< std::setw(10) << std::setprecision(4) << chi2[2] << " "
				<< std::setw(8) << std::setprecision(5) << dar_w << " "
				<< std::setw(8) << std::setprecision(5) << rms_r_w << " "
				<< std::setw(10) << std::setprecision(4) << chi2[1] << " "
				<< std::setw(8) << std::setprecision(5) << dar << " "
				<< std::setw(8) << std::setprecision(5) << rms_r << " "
				<< std::setw(10) << std::setprecision(4) << chi2[3] << " "
				<< std::endl;




		}
	}
}

void Calc_divide_angle_lateral(Momentum_recon::Mom_chain& c, int pl, double& dal, double& rms, double& dal_w, double& rms_w) {
	dal = NAN;
	dal_w = NAN;

	double ax, ay, dax, day, dang, angle, sum2 = 0;
	double sum2_w = 0;
	int count = 0;
	int count_w = 0;
	int dpl = 133, dpl_w = 133;

	for (auto itr = c.base_pair.begin(); itr != c.base_pair.end(); itr++) {
		ax = itr->first.ax;
		ay = itr->first.ay;
		angle = sqrt(ax * ax + ay * ay);
		dax = itr->second.ax - ax;
		day = itr->second.ay - ay;
		if (angle < 0.01) {
			dang = dax;
		}
		else {
			dang = (- dax * ay + day * ax) / angle;
		}

		// water
		if (itr->first.pl % 2 == 1 && itr->second.pl == itr->first.pl + 1) {
			if (dpl_w > abs(itr->first.pl - pl)) {
				dal_w = dang;
			}
			if (itr->first.pl != pl) {
				sum2_w += dang * dang;
				count_w++;
			}
			dpl_w = abs(itr->first.pl - pl);

		}
		// iron
		if (itr->first.pl % 2 == 0 && itr->second.pl == itr->first.pl + 1) {
			if (dpl > abs(itr->first.pl - pl)) {
				dal = dang;
			}
			if (itr->first.pl != pl) {
				sum2 += dang * dang;
				count++;
			}
			dpl = abs(itr->first.pl - pl);
		}
	}

	if (count <= 1) {
		rms = NAN;
	}
	else {
		rms = sqrt(sum2 / count);
	}

	if (count_w <= 1) {
		rms_w = NAN;
	}
	else {
		rms_w = sqrt(sum2_w / count_w);
	}



}
void Calc_divide_angle_radial(Momentum_recon::Mom_chain& c, int pl, double& dar, double& rms, double& dar_w, double& rms_w) {
	dar = NAN;
	dar_w = NAN;

	double ax, ay, dax, day, dang, angle, sum2 = 0;
	double sum2_w = 0;
	int count = 0;
	int count_w = 0;
	int dpl = 133, dpl_w = 133;
	for (auto itr = c.base_pair.begin(); itr != c.base_pair.end(); itr++) {
		ax = itr->first.ax;
		ay = itr->first.ay;
		angle = sqrt(ax * ax + ay * ay);
		dax = itr->second.ax - ax;
		day = itr->second.ay - ay;
		if (angle < 0.01) {
			dang = day;
		}
		else {
			dang = (dax * ax + day * ay) / angle;
		}

		if (itr->first.pl == pl) {
			dar = dang;
		}

		// water
		if (itr->first.pl % 2 == 1 && itr->second.pl == itr->first.pl + 1) {
			if (abs(itr->first.pl - pl) < dpl) {
				dar_w = dang;
			}
			dpl_w = abs(itr->first.pl - pl);
			if (itr->first.pl != pl) {
				sum2_w += dang * dang;
				count_w++;
			}
		}
		// iron
		if (itr->first.pl % 2 == 0 && itr->second.pl == itr->first.pl + 1) {
			if (abs(itr->first.pl - pl) < dpl) {
				dar = dang;
			}
			dpl = abs(itr->first.pl - pl);
			if (itr->first.pl != pl) {
				sum2 += dang * dang;
				count++;
			}
		}

	}


	if (count <= 1) {
		rms = NAN;
	}
	else {
		rms = sqrt(sum2 / count);
	}
	if (count_w <= 1) {
		rms_w = NAN;
	}
	else {
		rms_w = sqrt(sum2_w / count_w);
	}

}

double Calc_md(Momentum_recon::Mom_chain& c, int pl, double& dz) {
	double md = -1;
	for (auto itr = c.base_pair.begin(); itr != c.base_pair.end(); itr++) {
		if (pl == itr->first.pl) {
			matrix_3D::vector_3D pos0, pos1, dir0, dir1;
			pos0.x = itr->first.x;
			pos0.y = itr->first.y;
			pos0.z = itr->first.z;
			dir0.x = itr->first.ax;
			dir0.y = itr->first.ay;
			dir0.z = 1;
			pos1.x = itr->second.x;
			pos1.y = itr->second.y;
			pos1.z = itr->second.z;
			dir1.x = itr->second.ax;
			dir1.y = itr->second.ay;
			dir1.z = 1;
			double extra[2], z_range[2];
			z_range[0] = pos0.z;
			z_range[1] = pos1.z;
			md = matrix_3D::minimum_distance(pos0, pos1, dir0, dir1, z_range, extra);
			dz = (fabs(extra[0]) + fabs(extra[1])) / 2;
			return md;

		}
	}
	return md;
}

