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
	// •úŽË’·
	double Xwater = 3610000;//um
};


void Calc_angle_diff(double pbeta, double X, std::string out_txt);
void Calc_angle_diff_per_angle(double X, std::string out_txt);

int main(int argc, char** argv) {
	if (argc < 3) {
		fprintf(stderr, "usage: out.txt mom/mode");
		exit(1);
	}

	std::string out_txt = argv[1];//input momch
	double mom = std::stod(argv[2]);

	Thickness_of_Materials m;

	if (mom > 0) {
		Calc_angle_diff(mom, m.Xwater, out_txt);
	}
	else {
		Calc_angle_diff_per_angle(m.Xwater, out_txt);
	}

}
void Calc_angle_diff(double pbeta, double X, std::string out_txt) {
	std::ofstream ofs(out_txt);
	for (double dz = -3200; dz < 50;) {
		if (dz == 0) {
			dz -= 10;
		}
		for (double tan = 0.00001; tan < 8.0;) {
			double L, L2;
			L = fabs(dz) * tan;
			L2 = L / X;

			double dang = (1 / pow(3, 1 / 2)) * (13.6 / pbeta) * pow(L, 1 / 2) * (1 + 0.038 * log(L));// d(\theta_rms)


			double md_exp, L_scat, L_meas;

			L_scat = L * dang;
			L_meas = L * tan * pow(2, 1 / 2) * pow(0.45 * 0.45 + tan * tan * 16, 1 / 2)/210;//??

			md_exp = pow(L_scat * L_scat + L_meas * L_meas, 0.5);

			ofs << std::fixed << std::right
				<< std::setw(10) << std::setprecision(1) << pbeta << " "
				<< std::setw(10) << std::setprecision(1) << dz << " "
				<< std::setw(10) << std::setprecision(4) << tan << " "
				<< std::setw(10) << std::setprecision(1) << L << " "
				<< std::setw(10) << std::setprecision(4) << dang << " "
				<< std::setw(10) << std::setprecision(4) << tan * pow(2, 1 / 2) * pow(0.45 * 0.45 + tan * tan * 16, 1 / 2) /210<< " "
				<< std::setw(10) << std::setprecision(1) << L_scat << " "
				<< std::setw(10) << std::setprecision(1) << L_meas << " "
				<< std::setw(10) << std::setprecision(1) << md_exp << " "
				<< std::endl;
			tan += 0.2;
		}
		dz += 50;
	}

}

void Calc_angle_diff_per_angle(double X, std::string out_txt) {
	std::ofstream ofs(out_txt);

	double angle_min, angle_max, tan, pbeta;
	for (int i_ang = 0; i_ang < 20; i_ang++) {
		if (i_ang < 7) {
			// <=0.7
			angle_min = i_ang * 0.1;
			angle_max = (i_ang + 1) * 0.1;
		}
		else if (i_ang < 11) {
			// <=1.5
			angle_min = (i_ang - 7) * 0.2 + 0.7;
			angle_max = (i_ang - 7 + 1) * 0.2 + 0.7;
		}
		else if (i_ang < 15) {
			// <= 3.5
			angle_min = (i_ang - 11) * 0.4 + 1.5;
			angle_max = (i_ang - 11 + 1) * 0.4 + 1.5;
		}
		else if (i_ang < 19) {
			// <= 5.5
			angle_min = (i_ang - 15) * 0.6 + 3.1;
			angle_max = (i_ang - 15 + 1) * 0.6 + 3.1;
		}
		else {
			// 5.5 <= tan
			angle_min = (i_ang - 19) * 0.6 + 5.5;
			angle_max = 15.0;
		}


		tan = 0.5 * (angle_min + angle_max);
		if (i_ang==19) {
			tan = 6.5;
		}


		for (double dz = -3200; dz < 50;) {
			if (dz == 0) {
				dz -= 10;
			}
			double L, L2;
			L = fabs(dz) * tan;
			L2 = L / X;

			for (int pbeta = 50; pbeta < 850;) {
				double dang = (1 / pow(3, 1 / 2)) * (13.6 / pbeta) * pow(L, 1 / 2) * (1 + 0.038 * log(L));// d(\theta_rms)


				double md_exp, L_scat, L_meas;

				L_scat = L * dang;
				L_meas = L * tan * pow(2, 1 / 2) * pow(0.45 * 0.45 + tan * tan * 16, 1 / 2) / 210;//??

				md_exp = pow(L_scat * L_scat + L_meas * L_meas, 0.5);

				ofs << std::fixed << std::right
					<< std::setw(10) << std::setprecision(1) << pbeta << " "
					<< std::setw(10) << std::setprecision(1) << dz << " "
					<< std::setw(10) << std::setprecision(4) << tan << " "
					<< std::setw(10) << std::setprecision(1) << L << " "
					<< std::setw(10) << std::setprecision(4) << dang << " "
					<< std::setw(10) << std::setprecision(4) << tan * pow(2, 1 / 2) * pow(0.45 * 0.45 + tan * tan * 16, 1 / 2) / 210 << " "
					<< std::setw(10) << std::setprecision(1) << L_scat << " "
					<< std::setw(10) << std::setprecision(1) << L_meas << " "
					<< std::setw(10) << std::setprecision(1) << md_exp << " "
					<< std::endl;
				pbeta = pbeta + 50;
			}
			dz += 50;
		}
	}

}
