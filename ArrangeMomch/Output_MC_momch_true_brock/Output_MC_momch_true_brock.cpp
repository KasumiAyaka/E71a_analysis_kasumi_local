#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
# include <iostream>
# include <random>

#include <list>

struct tkey {
	int eid;
	double ip;
};
bool operator<(const tkey& lhs, const tkey& rhs) {
	return std::tie(lhs.eid, lhs.ip) < std::tie(rhs.eid, rhs.ip);
}

class stop_track {
public:
	int64_t chainid, groupid;
	int  nseg, npl, pl0, pl1, vph, rawid;
	double ax, ay, x, y, z;
	// ph-->pid
	int stoppl;
	int unixtime;
	double ip;
	int pid;
	double mom, rng;
	double mu_md, mu_dz;
	int stop_flg;
	double dl, dal, dr, dar;
	int d_pl;
	int vph2;//íºãﬂÇÃéüÇÃbasetrack
	double pb;
};
class track_pair {
public:
	int eventid;
	double x, y, z, md, oa;
	double dz;
	stop_track t[2];
};
class track_multi {
public:
	int eventid, pl;
	double x, y, z;
	std::vector< std::pair<int, stop_track>> trk;
	std::vector<track_pair>pair;
	int unixtime;
	double dz;

};

struct VE_flg {
	tkey k;
	int ecc;
	int ve1, ve2, ve3;
};

double minimum_distance_fixed(matrix_3D::vector_3D pos0, matrix_3D::vector_3D pos1, matrix_3D::vector_3D dir0, matrix_3D::vector_3D dir1, double z_range[2], double extra[2], double refz) {
	double extra0_distance, extra1_distance, delta;
	matrix_3D::vector_3D pos;
	pos.x = pos1.x - pos0.x;
	pos.y = pos1.y - pos0.y;
	pos.z = pos1.z - pos0.z;
	//ÇŸÇ⁄ïΩçsÇ»èÍçá
	if (opening_angle(dir0, dir1) < 0.0001) {
		extra0_distance = (pos1.z + pos0.z) / 2 - pos0.z;
		extra1_distance = (pos1.z + pos0.z) / 2 - pos1.z;
	}
	else {
		delta = dot(dir0, dir0) * dot(dir1, dir1) - pow(dot(dir0, dir1), 2.);
		extra0_distance = (+1 * dot(pos, dir0) * dot(dir1, dir1) - dot(dir0, dir1) * dot(pos, dir1)) / delta;
		extra1_distance = (-1 * dot(pos, dir1) * dot(dir0, dir0) + dot(dir0, dir1) * dot(pos, dir0)) / delta;
	}
	//range[0]:è¨,range[1]:ëÂ
	if (z_range[0] > z_range[1]) {
		double tmp_d = z_range[0];
		z_range[0] = z_range[1];
		z_range[1] = tmp_d;
	}

	matrix_3D::vector_3D extra0 = addition(pos0, const_multiple(dir0, extra0_distance));
	matrix_3D::vector_3D extra1 = addition(pos1, const_multiple(dir1, extra1_distance));

	if (extra0.z < refz + z_range[0] || extra1.z < refz + z_range[0]) {//2025/8/20 fixed
		extra0_distance = refz - pos0.z + z_range[0];
		extra1_distance = refz - pos1.z + z_range[0];
	}
	else if (extra0.z > refz + z_range[1] || extra1.z > refz + z_range[1]) {//2025/8/20 fixed
		extra0_distance = refz - pos0.z + z_range[1];
		extra1_distance = refz - pos1.z + z_range[1];
	}

	extra[0] = extra0_distance;
	extra[1] = extra1_distance;
	extra0 = addition(pos0, const_multiple(dir0, extra0_distance));
	extra1 = addition(pos1, const_multiple(dir1, extra1_distance));

	return distance(extra0, extra1);

}


void Output_mc_true_brock(std::vector<Momentum_recon::Event_information>& momch, std::string output);


int main(int argc, char** argv) {
	if (argc < 3) {
		fprintf(stderr, "usage:prg in-mfile.momch out.txt \n");
		exit(1);
	}
	std::string in_momch = argv[1];
	std::string output = argv[2];
	std::cout << argc <<  std::endl;

	//read momch
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);
	std::cout << "Finish reading momch. (" << momch.size() << ")" << std::endl;

	Output_mc_true_brock(momch, output);
	std::cout << " Fin." << std::endl;

}

void Output_mc_true_brock(std::vector<Momentum_recon::Event_information>& momch, std::string output) {
	stop_track stop_tmp;

	std::ofstream ofs(output);
	if (!ofs) {
		std::cerr << "\tFile Open error...\n" << output << std::endl;
	}

	for (auto& ev : momch) {
		ofs << std::right << std::fixed
			<< std::setw(12) << "groupid" << " "
			<< std::setw(12) << "unix_time" << " "
			<< std::setw(4) << "vpl" << " "
			<< std::setw(10) << "v_material" << " "
			<< std::setw(10) << "vx " << " "
			<< std::setw(10) << "vy" << " "
			<< std::setw(10) << "vz" << " "
			<< std::setw(5) << "recon" << " "
			<< std::setw(5) << "true" << " "
			<< std::setw(8) << "weight" << std::endl;

		ofs << std::right << std::fixed
			<< std::setw(12) << std::setprecision(0) << ev.groupid << " "
			<< std::setw(12) << std::setprecision(0) << ev.unix_time << " "
			<< std::setw(4) << std::setprecision(0) << ev.vertex_pl << " "
			<< std::setw(10) << std::setprecision(0) << ev.vertex_material << " "
			<< std::setw(10) << std::setprecision(1) << ev.vertex_position[0] << " "
			<< std::setw(10) << std::setprecision(1) << ev.vertex_position[1] << " "
			<< std::setw(10) << std::setprecision(1) << ev.vertex_position[2] << " "
			<< std::setw(5) << std::setprecision(0) << ev.chains.size() << " "
			<< std::setw(5) << std::setprecision(0) << ev.true_chains.size() << " "
			<< std::setw(8) << std::setprecision(1) << ev.weight << std::endl;

		ofs << std::right << std::fixed
			<< std::setw(12) << "chainid" << " "
			<< std::setw(12) << "particle_flg" << " "
			<< std::setw(12) << "mu rng" << " "
			<< std::setw(12) << "p rng" << " "
			<< std::setw(12) << "mu mcs" << " "
			<< std::setw(12) << "p mcs" << " "
			<< std::setw(12) << "bm rng" << " "
			<< std::setw(12) << "bm pi" << " "
			<< std::endl;
		for (auto& c : ev.true_chains) {

				ofs << std::right << std::fixed
					<< std::setw(12) << std::setprecision(0) << c.chainid << " "
					<< std::setw(12) << std::setprecision(0) << c.particle_flg << " "
					<< std::setw(12) << std::setprecision(1) << c.ecc_range_mom[0] << " "
					<< std::setw(12) << std::setprecision(1) << c.ecc_range_mom[1] << " "
					<< std::setw(12) << std::setprecision(1) << c.ecc_mcs_mom[0] << " "
					<< std::setw(12) << std::setprecision(1) << c.ecc_mcs_mom[1] << " "
					<< std::setw(12) << std::setprecision(1) << c.bm_range_mom << " "
					<< std::setw(12) << std::setprecision(1) << c.bm_curvature_mom << " "
					<< std::endl;

			}
	}


	//printf("\t* input fin.\n\t #of track : %d --> %d\n", cnt, use);

}
