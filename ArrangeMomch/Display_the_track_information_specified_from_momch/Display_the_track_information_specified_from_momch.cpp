// 2026/4/24
// kasumi
// Display_the_track_information_specified_from_momch
// To check track information

#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
# include <iostream>
# include <random>

#include <list>
#include <cassert>
#include <filesystem>


class stop_track {
public:
	int64_t chainid, groupid;
	int  nseg, npl, pl0, pl1, ph, rawid, pid;
	double ax, ay, x, y, z;
	// ph-->pid
	int stoppl;
	int unixtime;

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
	std::vector< std::pair<double, stop_track>> trk;
	std::vector<track_pair>pair;
	int unixtime;
	double dz;
};
struct vtx_point {
	double x, y, z;
};
struct DivideParam {
	int md, pnum, tnum;
};
struct Key {
	int gid, cid, pl;
};
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.gid, lhs.cid, lhs.pl) < std::tie(rhs.gid, rhs.cid, rhs.pl);
}
auto start = std::chrono::system_clock::now();//for measure working time

void MeasureProcessingTime(std::chrono::system_clock::time_point& start, std::chrono::system_clock::time_point& end) {
	auto dur = end - start;        // 要した時間を計算
	auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
	// 要した時間をミリ秒（1/1000秒）に変換して表示
	std::cout << msec << " milli sec \n";
	if (msec / 1000 < 60) {
		std::cout << msec / 1000 << "sec\n";
	}
	else if (msec / 1000 / 60 < 60) {
		std::cout << msec / 1000 / 60 << "min\n";
	}
	else if (msec / 1000 / 3600 < 24) {
		std::cout << msec / 1000 / 3600 << "h\n";
	}
	else {
		std::cout << (msec / 1000 / 3600) / 24 << "day" << (msec / 1000 / 3600) % 24 << "h\n";
	}
};

double minimum_distance_fixed(matrix_3D::vector_3D pos0, matrix_3D::vector_3D pos1, matrix_3D::vector_3D dir0, matrix_3D::vector_3D dir1, double z_range[2], double extra[2], double refz) {
	double extra0_distance, extra1_distance, delta;
	matrix_3D::vector_3D pos;
	pos.x = pos1.x - pos0.x;
	pos.y = pos1.y - pos0.y;
	pos.z = pos1.z - pos0.z;
	//ほぼ平行な場合
	if (opening_angle(dir0, dir1) < 0.0001) {
		extra0_distance = (pos1.z + pos0.z) / 2 - pos0.z;
		extra1_distance = (pos1.z + pos0.z) / 2 - pos1.z;
	}
	else {
		delta = dot(dir0, dir0) * dot(dir1, dir1) - pow(dot(dir0, dir1), 2.);
		extra0_distance = (+1 * dot(pos, dir0) * dot(dir1, dir1) - dot(dir0, dir1) * dot(pos, dir1)) / delta;
		extra1_distance = (-1 * dot(pos, dir1) * dot(dir0, dir0) + dot(dir0, dir1) * dot(pos, dir0)) / delta;
	}
	//range[0]:小,range[1]:大
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

void Display_recon_basetrack(std::vector<Momentum_recon::Event_information>& momch, Key k);
void Display_true_basetrack(std::vector<Momentum_recon::Event_information>& momch, Key k);
void Display_vertex_information(std::vector<Momentum_recon::Event_information>& momch, Key k);

int main(int argc, char** argv) {
	if (argc <2) {
		fprintf(stderr, "===============================================================================\n");
		fprintf(stderr, " usage:prg in.momch\n");
		fprintf(stderr, "===============================================================================\n");
		exit(1);
	}
	std::string in_momch = argv[1];// input momch
	//bool result = std::filesystem::create_directories(out_momch);
	//assert(result);
	//assert(std::filesystem::exists(out_momch));//directryの存在確認
	//assert(std::filesystem::is_directory(out_momch));//指定されたパスがディレクトリを指しているかを確認する

	auto start = std::chrono::system_clock::now();//for measure working time
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);

	Key k;
	int mode = 5;
	do {
		std::cout << "Do you want to display vertex information? Or do you want to display basetrack information?" << std::endl;
		std::cout << " vertex --> input 1\n basetrack --> input 2\n exit --> -1" << std::endl;
		std::cin >> mode;
		if (mode > 2)continue;

		k = { 0 };
		if (mode == 1) {// vertex
			std::cout << "Please input groupid" << std::endl;
			std::cin >> k.gid;
			std::cout << std::endl;

			Display_vertex_information(momch, k);

		}
		else if (mode == 2) {// btrk
			std::cout << " Which information should be displayed: true or recon?" << std::endl;
			std::cout << "   true  --> input 1\n   recon --> input 2\n   exit  --> -1" << std::endl;
			std::cin >> mode;
			if (mode < 0)break;
			if (mode > 2)continue;

			std::cout << "Please input groupid, chainid, pl. \nex) 912 0 72" << std::endl;
			std::cin >> k.gid >> k.cid >> k.pl;
			std::cout << std::endl;
			if (mode == 1) {
				Display_true_basetrack(momch, k);
			}
			else if (mode == 2) {
				Display_recon_basetrack(momch, k);
			}

		}
		std::cout << std::endl;

	} while (mode > 0);




	auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
	MeasureProcessingTime(start, end);
}


void Display_recon_basetrack(std::vector<Momentum_recon::Event_information>& momch, Key k) {

	int vpl;
	double mom;
	for (auto& ev : momch) {
		if(ev.groupid==k.gid){
			for (auto& c : ev.chains) {
				if (c.chainid == k.cid) {
					for (auto& b : c.base) {
						if (b.pl == k.pl) {
							mom = c.ecc_mcs_mom[0];
							if (c.particle_flg % 10000 == 2212) {
								mom = c.ecc_mcs_mom[1];
							}
							std::cout << std::left << std::fixed
								<< " groupid = " << std::setw(5) << ev.groupid << "\n"
								<< " vertex material = " << std::setw(3) << ev.vertex_material << " (0:wt, 1:base, 2:iron,5:emul, -2:pene,-5:dup)\n"
								<< " chainid = " << std::setw(3) << c.chainid << "\n"
								<< " pl      = " << std::setw(3) << b.pl << "\n"
								<< " ax      = " << std::setw(8) << std::setprecision(4) << b.ax << "\n"
								<< " ay      = " << std::setw(8) << std::setprecision(4) << b.ay << "\n"
								<< " x       = " << std::setw(8) << std::setprecision(1) << b.x << "\n"
								<< " y       = " << std::setw(8) << std::setprecision(1) << b.y << "\n"
								<< " z       = " << std::setw(8) << std::setprecision(1) << b.z << "\n"
								<< " VPH     = " << std::setw(4) << std::setprecision(0) << (b.m[0].ph + b.m[1].ph) % 10000 << "\n"
								<< " PID     = " << std::setw(4) << c.particle_flg % 10000 << "\n"
								<< " MCS     = " << std::setw(8) << std::setprecision(1) << mom << "\n"
								<< " muon pb = " << std::setw(8) << std::setprecision(1) << c.Get_muon_mcs_pb() << "\n"
								<< " bm rng  = " << std::setw(8) << std::setprecision(1) << c.bm_range_mom << "\n"
								<< " bm curv = " << std::setw(8) << std::setprecision(1) << c.bm_curvature_mom << "\n"
								<< " rawid  = " << std::setw(8) << std::setprecision(0) << b.rawid << "\n"
								<< std::endl;
						}
					}
				}
			}
		}
	}
}
void Display_true_basetrack(std::vector<Momentum_recon::Event_information>& momch, Key k) {

	double mom;
	for (auto& ev : momch) {
		if (ev.groupid == k.gid) {
			for (auto& c : ev.true_chains) {
				if (c.chainid == k.cid) {
					for (auto& b : c.base) {
						if (b.pl == k.pl) {
							mom = c.ecc_mcs_mom[0];
							if (c.particle_flg % 10000 == 2212) {
								mom = c.ecc_mcs_mom[1];
							}
							std::cout << std::left << std::fixed
								<< " groupid = " << std::setw(5) << ev.groupid << "\n"
								<< " vertex material = " << std::setw(3) << ev.vertex_material << " (0:wt, 1:base, 2:iron,5:emul, -2:pene,-5:dup)\n"
								<< " chainid = " << std::setw(3) << c.chainid << "\n"
								<< " pl      = " << std::setw(3) << b.pl << "\n"
								<< " ax      = " << std::setw(8) << std::setprecision(4) << b.ax << "\n"
								<< " ay      = " << std::setw(8) << std::setprecision(4) << b.ay << "\n"
								<< " x       = " << std::setw(8) << std::setprecision(1) << b.x << "\n"
								<< " y       = " << std::setw(8) << std::setprecision(1) << b.y << "\n"
								<< " z       = " << std::setw(8) << std::setprecision(1) << b.z << "\n"
								<< " VPH     = " << std::setw(4) << std::setprecision(0) << (b.m[0].ph + b.m[1].ph) % 10000 << "\n"
								<< " PID     = " << std::setw(4) << c.particle_flg % 10000 << "\n"
								<< " MCS     = " << std::setw(8) << std::setprecision(1) << mom << "\n"
								<< " muon pb = " << std::setw(8) << std::setprecision(1) << c.Get_muon_mcs_pb() << "\n"
								<< " bm rng  = " << std::setw(8) << std::setprecision(1) << c.bm_range_mom << "\n"
								<< " bm curv = " << std::setw(8) << std::setprecision(1) << c.bm_curvature_mom << "\n"
								<< " rawid   = " << std::setw(8) << std::setprecision(0) << b.rawid << "\n"
								<< std::endl;
						}
					}
				}
			}
		}
	}
}
void Display_vertex_information(std::vector<Momentum_recon::Event_information>& momch, Key k) {

	int vpl;
	for (auto& ev : momch) {
		if (ev.groupid == k.gid) {
			std::cout << std::left << std::fixed
				<< " groupid = "<< std::setw(5) << ev.groupid << "\n"
				<< " Timestamp = " << std::setw(3) << ev.unix_time << "\n"
				<< " Vertex material = " << std::setw(3) << ev.vertex_material << " (0:wt, 1:base, 2:iron,5:emul, -2:pene,-5:dup)\n"
				<< " vertex pl      = " << std::setw(3) << ev.vertex_pl << "\n"
				<< " # of true trk  = " << std::setw(3) << ev.true_chains.size() << "\n"
				<< " # of recon trk = " << std::setw(3) << ev.chains.size() << "\n"
				<< " ax of neutrino = " << std::setw(8) << std::setprecision(4) <<ev.nu_ax << "\n"
				<< " ay of neutrino = " << std::setw(8) << std::setprecision(4) << ev.nu_ay << "\n"
				<< " true  vx = " << std::setw(8) << std::setprecision(1) << ev.true_vertex_position[0] << "\n"
				<< " true  vy = " << std::setw(8) << std::setprecision(1) << ev.true_vertex_position[1] << "\n"
				<< " true  vz = " << std::setw(8) << std::setprecision(1) << ev.true_vertex_position[2] << "\n"
				<< " recon vx = " << std::setw(8) << std::setprecision(1) << ev.vertex_position[0] << "\n"
				<< " recon vy = " << std::setw(8) << std::setprecision(1) << ev.vertex_position[1] << "\n"
				<< " recon vz = " << std::setw(8) << std::setprecision(1) << ev.vertex_position[2] << "\n"
				<< " weight   = " << std::setw(8) << std::setprecision(4) <<ev.weight << "\n"
				<< std::endl;
		}
	}
}

