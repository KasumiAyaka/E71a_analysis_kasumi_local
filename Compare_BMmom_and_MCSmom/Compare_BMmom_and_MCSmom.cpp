// 2026/02/17
// kasumi
// Compare_BMmom_and_MCSmom
// Input : Momch, output : txt
// OutputFileFormat :
// timestamp bunch eventid BMmom mcserr- + MCSmom mcserr- + nseg npl stoppl ax ay x y charge

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

struct Lst {
	int timestamp;
	double BMmom, charge;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.timestamp, lhs.BMmom, lhs.charge) < std::tie(rhs.timestamp, rhs.BMmom, rhs.charge);
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

void Output(std::vector<Momentum_recon::Event_information>& momch, std::string output);


int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "usage:prg in.momch output.txt\n");
		exit(1);
	}
	std::string in_momch = argv[1];// input momch
	std::string out_txt = argv[2];

	//bool result = std::filesystem::create_directories(out_momch);
	//assert(result);
	//assert(std::filesystem::exists(out_momch));//directryの存在確認
	//assert(std::filesystem::is_directory(out_momch));//指定されたパスがディレクトリを指しているかを確認する

	auto start = std::chrono::system_clock::now();//for measure working time


	// read momch
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);

	// read sharingfile
	Output(momch, out_txt);

	
	auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
	MeasureProcessingTime(start, end);
}

void Output(std::vector<Momentum_recon::Event_information>& momch, std::string output) {
	std::ofstream ofs(output);

	int cnt = 0;
	for (auto& ev : momch) {

		for (auto& c : ev.chains) {
			//if (c.particle_flg == 13) {
			if (c.chainid == 0) {
				ofs << std::right << std::fixed
					<< std::setw(12) << std::setprecision(0) << ev.unix_time << " "
					<< std::setw(1) << std::setprecision(0) << ev.weight << " "//bunch
					<< std::setw(5) << std::setprecision(0) << ev.groupid << " "
					<< std::setw(3) << std::setprecision(0) << ev.vertex_pl << " "
					//<< std::setw(3) << std::setprecision(0) << ev.chains.size() << " "//#trk
					<< std::setw(3) << std::setprecision(0) << c.base.size() << " "//nseg
					<< std::setw(3) << std::setprecision(0) << c.base.rbegin()->pl - c.base.begin()->pl + 1 << " "//npl
					<< std::setw(12) << std::setprecision(4) << c.base.rbegin()->ax << " "
					<< std::setw(12) << std::setprecision(4) << c.base.rbegin()->ay << " "
					<< std::setw(12) << std::setprecision(1) << c.base.rbegin()->x << " "
					<< std::setw(12) << std::setprecision(1) << c.base.rbegin()->y << " "
					<< std::setw(3) << std::setprecision(0) << c.charge_sign << " "
					<< std::setw(12) << std::setprecision(1) << c.bm_range_mom << " "
					<< std::setw(12) << std::setprecision(4) << c.bm_range_mom_error[0] << " "
					<< std::setw(12) << std::setprecision(4) << c.bm_range_mom_error[1] << " "
					<< std::setw(12) << std::setprecision(1) << c.ecc_mcs_mom[0] << " "
					<< std::setw(12) << std::setprecision(1) << c.ecc_mcs_mom_error[0][0] << " "
					<< std::setw(12) << std::setprecision(1) << c.ecc_mcs_mom_error[0][1] << " "
					<< std::setw(12) << std::setprecision(0) << ev.unix_time % 1500000000 << 0 << ev.weight
					<< std::endl;

					// timestamp bunch eventid BMmom mcserr- + MCSmom mcserr- + nseg npl stoppl ax ay x y charge

				cnt++;
			}
		}
	}
	std::cout << "\t* The Number of add info events " << cnt << std::endl;

}
