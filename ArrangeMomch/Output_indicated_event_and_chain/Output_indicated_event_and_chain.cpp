// 2026/03/26
// flgを指定して任意のpartnerだけoutputする。

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
	int gid;
	double w;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.w) < std::tie(rhs.gid, rhs.w);
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


void SelectChains(std::vector<Momentum_recon::Event_information>& momch, int chaflg);
void SelectEvents(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, int matflg);
int main(int argc, char** argv) {
	if (argc != 5) {
		fprintf(stderr, "usage:prg in.momch output.momch material-flg charge-flg\n");
		exit(1);
	}
	std::string in_momch = argv[1];// input momch
	std::string out_momch = argv[2];// output momch
	int matflg = std::stoi(argv[3]);
	int chaflg = std::stoi(argv[4]);

	std::cout << "\n * Exclude vertex_material != " << matflg << "." << std::endl;
	std::cout << " * Exclude charge_sign     <  " << chaflg << ", except for muon." << std::endl;

	//bool result = std::filesystem::create_directories(out_momch);
	//assert(result);
	//assert(std::filesystem::exists(out_momch));//directryの存在確認
	//assert(std::filesystem::is_directory(out_momch));//指定されたパスがディレクトリを指しているかを確認する

	auto start = std::chrono::system_clock::now();//for measure working time


	// read momch
	std::vector<Momentum_recon::Event_information> momch0 = Momentum_recon::Read_Event_information_extension(in_momch);
	std::vector<Momentum_recon::Event_information> momch;
	SelectEvents(momch0, momch, matflg);
	// erasechain
	SelectChains(momch, chaflg);

	// write out
	Momentum_recon::Write_Event_information_extension(out_momch, momch);

	auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
	MeasureProcessingTime(start, end);
}

void SelectChains(std::vector<Momentum_recon::Event_information>& momch,int chaflg) {
	std::cout << "\n Select Chains" << std::endl;
	int count = 0;
	int remain = 0;


	for (auto& ev : momch) {
		for (auto c = ev.chains.begin(); c != ev.chains.end();) {
			// 条件一致した要素を削除する
			if (c->particle_flg != 13 && c->charge_sign < chaflg) {
				std::cout << "    erase : " << ev.groupid << ", " << c->chainid << " " << c->charge_sign << std::endl;;
				c = ev.chains.erase(c);
				count++;

			}
			// 要素削除をしない場合に、イテレータを進める
			else {
				c++;
				remain++;
			}
		}
	}
	std::cout << "Remain chain num : " << remain << std::endl;
	std::cout << "Erase  chain num : " << count << std::endl;
}
void SelectEvents(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, int matflg) {
	std::cout << "\n Select Events" << std::endl;

	int count = 0;
	int remain = 0;

	std::cout << " * Before erase : " << momch0.size() << std::endl;
	for (auto& ev : momch0) {
		if (ev.vertex_material != matflg) {// found erase event
			std::cout << "    erase : " << ev.groupid << ", " << ev.vertex_material << std::endl;
			count++;
		}
		else {//remain
			momch.push_back(ev);
			remain++;
		}
	}
	std::cout << " * After  erase : " << momch.size() << std::endl;

	std::cout << "Remain event num : " << remain << std::endl;
	std::cout << "Erase  event num : " << count << std::endl;
}
