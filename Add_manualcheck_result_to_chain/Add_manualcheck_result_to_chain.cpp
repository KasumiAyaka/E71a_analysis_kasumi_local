// 2026/02/09
// kasumi
// Add_manualcheck_result_to_chain
// input momch, groupid,chainid list
// Add flg to chains which are judged to erase by manual check.
// Chain header:charge sign
// -2:erased by manual check, -3:erased by final check

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
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
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


void read_chain_list(std::string in, std::set<std::pair<int, int>>& set);
void add_chain_flg_to_charge_sign(std::vector<Momentum_recon::Event_information>& momch, std::set<std::pair<int, int>>& set, int flg);
void add_chain_flg_to_charge_sign_2(std::vector<Momentum_recon::Event_information>& momch, std::set<std::pair<int, int>>& set, int flg);

int main(int argc, char** argv) {
	if (argc < 5) {
		fprintf(stderr, "usage:prg in.momch in.txt flg output.momch [0]\n");
		fprintf(stderr, "=============================================================\n");
		fprintf(stderr, " * flg\n");
		fprintf(stderr, " -2:erased by manual check, -3:erased by final check\n");
		fprintf(stderr, "=============================================================\n");
		exit(1);
	}
	std::string in_momch = argv[1];// input momch
	std::string in_list = argv[2];
	int flgnum = std::stoi(argv[3]);
	std::string out_momch = argv[4];// output momch
	int mode = -1;
	if (argc == 6) {
		mode = 1;
	}
	//bool result = std::filesystem::create_directories(out_momch);
	//assert(result);
	//assert(std::filesystem::exists(out_momch));//directryの存在確認
	//assert(std::filesystem::is_directory(out_momch));//指定されたパスがディレクトリを指しているかを確認する

	auto start = std::chrono::system_clock::now();//for measure working time


	// read momch
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);
	
	// read erase chain list
	std::set<std::pair<int, int>> set;
	read_chain_list(in_list, set);
	
	// add flg
	if (mode < 0) {
		add_chain_flg_to_charge_sign(momch, set, flgnum);
	}
	else {
		add_chain_flg_to_charge_sign_2(momch, set, flgnum);

	}
	// write out
	Momentum_recon::Write_Event_information_extension(out_momch, momch);

	auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
	MeasureProcessingTime(start, end);
}

void read_chain_list(std::string in,std::set<std::pair<int,int>>&set) {

	std::ifstream ifs(in);
	if (!ifs) {
		std::cerr << " File open error ! " << std::endl;
		std::cerr << in << std::endl;
		return;
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;
	
	Lst l;
	while (std::getline(ifs, str)) {
		str_v = StringSplit_with_tab(str);
		if (str_v.size() != 2) continue;
		set.insert(std::make_pair(std::stoi(str_v[0]), std::stoi(str_v[1])));// groupid chainid
	}

	std::cout << "\t* The Number of listed chain is " << set.size() << std::endl;


}

void add_chain_flg_to_charge_sign(std::vector<Momentum_recon::Event_information>& momch,std::set<std::pair<int,int>>&set,int flg) {
	std::cout << "\t  Add flg = " << flg << " to listed chain." << std::endl;

	int cnt = 0;
	int all = 0;
	for (auto& ev : momch) {
		
		for (auto& c : ev.chains) {
			all++;
			if (c.particle_flg == 13) continue;
			auto itr = set.find(std::make_pair(ev.groupid, c.chainid));
			if (itr != set.end()) {//find
				c.charge_sign = flg;
				cnt++;
			}
		}
	}
		
	std::cout << "\t* The Number of add flg chain is " << cnt << "(/" << all << ")" << std::endl;

}

void add_chain_flg_to_charge_sign_2(std::vector<Momentum_recon::Event_information>& momch, std::set<std::pair<int, int>>& set, int flg) {
	std::cout << "\t  Add flg = " << flg << " to listed chain." << std::endl;

	int cnt = 0;
	int all = 0;
	for (auto& ev : momch) {

		for (auto& c : ev.chains) {
			all++;
			if (c.particle_flg == 13) continue;
			auto itr = set.find(std::make_pair(ev.groupid, c.chainid));
			if (itr != set.end()) {//find
				//remain
			}
			else {
				//erase
				c.charge_sign = flg;
				cnt++;

			}
		}
	}

	std::cout << "\t* The Number of add flg chain is " << cnt << "(/" << all << ")" << std::endl;

}
