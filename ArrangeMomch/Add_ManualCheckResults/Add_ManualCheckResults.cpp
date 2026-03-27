// 2026/01/21
// Add_ManualCheckResults
// kasumi

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


void add_vertex_material_flg(std::vector<Momentum_recon::Event_information>& momch, std::set<int>& list, int mat);
void Set_manual_check_result(std::string input, std::set<int> &list);
void add_vertex_material_flg_all(std::vector<Momentum_recon::Event_information>& momch, int mat);
void add_sandmuon_flg(std::vector<Momentum_recon::Event_information>& momch, std::set<int>& list);


int main(int argc, char** argv) {
	if (argc < 4||argc>6) {
		fprintf(stderr, "===============================================================================\n");
		fprintf(stderr, " usage:prg in.momch output.momch materialflg [in_list.txt] [0:add sandmuon flg]\n\n * materialflg\n");
		fprintf(stderr, "\t 0 : water\t\n 1 : carbon(base)\t\n 2 : iron\t\n 5 : emulsion\t\n 6 : envelope\n");
		fprintf(stderr, "\t-2 : sand muon(penetrate)\t -3 :sideout\n");
		fprintf(stderr, "===============================================================================\n");
		exit(1);
	}
	std::string in_momch = argv[1];// input momch
	std::string out_momch = argv[2];// output momch
	int material = std::stoi(argv[3]);// material

	std::string in_list = "";
	if (argc >4) {
		in_list = argv[4];
	}

	//bool result = std::filesystem::create_directories(out_momch);
	//assert(result);
	//assert(std::filesystem::exists(out_momch));//directryの存在確認
	//assert(std::filesystem::is_directory(out_momch));//指定されたパスがディレクトリを指しているかを確認する

	auto start = std::chrono::system_clock::now();//for measure working time


	// read momch
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);

	if (argc == 5) {
		// to add material-flg to listed event
		std::cout << " * add material flg ( " << material << " ) to indicated event" << std::endl;
		// read flg
		std::set<int>ml;
		Set_manual_check_result(in_list, ml);

		// add vertex material
		add_vertex_material_flg(momch, ml,material);

	}
	else if (argc == 6) {
		// for sand muon
		material = -2;
		std::cout << " * add material flg ( " << material << " ) to indicated event" << std::endl;
		// read flg
		std::set<int>ml;
		Set_manual_check_result(in_list, ml);

		// add vertex material
		add_sandmuon_flg(momch, ml);

	}
	else {
		// to add vertex-material-flg to all event
		std::cout << " * add material flg ( " << material << " ) to all event" << std::endl;
		// add vertex material
		add_vertex_material_flg_all(momch, material);

	}



	// write out
	Momentum_recon::Write_Event_information_extension(out_momch, momch);

	auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
	MeasureProcessingTime(start, end);
}

void Set_manual_check_result(std::string input, std::set<int>& list) {
	
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "Failed to set manual check results." << std::endl;
		exit(0);
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward

	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);
		list.insert(std::stoi(str_v[0]));// group id

	}
	std::cout << "list size = " << list.size() << std::endl;
}

void add_vertex_material_flg_all(std::vector<Momentum_recon::Event_information>& momch, int mat) {

	int vpl;
	for (auto& ev : momch) {
		ev.vertex_material = mat;
		// 0 : water, 1: carbon(base), 2: iron, 5: emulsion, 6: envelope, 
		// -2 : penetrate, -3 : side out
		std::cout << ev.groupid << " " << ev.vertex_material << std::endl;
	}
}
void add_vertex_material_flg(std::vector<Momentum_recon::Event_information>& momch, std::set<int>& list, int mat) {

	int vpl;
	for (auto& ev : momch) {
		// 0 : water, 1: carbon(base), 2: iron, 5: emulsion, 6: envelope, 
		// -2 : penetrate, -3 : side out
		if (list.count(ev.groupid) > 0) {
			ev.vertex_material = mat;
		}
		std::cout << ev.groupid << " " << ev.vertex_material << std::endl;

	}
}

void add_sandmuon_flg(std::vector<Momentum_recon::Event_information>& momch, std::set<int>& list) {

	int cnt=0;
	for (auto& ev : momch) {
		// 0 : water, 1: carbon(base), 2: iron, 5: emulsion, 6: envelope, 
		// -2 : penetrate, -3 : side out
		if (list.count(ev.unix_time) > 0) {
			ev.vertex_material = -2;
			cnt++;
		}
		std::cout << ev.groupid << " " << ev.unix_time << " " << ev.vertex_material << std::endl;

	}

	std::cout << " * add flg num : " << cnt << std::endl;
}
