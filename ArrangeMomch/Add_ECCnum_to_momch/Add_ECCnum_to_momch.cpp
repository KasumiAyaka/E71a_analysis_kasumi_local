// 2026/04/18
// kasumi
// Add_ECCnum_to_momch.cpp
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


void Set_manual_check_result(std::string input, std::map<int, int>& list);
void add_eccnum(std::vector<Momentum_recon::Event_information>& momch, int eccnum);


int main(int argc, char** argv) {
	if (argc < 3 || argc>4) {
		fprintf(stderr, "===============================================================================\n");
		fprintf(stderr, " usage:prg in.momch #ECC [output.momch]\n\n");
		fprintf(stderr, " If you specify an output filename, a new output file will be created.\n If you do not specify a filename, the original file will be overwritten with the new ECC number.\n");
		fprintf(stderr, "===============================================================================\n");
		exit(1);
	}
	std::string in_momch = argv[1];// input momch
	int eccnum = std::stoi(argv[2]);
	std::string out_momch;
	if (argc > 3) {
		out_momch = argv[3];// output momch
	}
	else {
		out_momch = in_momch;
	}

	if (eccnum < 1) {
		std::cout << "Warning!\nPlease enter the appropriate ECC numbers." << std::endl;
		std::cin >> eccnum;
	}
	else if (eccnum >9) {
		std::cout << "Warning!\nIs this ECC number valid?\nNote that in the E71a experiment, ECC numbers range from 1 to 9." << std::endl;
		std::cout << "Please enter \"y\" if this is what you expected, or \"n\" if not." << std::endl;
		char a;
		std::cin >> a;
		if (a == 'n') {
			std::cout << "Please input correct ECC number." << std::endl;
			std::cin >> eccnum;
		}
		else if (a == 'y') {
			// no problem
		}
		else {
			std::cout << "?" << std::endl;
			exit(1);
		}
	}

									
	//bool result = std::filesystem::create_directories(out_momch);
	//assert(result);
	//assert(std::filesystem::exists(out_momch));//directryの存在確認
	//assert(std::filesystem::is_directory(out_momch));//指定されたパスがディレクトリを指しているかを確認する

	auto start = std::chrono::system_clock::now();//for measure working time


	// read momch
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);

	add_eccnum(momch, eccnum);

	// write out
	Momentum_recon::Write_Event_information_extension(out_momch, momch);

	auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
	MeasureProcessingTime(start, end);
}


void add_eccnum(std::vector<Momentum_recon::Event_information>& momch, int eccnum) {

	int cnt = 0;
	for (auto& ev : momch) {
		ev.ECC_id = eccnum;
		cnt++;
	}

	std::cout << " * ECC number ( " << eccnum << " ) have been assigned to " << cnt << " events. " << std::endl;
}
