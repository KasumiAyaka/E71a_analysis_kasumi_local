//　2026/08/24
// 各エリアから濃い飛跡を選んでくる
// kasumi
#pragma comment(lib, "VxxReader.lib")
#include "VxxReader.h"

#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
#include <set>
#include <random>
#include <filesystem>

void Read_btrk(std::string filename, int pl, std::ofstream& ofs, int ecc, int zon, int tnum);

int main(int argc, char** argv) {
	if (argc < 3) {
		//fprintf(stderr, "usage : prg_name [input m-file-bin] groupnum gmin gmax [output m-file-bin]\n");
		fprintf(stderr, "usage : prg_name [input m-file-bin] [output m-file-txt] [pickup num(default==3)]\n");
		exit(1);
	}
	int eccnum = std::stoi(argv[1]);
	//std::string file_in_ECC;// = argv[1];
	//uint64_t gnum = std::stoll(argv[2]);
	//uint64_t gmin = std::stoll(argv[3]);
	//uint64_t gmax = std::stoll(argv[4]);
	std::string file_out_path = argv[2];
	int tnum = 3;
	if (argc > 3)tnum = std::stoi(argv[3]);

	// input path
	std::stringstream file_in_ECC;
	if (eccnum < 4) {
		file_in_ECC << "T:\\NINJA\\E71a\\ECC" << eccnum;
	}
	else if (eccnum == 4) {
		file_in_ECC << "K:\\NINJA\\E71a\\ECC" << eccnum;
	}
	else if (eccnum < 7) {
		file_in_ECC << "I:\\NINJA\\E71a\\ECC" << eccnum;
	}
	else {
		file_in_ECC << "K:\\NINJA\\E71a\\ECC" << eccnum;
	}

	// output path
	std::filesystem::path dir_path = file_out_path;
	if (!std::filesystem::exists(dir_path)) {
		// 2. 存在しない場合は作成
		if (std::filesystem::create_directories(dir_path)) {
			std::cout << "ディレクトリを作成しました: " << dir_path << std::endl;
		}
		else {
			std::cout << "ディレクトリの作成に失敗しました。" << std::endl;
		}
	}
	else {
		std::cout << "ディレクトリは既に存在します: " << dir_path << std::endl;
	}


	for (int pl = 3; pl <= 133; pl++) {
		std::stringstream file_output;
		file_output << file_out_path << "\\PL" << std::setw(3) << std::setfill('0') << pl << ".yaml";
		std::ofstream ofs(file_output.str());
		printf("PL%03d start\n", pl);

		for (int area = 1; area <= 6; area++) {
			std::stringstream file_in_base;
			file_in_base << file_in_ECC.str() << "\\Area" << area << "\\PL" << std::setw(3) << std::setfill('0') << pl << "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
			if (!std::filesystem::exists(file_in_base.str())) {
				std::cout << file_in_base.str() << " doesn't exist." << std::endl;
				continue;
			}
			//printf("input file  [%s]\n", file_in_base.str().c_str());
			//printf("output file [%s]\n", file_output.str().c_str());

			Read_btrk(file_in_base.str(), pl, ofs, eccnum, area, tnum);
		}
	}
}

void Read_btrk(std::string filename, int pl, std::ofstream& ofs, int ecc, int zon, int tnum) {
	//3. オプションを与えたい場合。
	//dump_bvxxの--index、--a、--c、--phに相当するオプションをBvxxReaderでも使うことができる。
	//それぞれ次の型の引数を要求する。
	//index : const std::array<int, 2>&
	//a     : const std::vector<CutArea>&
	//c     : const std::string&
	//ph    : int
	vxx::BvxxReader br;
	std::vector<vxx::CutArea> area;
	area.push_back(vxx::CutArea(70000, 80000, 50000, 60000, -1.0, 1.0, -1.0, 1.0));//xmin, xmax, ymin, ymax, axmin, axmax, aymin, aymax
	std::vector<vxx::base_track_t> base = br.ReadAll(filename, pl, 0, vxx::opt::ph = 32, vxx::opt::a = area);
	std::vector<vxx::base_track_t> base_selected;
	const int NumberOfImager = 72;
	uint32_t ShotID;
	int view, imager;
	for (auto& b : base) {
		if (b.ax * b.ax + b.ay * b.ay < 0.7) {
			//fprintf(stderr, "%d %d %4d %7.4lf %7.4lf %9.1lf %9.1lf\n", b.pl, b.rawid, b.m[0].ph % 10000 + b.m[1].ph % 10000, b.ax, b.ay, b.x, b.y);
			base_selected.push_back(b);
		}
	}
	int n = tnum;
	// 要素が3つ未満の場合は選択できないのでチェック
	if (base_selected.size() < n) {
		std::cout << "条件を満たすtrackが" << tnum << "つ未満です。" << std::endl;
		n = base_selected.size();
	}

	// 3. 乱数エンジンを用意してシャッフル
	std::random_device seed;
	std::mt19937 engine(seed());
	std::shuffle(base_selected.begin(), base_selected.end(), engine);

	// 4. 先頭の3つを出力（取り出し）
	//std::cout << "Selected track" << std::endl;
	for (int i = 0; i < n; ++i) {
		//fprintf(stderr, "%d %d %4d %7.4lf %7.4lf %9.1lf %9.1lf\n", base_selected[i].pl, base_selected[i].rawid, base_selected[i].m[0].ph % 10000 + base_selected[i].m[1].ph % 10000, base_selected[i].ax, base_selected[i].ay, base_selected[i].x, base_selected[i].y);


		ofs << "  - Plate: " << pl << std::endl;
		ofs << "    RawID: " << base_selected[i].rawid << std::endl;
		ofs << "#    Zone: " << zon << std::endl;
		ofs << "    Surface: Both" << std::endl;
		ofs << "    SaveTo: \"K:/NINJA/E71a/ManualCheck/ECC" << ecc << "/IMG/Ref/PL" << std::setw(3) << std::setfill('0') << pl << "\"" << std::endl;

	}
}