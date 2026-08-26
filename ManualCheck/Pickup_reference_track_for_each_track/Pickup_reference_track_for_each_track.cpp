// 2026/08/25
// �}�j���A���`�F�b�N�����Ղ��ƂɃ��t�@�����X�g���b�N����{�I��ł���
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


// TrackList/PLxxx.yaml の各エントリ。
// RawID が既知の場合: "- Plate / RawID / Surface / SaveTo"
// RawID が未確定で位置から探索する場合: "- Plate / Track: [ax,ay,x,y] / Zone / Surface / SaveTo"
struct TrackEntry {
	int Plate = 0;
	long long RawID = -1;   // -1 : 未指定（Track から探索する）
	bool HasTrack = false;
	double TrackAx = 0, TrackAy = 0, TrackX = 0, TrackY = 0;
	int Zone = -1;          // -1 : 未指定
	std::string Surface;
	std::string SaveTo;
};
std::vector<TrackEntry> Read_yaml(const std::string& filename);
void Read_btrk(std::string filename, std::string cm_path, int pl, std::ofstream& ofs, int ecc, TrackEntry trk);
void Search_ref_btrk(std::string filename, int pl, std::ofstream& ofs, int ecc, int zone, TrackEntry trk);
static std::string trim(const std::string& s);

int main(int argc, char** argv) {
	if (argc < 3) {
		//fprintf(stderr, "usage : prg_name [input m-file-bin] groupnum gmin gmax [output m-file-bin]\n");
		fprintf(stderr, "usage : prg_name #ECC [in-YAML-file-path] [output_dir_path] [pickup num(default==3)]\n");
		exit(1);
	}
	int eccnum = std::stoi(argv[1]);
	std::string in_dir_path = argv[2];
	std::string out_dir_path = argv[3];

	// ECC path
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
	std::filesystem::path dir_path = out_dir_path;
	if (!std::filesystem::exists(dir_path)) {
		// 2. ���݂��Ȃ��ꍇ�͍쐬
		if (std::filesystem::create_directories(dir_path)) {
			std::cout << "Succeed to create dirctory: " << dir_path << std::endl;
		}
		else {
			std::cout << "Failed to create dirctory" << std::endl;
		}
	}
	else {
		std::cout << "This directory already exists.: " << dir_path << std::endl;
		exit;
	}

	std::filesystem::path file_path;
	for (int pl = 3; pl <= 133; pl++) {
		std::stringstream file_input,file_output;
		file_input << in_dir_path << "\\PL" << std::setw(3) << std::setfill('0') << pl << ".yaml";
		file_output << out_dir_path << "\\PL" << std::setw(3) << std::setfill('0') << pl << ".yaml";
		file_path = file_input.str();
		if (!std::filesystem::exists(file_path)) continue;

		printf("PL%03d start\n", pl);
		std::ofstream ofs(file_output.str());
		// set track list
		std::vector<TrackEntry>tlist = Read_yaml(file_input.str());

		for (auto itr = tlist.begin(); itr != tlist.end();itr++) {
			if (itr->RawID == -1) {// prediction track
				std::stringstream file_in_base;
				file_in_base << file_in_ECC.str() << "\\Area" << itr->Zone << "\\PL" << std::setw(3) << std::setfill('0') << pl << "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
				if (!std::filesystem::exists(file_in_base.str())) {
					std::cout << file_in_base.str() << " doesn't exist." << std::endl;
					continue;
				}
				//printf("input file  [%s]\n", file_in_base.str().c_str());
				//printf("output file [%s]\n", file_output.str().c_str());

				Search_ref_btrk(file_in_base.str(), pl, ofs, eccnum, itr->Zone, *itr);

			}
			else {// real track
				std::stringstream file_in_base,file_in_cm;
				file_in_base << file_in_ECC.str() << "\\Area0\\PL" << std::setw(3) << std::setfill('0') << pl << "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
				file_in_cm << file_in_ECC.str() << "\\Area0\\0\align\corrmap-abs.lst";
				if (!std::filesystem::exists(file_in_base.str())) {
					std::cout << file_in_base.str() << " doesn't exist." << std::endl;
					continue;
				}
				//printf("input file  [%s]\n", file_in_base.str().c_str());
				//printf("output file [%s]\n", file_output.str().c_str());

				Read_btrk(file_in_base.str(),file_in_cm.str(), pl, ofs, eccnum,*itr);

			}
		}
	}
}

int CheckArea(double x, double y) {
	int area;
	if (x < 131000) {//area 1,3,6
		if (y < 88000) {// are1
			area = 1;
		}
		else if (y < 171000) {
			area = 3;
		}
		else {
			area = 5;
		}
	}
	else {
		if (y < 88000) {// are1
			area = 2;
		}
		else if (y < 171000) {
			area = 4;
		}
		else {
			area = 6;
		}
	}
	return area;
}
void Read_btrk(std::string filename, std::string cm_path,int pl, std::ofstream& ofs, int ecc, TrackEntry trk) {

	// set track
	int zon;
	vxx::BvxxReader br;
	std::string cm = "path to corrmap-abs";
	std::array<int, 2> index = { trk.RawID,trk.RawID + 1 };//1234<=rawid<5678であるようなものだけを読む。
	if (br.Begin(filename, pl, 0, vxx::opt::c = cm, vxx::opt::index = index))
	{
		vxx::HashEntry h;
		vxx::base_track_t b;
		while (br.NextHashEntry(h))
		{
			while (br.NextBaseTrack(b))
			{
				fprintf(stderr, "%d %d %7.4lf %7.4lf %9.1lf %9.1lf\n", b.pl, b.rawid, b.ax, b.ay, b.x, b.y);
				trk.TrackX = b.x;
				trk.TrackY = b.y;
			}
		}
		br.End();
	}

	ofs << "  - Plate: " << trk.Plate << std::endl;
	ofs << "    RawID: " << trk.RawID << std::endl;
	ofs << "    Surface: Both" << std::endl;
	ofs << "    SaveTo: " << trk.SaveTo << std::endl;

	// search ref
	vxx::BvxxReader br2;
	std::vector<vxx::CutArea> area;
	std::vector<vxx::base_track_t> base_sel;
	int refnum = 0;
	double d = 5000;
	do {
		// target positionから+- d[um]の範囲から角度の立ったproton likeな飛跡を探索
		area.push_back(vxx::CutArea(trk.TrackX - d, trk.TrackX + d, trk.TrackY - d, trk.TrackY + d, -1.0, 1.0, -1.0, 1.0));//xmin, xmax, ymin, ymax, axmin, axmax, aymin, aymax
		std::vector<vxx::base_track_t> base = br.ReadAll(filename, pl, 0, vxx::opt::ph = 32, vxx::opt::a = area);//pl,pos
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

		refnum = base_selected.size();
		// 要素が0 --> 探索領域を広げて再探索
		if (refnum < 1) {
			std::cout << "Cannot find reference trk. ( " << d << " )" << std::endl;
			d = d + 2000;
		}
		else {
			base_sel = base_selected;
		}
	} while (refnum == 0);

	// 3. 乱数エンジンを用意してシャッフル
	std::random_device seed;
	std::mt19937 engine(seed());
	std::shuffle(base_sel.begin(), base_sel.end(), engine);

	// 4. 先頭の1つを出力（取り出し）
	ofs << "  - Plate: " << pl << std::endl;
	ofs << "    RawID: " << base_sel[0].rawid << std::endl;
	//ofs << "    #Zone: " << zone << std::endl;
	ofs << "    Surface: Both" << std::endl;
	ofs << "    SaveTo: \"G:/ECC" << ecc << "/IMG/Ref/PL" << std::setw(3) << std::setfill('0') << pl << "\\raw" << std::setw(5) << std::setfill('0') << base_sel[0].rawid << "\"" << std::endl;
}

static std::string trim(const std::string& s) {
	size_t b = s.find_first_not_of(" \t\r\n");
	if (b == std::string::npos) return "";
	size_t e = s.find_last_not_of(" \t\r\n");
	return s.substr(b, e - b + 1);
}

// �{�v���O�������o�͂��邠�PLxxx.yaml�i- Plate / RawID / Surface / SaveTo �̃��X�g�j��ǂݍ���
std::vector<TrackEntry> Read_yaml(const std::string& filename) {
	std::vector<TrackEntry> entries;
	std::ifstream ifs(filename);
	if (!ifs) {
		std::cout << filename << " could not be opened." << std::endl;
		return entries;
	}

	std::string line;
	while (std::getline(ifs, line)) {
		std::string trimmed = trim(line);
		if (trimmed.empty() || trimmed[0] == '#') continue;

		bool isNewEntry = (trimmed.rfind("- ", 0) == 0);
		if (isNewEntry) trimmed = trim(trimmed.substr(2));

		size_t colon = trimmed.find(':');
		if (colon == std::string::npos) continue;
		std::string key = trim(trimmed.substr(0, colon));
		std::string value = trim(trimmed.substr(colon + 1));
		if (value.size() >= 2 && value.front() == '"' && value.back() == '"') {
			value = value.substr(1, value.size() - 2);
		}

		if (isNewEntry) entries.emplace_back();
		if (entries.empty()) continue;

		TrackEntry& e = entries.back();
		if (key == "Plate") e.Plate = std::stoi(value);
		else if (key == "RawID") e.RawID = std::stoll(value);
		else if (key == "Zone") e.Zone = std::stoi(value);
		else if (key == "Surface") e.Surface = value;
		else if (key == "SaveTo") e.SaveTo = value;
		else if (key == "Track") {
			// value : "[ax,ay,x,y]"
			if (!value.empty() && value.front() == '[') value.erase(0, 1);
			if (!value.empty() && value.back() == ']') value.pop_back();
			std::vector<double> nums;
			std::stringstream ss(value);
			std::string tok;
			while (std::getline(ss, tok, ',')) {
				nums.push_back(std::stod(trim(tok)));
			}
			if (nums.size() == 4) {
				e.HasTrack = true;
				e.TrackAx = nums[0];
				e.TrackAy = nums[1];
				e.TrackX = nums[2];
				e.TrackY = nums[3];
			}
		}
	}
	return entries;
}
void Search_ref_btrk(std::string filename, int pl, std::ofstream& ofs, int ecc, int zone, TrackEntry trk) {
	ofs << "  - Plate: " << trk.Plate << std::endl;
	ofs << "    Track: [" << trk.TrackAx<<","<<trk.TrackAy<<","<<trk.TrackX<<","<<trk.TrackY<<"]" << std::endl;
	ofs << "    Zone: " << trk.Zone << std::endl;
	ofs << "    Surface: Both" << std::endl;
	ofs << "    SaveTo: " << trk.SaveTo << std::endl;

	double d = 5000;//um
	//3. オプションを与えたい場合。
	//dump_bvxxの--index、--a、--c、--phに相当するオプションをBvxxReaderでも使うことができる。
	//それぞれ次の型の引数を要求する。
	//index : const std::array<int, 2>&
	//a     : const std::vector<CutArea>&
	//c     : const std::string&
	//ph    : int
	vxx::BvxxReader br;
	std::vector<vxx::CutArea> area;
	std::vector<vxx::base_track_t> base_sel;
	int refnum = 0;
	do {
		// target positionから+- d[um]の範囲から角度の立ったproton likeな飛跡を探索
		area.push_back(vxx::CutArea(trk.TrackX - d, trk.TrackX + d, trk.TrackY - d, trk.TrackY + d, -1.0, 1.0, -1.0, 1.0));//xmin, xmax, ymin, ymax, axmin, axmax, aymin, aymax
		std::vector<vxx::base_track_t> base = br.ReadAll(filename, pl, 0, vxx::opt::ph = 32, vxx::opt::a = area);//pl,pos
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

		refnum = base_selected.size();
		// 要素が0 --> 探索領域を広げて再探索
		if (refnum < 1) {
			std::cout << "Cannot find reference trk. ( " << d << " )" << std::endl;
			d = d + 2000;
		}
		else {
			base_sel = base_selected;
		}
	} while (refnum == 0);

	// 3. 乱数エンジンを用意してシャッフル
	std::random_device seed;
	std::mt19937 engine(seed());
	std::shuffle(base_sel.begin(), base_sel.end(), engine);

	// 4. 先頭の1つを出力（取り出し）
	ofs << "  - Plate: " << pl << std::endl;
	ofs << "    RawID: " << base_sel[0].rawid << std::endl;
	ofs << "    Zone: " << zone << std::endl;
	ofs << "    Surface: Both" << std::endl;
	ofs << "    SaveTo: \"G:\ECC" << ecc << "/IMG/Ref/PL" << std::setw(3) << std::setfill('0') << pl << "\\raw"<< std::setw(5) << std::setfill('0') << base_sel[0].rawid <<"\"" << std::endl;

}