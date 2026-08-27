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
struct Key {
	double x, y;
	int area;
};
std::vector<TrackEntry> Read_yaml(const std::string& filename);
void Search_ref_btrk(std::string filename, int pl, std::ofstream& ofs, int ecc, int zone, TrackEntry& trk);
static std::string trim(const std::string& s);
std::string InsertRefFolder(const std::string& saveTo);
int ExtractEventNumber(const std::string& saveTo);
void Search_basetrack(std::string ECC_path, int pl, std::string cm_path, Key k, int64_t m[2], TrackEntry &trk);
void Large_area2Small_area(std::string ECC_path, std::string cm_path, TrackEntry &trk);

int main(int argc, char** argv) {
	if (argc < 3) {
		//fprintf(stderr, "usage : prg_name [input m-file-bin] groupnum gmin gmax [output m-file-bin]\n");
		fprintf(stderr, "usage : prg_name #ECC [in-YAML-file-path] [output_dir_path] [aera,Basetracklist]\n");
		exit(1);
	}
	int eccnum = std::stoi(argv[1]);
	std::string in_dir_path = argv[2];
	std::string out_dir_path = argv[3];
	std::string output_log = argv[4];

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
			std::cout << "\tSucceed to create dirctory: " << dir_path << std::endl;
		}
		else {
			std::cout << "\tFailed to create dirctory" << std::endl;
		}
	}
	else {
		std::cout << "\tThis directory already exists.: " << dir_path << std::endl;
		exit;
	}

	std::ofstream ofs_log(output_log);
	std::filesystem::path file_path;
	for (int pl = 3; pl <= 133; pl++) {
		std::stringstream file_input,file_output;
		file_input << in_dir_path << "\\PL" << std::setw(3) << std::setfill('0') << pl << ".yaml";
		file_output << out_dir_path << "\\PL" << std::setw(3) << std::setfill('0') << pl << ".yaml";
		file_path = file_input.str();
		if (!std::filesystem::exists(file_path)) continue;

		printf("\n\n * PL%03d start\n", pl);
		std::ofstream ofs(file_output.str());
		// set track list
		std::vector<TrackEntry>tlist = Read_yaml(file_input.str());
		ofs_log << pl << std::endl;
		for (auto itr = tlist.begin(); itr != tlist.end();itr++) {
			if (itr->RawID == -1) {
				// prediction track
				ofs << "  - Plate: " << pl << std::endl;
				ofs << "    Track: [" << itr->TrackAx << "," << itr->TrackAy << "," << itr->TrackX << "," << itr->TrackY << "]" << std::endl;
				ofs << "    Zone: " << itr->Zone << std::endl;
				ofs << "    Surface: Both" << std::endl;
				ofs << "    SaveTo: " << itr->SaveTo << std::endl;
				
				ofs_log << itr->Zone << " " << ExtractEventNumber(itr->SaveTo) << " " << itr->RawID << " " << 0 << std::endl;
				// search reference track
				std::stringstream file_in_base;
				file_in_base << file_in_ECC.str() << "\\Area" << itr->Zone << "\\PL" << std::setw(3) << std::setfill('0') << pl << "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
				if (!std::filesystem::exists(file_in_base.str())) {
					std::cout << file_in_base.str() << " doesn't exist." << std::endl;
					continue;
				}
				Search_ref_btrk(file_in_base.str(), pl, ofs, eccnum, itr->Zone, *itr);
				ofs_log << itr->Zone << " " << ExtractEventNumber(itr->SaveTo) << " " << itr->RawID << " " << -1 << std::endl;

			}
			else {
				// real track
				ofs << "  - Plate: " << itr->Plate << std::endl;
				ofs << "    RawID: " << itr->RawID << std::endl;
				ofs << "    Surface: Both" << std::endl;
				ofs << "    SaveTo: " << itr->SaveTo << std::endl;
				int tmp = itr->RawID;
				// search reference track
				std::stringstream file_in_base,file_in_cm;
				file_in_base << file_in_ECC.str() << "\\Area0\\PL" << std::setw(3) << std::setfill('0') << pl << "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
				file_in_cm << file_in_ECC.str() << "\\Area0\\0\\align\\corrmap-abs.lst";
				if (!std::filesystem::exists(file_in_base.str())) {
					std::cout << file_in_base.str() << " doesn't exist." << std::endl;
					continue;
				}
				Large_area2Small_area(file_in_ECC.str(), file_in_cm.str(), *itr);
				ofs_log << itr->Zone << " " << ExtractEventNumber(itr->SaveTo) << " " << tmp << " " << itr->RawID << std::endl;// area,eventid,rawid,original rawid
				Search_ref_btrk(file_in_base.str(), pl, ofs, eccnum, itr->Zone, *itr);
				ofs_log << itr->Zone << " " << ExtractEventNumber(itr->SaveTo) << " " << itr->RawID << " " << -1 << std::endl;

				//printf("input file  [%s]\n", file_in_base.str().c_str());
				//printf("output file [%s]\n", file_output.str().c_str());

			}
		}
	}
}

void CheckArea(double x, double y,Key& ret) {
	if (x < 131000) {//area 1,3,6
		if (y < 88000) {// are1
			ret.area = 1;
			ret.x = x;
			ret.y = y;
		}
		else if (y < 171000) {
			ret.area = 3;
			ret.x = x;
			ret.y = y - 75000;
		}
		else {
			ret.area = 5;
			ret.x = x;
			ret.y = y - 155000;
		}
	}
	else {
		if (y < 88000) {// are2
			ret.area = 2;
			ret.x = x - 120000;
			ret.y = y;
		}
		else if (y < 171000) {
			ret.area = 4;
			ret.x = x - 120000;
			ret.y = y - 75000;
		}
		else {
			ret.area = 6;
			ret.x = x - 120000;
			ret.y = y - 155000;
		}
	}
}
void Large_area2Small_area(std::string ECC_path, std::string cm_path, TrackEntry& trk) {
	std::stringstream Marged;
	Marged << ECC_path << "\\Area0\\PL" << std::setw(3) << std::setfill('0') << trk.Plate << "\\b" << std::setw(3) << std::setfill('0') << trk.Plate << ".sel.cor.vxx";
	// set track
	int64_t m[2]; 
	Key k = { 0 };
	vxx::BvxxReader br;
	std::array<int, 2> index = { trk.RawID, trk.RawID + 1 };//1234<=rawid<5678であるようなものだけを読む。
	if (br.Begin(Marged.str(), trk.Plate, 0, vxx::opt::c = cm_path, vxx::opt::index = index))
	{
		vxx::HashEntry h;
		vxx::base_track_t b;
		while (br.NextHashEntry(h))
		{
			while (br.NextBaseTrack(b))
			{
				fprintf(stderr, "%d %d %7.4lf %7.4lf %9.1lf %9.1lf\n", b.pl, b.rawid, b.ax, b.ay, b.x, b.y);
				m[0] = b.m[0].rawid;
				m[1] = b.m[1].rawid;
				CheckArea(b.x, b.y, k);
			}
		}
		br.End();
	}

	std::cout << "\t=== Area" << k.area << " ===" << std::endl;
	Search_basetrack(ECC_path, trk.Plate, cm_path, k, m, trk);//set area & pos
	int tmp = k.area;
	if (trk.TrackAx==0 && trk.TrackAy==0) {
		for (int i = 1; i < 7; i++) {
			if (i == tmp)continue;
			std::cout << "\t=== Area" << i << " ===" << std::endl;
			k.area = i;
			Search_basetrack(ECC_path, trk.Plate, cm_path, k, m, trk);//set area & pos
			if (trk.TrackAx != 0 || trk.TrackAy != 0) 	break;
		}
	}
	if (trk.TrackAx == 0 && trk.TrackAy == 0) {
		std::cout << "\tBasetrack " << trk.RawID << "was not found in small area.\n\n" << std::endl;
		std::exit(EXIT_SUCCESS);
	}
}
void Search_basetrack(std::string ECC_path,int pl, std::string cm_path, Key k,int64_t m[2], TrackEntry& trk) {
	std::stringstream Origin;

	Origin << ECC_path << "\\Area" << k.area << "\\PL" << std::setw(3) << std::setfill('0') << trk.Plate << "\\b" << std::setw(3) << std::setfill('0') << trk.Plate << ".sel.cor.vxx";


	double d =50000;//um
	vxx::BvxxReader br;
	std::vector<vxx::CutArea> area;
	std::vector<vxx::base_track_t> base_sel;
	int refnum = 0;
		// target positionから+- d[um]の範囲から角度の立ったproton likeな飛跡を探索
		area.push_back(vxx::CutArea(k.x - d, k.x + d, k.y - d, k.y + d));
		//std::vector<vxx::base_track_t> base = br.ReadAll(Origin.str(), pl, k.area, vxx::opt::a = area);//pl,pos
		std::vector<vxx::base_track_t> base = br.ReadAll(Origin.str(), pl, k.area);//pl,pos
		std::vector<vxx::base_track_t> base_selected;

		for (auto& b : base) {

			if (b.m[0].rawid == m[0] && b.m[1].rawid == m[1]) {
				std::cout << "\tfound\n\n" << std::endl;
				trk.RawID = b.rawid;//small area
				trk.Zone = k.area;
				trk.TrackAx = b.ax;
				trk.TrackAy = b.ay;
				trk.TrackX = b.x;
				trk.TrackY = b.y;
			}
		}
}

static std::string trim(const std::string& s) {
	size_t b = s.find_first_not_of(" \t\r\n");
	if (b == std::string::npos) return "";
	size_t e = s.find_last_not_of(" \t\r\n");
	return s.substr(b, e - b + 1);
}

// "G:ECC6/IMG/Event01260/PL019" -> "G:ECC6/Ref/IMG/Event01260/PL019"
std::string InsertRefFolder(const std::string& saveTo) {
	size_t pos = saveTo.find('/');
	if (pos == std::string::npos) return saveTo;
	return saveTo.substr(0, pos + 1) + "Ref/" + saveTo.substr(pos + 1);
}

// "G:ECC6/IMG/Event01260/PL019" -> 1260
int ExtractEventNumber(const std::string& saveTo) {
	size_t pos = saveTo.find("Event");
	if (pos == std::string::npos) return -1;
	return std::stoi(saveTo.substr(pos + 5, 5));
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
void Search_ref_btrk(std::string filename, int pl, std::ofstream& ofs, int ecc, int zone, TrackEntry& trk) {

	double d = 1000;//um
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
			std::cout << "\tCannot find reference trk. ( " << d << " )" << std::endl;
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
	ofs << "    SaveTo: " << InsertRefFolder(trk.SaveTo) << std::endl;

	trk.RawID = base_sel[0].rawid;
	trk.Zone = zone;

}