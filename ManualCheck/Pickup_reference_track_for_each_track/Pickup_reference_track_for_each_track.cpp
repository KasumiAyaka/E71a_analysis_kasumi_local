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
#include <algorithm>
#include <sstream>
#include <cmath>
#include <iomanip>
#define NOMINMAX
#include <windows.h>


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

	// Large_area2Small_area でのrawidマッピング記録用
	bool Mapped = false;               // 小エリアでのbasetrackマッピングに成功したか
	int64_t MergedPh = -1;             // merged(Area0)ファイル上のbasetrackのph%10000
	double MergedAx = 0, MergedAy = 0, MergedX = 0, MergedY = 0; // mergedファイル上のax,ay,x,y（マッピング失敗時のフォールバック用）
	int64_t SmallAreaPh = -1;          // 小エリアで見つかったbasetrackのph%10000
};
struct Key {
	double x, y;
	int area;
};
// Area1~6のいずれにも一致するbasetrackが見つからなかった場合の記録用
// (プログラム最後にまとめて出力する)
struct SmallAreaNotFoundEntry {
	int pl;
	int event;
	long long rawidMerged;
	double ax, ay, x, y;
};
struct OutputBlock {
	int zone;
	std::string yamlText;
	std::string logText;
	std::string mapText;
	std::string checkText; // クロスチェック用ファイルに書き出すブロック(小エリア側の情報を含む)
};

// SafeReadAreaBvxxが子プロセスから受け取るbasetrackの必要最小限の情報。
struct SimpleBaseTrack {
	int64_t rawid = 0;
	double ax = 0, ay = 0, x = 0, y = 0;
	int64_t m0rawid = 0, m1rawid = 0;
	int m0ph = 0, m1ph = 0;
};
// 小エリアbvxxファイルの読み込みは、原因不明のクラッシュ(0xC0000409等)を起こすことがあるため、
// 自分自身を"--readarea"モードで子プロセスとして起動して読み込みを行い、
// 子プロセスがクラッシュしても親プロセスは影響を受けず、falseを返して呼び出し元にフォールバックさせる。
// filtered=falseの場合、xmin~ymaxは無視され、CutAreaフィルタなしで全件読み込む(従来の挙動と同じ)。
static bool SafeReadAreaBvxx(const std::string& filename, int pl, int area, bool filtered,
	double xmin, double xmax, double ymin, double ymax, std::vector<SimpleBaseTrack>& out);

std::vector<TrackEntry> Read_yaml(const std::string& filename);
// refSelectMode: 0 = 候補からランダムに1つ選ぶ(従来通り), 1 = basetrackを構成するmicrotrackのVPHの和(m[0].ph%10000+m[1].ph%10000)が最も大きいものを選ぶ
// checkOfs: クロスチェック用ファイルへの出力先(小エリア側のTrackを含むReferenceブロックを書き出す)
void Search_ref_btrk(std::string filename, int pl, std::ostream& ofs, std::ostream& checkOfs, int ecc, int zone, TrackEntry& trk, int refSelectMode);
static std::string trim(const std::string& s);
std::string InsertRefFolder(const std::string& saveTo);
int ExtractEventNumber(const std::string& saveTo);
void Search_basetrack(std::string ECC_path, int pl, std::string cm_path, Key k, int64_t m[2], TrackEntry &trk, int eccnum,int mode);
bool Search_nearest_basetrack(std::string ECC_path, int pl, TrackEntry& trk, Key k, int eccnum, vxx::base_track_t& best, double& bestDist);
void Large_area2Small_area(std::string ECC_path, std::string cm_path, TrackEntry &trk, int eccnum, std::vector<SmallAreaNotFoundEntry>& notFoundLog);

int main(int argc, char** argv) {
	// --readarea モード: 危険な小エリアbvxx読み込みを単発で行うワーカーとして起動された場合。
	// このプロセス自体がクラッシュしても、起動元の親プロセスには影響しない。
	if (argc >= 11 && std::string(argv[1]) == "--readarea") {
		std::string file = argv[2];
		int pl = std::stoi(argv[3]);
		int area = std::stoi(argv[4]);
		bool filtered = (std::stoi(argv[5]) != 0);
		double xmin = std::stod(argv[6]);
		double xmax = std::stod(argv[7]);
		double ymin = std::stod(argv[8]);
		double ymax = std::stod(argv[9]);
		std::string outFile = argv[10];

		vxx::BvxxReader br;
		std::vector<vxx::base_track_t> base;
		if (filtered) {
			std::vector<vxx::CutArea> cut;
			cut.push_back(vxx::CutArea(xmin, xmax, ymin, ymax));
			base = br.ReadAll(file, pl, area, vxx::opt::a = cut);
		}
		else {
			base = br.ReadAll(file, pl, area);
		}

		std::ofstream ofs(outFile, std::ios::binary);
		int64_t n = (int64_t)base.size();
		ofs.write((char*)&n, sizeof(n));
		for (auto& b : base) {
			ofs.write((char*)&b.rawid, sizeof(b.rawid));
			ofs.write((char*)&b.ax, sizeof(b.ax));
			ofs.write((char*)&b.ay, sizeof(b.ay));
			ofs.write((char*)&b.x, sizeof(b.x));
			ofs.write((char*)&b.y, sizeof(b.y));
			ofs.write((char*)&b.m[0].rawid, sizeof(b.m[0].rawid));
			ofs.write((char*)&b.m[1].rawid, sizeof(b.m[1].rawid));
			ofs.write((char*)&b.m[0].ph, sizeof(b.m[0].ph));
			ofs.write((char*)&b.m[1].ph, sizeof(b.m[1].ph));
		}
		ofs.close();
		return 0;
	}

	if (argc < 5) {
		//fprintf(stderr, "usage : prg_name [input m-file-bin] groupnum gmin gmax [output m-file-bin]\n");
		fprintf(stderr, "usage : prg_name #ECC [in-YAML-file-path] [output_dir_path] [output_log_filename:aera,Basetracklist] [refSelectMode(0=random,1=maxVPH)]\n");
		exit(1);
	}
	int eccnum = std::stoi(argv[1]);
	std::string in_dir_path = argv[2];
	std::string out_dir_path = argv[3];
	std::string output_log = argv[4];
	int refSelectMode = 0; // 0: ランダムに1つ選ぶ(既定), 1: VPHの和が最大のものを選ぶ
	if (argc >= 6) refSelectMode = std::stoi(argv[5]);

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
		exit(1);
	}

	//
	std::vector<std::pair<int,int>> errlog;
	// Area1~6のいずれにも一致するbasetrackが見つからなかった場合の記録
	std::vector<SmallAreaNotFoundEntry> notFoundInSmallAreaLog;


	std::ofstream ofs_log(output_log);
	// [PL eventid rawid1 zone rawid2 ax ay x y ph%10000] を出力する
	std::filesystem::path logPath(output_log);
	std::filesystem::path mapPath = logPath.parent_path() / (logPath.stem().string() + "_rawidmap" + logPath.extension().string());
	std::ofstream ofs_map(mapPath.string());
	// クロスチェック用ファイル(出力YAMLと同じディレクトリに1つ作成し、小エリア側のRawID/Track/Zoneを含めて出力する)
	std::filesystem::path checkPath = std::filesystem::path(out_dir_path) / "crosscheck.yaml";
	std::ofstream ofs_check(checkPath.string());
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
		std::vector<OutputBlock> blocks;
		for (auto itr = tlist.begin(); itr != tlist.end();itr++) {
			std::ostringstream yamlSs, logSs, mapSs, checkSs;
			if (itr->RawID == -1) {
				printf("\n   --- PREDICTION TRACK ---\n", pl);
				// prediction track
				long long origRawID = itr->RawID;
				yamlSs << "  - Plate: " << pl << std::endl;
				yamlSs << "    Track: [" << itr->TrackAx << "," << itr->TrackAy << "," << itr->TrackX << "," << itr->TrackY << "]" << std::endl;
				yamlSs << "    Zone: " << itr->Zone << std::endl;
				yamlSs << "    Surface: Both" << std::endl;
				yamlSs << "    SaveTo: " << itr->SaveTo << std::endl;
				// 外挿トラック自体には小エリア側の追加情報が無いため、クロスチェック用にも同じ内容を出力する
				checkSs << yamlSs.str();

				logSs << itr->Zone << " " << ExtractEventNumber(itr->SaveTo) << " " << itr->RawID << " " << 0 << std::endl;
				// prediction trackはph%10000を-1として出力する
				mapSs << pl << " " << ExtractEventNumber(itr->SaveTo) << " " << origRawID << " " << itr->Zone << " " << 0
					<< " " << itr->TrackAx << " " << itr->TrackAy << " " << itr->TrackX << " " << itr->TrackY << " " << -1 << std::endl;
				// search reference track
				std::stringstream file_in_base;
				if (eccnum == 4) {
					file_in_ECC.clear();
					file_in_ECC.str("");
					if (itr->Zone == 1) {
						file_in_ECC << "K:\\NINJA\\E71a\\ECC" << eccnum;
					}
					else {
						file_in_ECC << "I:\\NINJA\\E71a\\ECC" << eccnum;
					}
				}
				file_in_base << file_in_ECC.str() << "\\Area" << itr->Zone << "\\PL" << std::setw(3) << std::setfill('0') << pl << "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
				if (!std::filesystem::exists(file_in_base.str())) {
					std::cout << file_in_base.str() << " doesn't exist." << std::endl;
					blocks.push_back({ itr->Zone, yamlSs.str(), logSs.str(), mapSs.str(), checkSs.str() });
					continue;
				}
				Search_ref_btrk(file_in_base.str(), pl, yamlSs, checkSs, eccnum, itr->Zone, *itr, refSelectMode);
				logSs << itr->Zone << " " << ExtractEventNumber(itr->SaveTo) << " " << itr->RawID << " " << -1 << std::endl;

			}
			else {
				printf("\n   --- REAL TRACK ---\n", pl);
				if (eccnum == 4) {
					file_in_ECC.clear();
					file_in_ECC.str("");
					file_in_ECC << "K:\\NINJA\\E71a\\ECC" << eccnum;
				}
				// real track
				yamlSs << "  - Plate: " << itr->Plate << std::endl;
				yamlSs << "    RawID: " << itr->RawID << std::endl;
				int tmp = itr->RawID;
				// matching real track
				std::stringstream file_in_base,file_in_cm;
				file_in_base << file_in_ECC.str() << "\\Area0\\PL" << std::setw(3) << std::setfill('0') << pl << "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
				file_in_cm << file_in_ECC.str() << "\\Area0\\0\\align\\corrmap-abs.lst";
				if (!std::filesystem::exists(file_in_base.str())) {
					std::cout << file_in_base.str() << " doesn't exist." << std::endl;
					yamlSs << "    Surface: Both" << std::endl;
					yamlSs << "    SaveTo: " << itr->SaveTo << std::endl;
					// 小エリア側の情報が無いため、クロスチェック用にも同じ内容を出力する
					checkSs << yamlSs.str();
					blocks.push_back({ itr->Zone, yamlSs.str(), logSs.str(), mapSs.str(), checkSs.str() });
					continue;
				}
				Large_area2Small_area(file_in_ECC.str(), file_in_cm.str(), *itr, eccnum, notFoundInSmallAreaLog);
				yamlSs << "    Surface: Both" << std::endl;
				yamlSs << "    SaveTo: " << itr->SaveTo << std::endl;
				// 実在トラックの小エリア側RawID/Track/Zone(マッピング後の値)はクロスチェック用ファイルに出力する
				checkSs << "  - Plate: " << itr->Plate << std::endl;
				checkSs << "    RawID: " << tmp << std::endl;
				checkSs << "    RawID_small: " << itr->RawID << std::endl;
				checkSs << "    Track_small: [" << itr->TrackAx << "," << itr->TrackAy << "," << itr->TrackX << "," << itr->TrackY << "]" << std::endl;
				checkSs << "    Zone_small: " << itr->Zone << std::endl;
				checkSs << "    Surface: Both" << std::endl;
				checkSs << "    SaveTo: " << itr->SaveTo << std::endl;
				logSs << itr->Zone << " " << ExtractEventNumber(itr->SaveTo) << " " << tmp << " " << itr->RawID << std::endl;// area,eventid,rawid,original rawid
				if (itr->RawID == 10) {
					errlog.push_back(std::make_pair(pl,tmp));
				}
				// rawidマッピングに成功していればrawid1,rawid2のph%10000の和を、
				// 失敗していればmergedファイル上のax,ay,x,yと-10を出力する
				if (itr->Mapped) {
					long long phSum = itr->MergedPh + itr->SmallAreaPh;
					mapSs << pl << " " << ExtractEventNumber(itr->SaveTo) << " " << tmp << " " << itr->Zone << " " << itr->RawID
						<< " " << itr->TrackAx << " " << itr->TrackAy << " " << itr->TrackX << " " << itr->TrackY << " " << phSum << std::endl;
				}
				else {
					mapSs << pl << " " << ExtractEventNumber(itr->SaveTo) << " " << tmp << " " << itr->Zone << " " << itr->RawID
						<< " " << itr->MergedAx << " " << itr->MergedAy << " " << itr->MergedX << " " << itr->MergedY << " " << -10 << std::endl;
				}
																															 // search reference track
				// Referenceは実在トラックと同じ小エリア(itr->Zone、マッピング後の値)から探すため、
				// file_in_base(Area0のmergedファイル)ではなく、Area{itr->Zone}の小エリアファイルを使う
				std::stringstream file_in_ref_area;
				if (eccnum == 4) {
					if (itr->Zone == 1) file_in_ref_area << "K:\\NINJA\\E71a\\ECC" << eccnum;
					else file_in_ref_area << "I:\\NINJA\\E71a\\ECC" << eccnum;
				}
				else {
					file_in_ref_area << file_in_ECC.str();
				}
				file_in_ref_area << "\\Area" << itr->Zone << "\\PL" << std::setw(3) << std::setfill('0') << pl << "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
				if (!std::filesystem::exists(file_in_ref_area.str())) {
					std::cout << file_in_ref_area.str() << " doesn't exist." << std::endl;
				}
				else {
					Search_ref_btrk(file_in_ref_area.str(), pl, yamlSs, checkSs, eccnum, itr->Zone, *itr, refSelectMode);
				}
				logSs << itr->Zone << " " << ExtractEventNumber(itr->SaveTo) << " " << itr->RawID << " " << -1 << std::endl;

				//printf("input file  [%s]\n", file_in_base.str().c_str());
				//printf("output file [%s]\n", file_output.str().c_str());

			}
			blocks.push_back({ itr->Zone, yamlSs.str(), logSs.str(), mapSs.str(), checkSs.str() });
		}
		// zoneの値の昇順にソートしてから出力する
		std::stable_sort(blocks.begin(), blocks.end(), [](const OutputBlock& a, const OutputBlock& b) {
			return a.zone < b.zone;
			});
		for (const auto& blk : blocks) {
			ofs << blk.yamlText;
			ofs_log << blk.logText;
			ofs_map << blk.mapText;
			ofs_check << blk.checkText;
		}
	}

	if (errlog.size() > 0) {
		std::cout << "\n\n ERRORs:Cannot find this basetrack in a small area" << std::endl;
		std::cout << "PL , rawid(merged) = " << std::endl;
		for (auto itr = errlog.begin(); itr != errlog.end(); itr++) {
			std::cout << itr->first << ", " << itr->second << std::endl;
		}
	}

	if (notFoundInSmallAreaLog.size() > 0) {
		std::cout << "\n\n WARNINGs: Area1~6 has no basetrack matching the Area0 basetrack's microtracks" << std::endl;
		std::cout << "PL , event , rawid(merged) , ax , ay , x , y = " << std::endl;
		for (auto itr = notFoundInSmallAreaLog.begin(); itr != notFoundInSmallAreaLog.end(); itr++) {
			std::cout << itr->pl << ", " << itr->event << ", " << itr->rawidMerged
				<< ", " << itr->ax << ", " << itr->ay << ", " << itr->x << ", " << itr->y << std::endl;
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
// 見つかった小エリア側basetrackの位置(trk.TrackX,trk.TrackY)が、
// 目標位置(k.x,k.y)からSearch_ref_btrkの初期探索半径(d=1000um)±500um以内にあるかを確認する。
static void CheckMatchDistance(const TrackEntry& trk, const Key& k) {
	double dist = std::sqrt(std::pow(trk.TrackX - k.x, 2) + std::pow(trk.TrackY - k.y, 2));
	const double refD = 1000; // Search_ref_btrkの初期d(um)と同じ基準値
	const double tolerance = 500;
	if (dist > refD + tolerance) {
		std::cout << "\t[WARN] matched basetrack (rawid=" << trk.RawID << ") is far from target position: dist="
			<< dist << "um (target: d=" << refD << "+-" << tolerance << "um)\n" << std::endl;
	}
}

void Large_area2Small_area(std::string ECC_path, std::string cm_path, TrackEntry& trk, int eccnum, std::vector<SmallAreaNotFoundEntry>& notFoundLog) {
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
				trk.MergedAx = b.ax;
				trk.MergedAy = b.ay;
				trk.MergedX = b.x;
				trk.MergedY = b.y;
				trk.MergedPh = b.m[0].ph % 10000 + b.m[1].ph % 10000;
			}
		}
		br.End();
	}

	std::cout << "\t=== Area" << k.area << " ===" << std::endl;
	Search_basetrack(ECC_path, trk.Plate, cm_path, k, m, trk,eccnum,1);//set area & pos
	int tmp = k.area;
	if (trk.TrackAx==0 && trk.TrackAy==0) {
		for (int i = 1; i < 7; i++) {
			if (i == tmp)continue;
			std::cout << "\t=== Area" << i << " ===" << std::endl;
			k.area = i;
			Search_basetrack(ECC_path, trk.Plate, cm_path, k, m, trk, eccnum,-1);//set area & pos
			if (trk.TrackAx != 0 || trk.TrackAy != 0) 	break;
		}
	}
	if (trk.TrackAx != 0 || trk.TrackAy != 0) {
		CheckMatchDistance(trk, k);
	}
	if (trk.TrackAx == 0 && trk.TrackAy == 0) {
		std::cout << "\tBasetrack " << trk.RawID << "was not found in small area.\n\n" << std::endl;
		notFoundLog.push_back({ trk.Plate, ExtractEventNumber(trk.SaveTo), trk.RawID, trk.MergedAx, trk.MergedAy, trk.MergedX, trk.MergedY });

		// 完全一致するbasetrackが見つからなかった場合、
		// mergedファイルで見つけたbasetrackの位置に最も近いbasetrackを各Areaから探索して代わりに使う
		vxx::base_track_t nearest;
		double nearestDist = 0;
		bool found = false;
		int foundArea = -1;
		for (int i = 1; i < 7; i++) {
			vxx::base_track_t candidate;
			double dist = 0;
			Key ki = k;
			ki.area = i;
			if (Search_nearest_basetrack(ECC_path, trk.Plate, trk, ki, eccnum, candidate, dist)) {
				if (!found || dist < nearestDist) {
					nearest = candidate;
					nearestDist = dist;
					found = true;
					foundArea = i;
				}
			}
		}

		if (found) {
			std::cout << "\tnearest basetrack found. rawid=" << nearest.rawid << " area=" << foundArea << " dist=" << nearestDist << " um" << std::endl;
			trk.RawID = nearest.rawid;
			trk.Zone = foundArea;
			trk.TrackAx = nearest.ax;
			trk.TrackAy = nearest.ay;
			trk.TrackX = nearest.x;
			trk.TrackY = nearest.y;
			trk.SmallAreaPh = nearest.m[0].ph % 10000 + nearest.m[1].ph % 10000;
			trk.Mapped = true;
			CheckMatchDistance(trk, k);
		}
		else {
			std::cout << "\tNearest basetrack search also failed.\n\n" << std::endl;
			//std::exit(EXIT_SUCCESS);
			trk.RawID = 10;//instead of small area rawid
		}
	}
}
static bool SafeReadAreaBvxx(const std::string& filename, int pl, int area, bool filtered,
	double xmin, double xmax, double ymin, double ymax, std::vector<SimpleBaseTrack>& out) {
	char exePath[MAX_PATH];
	if (GetModuleFileNameA(NULL, exePath, MAX_PATH) == 0) {
		std::cout << "\t[WARN] GetModuleFileNameA failed" << std::endl;
		return false;
	}

	char tempDir[MAX_PATH];
	GetTempPathA(MAX_PATH, tempDir);
	char tempFile[MAX_PATH];
	if (GetTempFileNameA(tempDir, "brk", 0, tempFile) == 0) {
		std::cout << "\t[WARN] GetTempFileNameA failed" << std::endl;
		return false;
	}

	std::ostringstream cmd;
	cmd << "\"" << exePath << "\" --readarea \"" << filename << "\" " << pl << " " << area
		<< " " << (filtered ? 1 : 0) << " " << std::setprecision(17)
		<< xmin << " " << xmax << " " << ymin << " " << ymax << " \"" << tempFile << "\"";

	STARTUPINFOA si; ZeroMemory(&si, sizeof(si)); si.cb = sizeof(si);
	PROCESS_INFORMATION pi; ZeroMemory(&pi, sizeof(pi));
	std::string cmdStr = cmd.str();
	std::vector<char> cmdBuf(cmdStr.begin(), cmdStr.end());
	cmdBuf.push_back('\0');

	if (!CreateProcessA(NULL, cmdBuf.data(), NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi)) {
		std::cout << "\t[WARN] failed to launch worker process for " << filename << std::endl;
		DeleteFileA(tempFile);
		return false;
	}
	WaitForSingleObject(pi.hProcess, INFINITE);
	DWORD exitCode = 1;
	GetExitCodeProcess(pi.hProcess, &exitCode);
	CloseHandle(pi.hProcess);
	CloseHandle(pi.hThread);

	if (exitCode != 0) {
		std::cout << "\t[WARN] worker process crashed/failed (exit=0x" << std::hex << exitCode << std::dec
			<< ") while reading " << filename << " -- treating as not found" << std::endl;
		DeleteFileA(tempFile);
		return false;
	}

	std::ifstream ifs(tempFile, std::ios::binary);
	if (!ifs) {
		DeleteFileA(tempFile);
		return false;
	}
	int64_t n = 0;
	ifs.read((char*)&n, sizeof(n));
	out.clear();
	out.resize((size_t)n);
	for (int64_t i = 0; i < n; i++) {
		ifs.read((char*)&out[i].rawid, sizeof(out[i].rawid));
		ifs.read((char*)&out[i].ax, sizeof(out[i].ax));
		ifs.read((char*)&out[i].ay, sizeof(out[i].ay));
		ifs.read((char*)&out[i].x, sizeof(out[i].x));
		ifs.read((char*)&out[i].y, sizeof(out[i].y));
		ifs.read((char*)&out[i].m0rawid, sizeof(out[i].m0rawid));
		ifs.read((char*)&out[i].m1rawid, sizeof(out[i].m1rawid));
		ifs.read((char*)&out[i].m0ph, sizeof(out[i].m0ph));
		ifs.read((char*)&out[i].m1ph, sizeof(out[i].m1ph));
	}
	ifs.close();
	DeleteFileA(tempFile);
	return true;
}

bool Search_nearest_basetrack(std::string ECC_path, int pl, TrackEntry& trk, Key k, int eccnum, vxx::base_track_t& best, double& bestDist) {
	std::stringstream Origin, str_tmp;
	std::string ECC_Area = ECC_path;
	if (eccnum == 4) {
		if (k.area == 1) {
			str_tmp << "K:\\NINJA\\E71a\\ECC" << eccnum;
		}
		else {
			str_tmp << "I:\\NINJA\\E71a\\ECC" << eccnum;
		}
		ECC_Area = str_tmp.str();
	}

	Origin << ECC_Area << "\\Area" << k.area << "\\PL" << std::setw(3) << std::setfill('0') << trk.Plate << "\\b" << std::setw(3) << std::setfill('0') << trk.Plate << ".sel.cor.vxx";

	double d = 50000;//um : Search_basetrackと同じ探索範囲
	// フィルタ(CutArea)を使うと見つからないことがあったため、従来通りフィルタなしで全件読み込む。
	// ただし読み込み自体は子プロセスに分離し、クラッシュしても親プロセスは継続できるようにする。
	std::vector<SimpleBaseTrack> base;
	if (!SafeReadAreaBvxx(Origin.str(), pl, k.area, false, 0, 0, 0, 0, base)) {
		std::cout << "\t[WARN] could not read " << Origin.str() << " (worker process failed) -- treating as not found" << std::endl;
		return false;
	}

	bool found = false;
	for (auto& b : base) {
		double dx = b.x - k.x;
		double dy = b.y - k.y;
		if (std::abs(dx) > d || std::abs(dy) > d) continue;
		double dist = std::sqrt(dx * dx + dy * dy);
		if (!found || dist < bestDist) {
			best.rawid = b.rawid;
			best.ax = b.ax; best.ay = b.ay; best.x = b.x; best.y = b.y;
			best.m[0].rawid = b.m0rawid; best.m[1].rawid = b.m1rawid;
			best.m[0].ph = b.m0ph; best.m[1].ph = b.m1ph;
			bestDist = dist;
			found = true;
		}
	}
	return found;
}
void Search_basetrack(std::string ECC_path, int pl, std::string cm_path, Key k, int64_t m[2], TrackEntry& trk, int eccnum, int mode) {
	std::stringstream Origin, str_tmp;
	std::string ECC_Area = ECC_path;
	if (eccnum == 4) {
		str_tmp.clear();
		str_tmp.str("");
		if (k.area == 1) {
			str_tmp << "K:\\NINJA\\E71a\\ECC" << eccnum;
		}
		else {
			str_tmp << "I:\\NINJA\\E71a\\ECC" << eccnum;
		}
		ECC_Area = str_tmp.str();
	}
	Origin << ECC_Area << "\\Area" << k.area << "\\PL" << std::setw(3) << std::setfill('0') << trk.Plate << "\\b" << std::setw(3) << std::setfill('0') << trk.Plate << ".sel.cor.vxx";
	std::cout << Origin.str() << std::endl;

	double d = 50000;//um
	// 読み込み自体は子プロセスに分離し、クラッシュしても親プロセスは継続できるようにする。
	std::vector<SimpleBaseTrack> base;
	if (!SafeReadAreaBvxx(Origin.str(), pl, k.area, false, 0, 0, 0, 0, base)) {
		std::cout << "\t[WARN] could not read " << Origin.str() << " (worker process failed) -- treating as not found" << std::endl;
		return;
	}

	for (auto& b : base) {

		if (b.m0rawid == m[0] && b.m1rawid == m[1]) {
			std::cout << "\tFound!\n\n" << std::endl;
			trk.RawID = b.rawid;//small area
			trk.Zone = k.area;
			trk.TrackAx = b.ax;
			trk.TrackAy = b.ay;
			trk.TrackX = b.x;
			trk.TrackY = b.y;
			trk.SmallAreaPh = b.m0ph % 10000 + b.m1ph % 10000;
			trk.Mapped = true;
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
void Search_ref_btrk(std::string filename, int pl, std::ostream& ofs, std::ostream& checkOfs, int ecc, int zone, TrackEntry& trk, int refSelectMode) {

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
			std::cout << "\tFound. ( " << d << " )" << std::endl;

			base_sel = base_selected;
		}
	} while (refnum == 0);

	// 3. 候補から1つ選ぶ
	vxx::base_track_t chosen;
	if (refSelectMode == 1) {
		// basetrackを構成するmicrotrackのVPHの和(m[0].ph%10000+m[1].ph%10000)が最も大きいものを選ぶ
		long long bestVph = -1;
		for (auto& b : base_sel) {
			long long vph = (long long)(b.m[0].ph % 10000) + (long long)(b.m[1].ph % 10000);
			if (vph > bestVph) {
				bestVph = vph;
				chosen = b;
			}
		}
	}
	else {
		// 乱数エンジンを用意してシャッフルし、先頭の1つを取り出す(従来通り)
		std::random_device seed;
		std::mt19937 engine(seed());
		std::shuffle(base_sel.begin(), base_sel.end(), engine);
		chosen = base_sel[0];
	}

	// 4. 選んだ1つを出力(本体のYAMLには元の形式のまま出力し、Track_smallはクロスチェック用ファイルにのみ出力する)
	ofs << "  - Plate: " << pl << std::endl;
	ofs << "    RawID: " << chosen.rawid << std::endl;
	ofs << "    Zone: " << zone << std::endl;
	ofs << "    Surface: Both" << std::endl;
	ofs << "    SaveTo: " << InsertRefFolder(trk.SaveTo) << std::endl;

	checkOfs << "  - Plate: " << pl << std::endl;
	checkOfs << "    RawID: " << chosen.rawid << std::endl;
	checkOfs << "    Track_small: [" << chosen.ax << "," << chosen.ay << "," << chosen.x << "," << chosen.y << "]" << std::endl;
	checkOfs << "    Zone: " << zone << std::endl;
	checkOfs << "    Surface: Both" << std::endl;
	checkOfs << "    SaveTo: " << InsertRefFolder(trk.SaveTo) << std::endl;

	trk.RawID = chosen.rawid;
	trk.Zone = zone;

}