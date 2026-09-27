// 2026/09/25
// Pickup_reference_track_for_each_track.cppが出力したTrackList(PLxxx.yaml)を読み、
// 実在(Real)/Referenceの各飛跡が下記3つのvxxファイルに共通して存在するか確認する。
//   - Area{zone}\PLxxx\bxxx.sel.cor.vxx   (小エリア側 basetrack)
//   - Area{zone}\PLxxx\f{PL}{surface}_thick_0.vxx (小エリア側 microtrack, surface=1,2)
//   - Area0\PLxxx\bxxx.sel.cor.vxx        (mergeエリア側 basetrack)
// basetrackを構成するmicrotrackのzone/rawidが両basetrackファイル間・fvxxファイルとで
// 一致していることを確認し、一致しない(=飛跡を辿れない)場合はPL,event,basetrackのrawid,
// area(zone),area(zone)のrawid,Trackを出力する。
//
// 入力YAML(Pickup_reference_track_for_each_track.cppが出力する各PLxxx.yaml)は、
// 実在: Plate / RawID(Area0側) / Surface / SaveTo
// 外挿: Plate / Track / Zone / Surface / SaveTo
// Reference: Plate / RawID(小エリア側) / Zone / Surface / SaveTo
// という素の形式で、小エリア側の情報は含まれていない。実在トラックがどの小エリア(1~6)に
// 属するか、その小エリア側のrawid/Trackは何かは、Pickup_reference_track_for_each_track.cppの
// 出力に依存せず、このプログラム自身がArea0とArea1~6のbvxxを読み比べて独立に導出する。
// kasumi
#pragma comment(lib, "VxxReader.lib")
#include "VxxReader.h"

#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
#include <string>
#include <vector>
#include <map>
#include <array>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iostream>

// TrackListの1エントリ("- Plate: ..."で始まる1ブロック)から読み取った内容(素の形式のまま)。
struct Entry {
	int Plate = 0;
	std::string SaveTo;

	bool HasRawID = false; // 実在:Area0側rawid, Reference:小エリア側rawid
	long long RawID = -1;

	bool HasZone = false;  // 外挿/Referenceで使用(小エリアのzone)
	int Zone = -1;

	bool HasTrack = false; // 外挿で使用
	double TrackAx = 0, TrackAy = 0, TrackX = 0, TrackY = 0;
};

static std::string trim(const std::string& s) {
	size_t b = s.find_first_not_of(" \t\r\n");
	if (b == std::string::npos) return "";
	size_t e = s.find_last_not_of(" \t\r\n");
	return s.substr(b, e - b + 1);
}

// "G:ECC6/IMG/Event01260/PL019" (Referenceの場合は ".../Ref/IMG/...") -> 1260
static int ExtractEventNumber(const std::string& saveTo) {
	size_t pos = saveTo.find("Event");
	if (pos == std::string::npos) return -1;
	return std::stoi(saveTo.substr(pos + 5, 5));
}

// Pickup_reference_track_for_each_track.cppの出力YAML(PLxxx.yaml)を読む。
// "#"で始まる行は他のツールと同様にコメントとして無視する(素の形式のみを扱う)。
static std::vector<Entry> ReadOutputYaml(const std::string& filename) {
	std::vector<Entry> entries;
	std::ifstream ifs(filename);
	if (!ifs) {
		std::cout << filename << " could not be opened." << std::endl;
		return entries;
	}

	std::string line;
	while (std::getline(ifs, line)) {
		std::string body = trim(line);
		if (body.empty() || body[0] == '#') continue;

		bool isNewEntry = (body.rfind("- ", 0) == 0);
		if (isNewEntry) body = trim(body.substr(2));

		size_t colon = body.find(':');
		if (colon == std::string::npos) continue;
		std::string key = trim(body.substr(0, colon));
		std::string value = trim(body.substr(colon + 1));

		if (isNewEntry) entries.emplace_back();
		if (entries.empty()) continue;
		Entry& e = entries.back();

		if (key == "Plate") {
			e.Plate = std::stoi(value);
		}
		else if (key == "RawID") {
			e.HasRawID = true; e.RawID = std::stoll(value);
		}
		else if (key == "Zone") {
			e.HasZone = true; e.Zone = std::stoi(value);
		}
		else if (key == "Track") {
			std::string v = value;
			if (!v.empty() && v.front() == '[') v.erase(0, 1);
			if (!v.empty() && v.back() == ']') v.pop_back();
			std::vector<double> nums;
			std::stringstream ss(v);
			std::string tok;
			while (std::getline(ss, tok, ',')) nums.push_back(std::stod(trim(tok)));
			if (nums.size() == 4) {
				e.HasTrack = true;
				e.TrackAx = nums[0]; e.TrackAy = nums[1]; e.TrackX = nums[2]; e.TrackY = nums[3];
			}
		}
		else if (key == "SaveTo") {
			e.SaveTo = value;
		}
	}
	return entries;
}

// ECC4は小エリアがArea0,1=K:,Area2~6=I:に分かれて保存されている。
static std::string ECCAreaPath(int eccnum, int zone) {
	std::stringstream s;
	if (eccnum < 4) {
		s << "T:\\NINJA\\E71a\\ECC" << eccnum;
	}
	else if (eccnum == 4) {
		if (zone == 1 || zone == 0) s << "K:\\NINJA\\E71a\\ECC" << eccnum;
		else s << "I:\\NINJA\\E71a\\ECC" << eccnum;
	}
	else if (eccnum < 7) {
		s << "I:\\NINJA\\E71a\\ECC" << eccnum;
	}
	else {
		s << "K:\\NINJA\\E71a\\ECC" << eccnum;
	}
	return s.str();
}

// Area0(merge)側の座標から、小エリア(1~6)のうちどこに属するかを推定し、
// その小エリア内でのローカル座標も返す(Pickup_reference_track_for_each_track.cppのCheckAreaと同じ規則)。
static void GuessArea(double x, double y, int& area, double& lx, double& ly) {
	if (x < 131000) {
		if (y < 88000) { area = 1; lx = x; ly = y; }
		else if (y < 171000) { area = 3; lx = x; ly = y - 75000; }
		else { area = 5; lx = x; ly = y - 155000; }
	}
	else {
		if (y < 88000) { area = 2; lx = x - 120000; ly = y; }
		else if (y < 171000) { area = 4; lx = x - 120000; ly = y - 75000; }
		else { area = 6; lx = x - 120000; ly = y - 155000; }
	}
}

// 指定したareaに対する、Area0(merge)座標から小エリアのローカル座標へのオフセット。
// GuessAreaの範囲判定とは独立に、任意のareaについてオフセットだけを計算する
// (GuessAreaで当たらなかった場合に他areaも試すためのフォールバック用)。
static void AreaOffset(int area, double& dx, double& dy) {
	switch (area) {
	case 1: dx = 0; dy = 0; break;
	case 2: dx = -120000; dy = 0; break;
	case 3: dx = 0; dy = -75000; break;
	case 4: dx = -120000; dy = -75000; break;
	case 5: dx = 0; dy = -155000; break;
	case 6: dx = -120000; dy = -155000; break;
	default: dx = 0; dy = 0; break;
	}
}

// FvxxReaderのReadAll/Beginの第2引数はPL番号ではなく"POS"であり、
// POS = PL番号*10 + face(1 or 2。fvxxファイル名の"f{PL}{face}_thick_0.vxx"のfaceに対応)
// で計算する必要がある(VxxReaderライブラリのSample.cppとVxxReader.cppのMakePosより)。
static int MakeFvxxPos(int pl, int face) { return pl * 10 + face; }

// 小エリア側fvxxファイルのうちsurface=2側("f{PL}2_thick_0.vxx")のみを確認する。
// basetrackを構成する2つのmicrotrack(m[0],m[1])のうち、片方(おそらくm[1])の
// rawidがこのファイルに存在すれば良い(両方存在することは要求しない)。
static bool MicroTrackPairExistsInFvxxSurface2(int eccnum, int pl, int zone, const vxx::base_track_t& b) {
	std::stringstream f;
	f << ECCAreaPath(eccnum, zone) << "\\Area" << zone << "\\PL" << std::setw(3) << std::setfill('0') << pl
		<< "\\f" << std::setw(3) << std::setfill('0') << pl << "2_thick_0.vxx";
	if (!std::filesystem::exists(f.str())) {
		std::cout << "\t[DEBUG] fvxx file not found: " << f.str() << std::endl;
		return false;
	}

	int pos = MakeFvxxPos(pl, 2);
	vxx::FvxxReader fr;
	std::vector<vxx::micro_track_t> all = fr.ReadAll(f.str(), pos, zone);
	for (auto& mt : all) {
		if (mt.rawid == b.m[0].rawid || mt.rawid == b.m[1].rawid) return true;
	}

	std::cout << "\t[DEBUG] " << f.str() << " (pos=" << pos << ", zone=" << zone << ") : " << all.size() << " records, target rawid m[0]=" << b.m[0].rawid
		<< " m[1]=" << b.m[1].rawid;
	if (!all.empty()) {
		std::cout << " (sample: rawid=" << all[0].rawid << " zone=" << all[0].zone << ")";
	}
	std::cout << std::endl;
	return false;
}

static bool MicrotrackPairMatches(const vxx::base_track_t& a, const vxx::base_track_t& b) {
	bool matchDirect = (a.m[0].zone == b.m[0].zone && a.m[0].rawid == b.m[0].rawid && a.m[1].zone == b.m[1].zone && a.m[1].rawid == b.m[1].rawid);
	bool matchSwap = (a.m[0].zone == b.m[1].zone && a.m[0].rawid == b.m[1].rawid && a.m[1].zone == b.m[0].zone && a.m[1].rawid == b.m[0].rawid);
	return matchDirect || matchSwap;
}

// Area{area}\PL{pl}\b{pl}.sel.cor.vxx から、mergeBの持つmicrotrackペア(m[0],m[1])と
// 一致するbasetrackを探す(位置(lx,ly)の周辺±dのみを読み、負荷を抑える)。
static bool FindMatchingBasetrackInArea(int eccnum, int pl, int area, double lx, double ly,
	const vxx::base_track_t& mergeB, vxx::base_track_t& out) {
	std::stringstream f;
	f << ECCAreaPath(eccnum, area) << "\\Area" << area << "\\PL" << std::setw(3) << std::setfill('0') << pl
		<< "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
	if (!std::filesystem::exists(f.str())) return false;

	double d = 50000;//um
	vxx::BvxxReader br;
	std::vector<vxx::CutArea> cut;
	cut.push_back(vxx::CutArea(lx - d, lx + d, ly - d, ly + d));
	std::vector<vxx::base_track_t> base = br.ReadAll(f.str(), pl, area, vxx::opt::a = cut);
	for (auto& b : base) {
		if (MicrotrackPairMatches(b, mergeB)) { out = b; return true; }
	}
	return false;
}

struct VerifyResult {
	bool ok = false;
	long long MergedRawID = -1; // basetrackのrawid(Area0側)
	int Zone = -1;              // basetrackが属する小エリア(Area1~6)
	long long RawIDSmall = -1;  // 小エリア側のrawid
	double Ax = 0, Ay = 0, X = 0, Y = 0; // 小エリア側のTrack
	std::string reason;
};

// 実在トラック: Area0側のrawid(known)から出発し、Area0のbasetrackが持つmicrotrackペアと
// 一致するbasetrackをArea1~6から探して小エリア(zone)とrawidを独立に導出する。
// 見つかったbasetrackについて、fvxx(surface=2)にそのmicrotrackが存在するかも確認する。
static VerifyResult VerifyReal(int eccnum, int pl, long long mergedRawid) {
	VerifyResult r;

	std::stringstream area0File, cmFile;
	area0File << ECCAreaPath(eccnum, 0) << "\\Area0\\PL" << std::setw(3) << std::setfill('0') << pl
		<< "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
	cmFile << ECCAreaPath(eccnum, 0) << "\\Area0\\0\\align\\corrmap-abs.lst";
	if (!std::filesystem::exists(area0File.str())) {
		r.reason = "Area0 bvxx file not found";
		return r;
	}

	vxx::BvxxReader br0;
	std::array<int, 2> idx0 = { (int)mergedRawid, (int)mergedRawid + 1 };
	std::vector<vxx::base_track_t> res0 = br0.ReadAll(area0File.str(), pl, 0, vxx::opt::c = cmFile.str(), vxx::opt::index = idx0);
	vxx::base_track_t mergeB;
	bool foundMerged = false;
	for (auto& t : res0) {
		if (t.rawid == mergedRawid) { mergeB = t; foundMerged = true; break; }
	}
	if (!foundMerged) {
		r.reason = "basetrack not found in Area0 bvxx";
		return r;
	}
	r.MergedRawID = mergedRawid;

	// どの小エリアに属するかを推定し、まずそこから探す
	int guessArea; double lx, ly;
	GuessArea(mergeB.x, mergeB.y, guessArea, lx, ly);

	vxx::base_track_t smallB;
	bool found = FindMatchingBasetrackInArea(eccnum, pl, guessArea, lx, ly, mergeB, smallB);
	int foundArea = guessArea;
	if (!found) {
		for (int i = 1; i <= 6 && !found; i++) {
			if (i == guessArea) continue;
			double dx, dy;
			AreaOffset(i, dx, dy);
			double ax = mergeB.x + dx, ay = mergeB.y + dy;
			if (FindMatchingBasetrackInArea(eccnum, pl, i, ax, ay, mergeB, smallB)) {
				found = true;
				foundArea = i;
			}
		}
	}
	if (!found) {
		r.reason = "basetrack not found in any small area(1~6)";
		return r;
	}

	r.Zone = foundArea;
	r.RawIDSmall = smallB.rawid;
	r.Ax = smallB.ax; r.Ay = smallB.ay; r.X = smallB.x; r.Y = smallB.y;

	if (!MicroTrackPairExistsInFvxxSurface2(eccnum, pl, foundArea, smallB)) {
		r.reason = "microtrack not found in Area" + std::to_string(foundArea) + " fvxx(surface2)";
		return r;
	}

	r.ok = true;
	return r;
}

// Referenceトラック: 小エリア側のzone/rawidは既知(YAMLに直接記録されている)ので、
// 1) その小エリアbvxxにbasetrackが存在するか
// 2) そのmicrotrackがfvxx(surface=2)に存在するか
// 3) 同じmicrotrackペアを持つbasetrackがArea0(merge)bvxxに存在するか
// を確認する。Area0側のrawidは不明なので、PLごとに一度だけ全件読み込みキャッシュして検索する。
static VerifyResult VerifyReference(int eccnum, int pl, int zone, long long rawidSmall,
	std::map<int, std::vector<vxx::base_track_t>>& area0Cache) {
	VerifyResult r;
	r.Zone = zone;
	r.RawIDSmall = rawidSmall;

	std::stringstream smallFile;
	smallFile << ECCAreaPath(eccnum, zone) << "\\Area" << zone << "\\PL" << std::setw(3) << std::setfill('0') << pl
		<< "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
	if (!std::filesystem::exists(smallFile.str())) {
		r.reason = "Area" + std::to_string(zone) + " bvxx file not found";
		return r;
	}
	vxx::BvxxReader br;
	std::array<int, 2> idx = { (int)rawidSmall, (int)rawidSmall + 1 };
	std::vector<vxx::base_track_t> res = br.ReadAll(smallFile.str(), pl, zone, vxx::opt::index = idx);
	vxx::base_track_t b;
	bool found = false;
	for (auto& t : res) {
		if (t.rawid == rawidSmall) { b = t; found = true; break; }
	}
	if (!found) {
		r.reason = "basetrack not found in Area" + std::to_string(zone) + " bvxx";
		return r;
	}
	r.Ax = b.ax; r.Ay = b.ay; r.X = b.x; r.Y = b.y;

	if (!MicroTrackPairExistsInFvxxSurface2(eccnum, pl, zone, b)) {
		r.reason = "microtrack not found in Area" + std::to_string(zone) + " fvxx(surface2)";
		return r;
	}

	std::stringstream area0File, cmFile;
	area0File << ECCAreaPath(eccnum, 0) << "\\Area0\\PL" << std::setw(3) << std::setfill('0') << pl
		<< "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
	cmFile << ECCAreaPath(eccnum, 0) << "\\Area0\\0\\align\\corrmap-abs.lst";
	if (!std::filesystem::exists(area0File.str())) {
		r.reason = "Area0 bvxx file not found";
		return r;
	}
	if (area0Cache.find(pl) == area0Cache.end()) {
		vxx::BvxxReader br0;
		area0Cache[pl] = br0.ReadAll(area0File.str(), pl, 0, vxx::opt::c = cmFile.str());
	}
	bool foundMerged = false;
	for (auto& mb : area0Cache[pl]) {
		if (MicrotrackPairMatches(mb, b)) { r.MergedRawID = mb.rawid; foundMerged = true; break; }
	}
	if (!foundMerged) {
		r.reason = "no matching basetrack found in Area0 bvxx";
		return r;
	}

	r.ok = true;
	return r;
}

int main(int argc, char** argv) {
	if (argc < 4) {
		fprintf(stderr, "usage : prg_name #ECC [TrackList-YAML-dir(output of Pickup_reference_track_for_each_track)] [output-log-file-path]\n");
		exit(1);
	}
	int eccnum = std::stoi(argv[1]);
	std::string in_dir_path = argv[2];
	std::string output_log = argv[3];

	// 出力先ディレクトリが存在しない場合は作成する
	std::filesystem::path outLogPath(output_log);
	std::filesystem::path outLogDir = outLogPath.parent_path();
	if (!outLogDir.empty() && !std::filesystem::exists(outLogDir)) {
		if (std::filesystem::create_directories(outLogDir)) {
			std::cout << "\tSucceed to create directory: " << outLogDir << std::endl;
		}
		else {
			std::cout << "\tFailed to create directory: " << outLogDir << std::endl;
		}
	}

	std::ofstream ofs_log(output_log);
	ofs_log << "PL event basetrack_rawid(merged) area(zone) area_rawid ax ay x y type reason" << std::endl;

	std::map<int, std::vector<vxx::base_track_t>> area0Cache;
	int totalChecked = 0, totalFailed = 0;

	for (int pl = 3; pl <= 133; pl++) {
		std::stringstream file_in;
		file_in << in_dir_path << "\\PL" << std::setw(3) << std::setfill('0') << pl << ".yaml";
		if (!std::filesystem::exists(file_in.str())) continue;

		std::vector<Entry> entries = ReadOutputYaml(file_in.str());
		printf("\n * PL%03d : %zu entries\n", pl, entries.size());

		for (auto& e : entries) {
			bool isReference = (e.SaveTo.find("/Ref/") != std::string::npos);

			VerifyResult vr;
			std::string type;

			if (isReference) {
				if (!e.HasRawID || !e.HasZone) continue; // 情報が揃っていないものはスキップ
				type = "Reference";
				totalChecked++;
				vr = VerifyReference(eccnum, pl, e.Zone, e.RawID, area0Cache);
			}
			else if (e.HasRawID) {
				// 実在トラック(RawIDがあるもの。外挿はRawIDを持たないので対象外)
				type = "Real";
				totalChecked++;
				vr = VerifyReal(eccnum, pl, e.RawID);
			}
			else {
				continue; // 外挿トラックは対象外
			}

			if (!vr.ok) {
				totalFailed++;
				int event = ExtractEventNumber(e.SaveTo);
				std::cout << "\t[NG] PL" << pl << " event" << event << " (" << type << ") : " << vr.reason << std::endl;
				ofs_log << pl << " " << event << " " << vr.MergedRawID << " " << vr.Zone << " " << vr.RawIDSmall
					<< " " << vr.Ax << " " << vr.Ay << " " << vr.X << " " << vr.Y
					<< " " << type << " \"" << vr.reason << "\"" << std::endl;
			}
		}
	}

	std::cout << "\n\n === summary ===" << std::endl;
	std::cout << "checked : " << totalChecked << std::endl;
	std::cout << "failed  : " << totalFailed << std::endl;

	return 0;
}
