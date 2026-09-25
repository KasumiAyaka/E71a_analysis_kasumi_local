// 2026/09/25
// Pickup_reference_track_for_each_track.cppが出力したTrackList(PLxxx.yaml)を読み、
// 実在(Real)/Referenceの各飛跡が下記3つのvxxファイルに共通して存在するか確認する。
//   - Area{zone}\PLxxx\bxxx.sel.cor.vxx   (小エリア側 basetrack)
//   - Area{zone}\PLxxx\f{PL}{surface}_thick_0.vxx (小エリア側 microtrack, surface=1,2)
//   - Area0\PLxxx\bxxx.sel.cor.vxx        (mergeエリア側 basetrack)
// basetrackを構成するmicrotrackのzone/rawidが両basetrackファイル間・fvxxファイルとで
// 一致していることを確認し、一致しない(=飛跡を辿れない)場合はPL,event,basetrackのrawid,
// area(zone),area(zone)のrawid,Trackを出力する。
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

// TrackListの1エントリ("- Plate: ..."で始まる1ブロック)から読み取った内容。
// 出力YAMLは、コメントアウトされた行(# RawID / #Track / #Zone)にも意味のある値が
// 入っているため、コメントの有無を区別して両方保持する。
struct Entry {
	int Plate = 0;
	std::string SaveTo;

	bool HasRawID = false;          // "RawID:" (コメントなし)
	long long RawID = -1;
	bool HasRawIDComment = false;   // "#RawID:" (実在の小エリア側rawid)
	long long RawIDComment = -1;

	bool HasZone = false;           // "Zone:" (コメントなし。外挿/Referenceで使用)
	int Zone = -1;
	bool HasZoneComment = false;    // "#Zone:" (実在の小エリア側zone)
	int ZoneComment = -1;

	bool HasTrack = false;          // "Track:" (コメントなし。外挿で使用)
	double TrackAx = 0, TrackAy = 0, TrackX = 0, TrackY = 0;
	bool HasTrackComment = false;   // "#Track:" (実在/Referenceの小エリア側Track)
	double CTrackAx = 0, CTrackAy = 0, CTrackX = 0, CTrackY = 0;
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

// Pickup_reference_track_for_each_track.cppの出力YAMLを読む。
// 他のツールと異なり、"#"で始まるコメント行も(コメントであることを記録した上で)値として読み取る。
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
		if (body.empty()) continue;

		bool commented = false;
		if (body[0] == '#') {
			commented = true;
			body = trim(body.substr(1));
			if (body.empty()) continue;
		}

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
			if (commented) { e.HasRawIDComment = true; e.RawIDComment = std::stoll(value); }
			else { e.HasRawID = true; e.RawID = std::stoll(value); }
		}
		else if (key == "Zone") {
			if (commented) { e.HasZoneComment = true; e.ZoneComment = std::stoi(value); }
			else { e.HasZone = true; e.Zone = std::stoi(value); }
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
				if (commented) {
					e.HasTrackComment = true;
					e.CTrackAx = nums[0]; e.CTrackAy = nums[1]; e.CTrackX = nums[2]; e.CTrackY = nums[3];
				}
				else {
					e.HasTrack = true;
					e.TrackAx = nums[0]; e.TrackAy = nums[1]; e.TrackX = nums[2]; e.TrackY = nums[3];
				}
			}
		}
		else if (key == "SaveTo") {
			e.SaveTo = value; // SaveToは常にコメントなしで出力される
		}
	}
	return entries;
}

// ECC4は小エリアがArea1=K:,Area2~6=I:に分かれて保存されている。
// zone=0はArea0(mergeエリア)を意味し、K:に保存されている。
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

// 小エリア側fvxxファイル(surface=1,2)のどちらかに、rawidが一致するmicrotrackが
// 存在するかを確認する。Area/PL/surfaceでファイルを絞り込み済みのため、
// キーはrawidのみで判定する(basetrack内のm[].zoneフィールドは、fvxxファイル自体の
// microtrackレコードが持つzoneフィールドと同じ意味とは限らないため使用しない)。
static bool MicroTrackExistsInFvxx(int eccnum, int pl, int zone, const vxx::micro_track_subset_t& m) {
	std::string base = ECCAreaPath(eccnum, zone);
	for (int surface = 1; surface <= 2; surface++) {
		std::stringstream f;
		f << base << "\\Area" << zone << "\\PL" << std::setw(3) << std::setfill('0') << pl
			<< "\\f" << std::setw(3) << std::setfill('0') << pl << surface << "_thick_0.vxx";
		if (!std::filesystem::exists(f.str())) continue;

		vxx::FvxxReader fr;
		std::array<int, 2> idx = { (int)m.rawid, (int)m.rawid + 1 };
		std::vector<vxx::micro_track_t> res = fr.ReadAll(f.str(), pl, zone, vxx::opt::index = idx);
		for (auto& mt : res) {
			if (mt.rawid == m.rawid) return true;
		}
	}
	return false;
}

struct VerifyResult {
	bool ok = false;
	long long MergedRawID = -1; // basetrackのrawid(Area0側)。不明な場合は-1のまま。
	std::string reason;
};

// zone(小エリア,1~6)側のbasetrack(rawidSmall)を起点に、
// 1) 小エリアbvxxにbasetrackが存在するか
// 2) そのbasetrackを構成するmicrotrack(m[0],m[1])が小エリアfvxxに存在するか
// 3) 同じmicrotrackペアを持つbasetrackがArea0(merge)bvxxに存在するか
// を確認する。knownMergedRawidが既知(実在トラック)の場合は、その値でArea0側を直接検索する。
// 不明(Referenceトラック)の場合は、Area0側の同一PL分を一度だけ読み込みキャッシュして検索する。
static VerifyResult VerifyTrack(int eccnum, int pl, int zone, long long rawidSmall, long long knownMergedRawid,
	std::map<int, std::vector<vxx::base_track_t>>& area0Cache) {
	VerifyResult r;

	// 1. 小エリアbvxxでbasetrackを検索
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

	// 2. 小エリアfvxxにmicrotrackが存在するか
	if (!MicroTrackExistsInFvxx(eccnum, pl, zone, b.m[0]) || !MicroTrackExistsInFvxx(eccnum, pl, zone, b.m[1])) {
		r.reason = "microtrack not found in Area" + std::to_string(zone) + " fvxx";
		return r;
	}

	// 3. Area0(merge)側の確認
	std::stringstream area0File, cmFile;
	area0File << ECCAreaPath(eccnum, 0) << "\\Area0\\PL" << std::setw(3) << std::setfill('0') << pl
		<< "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
	cmFile << ECCAreaPath(eccnum, 0) << "\\Area0\\0\\align\\corrmap-abs.lst";
	if (!std::filesystem::exists(area0File.str())) {
		r.reason = "Area0 bvxx file not found";
		return r;
	}

	if (knownMergedRawid != -1) {
		// 実在トラック: mergeエリア側のrawidが既知なので直接検索する
		vxx::BvxxReader br0;
		std::array<int, 2> idx0 = { (int)knownMergedRawid, (int)knownMergedRawid + 1 };
		std::vector<vxx::base_track_t> res0 = br0.ReadAll(area0File.str(), pl, 0, vxx::opt::c = cmFile.str(), vxx::opt::index = idx0);
		vxx::base_track_t mb;
		bool foundMerged = false;
		for (auto& t : res0) {
			if (t.rawid == knownMergedRawid) { mb = t; foundMerged = true; break; }
		}
		if (!foundMerged) {
			r.reason = "basetrack not found in Area0 bvxx";
			return r;
		}
		bool matchDirect = (mb.m[0].zone == b.m[0].zone && mb.m[0].rawid == b.m[0].rawid && mb.m[1].zone == b.m[1].zone && mb.m[1].rawid == b.m[1].rawid);
		bool matchSwap = (mb.m[0].zone == b.m[1].zone && mb.m[0].rawid == b.m[1].rawid && mb.m[1].zone == b.m[0].zone && mb.m[1].rawid == b.m[0].rawid);
		if (!matchDirect && !matchSwap) {
			r.reason = "Area0 basetrack's microtracks don't match Area" + std::to_string(zone) + " basetrack's microtracks";
			return r;
		}
		r.MergedRawID = mb.rawid;
	}
	else {
		// Referenceトラック: mergeエリア側のrawidが不明なので、同一microtrackペアを持つbasetrackを探す
		if (area0Cache.find(pl) == area0Cache.end()) {
			vxx::BvxxReader br0;
			area0Cache[pl] = br0.ReadAll(area0File.str(), pl, 0, vxx::opt::c = cmFile.str());
		}
		bool foundMerged = false;
		for (auto& mb : area0Cache[pl]) {
			bool matchDirect = (mb.m[0].zone == b.m[0].zone && mb.m[0].rawid == b.m[0].rawid && mb.m[1].zone == b.m[1].zone && mb.m[1].rawid == b.m[1].rawid);
			bool matchSwap = (mb.m[0].zone == b.m[1].zone && mb.m[0].rawid == b.m[1].rawid && mb.m[1].zone == b.m[0].zone && mb.m[1].rawid == b.m[0].rawid);
			if (matchDirect || matchSwap) { r.MergedRawID = mb.rawid; foundMerged = true; break; }
		}
		if (!foundMerged) {
			r.reason = "no matching basetrack found in Area0 bvxx";
			return r;
		}
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

			int zone = -1;
			long long rawidSmall = -1;
			long long knownMergedRawid = -1;
			double tax = 0, tay = 0, tx = 0, ty = 0;
			std::string type;

			if (isReference) {
				if (!e.HasRawID || !e.HasZone) continue; // 情報が揃っていないものはスキップ
				zone = e.Zone;
				rawidSmall = e.RawID;
				knownMergedRawid = -1;
				if (e.HasTrackComment) { tax = e.CTrackAx; tay = e.CTrackAy; tx = e.CTrackX; ty = e.CTrackY; }
				type = "Reference";
			}
			else if (e.HasRawID && e.HasRawIDComment && e.HasZoneComment) {
				// 実在トラック(小エリアへのマッピングに成功しているもののみ確認対象)
				zone = e.ZoneComment;
				rawidSmall = e.RawIDComment;
				knownMergedRawid = e.RawID;
				if (e.HasTrackComment) { tax = e.CTrackAx; tay = e.CTrackAy; tx = e.CTrackX; ty = e.CTrackY; }
				type = "Real";
			}
			else {
				continue; // 外挿トラック、またはマッピング情報が無い実在トラックは対象外
			}

			totalChecked++;
			VerifyResult vr = VerifyTrack(eccnum, pl, zone, rawidSmall, knownMergedRawid, area0Cache);
			if (!vr.ok) {
				totalFailed++;
				int event = ExtractEventNumber(e.SaveTo);
				long long mergedOut = (vr.MergedRawID != -1) ? vr.MergedRawID : knownMergedRawid;
				std::cout << "\t[NG] PL" << pl << " event" << event << " (" << type << ") : " << vr.reason << std::endl;
				ofs_log << pl << " " << event << " " << mergedOut << " " << zone << " " << rawidSmall
					<< " " << tax << " " << tay << " " << tx << " " << ty
					<< " " << type << " \"" << vr.reason << "\"" << std::endl;
			}
		}
	}

	std::cout << "\n\n === summary ===" << std::endl;
	std::cout << "checked : " << totalChecked << std::endl;
	std::cout << "failed  : " << totalFailed << std::endl;

	return 0;
}
