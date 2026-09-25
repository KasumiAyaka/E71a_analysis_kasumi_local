// 2026/09/10
// K:\NINJA\E71a\ManualCheck\TrackList\ECC* の中の PLxxx.yaml
// (Pickup_reference_track_for_each_track.cpp の出力: 元のトラック + 見つかった近くのトラックのペア)
// を読み込み、btrklist{ECC番号}_rawidmap.txt のフォーマットで出力する。
// 出力フォーマット: PL eventid rawid1 zone rawid2 ax ay x y ph1 ph2 trk_type
// 実在・外挿・近くのトラック(Reference)それぞれに対して1行ずつ出力する(1ペア=2行)。
//   trk_type = 0 (実在)   : rawid1=自身のRawID, rawid2=-1
//                            ax,ay,x,y,ph1,ph2 は merged(Area0)側 basetrack の値(vxxファイルから、%10000しない生値)
//                            zone は Area1-6をクロスチェックし実際に見つかったエリア(フォールバック時は近くのトラックのzone)
//   trk_type = 1 (外挿)   : rawid1=-1, rawid2=-1, zone=YAML記載のZone
//                            ax,ay,x,y は YAML記載のTrack(予測位置)。ph1,ph2はvxx情報を持たないため-1
//   trk_type = -1 (近くのトラック/Reference) : rawid1=-1, rawid2=自身のRawID, zone=YAML記載のZone
//                            ax,ay,x,y,ph1,ph2 は自身のzoneの basetrack の値(vxxファイルから、%10000しない生値)
// kasumi
#pragma comment(lib, "VxxReader.lib")
#include "VxxReader.h"

#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <utility>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <regex>

// TrackList\ECC*\PLxxx.yaml の1エントリ。
// 実在するトラック : "- Plate / RawID / Surface / SaveTo"              (RawIDあり, Zoneなし, Trackなし)
// 外挿したトラック : "- Plate / Track: [ax,ay,x,y] / Zone / Surface / SaveTo" (RawIDなし, Zoneあり, Trackあり)
// 近くのトラック   : "- Plate / RawID / Zone / Surface / SaveTo"        (RawIDあり, Zoneあり, Trackなし)
struct TrackEntry {
	int Plate = 0;
	long long RawID = -1;
	bool HasTrack = false;
	double TrackAx = 0, TrackAy = 0, TrackX = 0, TrackY = 0;
	int Zone = -1;
	std::string Surface;
	std::string SaveTo;
};

enum class EntryKind { Real, Extrapolated, Nearby, Unknown };

// クロスチェックでArea1-6のどこにも見つからなかった実在トラック
struct NotFoundEntry {
	int pl;
	int eventid;
	long long rawid;
	double ax, ay, x, y;
};
// クロスチェックで見つかったが、近くのトラックのzoneと異なっていた実在トラック
struct ZoneMismatchEntry {
	int pl;
	int eventid;
	long long rawid;
};

static std::string trim(const std::string& s) {
	size_t b = s.find_first_not_of(" \t\r\n");
	if (b == std::string::npos) return "";
	size_t e = s.find_last_not_of(" \t\r\n");
	return s.substr(b, e - b + 1);
}

// "G:ECC3/IMG/Event14570/PL016" -> 14570
static int ExtractEventNumber(const std::string& saveTo) {
	size_t pos = saveTo.find("Event");
	if (pos == std::string::npos) return -1;
	return std::stoi(saveTo.substr(pos + 5, 5));
}

static EntryKind ClassifyEntry(const TrackEntry& e) {
	bool hasRawID = (e.RawID != -1);
	bool hasZone = (e.Zone != -1);
	if (hasRawID && !e.HasTrack && !hasZone) return EntryKind::Real;
	if (!hasRawID && e.HasTrack && hasZone) return EntryKind::Extrapolated;
	if (hasRawID && !e.HasTrack && hasZone) return EntryKind::Nearby;
	return EntryKind::Unknown;
}

// "G:ECC6/IMG/Event01260/PL019" -> "G:ECC6/Ref/IMG/Event01260/PL019"
// "G:ECC2/Ref/Event08521/PL020" -> "G:ECC2/Ref/Ref/Event08521/PL020"
// Pickup_reference_track_for_each_track.cpp の InsertRefFolder と同じ変換。
// 近くのトラックのSaveToは、元のトラックのSaveToをこの変換で得られるものと一致するはず。
static std::string InsertRefFolder(const std::string& saveTo) {
	size_t pos = saveTo.find('/');
	if (pos == std::string::npos) return saveTo;
	return saveTo.substr(0, pos + 1) + "Ref/" + saveTo.substr(pos + 1);
}

static std::vector<TrackEntry> Read_yaml(const std::string& filename) {
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

// eccnum, zone (Area0=merged の場合は IsMerged=true) から basetrackファイルの置かれているドライブ/ルートを決める。
// Pickup_reference_track_for_each_track.cpp の file_in_ECC / Search_basetrack の分岐をそのまま踏襲している。
static std::string ResolveEccRoot(int eccnum, int zone, bool isMerged) {
	std::stringstream root;
	if (eccnum == 4) {
		if (isMerged) {
			// 実トラックのArea0(merged)参照は常にK: (main()の実トラック分岐と同じ)
			root << "K:\\NINJA\\E71a\\ECC" << eccnum;
		}
		else if (zone == 1) {
			root << "K:\\NINJA\\E71a\\ECC" << eccnum;
		}
		else {
			root << "I:\\NINJA\\E71a\\ECC" << eccnum;
		}
	}
	else if (eccnum < 4) {
		root << "T:\\NINJA\\E71a\\ECC" << eccnum;
	}
	else if (eccnum < 7) {
		root << "I:\\NINJA\\E71a\\ECC" << eccnum;
	}
	else {
		root << "K:\\NINJA\\E71a\\ECC" << eccnum;
	}
	return root.str();
}

static std::string BasetrackFilePath(int eccnum, int zone, int pl, bool isMerged) {
	std::stringstream path;
	path << ResolveEccRoot(eccnum, zone, isMerged)
		<< "\\Area" << (isMerged ? 0 : zone)
		<< "\\PL" << std::setw(3) << std::setfill('0') << pl
		<< "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
	return path.str();
}

// 1プレート分の処理中、同じvxxファイルを何度も読み直さないためのキャッシュ。
// (K:/I:/T:はネットワークドライブのため、同じファイルの再オープンが特に高コスト)
struct PlateCache {
	// zone -> (rawid -> basetrack)。「近くのトラック」参照用。zone単位で遅延ロードする。
	std::map<int, std::map<long long, vxx::base_track_t>> zoneTracks;
	// (m[0].rawid, m[1].rawid) -> zone。クロスチェック用。Area1-6を一度に読み込んで構築する。
	std::map<std::pair<long long, long long>, int> crossCheckIndex;
	bool allZonesLoaded = false;
	// Area0(merged)の rawid -> basetrack。実在トラックのph1/ph2/クロスチェック用。
	std::map<long long, vxx::base_track_t> mergedTracks;
	bool mergedLoaded = false;
};

// 指定zoneのvxxファイルをまだ読んでいなければ1回だけ読み込み、
// rawid検索用マップとクロスチェック用インデックスの両方に登録する。
static void LoadZoneFile(PlateCache& cache, int eccnum, int pl, int zone) {
	if (cache.zoneTracks.count(zone)) return; // 読み込み済み
	auto& m = cache.zoneTracks[zone]; // 存在しなくても空マップを登録し、二重読み込みを防ぐ
	std::string file = BasetrackFilePath(eccnum, zone, pl, false);
	if (!std::filesystem::exists(file)) return;
	vxx::BvxxReader br;
	std::vector<vxx::base_track_t> base = br.ReadAll(file, pl, zone);
	for (auto& b : base) {
		m[b.rawid] = b;
		cache.crossCheckIndex[{ b.m[0].rawid, b.m[1].rawid }] = zone;
	}
}

// クロスチェックはArea1-6のどこにあるか分からないため、全zoneをまとめてロードする。
static void EnsureAllZonesLoaded(PlateCache& cache, int eccnum, int pl) {
	if (cache.allZonesLoaded) return;
	for (int zone = 1; zone <= 6; zone++) {
		LoadZoneFile(cache, eccnum, pl, zone);
	}
	cache.allZonesLoaded = true;
}

// Area0(merged)のvxxファイルをまだ読んでいなければ1回だけ読み込む。
static void EnsureMergedLoaded(PlateCache& cache, int eccnum, int pl) {
	if (cache.mergedLoaded) return;
	cache.mergedLoaded = true;
	std::string file = BasetrackFilePath(eccnum, 0, pl, true);
	if (!std::filesystem::exists(file)) return;
	vxx::BvxxReader br;
	std::vector<vxx::base_track_t> base = br.ReadAll(file, pl, 0);
	for (auto& b : base) {
		cache.mergedTracks[b.rawid] = b;
	}
}

// TrackList\ECC{N} の1プレート分 (PLxxx.yaml) を処理し、rawidmapの行をofsに書き出す。
static void ProcessPlate(const std::string& yamlPath, int eccnum, int pl, std::ofstream& ofs, int& written, int& skipped,
	std::vector<NotFoundEntry>& notFoundLog, std::vector<ZoneMismatchEntry>& zoneMismatchLog, PlateCache& cache) {
	std::vector<TrackEntry> entries = Read_yaml(yamlPath);

	size_t i = 0;
	while (i < entries.size()) {
		EntryKind kind = ClassifyEntry(entries[i]);
		if (kind != EntryKind::Real && kind != EntryKind::Extrapolated) {
			std::cout << "\t[WARN] " << yamlPath << " : entry " << i << " is not a real/extrapolated track. skipped." << std::endl;
			++i;
			++skipped;
			continue;
		}
		if (i + 1 >= entries.size() || ClassifyEntry(entries[i + 1]) != EntryKind::Nearby) {
			std::cout << "\t[WARN] " << yamlPath << " : entry " << i << " has no paired nearby track. skipped." << std::endl;
			++i;
			++skipped;
			continue;
		}

		const TrackEntry& orig = entries[i];
		const TrackEntry& nearby = entries[i + 1];

		// 近くのトラックのSaveToが、元のトラックのSaveToにRef/を挿入したものと一致するか確認する。
		std::string expectedNearbySaveTo = InsertRefFolder(orig.SaveTo);
		if (nearby.SaveTo != expectedNearbySaveTo) {
			std::cout << "\t[WARN] " << yamlPath << " : entry " << i
				<< " nearby SaveTo mismatch. expected=" << expectedNearbySaveTo
				<< " actual=" << nearby.SaveTo << std::endl;
		}

		int eventid = ExtractEventNumber(orig.SaveTo);
		int nearbyZone = nearby.Zone;   // 近くのトラック(Reference)自身が記録されているエリア
		long long nearbyRawid = nearby.RawID;

		// 近くのトラックは常にvxxファイルから読む(キャッシュ経由。zoneファイルはプレート内で使い回す)
		LoadZoneFile(cache, eccnum, pl, nearbyZone);
		auto& nearbyZoneMap = cache.zoneTracks[nearbyZone];
		auto nearbyIt = nearbyZoneMap.find(nearbyRawid);
		bool nearbyOk = (nearbyIt != nearbyZoneMap.end());

		// --- 1行目: 元のトラック(実在 or 外挿) ---
		long long origRawid1 = -1;
		int origZone;
		double origAx = 0, origAy = 0, origX = 0, origY = 0;
		long long origPh1 = -1, origPh2 = -1;
		int origType;

		if (kind == EntryKind::Real) {
			origType = 0;
			origRawid1 = orig.RawID;
			origZone = nearbyZone; // クロスチェックできない場合のフォールバック

			EnsureMergedLoaded(cache, eccnum, pl);
			auto mergedIt = cache.mergedTracks.find(origRawid1);
			bool mergedOk = (mergedIt != cache.mergedTracks.end());
			if (mergedOk) {
				origAx = mergedIt->second.ax; origAy = mergedIt->second.ay;
				origX = mergedIt->second.x; origY = mergedIt->second.y;
				origPh1 = mergedIt->second.m[0].ph; // %10000しない生値
				origPh2 = mergedIt->second.m[1].ph;

				// クロスチェック: mergedのbasetrackが持つ子rawid(m[0].rawidとm[1].rawidのペアが一致)で
				// Area1〜6を総当たりし、実際にどのエリアにこの実在トラックがあるかを確認する。
				EnsureAllZonesLoaded(cache, eccnum, pl);
				auto czIt = cache.crossCheckIndex.find({ mergedIt->second.m[0].rawid, mergedIt->second.m[1].rawid });
				if (czIt == cache.crossCheckIndex.end()) {
					// Area1-6のどこにも見つからなかった -> 最後にまとめて出力する
					notFoundLog.push_back({ pl, eventid, origRawid1, mergedIt->second.ax, mergedIt->second.ay, mergedIt->second.x, mergedIt->second.y });
				}
				else {
					if (czIt->second != nearbyZone) {
						// 見つかったが、近くのトラックのzoneと異なる -> 最後にまとめて出力する
						zoneMismatchLog.push_back({ pl, eventid, origRawid1 });
					}
					origZone = czIt->second;
				}
			}
		}
		else { // Extrapolated
			origType = 1;
			origZone = orig.Zone;
			origAx = orig.TrackAx; origAy = orig.TrackAy;
			origX = orig.TrackX; origY = orig.TrackY;
			// ph1,ph2はvxx情報を持たないため-1のまま
		}

		ofs << pl << " " << eventid << " " << origRawid1 << " " << origZone << " " << -1 << " "
			<< origAx << " " << origAy << " " << origX << " " << origY << " " << origPh1 << " " << origPh2 << " " << origType << std::endl;

		// --- 2行目: 近くのトラック(Reference) ---
		double refAx = 0, refAy = 0, refX = 0, refY = 0;
		long long refPh1 = -1, refPh2 = -1;
		if (nearbyOk) {
			refAx = nearbyIt->second.ax; refAy = nearbyIt->second.ay;
			refX = nearbyIt->second.x; refY = nearbyIt->second.y;
			refPh1 = nearbyIt->second.m[0].ph; // %10000しない生値
			refPh2 = nearbyIt->second.m[1].ph;
		}

		ofs << pl << " " << eventid << " " << -1 << " " << nearbyZone << " " << nearbyRawid << " "
			<< refAx << " " << refAy << " " << refX << " " << refY << " " << refPh1 << " " << refPh2 << " " << -1 << std::endl;

		written += 2;

		i += 2;
	}
}

int main(int argc, char** argv) {
	// usage: Make_rawidmap_from_TrackList.exe [baseDir] [eccnum]
	//   baseDir 省略時: K:\NINJA\E71a\ManualCheck\TrackList
	//   eccnum  省略時: baseDir直下の ECC<数字> フォルダをすべて処理する
	//           指定時: そのECC番号のフォルダのみ処理する
	std::string baseDir = "K:\\NINJA\\E71a\\ManualCheck\\TrackList";
	int onlyEccnum = -1;
	if (argc >= 2) baseDir = argv[1];
	if (argc >= 3) onlyEccnum = std::stoi(argv[2]);

	if (!std::filesystem::exists(baseDir)) {
		std::cerr << baseDir << " doesn't exist." << std::endl;
		return 1;
	}

	std::regex eccDirPattern("ECC(\\d+)");
	std::vector<NotFoundEntry> notFoundLog;
	std::vector<ZoneMismatchEntry> zoneMismatchLog;
	int eccProcessed = 0;

	for (auto& entry : std::filesystem::directory_iterator(baseDir)) {
		if (!entry.is_directory()) continue;
		std::string dirName = entry.path().filename().string();
		std::smatch m;
		if (!std::regex_match(dirName, m, eccDirPattern)) continue;

		int eccnum = std::stoi(m[1].str());
		if (onlyEccnum != -1 && eccnum != onlyEccnum) continue;
		++eccProcessed;
		std::string eccDir = entry.path().string();
		std::string outputLog = baseDir + "\\btrklist" + std::to_string(eccnum) + "_rawidmap.txt";

		std::cout << "\n=== ECC" << eccnum << " (" << eccDir << ") ===" << std::endl;
		if (std::filesystem::exists(outputLog)) {
			std::cout << "\t" << outputLog << " already exists. it will be overwritten." << std::endl;
		}

		std::ofstream ofs(outputLog);
		int written = 0, skipped = 0;
		for (int pl = 3; pl <= 133; pl++) {
			std::stringstream yamlPath;
			yamlPath << eccDir << "\\PL" << std::setw(3) << std::setfill('0') << pl << ".yaml";
			if (!std::filesystem::exists(yamlPath.str())) continue;
			PlateCache cache; // vxxファイルのキャッシュはプレートごとに独立(ファイルが別物なので使い回さない)
			ProcessPlate(yamlPath.str(), eccnum, pl, ofs, written, skipped, notFoundLog, zoneMismatchLog, cache);
		}
		std::cout << "\t" << written << " rows written, " << skipped << " entries skipped -> " << outputLog << std::endl;
	}

	if (onlyEccnum != -1 && eccProcessed == 0) {
		std::cerr << "ECC" << onlyEccnum << " folder not found under " << baseDir << std::endl;
		return 1;
	}

	if (!notFoundLog.empty()) {
		std::cout << "\n\n=== Real tracks NOT found in Area1-6 (cross-check) ===" << std::endl;
		std::cout << "PL Eventid rawid ax ay x y" << std::endl;
		for (auto& e : notFoundLog) {
			std::cout << e.pl << " " << e.eventid << " " << e.rawid << " "
				<< e.ax << " " << e.ay << " " << e.x << " " << e.y << std::endl;
		}
	}
	if (!zoneMismatchLog.empty()) {
		std::cout << "\n\n=== Real tracks whose cross-checked zone differs from the nearby track's zone ===" << std::endl;
		std::cout << "PL Eventid rawid" << std::endl;
		for (auto& e : zoneMismatchLog) {
			std::cout << e.pl << " " << e.eventid << " " << e.rawid << std::endl;
		}
	}

	return 0;
}
