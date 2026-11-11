// 2026/11/10
// kasumi
// Check_BMmom_timestamp_and_base_consistency
// momch内のmuon chain(particle_flg==13,chainid==0)について、以下の4種類の重複/競合を検出する。
// event headerのweightにbunch情報が入っているため、従来timestamp単独でキーにしていた判定は
// 全てtimestampとbunchのペアをキーとする判定に置き換えている(bunch = (int)ev.weight)。
//   (1) 完全重複      : Key{timestamp, bunch, max_pl, max_pl_rawid, bm_range_mom, ecc_mcs_mom} が全て一致
//   (2) timestamp重複 : Key{timestamp, bunch} が一致 (運動量評価の母集団)
//   (3) ECC分岐       : Key{timestamp, bunch, bm_range_mom} が一致 (BMは同一だがECCの読み取り結果が分岐)
//   (4) timestamp競合 : Key{max_pl, max_pl_rawid} が一致するのに(timestamp,bunch)が異なる
//                        (同一のECC basetrackが異なる(timestamp,bunch)に重複して割り当てられている)

#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
# include <iostream>
# include <random>

#include <list>
#include <cassert>
#include <filesystem>

struct FullKey {
	int timestamp, bunch, max_pl, max_pl_rawid;
	double bm_mom, ecc_mcs_mom;
};
bool operator<(const FullKey& lhs, const FullKey& rhs) {
	return std::tie(lhs.timestamp, lhs.bunch, lhs.max_pl, lhs.max_pl_rawid, lhs.bm_mom, lhs.ecc_mcs_mom)
		< std::tie(rhs.timestamp, rhs.bunch, rhs.max_pl, rhs.max_pl_rawid, rhs.bm_mom, rhs.ecc_mcs_mom);
}

struct TsBunchKey {
	int timestamp, bunch;
};
bool operator<(const TsBunchKey& lhs, const TsBunchKey& rhs) {
	return std::tie(lhs.timestamp, lhs.bunch) < std::tie(rhs.timestamp, rhs.bunch);
}

struct TsBunchBmKey {
	int timestamp, bunch;
	double bm_mom;
};
bool operator<(const TsBunchBmKey& lhs, const TsBunchBmKey& rhs) {
	return std::tie(lhs.timestamp, lhs.bunch, lhs.bm_mom) < std::tie(rhs.timestamp, rhs.bunch, rhs.bm_mom);
}

struct PlKey {
	int max_pl, max_pl_rawid;
};
bool operator<(const PlKey& lhs, const PlKey& rhs) {
	return std::tie(lhs.max_pl, lhs.max_pl_rawid) < std::tie(rhs.max_pl, rhs.max_pl_rawid);
}

struct Entry {
	int groupid, timestamp, bunch;
	double bm_mom, ecc_mcs_mom;
	int max_pl, max_pl_rawid;// baseの中で最も大きいPLと、そのPLのrawid
};
bool sort_timestamp_bunch_groupid(const Entry& lhs, const Entry& rhs) {
	return std::tie(lhs.timestamp, lhs.bunch, lhs.groupid) < std::tie(rhs.timestamp, rhs.bunch, rhs.groupid);
}

void MeasureProcessingTime(std::chrono::system_clock::time_point& start, std::chrono::system_clock::time_point& end) {
	auto dur = end - start;        // 計測時間を計算
	auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
	// 計測時間をミリ秒（1/1000秒）に変換して表示
	if (msec < 1000) {
		std::cout << msec << " milli sec \n";
	}
	else if (msec / 1000 < 60) {
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

void collect_muon_chain_info(std::vector<Momentum_recon::Event_information>& momch, std::vector<Entry>& all_entries);
void pick_full_duplicates(std::vector<Entry>& all_entries, std::vector<Entry>& matched_list);
void pick_timestamp_duplicates(std::vector<Entry>& all_entries, std::vector<Entry>& ts_dup_list);
void pick_ecc_branch_duplicates(std::vector<Entry>& all_entries, std::vector<Entry>& ecc_branch_list);
void pick_timestamp_conflicts(std::vector<Entry>& all_entries, std::vector<Entry>& ts_conflict_list);
void write_matched_list(std::ostream& os, std::vector<Entry>& matched_list, int list_mode);
void write_entry_list(std::ostream& os, std::vector<Entry>& list);


int main(int argc, char** argv) {
	if (argc < 2 || argc > 7) {
		fprintf(stderr, "===============================================================================\n");
		fprintf(stderr, " usage : prg in.momch [out_full.txt] [out_timestamp_dup.txt] [out_ecc_branch.txt] [out_timestamp_conflict.txt] [list_mode:default==0]\n\n");
		fprintf(stderr, " muon chain(particle_flg==13,chainid==0)について、以下の4種類を検出する。\n");
		fprintf(stderr, " event headerのweightにbunch情報が入っているため、timestampのみでキーにしていた判定は\n");
		fprintf(stderr, " 全てtimestampとbunchのペアをキーとする判定に置き換えている。\n");
		fprintf(stderr, "  (1) 完全重複      : Key{timestamp,bunch,max_pl,rawid(at max_pl),bm_range_mom,ecc_mcs_mom} が全て一致\n");
		fprintf(stderr, "                      -> out_full.txt\n");
		fprintf(stderr, "  (2) timestamp重複 : Key{timestamp,bunch} が一致 (運動量評価の母集団)\n");
		fprintf(stderr, "                      -> out_timestamp_dup.txt\n");
		fprintf(stderr, "  (3) ECC分岐       : Key{timestamp,bunch,bm_range_mom} が一致 (BMは同一だがECCの読み取り結果が分岐)\n");
		fprintf(stderr, "                      -> out_ecc_branch.txt\n");
		fprintf(stderr, "  (4) timestamp競合 : Key{max_pl,rawid(at max_pl)} が一致するのに(timestamp,bunch)が異なる\n");
		fprintf(stderr, "                      -> out_timestamp_conflict.txt\n\n");
		fprintf(stderr, " 出力順はいずれも timestamp昇順, bunch昇順, groupid昇順。\n");
		fprintf(stderr, " * list_mode (out_full.txtの内容。画面表示や他の出力には影響しない)\n");
		fprintf(stderr, "   0:timestamp,bunch,groupid,max_pl,rawid(at max_pl) (default)\n");
		fprintf(stderr, "   1:timestamp,bunch,groupid\n");
		fprintf(stderr, "   2:groupid\n");
		fprintf(stderr, "===============================================================================\n");
		exit(1);
	}
	std::string in_momch = argv[1];// input momch
	std::string out_full, out_ts_dup, out_ecc_branch, out_ts_conflict;
	int list_mode = 0;
	if (argc >= 3) {
		out_full = argv[2];
	}
	if (argc >= 4) {
		out_ts_dup = argv[3];
	}
	if (argc >= 5) {
		out_ecc_branch = argv[4];
	}
	if (argc >= 6) {
		out_ts_conflict = argv[5];
	}
	if (argc == 7) {
		list_mode = std::stoi(argv[6]);
	}

	auto start = std::chrono::system_clock::now();//for measure working time

	// read momch
	std::cout << " Reading momch" << std::endl;
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);

	// collect muon chain(particle_flg==13,chainid==0) info
	std::vector<Entry> all_entries;
	collect_muon_chain_info(momch, all_entries);

	// (1) 完全重複 : timestamp,bunch,max_pl,rawid,bm_range_mom,ecc_mcs_momが全て一致
	std::vector<Entry> matched_list;
	pick_full_duplicates(all_entries, matched_list);

	// (2) timestamp重複 : timestamp,bunchが一致(運動量評価の母集団)
	std::vector<Entry> ts_dup_list;
	pick_timestamp_duplicates(all_entries, ts_dup_list);

	// (3) ECC分岐 : timestamp,bunch,bm_range_momが一致(BMは同一だがECCの結果が分岐)
	std::vector<Entry> ecc_branch_list;
	pick_ecc_branch_duplicates(all_entries, ecc_branch_list);

	// (4) timestamp競合 : max_pl,rawidが一致するのに(timestamp,bunch)が異なる
	std::vector<Entry> ts_conflict_list;
	pick_timestamp_conflicts(all_entries, ts_conflict_list);

	// timestamp昇順,bunch昇順,groupid昇順に並べ替える
	std::sort(matched_list.begin(), matched_list.end(), sort_timestamp_bunch_groupid);
	std::sort(ts_dup_list.begin(), ts_dup_list.end(), sort_timestamp_bunch_groupid);
	std::sort(ecc_branch_list.begin(), ecc_branch_list.end(), sort_timestamp_bunch_groupid);
	std::sort(ts_conflict_list.begin(), ts_conflict_list.end(), sort_timestamp_bunch_groupid);

	std::cout << "\n *** (1) Complete duplicate (timestamp,bunch,max_pl,rawid,bm_mom,ecc_mcs_mom) list [timestamp asc, bunch asc, groupid asc] ***" << std::endl;
	write_matched_list(std::cout, matched_list, list_mode);
	if (!out_full.empty()) {
		std::ofstream ofs(out_full);
		if (!ofs) {
			std::cerr << "\nFailed to create the output file. :" << out_full << std::endl;
		}
		else {
			write_matched_list(ofs, matched_list, list_mode);
			std::cout << "\nWriting the complete duplicate list to the output file is complete : " << out_full << std::endl;
		}
	}

	std::cout << "\n *** (2) Timestamp duplicate (for momentum check) list [timestamp asc, bunch asc, groupid asc] ***" << std::endl;
	write_entry_list(std::cout, ts_dup_list);
	if (!out_ts_dup.empty()) {
		std::ofstream ofs2(out_ts_dup);
		if (!ofs2) {
			std::cerr << "\nFailed to create the output file. :" << out_ts_dup << std::endl;
		}
		else {
			write_entry_list(ofs2, ts_dup_list);
			std::cout << "\nWriting the timestamp duplicate list to the output file is complete : " << out_ts_dup << std::endl;
		}
	}

	std::cout << "\n *** (3) ECC branch (timestamp,bunch,bm_mom match but ECC differs) list [timestamp asc, bunch asc, groupid asc] ***" << std::endl;
	write_entry_list(std::cout, ecc_branch_list);
	if (!out_ecc_branch.empty()) {
		std::ofstream ofs3(out_ecc_branch);
		if (!ofs3) {
			std::cerr << "\nFailed to create the output file. :" << out_ecc_branch << std::endl;
		}
		else {
			write_entry_list(ofs3, ecc_branch_list);
			std::cout << "\nWriting the ECC branch list to the output file is complete : " << out_ecc_branch << std::endl;
		}
	}

	std::cout << "\n *** (4) Timestamp conflict (max_pl,rawid match but (timestamp,bunch) differs) list [timestamp asc, bunch asc, groupid asc] ***" << std::endl;
	write_entry_list(std::cout, ts_conflict_list);
	if (!out_ts_conflict.empty()) {
		std::ofstream ofs4(out_ts_conflict);
		if (!ofs4) {
			std::cerr << "\nFailed to create the output file. :" << out_ts_conflict << std::endl;
		}
		else {
			write_entry_list(ofs4, ts_conflict_list);
			std::cout << "\nWriting the timestamp conflict list to the output file is complete : " << out_ts_conflict << std::endl;
		}
	}

	std::cout << "\n *** Summary ***" << std::endl;
	std::cout << "  # of muon chain (particle_flg==13,chainid==0) : " << all_entries.size() << std::endl;
	std::cout << "  # of (1) complete duplicate events             : " << matched_list.size() << std::endl;
	std::cout << "  # of (2) timestamp duplicate events            : " << ts_dup_list.size() << std::endl;
	std::cout << "  # of (3) ECC branch events                     : " << ecc_branch_list.size() << std::endl;
	std::cout << "  # of (4) timestamp conflict events             : " << ts_conflict_list.size() << std::endl;

	auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
	MeasureProcessingTime(start, end);
	std::cout << "Fin." << std::endl;
}

void collect_muon_chain_info(std::vector<Momentum_recon::Event_information>& momch, std::vector<Entry>& all_entries) {

	for (auto& ev : momch) {
		for (auto& c : ev.chains) {
			if (c.particle_flg != 13) continue;
			if (c.chainid != 0) continue;
			// bm_range_momとtimestampが無効値(未設定)のchainは対象から除外する
			if (c.bm_range_mom <= 0) continue;
			if (ev.unix_time <= 0) continue;
			if (c.base.empty()) continue;// baseが無いchainは最大PLが定義できないため除外する

			int max_pl = std::numeric_limits<int>::min();
			int max_pl_rawid = -1;
			for (auto& b : c.base) {
				if (b.pl > max_pl) {
					max_pl = b.pl;
					max_pl_rawid = b.rawid;
				}
			}

			Entry e;
			e.groupid = ev.groupid;
			e.timestamp = ev.unix_time;
			e.bunch = (int)ev.weight;// event headerのweightにbunch情報が入っている
			e.bm_mom = c.bm_range_mom;
			e.ecc_mcs_mom = c.ecc_mcs_mom[0];// [0]:muon
			e.max_pl = max_pl;
			e.max_pl_rawid = max_pl_rawid;

			all_entries.push_back(e);
		}
	}

	std::cout << "\t* The Number of muon chain (particle_flg==13,chainid==0) is " << all_entries.size() << std::endl;
}

void pick_full_duplicates(std::vector<Entry>& all_entries, std::vector<Entry>& matched_list) {

	std::map<FullKey, std::vector<int>> full_map;
	for (size_t i = 0; i < all_entries.size(); i++) {
		FullKey k{ all_entries[i].timestamp, all_entries[i].bunch, all_entries[i].max_pl, all_entries[i].max_pl_rawid,
			all_entries[i].bm_mom, all_entries[i].ecc_mcs_mom };
		full_map[k].push_back((int)i);
	}

	int dup_group_num = 0;
	for (auto& kv : full_map) {
		std::vector<int>& idx = kv.second;
		if (idx.size() < 2) continue;// 重複していないchainは対象外
		dup_group_num++;

		std::cout << "\n [FULL DUP] timestamp = " << kv.first.timestamp << " , bunch = " << kv.first.bunch
			<< " , max_pl = " << kv.first.max_pl << " , rawid(at max_pl) = " << kv.first.max_pl_rawid
			<< " , bm_range_mom = " << kv.first.bm_mom << " , ecc_mcs_mom = " << kv.first.ecc_mcs_mom
			<< " : " << idx.size() << " events." << std::endl;
		for (int i : idx) {
			Entry& e = all_entries[i];
			std::cout << "\t timestamp = " << e.timestamp << " , bunch = " << e.bunch << " , groupid = " << e.groupid << std::endl;
			matched_list.push_back(e);
		}
	}

	std::cout << "\t* The Number of distinct (timestamp,bunch,max_pl,rawid,bm_mom,ecc_mcs_mom) is " << full_map.size() << std::endl;
	std::cout << "\t* The Number of complete-duplicate groups (2+ events) is " << dup_group_num << std::endl;
}

void pick_timestamp_duplicates(std::vector<Entry>& all_entries, std::vector<Entry>& ts_dup_list) {

	std::map<TsBunchKey, std::vector<int>> ts_map;
	for (size_t i = 0; i < all_entries.size(); i++) {
		TsBunchKey k{ all_entries[i].timestamp, all_entries[i].bunch };
		ts_map[k].push_back((int)i);
	}

	for (auto& kv : ts_map) {
		std::vector<int>& idx = kv.second;
		if (idx.size() < 2) continue;// timestamp,bunchで重複していないchainは対象外

		for (int i : idx) {
			Entry& e = all_entries[i];
			std::cout << "\n [TIMESTAMP DUP] timestamp = " << e.timestamp << " , bunch = " << e.bunch << " , groupid = " << e.groupid
				<< " , max_pl = " << e.max_pl << " , rawid(at max_pl) = " << e.max_pl_rawid
				<< " , bm_range_mom = " << e.bm_mom << " , ecc_mcs_mom = " << e.ecc_mcs_mom << std::endl;
			ts_dup_list.push_back(e);
		}
	}
}

void pick_ecc_branch_duplicates(std::vector<Entry>& all_entries, std::vector<Entry>& ecc_branch_list) {

	std::map<TsBunchBmKey, std::vector<int>> branch_map;
	for (size_t i = 0; i < all_entries.size(); i++) {
		TsBunchBmKey k{ all_entries[i].timestamp, all_entries[i].bunch, all_entries[i].bm_mom };
		branch_map[k].push_back((int)i);
	}

	for (auto& kv : branch_map) {
		std::vector<int>& idx = kv.second;
		if (idx.size() < 2) continue;// timestamp,bunch,bm_momが一致していないchainは対象外

		for (int i : idx) {
			Entry& e = all_entries[i];
			std::cout << "\n [ECC BRANCH] timestamp = " << e.timestamp << " , bunch = " << e.bunch << " , bm_range_mom = " << e.bm_mom
				<< " , groupid = " << e.groupid << " , max_pl = " << e.max_pl
				<< " , rawid(at max_pl) = " << e.max_pl_rawid << " , ecc_mcs_mom = " << e.ecc_mcs_mom << std::endl;
			ecc_branch_list.push_back(e);
		}
	}
}

void pick_timestamp_conflicts(std::vector<Entry>& all_entries, std::vector<Entry>& ts_conflict_list) {

	std::map<PlKey, std::vector<int>> pl_map;
	for (size_t i = 0; i < all_entries.size(); i++) {
		PlKey k{ all_entries[i].max_pl, all_entries[i].max_pl_rawid };
		pl_map[k].push_back((int)i);
	}

	for (auto& kv : pl_map) {
		std::vector<int>& idx = kv.second;
		if (idx.size() < 2) continue;// max_pl,rawidが一致していないchainは対象外

		// 同一(timestamp,bunch)内での一致は「timestamp競合」ではないため、
		// 異なる(timestamp,bunch)が混在している場合のみ対象とする
		bool has_diff_ts = false;
		for (size_t a = 1; a < idx.size(); a++) {
			if (all_entries[idx[a]].timestamp != all_entries[idx[0]].timestamp ||
				all_entries[idx[a]].bunch != all_entries[idx[0]].bunch) {
				has_diff_ts = true;
				break;
			}
		}
		if (!has_diff_ts) continue;

		std::cout << "\n [TIMESTAMP CONFLICT] max_pl = " << kv.first.max_pl << " , rawid(at max_pl) = " << kv.first.max_pl_rawid
			<< " : " << idx.size() << " events." << std::endl;
		for (int i : idx) {
			Entry& e = all_entries[i];
			std::cout << "\t timestamp = " << e.timestamp << " , bunch = " << e.bunch << " , groupid = " << e.groupid
				<< " , bm_range_mom = " << e.bm_mom << " , ecc_mcs_mom = " << e.ecc_mcs_mom << std::endl;
			ts_conflict_list.push_back(e);
		}
	}
}

void write_matched_list(std::ostream& os, std::vector<Entry>& matched_list, int list_mode) {
	for (auto& e : matched_list) {
		if (list_mode == 2) {
			// groupidのみ
			os << e.groupid << std::endl;
		}
		else if (list_mode == 1) {
			// timestamp,bunch,groupid
			os << e.timestamp << "\t" << e.bunch << "\t" << e.groupid << std::endl;
		}
		else {
			// timestamp,bunch,groupid,max_pl,rawid(at max_pl)
			os << e.timestamp << "\t" << e.bunch << "\t" << e.groupid << "\t" << e.max_pl << "\t" << e.max_pl_rawid << std::endl;
		}
	}
}

void write_entry_list(std::ostream& os, std::vector<Entry>& list) {
	// timestamp,bunch,groupid,max_pl,rawid(at max_pl),bm_range_mom,ecc_mcs_mom
	for (auto& e : list) {
		os << e.timestamp << "\t" << e.bunch << "\t" << e.groupid << "\t" << e.max_pl << "\t" << e.max_pl_rawid
			<< "\t" << e.bm_mom << "\t" << e.ecc_mcs_mom << std::endl;
	}
}
