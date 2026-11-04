// 2026/11/04
// kasumi
// Check_BMmom_timestamp_and_base_consistency
// momch内のmuon chain(particle_flg==13,chainid==0)について、
// (timestamp,max_pl,rawid(at max_pl))が完全に一致する場合は「完全な重複」として一覧を表示・出力する。
// timestampが重複しているイベント(完全な重複か否かを問わない)については、運動量評価用に
// timestamp,eventid,max_pl,rawid,BM_range_mom,ECC_MCS_momを別のテキストファイルに出力する。

#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
# include <iostream>
# include <random>

#include <list>
#include <cassert>
#include <filesystem>

struct Key {
	int timestamp, max_pl, max_pl_rawid;
};
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.timestamp, lhs.max_pl, lhs.max_pl_rawid) < std::tie(rhs.timestamp, rhs.max_pl, rhs.max_pl_rawid);
}

struct Entry {
	int groupid, timestamp;
	double bm_mom, ecc_mcs_mom;
	int max_pl, max_pl_rawid;// baseの中で最も大きいPLと、そのPLのrawid
};
bool sort_timestamp_groupid(const Entry& lhs, const Entry& rhs) {
	return std::tie(lhs.timestamp, lhs.groupid) < std::tie(rhs.timestamp, rhs.groupid);
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
void write_matched_list(std::ostream& os, std::vector<Entry>& matched_list, int list_mode);
void write_timestamp_duplicate_list(std::ostream& os, std::vector<Entry>& ts_dup_list);


int main(int argc, char** argv) {
	if (argc < 2 || argc > 5) {
		fprintf(stderr, "===============================================================================\n");
		fprintf(stderr, " usage : prg in.momch [output.txt] [output_timestamp_dup.txt] [list_mode:default==0]\n\n");
		fprintf(stderr, " muon chain(particle_flg==13,chainid==0)について、\n");
		fprintf(stderr, " (timestamp,max_pl,rawid(at max_pl))が完全に一致するものを「完全な重複」として\n");
		fprintf(stderr, " 一覧を画面に表示し、output.txtを指定した場合はファイルにも出力する(timestamp昇順,groupid昇順)。\n");
		fprintf(stderr, " timestampが重複しているイベント(完全な重複か否かを問わない)は運動量評価のため、\n");
		fprintf(stderr, " timestamp,eventid,max_pl,rawid,BM_range_mom,ECC_MCS_momをoutput_timestamp_dup.txtに出力する。\n\n");
		fprintf(stderr, " * list_mode (output.txtの内容。画面表示には影響しない)\n");
		fprintf(stderr, "   0:timestamp,groupid,max_pl,rawid(at max_pl) (default)\n");
		fprintf(stderr, "   1:timestamp,groupid\n");
		fprintf(stderr, "   2:groupid\n");
		fprintf(stderr, "===============================================================================\n");
		exit(1);
	}
	std::string in_momch = argv[1];// input momch
	std::string out_txt, out_ts_dup_txt;
	int list_mode = 0;
	if (argc >= 3) {
		out_txt = argv[2];
	}
	if (argc >= 4) {
		out_ts_dup_txt = argv[3];
	}
	if (argc == 5) {
		list_mode = std::stoi(argv[4]);
	}

	auto start = std::chrono::system_clock::now();//for measure working time

	// read momch
	std::cout << " Reading momch" << std::endl;
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);

	// collect muon chain(particle_flg==13,chainid==0) info
	std::vector<Entry> all_entries;
	collect_muon_chain_info(momch, all_entries);

	// (timestamp,max_pl,rawid(at max_pl))が完全に一致するもの
	std::vector<Entry> matched_list;
	pick_full_duplicates(all_entries, matched_list);

	// timestampが重複しているもの(完全な重複も含む)
	std::vector<Entry> ts_dup_list;
	pick_timestamp_duplicates(all_entries, ts_dup_list);

	// timestamp昇順,groupid昇順に並べ替える
	std::sort(matched_list.begin(), matched_list.end(), sort_timestamp_groupid);
	std::sort(ts_dup_list.begin(), ts_dup_list.end(), sort_timestamp_groupid);

	std::cout << "\n *** Complete duplicate (timestamp,max_pl,rawid) list [timestamp asc, groupid asc] ***" << std::endl;
	write_matched_list(std::cout, matched_list, list_mode);
	if (!out_txt.empty()) {
		std::ofstream ofs(out_txt);
		if (!ofs) {
			std::cerr << "\nFailed to create the output file. :" << out_txt << std::endl;
		}
		else {
			write_matched_list(ofs, matched_list, list_mode);
			std::cout << "\nWriting the complete duplicate list to the output file is complete : " << out_txt << std::endl;
		}
	}

	std::cout << "\n *** Timestamp duplicate (for momentum check) list [timestamp asc, groupid asc] ***" << std::endl;
	write_timestamp_duplicate_list(std::cout, ts_dup_list);
	if (!out_ts_dup_txt.empty()) {
		std::ofstream ofs2(out_ts_dup_txt);
		if (!ofs2) {
			std::cerr << "\nFailed to create the output file. :" << out_ts_dup_txt << std::endl;
		}
		else {
			write_timestamp_duplicate_list(ofs2, ts_dup_list);
			std::cout << "\nWriting the timestamp duplicate list to the output file is complete : " << out_ts_dup_txt << std::endl;
		}
	}

	std::cout << "\n *** Summary ***" << std::endl;
	std::cout << "  # of muon chain (particle_flg==13,chainid==0)                 : " << all_entries.size() << std::endl;
	std::cout << "  # of complete duplicate events (timestamp,max_pl,rawid match) : " << matched_list.size() << std::endl;
	std::cout << "  # of timestamp duplicate events (for momentum check)          : " << ts_dup_list.size() << std::endl;

	auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
	MeasureProcessingTime(start, end);
	std::cout << "Fin." << std::endl;
}

void collect_muon_chain_info(std::vector<Momentum_recon::Event_information>& momch, std::vector<Entry>& all_entries) {

	for (auto& ev : momch) {
		for (auto& c : ev.chains) {
			if (c.particle_flg != 13) continue;
			if (c.chainid != 0) continue;
			// bm_range_momとtimestampが初期値(未設定)のchainは対象から除外する
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

	std::map<Key, std::vector<int>> full_map;
	for (size_t i = 0; i < all_entries.size(); i++) {
		Key k{ all_entries[i].timestamp, all_entries[i].max_pl, all_entries[i].max_pl_rawid };
		full_map[k].push_back((int)i);
	}

	int dup_group_num = 0;
	for (auto& kv : full_map) {
		std::vector<int>& idx = kv.second;
		if (idx.size() < 2) continue;// 重複していないchainは対象外
		dup_group_num++;

		std::cout << "\n [MATCH] timestamp = " << kv.first.timestamp << " , max_pl = " << kv.first.max_pl
			<< " , rawid(at max_pl) = " << kv.first.max_pl_rawid << " : " << idx.size() << " events." << std::endl;
		for (int i : idx) {
			Entry& e = all_entries[i];
			std::cout << "\t timestamp = " << e.timestamp << " , groupid = " << e.groupid
				<< " , bm_range_mom = " << e.bm_mom
				<< " , max_pl = " << e.max_pl << " , rawid(at max_pl) = " << e.max_pl_rawid << std::endl;
			matched_list.push_back(e);
		}
	}

	std::cout << "\t* The Number of distinct (timestamp,max_pl,rawid) is " << full_map.size() << std::endl;
	std::cout << "\t* The Number of (timestamp,max_pl,rawid) shared by 2+ events is " << dup_group_num << std::endl;
}

void pick_timestamp_duplicates(std::vector<Entry>& all_entries, std::vector<Entry>& ts_dup_list) {

	std::map<int, std::vector<int>> ts_map;
	for (size_t i = 0; i < all_entries.size(); i++) {
		ts_map[all_entries[i].timestamp].push_back((int)i);
	}

	for (auto& kv : ts_map) {
		std::vector<int>& idx = kv.second;
		if (idx.size() < 2) continue;// timestampが他と重複していないchainは対象外

		for (int i : idx) {
			Entry& e = all_entries[i];
			std::cout << "\n [TIMESTAMP DUP] timestamp = " << e.timestamp << " , groupid = " << e.groupid
				<< " , max_pl = " << e.max_pl << " , rawid(at max_pl) = " << e.max_pl_rawid
				<< " , bm_range_mom = " << e.bm_mom << " , ecc_mcs_mom = " << e.ecc_mcs_mom << std::endl;
			ts_dup_list.push_back(e);
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
			// timestamp,groupid
			os << e.timestamp << "\t" << e.groupid << std::endl;
		}
		else {
			// timestamp,groupid,max_pl,rawid(at max_pl)
			os << e.timestamp << "\t" << e.groupid << "\t" << e.max_pl << "\t" << e.max_pl_rawid << std::endl;
		}
	}
}

void write_timestamp_duplicate_list(std::ostream& os, std::vector<Entry>& ts_dup_list) {
	// timestamp,eventid,max_pl,rawid,BM_range_mom,ECC_MCS_mom
	for (auto& e : ts_dup_list) {
		os << e.timestamp << "\t" << e.groupid << "\t" << e.max_pl << "\t" << e.max_pl_rawid
			<< "\t" << e.bm_mom << "\t" << e.ecc_mcs_mom << std::endl;
	}
}
