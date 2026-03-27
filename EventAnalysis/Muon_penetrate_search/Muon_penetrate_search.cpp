// 2025/1/29
// kasumi
// This code is written kasum based on "Crazy_penetrate_search_for_muon_v2".

// 2025/01/25
// kasumi
// ma,od,the PL range of searchを引数で与えるように変更

//usage
//M:\data\NINJA\prg\Crazy_penetrate_search_for_muon_v2.exe muon_add.all T:\Udrive116\Linklet1\base_use.bin T:\NINJA\E71a\ECC1 out.txt

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <set>
#include <map>
#include <chrono>
#pragma comment(lib, "VxxReader.lib")
#include "VxxReader.h"
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.h>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <omp.h>
namespace mfile0 {
	bool operator<(const mfile0::M_Chain& left, const mfile0::M_Chain& right) {
		if (left.basetracks.begin()->group_id == right.basetracks.begin()->group_id)return left.chain_id < right.chain_id;
		return left.basetracks.begin()->group_id < right.basetracks.begin()->group_id;
	}
}

class basetrack_minimum {
public:
	int pl, rawid, ph;
	float ax, ay, x, y;
};
class basetrack_mfile :public basetrack_minimum {
public:
	float z;
};

class penetrate_check_area {
	int pl;
	double x_min, x_max, y_min, y_max;
	double ax_min, ax_max, ay_min, ay_max;
public:
	void Set_area_position(double in_xmin, double in_xmax, double in_ymin, double in_ymax);
	void Set_area_angle(double in_xmin, double in_xmax, double in_ymin, double in_ymax);
	void Set_pl(int in_pl);
	bool judage_inside(int pl, double x, double y, double ax, double ay);
};
void penetrate_check_area::Set_area_position(double in_xmin, double in_xmax, double in_ymin, double in_ymax) {
	x_min = in_xmin;
	x_max = in_xmax;
	y_min = in_ymin;
	y_max = in_ymax;
}
void penetrate_check_area::Set_area_angle(double in_xmin, double in_xmax, double in_ymin, double in_ymax) {
	ax_min = in_xmin;
	ax_max = in_xmax;
	ay_min = in_ymin;
	ay_max = in_ymax;
}
void penetrate_check_area::Set_pl(int in_pl) {
	pl = in_pl;
}
bool penetrate_check_area::judage_inside(int in_pl, double x, double y, double ax, double ay) {
	if (in_pl != pl)return false;

	if (x < x_min)return false;
	if (x_max < x)return false;
	if (y < y_min)return false;
	if (y_max < y)return false;

	if (ax < ax_min)return false;
	if (ax_max < ax)return false;
	if (ay < ay_min)return false;
	if (ay_max < ay)return false;

	return true;
}



bool sort_basetrack_pl(const basetrack_minimum& left, const basetrack_minimum& right) {
	if (left.pl == right.pl)return left.rawid < right.rawid;
	return left.pl < right.pl;
}

std::map<int, std::vector< mfile0::M_Base>> read_base(std::string filename);
std::map<mfile0::M_Chain, std::vector<mfile0::M_Base>>Search_penetrate_cand(std::vector<mfile0::M_Chain>& chain, std::map<int, std::vector< mfile0::M_Base>>& base, double threshold_oa, double threshold_md, int thr_upl, int thr_dpl);


void output(std::string filename, std::map<mfile0::M_Chain, std::vector<mfile0::M_Base>>& pene_cand);
void trans_local(std::map<int, std::vector<mfile0::M_Base>>& base_map_single, std::map<int, std::vector<corrmap_3d::align_param2>>& corr);
int use_thread(double ratio, bool output);

int main(int argc, char** argv) {

	if (argc != 5 && argc != 9) {
		fprintf(stderr, "usage1\n file_in_mfile file_in_base file_in_ECC file_out\n");
		fprintf(stderr, "usage2\n file_in_mfile file_in_base file_in_ECC file_out thr_oa thr_md thr_d(upl) thr_d(dpl)\n");
		fprintf(stderr, " file_in_mfile file_in_base file_in_ECC file_out 0.2    200     4         -1\n");
		fprintf(stderr, "I:\\NINJA\\E71a\\work\\kasumi\\ECC\\MuonAnalysis\\x64\\Release\\Crazy_penetrate_search_for_muon_v2.exe muon.all T:\\Udrive116\\Linklet1\\base_use.bin T:\\NINJA\\E71a\\ECC1 output.txt\n");
		exit(1);
	}


	std::string file_in_muon_mfile = argv[1];
	std::string file_in_all_base = argv[2];
	std::string file_in_ECC = argv[3];
	std::string file_out = argv[4];

	double thr_oa = 0.2;
	double thr_md = 200;
	int pl_ustream_thr = 4;
	int pl_dstream_thr = -1;

	if (argc == 9) {
		thr_oa = std::stod(argv[5]);//0.2
		thr_md = std::stod(argv[6]);//200
		pl_ustream_thr = std::stod(argv[7]);//4
		pl_dstream_thr = std::stod(argv[8]);//-1
	}

	//double hash_size_angle = 0.1;
	//double hash_size_position = 2000;


	mfile0::Mfile m;
	mfile0::read_mfile(file_in_muon_mfile, m);

	std::map<int, std::vector< mfile0::M_Base>>  base = read_base(file_in_all_base);

	std::string file_in_corrmap = file_in_ECC + "\\Area0\\0\\align\\fine\\local\\corrmap-local-abs.lst";
	std::map<int, std::vector<corrmap_3d::align_param>>corrmap = corrmap_3d::read_ali_param_abs(file_in_corrmap, 1);
	std::map<int, std::vector<corrmap_3d::align_param2>>corrmap_dd = corrmap_3d::DelaunayDivide_map(corrmap);

	//base trans
	trans_local(base, corrmap_dd);

	//oaとmdの閾値の決定
	//double thr_oa = 0.2, thr_md = 200;
	//muon飛跡から探索すべき領域の決定
	//std::map<mfile0::M_Chain, std::vector<penetrate_check_area >>muon_area = Set_check_area(m.chains, z_map, thr_oa, thr_md);
	//領域内の飛跡をpickup
	std::map<mfile0::M_Chain, std::vector<mfile0::M_Base >>pene_cand = Search_penetrate_cand(m.chains, base, thr_oa, thr_md, pl_dstream_thr, pl_ustream_thr);

	//md,oaで探索
	//for (auto itr = pene_cand.begin(); itr != pene_cand.end(); itr++) {
	//	itr->second = judege_penetrate(itr->first, itr->second, thr_oa, thr_md);
	//}

	output(file_out, pene_cand);

}

int64_t Basetrack_num(std::string filename) {
	std::ifstream ifs(filename, std::ios::binary);
	//filesize取得
	ifs.seekg(0, std::ios::end);
	int64_t eofpos = ifs.tellg();
	ifs.clear();
	ifs.seekg(0, std::ios::beg);
	int64_t begpos = ifs.tellg();
	int64_t nowpos = ifs.tellg();
	int64_t size2 = eofpos - begpos;
	return size2 / sizeof(basetrack_minimum);
}
std::map<int, std::vector< mfile0::M_Base>> read_base(std::string filename) {
	std::map<int, std::vector< mfile0::M_Base>> ret;
	std::multimap<int, basetrack_minimum> base_multimap;
	int64_t track_num = Basetrack_num(filename);
	//ret.reserve(track_num);

	std::ifstream ifs(filename, std::ios::binary);
	//filesize取得
	ifs.seekg(0, std::ios::end);
	int64_t eofpos = ifs.tellg();
	ifs.clear();
	ifs.seekg(0, std::ios::beg);
	int64_t begpos = ifs.tellg();
	int64_t nowpos = ifs.tellg();
	int64_t size2 = eofpos - begpos;
	int64_t GB = size2 / (1000 * 1000 * 1000);
	int64_t MB = (size2 - GB * 1000 * 1000 * 1000) / (1000 * 1000);
	int64_t KB = (size2 - GB * 1000 * 1000 * 1000 - MB * 1000 * 1000) / (1000);
	if (GB > 0) {
		std::cout << "FILE size :" << GB << "." << MB << " [GB]" << std::endl;
	}
	else {
		std::cout << "FILE size :" << MB << "." << KB << " [MB]" << std::endl;
	}
	int64_t count = 0;
	basetrack_minimum b;
	while (ifs.read((char*)&b, sizeof(basetrack_minimum))) {
		if (count % 10000 == 0) {
			nowpos = ifs.tellg();
			auto size1 = nowpos - begpos;
			std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
		}
		count++;
		//ここ
		base_multimap.insert(std::make_pair(b.pl, b));
		//ret.emplace_back(b);
	}
	auto size1 = eofpos - begpos;
	std::cerr << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%" << std::endl;;
	if (count == 0) {
		fprintf(stderr, "%s no basetrack!\n", filename.c_str());
		exit(1);
	}


	for (auto itr = base_multimap.begin(); itr != base_multimap.end(); itr++) {
		std::vector< mfile0::M_Base>base_v;
		int count = base_multimap.count(itr->first);
		base_v.reserve(count);
		auto range = base_multimap.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			mfile0::M_Base b;
			b.pos = res->second.pl * 10 + 1;
			b.rawid = res->second.rawid;
			b.group_id = 0;
			b.ph = res->second.ph;

			b.ax = res->second.ax;
			b.ay = res->second.ay;
			b.x = res->second.x;
			b.y = res->second.y;
			b.z = 0;

			b.flg_d[0] = 0;
			b.flg_d[1] = 0;
			b.flg_i[0] = 0;
			b.flg_i[1] = 0;
			b.flg_i[2] = 0;
			b.flg_i[3] = 0;

			base_v.push_back(b);
		}
		ret.insert(std::make_pair(itr->first, base_v));

		itr = std::next(itr, count - 1);
	}


	return ret;

}

std::map<mfile0::M_Chain, std::vector<mfile0::M_Base>>Search_penetrate_cand(std::vector<mfile0::M_Chain>& chain, std::map<int, std::vector<mfile0::M_Base>>& base, double threshold_oa, double threshold_md, int thr_upl, int thr_dpl) {
	std::map<mfile0::M_Chain, std::vector<mfile0::M_Base>> ret;
	matrix_3D::vector_3D pos0, pos1, dir0, dir1;

	std::chrono::system_clock::time_point  start, end; // 型は auto で可
	start = std::chrono::system_clock::now(); // 計測開始時間

	for (auto& c : chain) {
		std::vector<mfile0::M_Base> connect_cand;
		//最上流basetrack pickup
		int mu_pl = c.basetracks.rbegin()->pos / 10;
		int mu_rawid = c.basetracks.rbegin()->rawid;
		double mu_z = c.basetracks.rbegin()->z;
		double mu_x = c.basetracks.rbegin()->x;
		double mu_y = c.basetracks.rbegin()->y;

		//muonの角度は上流の平均
		double mu_ax = 0, mu_ay = 0;
		int count = 0;
		for (int i = c.basetracks.size() - 1; i >= std::max((int)c.basetracks.size() - 6, 0); i--) {
			mu_ax += c.basetracks[i].ax;
			mu_ay += c.basetracks[i].ay;
			count++;
		}
		mu_ax = mu_ax / count;
		mu_ay = mu_ay / count;

		pos0.x = mu_x;
		pos0.y = mu_y;
		pos0.z = mu_z;
		dir0.x = mu_ax;
		dir0.y = mu_ay;
		dir0.z = 1;

		//PL番号 -1〜4の間で探索
		for (int i_pl = thr_dpl; i_pl <= thr_upl; i_pl++) {
			int t_pl = mu_pl + i_pl;
			//PL3-PL133の範囲外
			if (base.count(t_pl) == 0)continue;
			auto base_v = base.at(t_pl);
			for (auto& b : base_v) {
				if (mu_pl == b.pos / 10 && mu_rawid == b.rawid)continue;
				pos1.x = b.x;
				pos1.y = b.y;
				pos1.z = b.z;
				dir1.x = b.ax;
				dir1.y = b.ay;
				dir1.z = 1;
				double extra[2], md, oa, z_range[2] = { pos1.z,pos0.z };
				oa = matrix_3D::opening_angle(dir0, dir1);
				if (oa > threshold_oa)continue;
				md = matrix_3D::minimum_distance(pos0, pos1, dir0, dir1, z_range, extra);
				if (md > threshold_md)continue;
				connect_cand.push_back(b);

			}




		}
		ret.insert(std::make_pair(c, connect_cand));
	}

	end = std::chrono::system_clock::now();  // 計測終了時間
	double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(); //処理に要した時間をミリ秒に変換
	printf("time %.1lf[s]\n", elapsed / 1000);

	return ret;
}
void output(std::string filename, std::map<mfile0::M_Chain, std::vector<mfile0::M_Base>>& pene_cand) {

	std::ofstream ofs(filename);
	int all = pene_cand.size(), count = 0;

	for (auto itr = pene_cand.begin(); itr != pene_cand.end(); itr++) {
		if (count % 100 == 0) {
			fprintf(stderr, "\r write penetrate candidate %d/%d(%4.1lf%%)", count, all, count * 100. / all);
		}
		count++;

		for (int i = 0; i < itr->second.size(); i++) {
			ofs << std::right << std::fixed
				<< std::setw(10) << std::setprecision(0) << itr->first.basetracks.begin()->group_id << " "
				<< std::setw(10) << std::setprecision(0) << i + 1 << " "
				<< std::setw(10) << std::setprecision(0) << itr->second[i].pos / 10 << " "
				<< std::setw(10) << std::setprecision(0) << itr->second[i].rawid << std::endl;
		}
	}
	fprintf(stderr, "\r write penetrate candidate %d/%d(%4.1lf%%)\n", count, all, count * 100. / all);

}

void trans_local(std::map<int, std::vector<mfile0::M_Base>>& base_map_single, std::map<int, std::vector<corrmap_3d::align_param2>>& corr) {
	std::vector<int> all_pl;
	for (auto itr = base_map_single.begin(); itr != base_map_single.end(); itr++) {
		all_pl.push_back(itr->first);
	}

	int all = all_pl.size(), count = 0;
#pragma omp parallel for num_threads(use_thread(1.0,true)) schedule(dynamic,1)
	for (int i = 0; i < all_pl.size(); i++) {
		//for (auto itr = base_map_single.begin(); itr != base_map_single.end(); itr++) {
			//count = base_map_single.count(itr->first);
		int pl = all_pl[i];
		
		if (pl > 47 && pl < 55) {
#pragma omp critical
			{
				printf("PL%03d basetrack tans %3d/%3d\n", pl, count, all);
				count++;
				if (corr.count(pl) == 0) {
					fprintf(stderr, "PL%03d corrmap not found\n", pl);
					exit(1);
				}
				if (base_map_single.count(pl) == 0) {
					fprintf(stderr, "PL%03d basetrack not found\n", pl);
					exit(1);
				}
			}
			std::vector<corrmap_3d::align_param2> param = corr.at(pl);
			auto res = base_map_single.find(pl);
			std::vector< mfile0::M_Base*> base_trans;
			base_trans.reserve(res->second.size());
			for (auto itr = res->second.begin(); itr != res->second.end(); itr++) {
				base_trans.push_back(&(*itr));
			}
			std::vector <std::pair<mfile0::M_Base*, corrmap_3d::align_param2*>> base_trans_map = corrmap_3d::track_affineparam_correspondence(base_trans, param);
			trans_base_all(base_trans_map);
		}
	}

}

int use_thread(double ratio, bool output) {
	if (ratio >= 1) {
		ratio = 1;
	}
	else if (ratio < 0.1) {
		ratio = 0.1;
	}
	int num_all_thread = omp_get_max_threads();
	if (output) {
		printf("max thread = %d\n", num_all_thread);
		printf("ratio      = %3.2lf\n", ratio);
		printf("using... %d thread\n", int(num_all_thread * ratio));
	}
	return (int)(num_all_thread * ratio);
}

