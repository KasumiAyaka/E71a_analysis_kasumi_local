// 2025/11/27
// Add_partner_stop_flg
// kasumi
// * stop_flg
// 0: penetrate, 2: ECC_stop

#define _CRT_SECURE_NO_WARNINGS
#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#pragma comment(lib, "VxxReader.lib")
#include "VxxReader.h"
#include <FILE_structure.hpp>
#include <functions.hpp>
# include <iostream>
# include <random>

#include <list>
#include <cassert>
#include <filesystem>

class Point {
public:
	double x, y, z;
};
class Fiducial_Area_Local {
public:
	int pl;
	Point p[2];
	//double x0, y0, z0, x1, y1, z1;
};


class stop_track {
public:
	int64_t chainid, groupid;
	int  nseg, npl, pl0, pl1, ph, rawid, pid;
	double ax, ay, x, y, z;
	// ph-->pid
	int stoppl;
	int unixtime;

};
class track_pair {
public:
	int eventid;
	double x, y, z, md, oa;
	double dz;
	stop_track t[2];
};
class track_multi {
public:
	int eventid, pl;
	double x, y, z;
	std::vector< std::pair<double, stop_track>> trk;
	std::vector<track_pair>pair;
	int unixtime;
	double dz;
};
struct vtx_point {
	double x, y, z;
};
struct DivideParam {
	int md, pnum, tnum;
};
struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}
auto start = std::chrono::system_clock::now();//for measure working time

void MeasureProcessingTime(std::chrono::system_clock::time_point& start, std::chrono::system_clock::time_point& end) {
	auto dur = end - start;        // 要した時間を計算
	auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
	// 要した時間をミリ秒（1/1000秒）に変換して表示
	std::cout << msec << " milli sec \n";
	if (msec / 1000 < 60) {
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

double minimum_distance_fixed(matrix_3D::vector_3D pos0, matrix_3D::vector_3D pos1, matrix_3D::vector_3D dir0, matrix_3D::vector_3D dir1, double z_range[2], double extra[2], double refz) {
	double extra0_distance, extra1_distance, delta;
	matrix_3D::vector_3D pos;
	pos.x = pos1.x - pos0.x;
	pos.y = pos1.y - pos0.y;
	pos.z = pos1.z - pos0.z;
	//ほぼ平行な場合
	if (opening_angle(dir0, dir1) < 0.0001) {
		extra0_distance = (pos1.z + pos0.z) / 2 - pos0.z;
		extra1_distance = (pos1.z + pos0.z) / 2 - pos1.z;
	}
	else {
		delta = dot(dir0, dir0) * dot(dir1, dir1) - pow(dot(dir0, dir1), 2.);
		extra0_distance = (+1 * dot(pos, dir0) * dot(dir1, dir1) - dot(dir0, dir1) * dot(pos, dir1)) / delta;
		extra1_distance = (-1 * dot(pos, dir1) * dot(dir0, dir0) + dot(dir0, dir1) * dot(pos, dir0)) / delta;
	}
	//range[0]:小,range[1]:大
	if (z_range[0] > z_range[1]) {
		double tmp_d = z_range[0];
		z_range[0] = z_range[1];
		z_range[1] = tmp_d;
	}

	matrix_3D::vector_3D extra0 = addition(pos0, const_multiple(dir0, extra0_distance));
	matrix_3D::vector_3D extra1 = addition(pos1, const_multiple(dir1, extra1_distance));

	if (extra0.z < refz + z_range[0] || extra1.z < refz + z_range[0]) {//2025/8/20 fixed
		extra0_distance = refz - pos0.z + z_range[0];
		extra1_distance = refz - pos1.z + z_range[0];
	}
	else if (extra0.z > refz + z_range[1] || extra1.z > refz + z_range[1]) {//2025/8/20 fixed
		extra0_distance = refz - pos0.z + z_range[1];
		extra1_distance = refz - pos1.z + z_range[1];
	}

	extra[0] = extra0_distance;
	extra[1] = extra1_distance;
	extra0 = addition(pos0, const_multiple(dir0, extra0_distance));
	extra1 = addition(pos1, const_multiple(dir1, extra1_distance));

	return distance(extra0, extra1);

}


void add_stop_flg(std::vector<Momentum_recon::Event_information>& momch, std::map<int, std::vector<Fiducial_Area_Local>>& area, double edge_cut, int ex_pl_max);
int divide_edge_out(Momentum_recon::Mom_basetrack& btrk, std::map<int, std::vector<Fiducial_Area_Local>>& area, double edge_cut, int ex_pl_max);
int divide_edge_out_fwd(Momentum_recon::Mom_basetrack& btrk, std::map<int, std::vector<Fiducial_Area_Local>>& area, double edge_cut, int ex_pl_max);
bool judge_fiducial_area(std::vector<Fiducial_Area_Local>& area, Momentum_recon::Mom_basetrack& b);
bool judge_fiducial_area_fwd(std::vector<Fiducial_Area_Local>& area, Momentum_recon::Mom_basetrack& b);


std::map<int, std::vector<Fiducial_Area_Local>> read_fiducial_Area(std::string filename);
void trans_mfile_cordinate(std::vector<corrmap_3d::align_param2>& param, std::vector<Fiducial_Area_Local>& area);
void trans_base(std::vector<Point*>& p, corrmap_3d::align_param2* param);
void trans_base_all(std::vector < std::pair<Point*, corrmap_3d::align_param2*>>& track_pair);
std::vector <std::pair<Point*, corrmap_3d::align_param2*>>track_affineparam_correspondence(std::vector<Point*>& p, std::vector <corrmap_3d::align_param2>& param);
bool sort_M_base(const mfile0::M_Base& left, const mfile0::M_Base& right) {
	return left.pos < right.pos;
}

int main(int argc, char** argv) {
	if (argc != 5) {
		fprintf(stderr, "usage:prg in.momch I:\\NINJA\\E71a\\ECC5 I:\\E71a\\ECC5\Area0\0 output.momch\n");
		exit(1);
	}
	std::string in_momch = argv[1];// input momch
	std::string file_in_ECC = argv[2];// ECC path
	std::string file_in_area = argv[3];
	std::string out_momch = argv[4];// output momch

	//bool result = std::filesystem::create_directories(out_momch);
	//assert(result);
	//assert(std::filesystem::exists(out_momch));//directryの存在確認
	//assert(std::filesystem::is_directory(out_momch));//指定されたパスがディレクトリを指しているかを確認する

	auto start = std::chrono::system_clock::now();//for measure working time


	// read momch
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);
	
	//corrmap absの読み込み
	std::string file_in_corrmap = file_in_ECC + "\\Area0\\0\\align\\fine\\local\\corrmap-local-abs.lst";
	std::map<int, std::vector<corrmap_3d::align_param>>corrmap = corrmap_3d::read_ali_param_abs(file_in_corrmap, 1);
	std::map<int, std::vector<corrmap_3d::align_param2>>corrmap_dd = corrmap_3d::DelaunayDivide_map(corrmap);
	std::cout << " corrmap size = " << corrmap.size() << std::endl;
	std::cout << " corrmap_dd size = " << corrmap_dd.size() << std::endl;

	//fiducial areaの読み込み
	std::map<int, std::vector<Fiducial_Area_Local>> area = read_fiducial_Area(file_in_area);
	std::cout << " Area size = " << area.size() << std::endl;
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		if (corrmap_dd.count(itr->first) == 0) {
			fprintf(stderr, "corrmap local abs PL%03d not found\n", itr->first);
			exit(1);
		}
		std::vector<corrmap_3d::align_param2> param = corrmap_dd.at(itr->first);
		trans_mfile_cordinate(param, itr->second);
	}


	// add stop flg
	add_stop_flg(momch, area, 0, 4);


	// write out
	Momentum_recon::Write_Event_information_extension(out_momch, momch);

	auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
	MeasureProcessingTime(start, end);

	std::cout << "Finished." << std::endl;
}

void add_stop_flg(std::vector<Momentum_recon::Event_information>& momch, std::map<int, std::vector<Fiducial_Area_Local>>& area, double edge_cut, int ex_pl_max) {

	int vpl;
	int cnt = 0; int stop = 0;
	for (auto& ev : momch) {
		vpl = ev.vertex_pl;
		for (auto& c : ev.chains) {
			if (c.particle_flg == 13) {
				continue;
			}
			c.stop_flg = 2;// ECC stop
			cnt++;

			if (c.base.begin()->pl > vpl) {	// bwd
				if (c.base.rbegin()->pl > 130) {//penetrate
					c.stop_flg = 0;
				}
				else {
					//最上流から4PL外挿,edgeから5mm以内に入ったらedge out
					c.stop_flg = divide_edge_out(*c.base.rbegin(), area, edge_cut, ex_pl_max);
				}

			}
			else { //fwd
				if (c.base.begin()->pl < 4) {//penetrate
					c.stop_flg = 0;
				}
				else {
					//最上流から4PL外挿,edgeから5mm以内に入ったらedge out
					c.stop_flg = divide_edge_out_fwd(*c.base.begin(), area, edge_cut, ex_pl_max);
				}
			}
			
			if (c.stop_flg == 0)stop++;
		}
	}

	std::cout << " # of Stop partner track = " << stop << " / " << cnt << std::endl;

}

std::map<int, std::vector<Fiducial_Area_Local>> read_fiducial_Area(std::string filename) {

	std::ifstream ifs(filename);
	std::multimap<int, Fiducial_Area_Local> fa_multi;
	std::map<int, std::vector<Fiducial_Area_Local>> ret;
	Fiducial_Area_Local fa;
	while (ifs >> fa.pl >> fa.p[0].x >> fa.p[0].y >> fa.p[0].z >> fa.p[1].x >> fa.p[1].y >> fa.p[1].z) {
		fa_multi.insert(std::make_pair(fa.pl, fa));
	}

	int count = 0;
	for (auto itr = fa_multi.begin(); itr != fa_multi.end(); itr++) {
		count = fa_multi.count(itr->first);
		auto range = fa_multi.equal_range(itr->first);
		std::vector<Fiducial_Area_Local> fa_vec;
		for (auto itr2 = range.first; itr2 != range.second; itr2++) {
			fa_vec.push_back(itr2->second);
		}
		ret.insert(std::make_pair(itr->first, fa_vec));
	}

	return ret;

}
void trans_mfile_cordinate(std::vector<corrmap_3d::align_param2>& param, std::vector<Fiducial_Area_Local>& area) {

	std::vector< Point*> p_trans;
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		p_trans.push_back(&(itr->p[0]));
		p_trans.push_back(&(itr->p[1]));
	}
	std::vector <std::pair<Point*, corrmap_3d::align_param2*>> p_trans_map = track_affineparam_correspondence(p_trans, param);
	trans_base_all(p_trans_map);
}

int divide_edge_out(Momentum_recon::Mom_basetrack& btrk, std::map<int, std::vector<Fiducial_Area_Local>>& area, double edge_cut, int ex_pl_max) {


	double xmin, xmax, ymin, ymax;


	int up_pl;
	double up_z, up_x, up_y, up_ax, up_ay, ex_z, ex_x, ex_y;
	int flg = 2;

	up_pl = btrk.pl;// most upstream pl number
	//up_z = z_map.at(up_pl);
	//up_x = itr->basetracks.rbegin()->x;
	//up_y = itr->basetracks.rbegin()->y;
	//up_ax = itr->basetracks.rbegin()->ax;
	//up_ay = itr->basetracks.rbegin()->ay;
	for (int ex_pl = 0; ex_pl <= ex_pl_max; ex_pl++) {
		if (area.count(up_pl + ex_pl) == 0)continue;
		//ex_z = z_map.at(up_pl + ex_pl);
		//ex_x = up_x + up_ax * (ex_z - up_z);
		//ex_y = up_y + up_ay * (ex_z - up_z);
		//std::cout << "ok?" << std::endl;

		if (!judge_fiducial_area(area.at(up_pl + ex_pl), btrk)) {
			//true でArea内　falseでarea外,if (!条件) : 条件が「偽(false)」なら実行

			flg = 0;//sideout
		}
	}
	return flg;
}
bool judge_fiducial_area(std::vector<Fiducial_Area_Local>& area, Momentum_recon::Mom_basetrack &b) {

	std::map<double, Point> point_map;
	double ex_x, ex_y, dist;
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		// extraporate upstram direction
		ex_x = b.x + b.ax * (itr->p[0].z - b.z);
		ex_y = b.y + b.ay * (itr->p[0].z - b.z);
		dist = pow(ex_x - itr->p[0].x, 2) + pow(ex_y - itr->p[0].y, 2);
		point_map.insert(std::make_pair(dist, itr->p[0]));
	}

	//std::cout << "(" << ex_x << ", " << ex_y << ") " << std::endl;
	//if (ex_x < 5000 || ex_x>250000 - 5000) return false;
	//if (ex_y < 5000 || ex_y>250000 - 5000) return false;
	//return true;

	/**/
	//外挿先から距離の一番近い点のz座標を使用
	double z = point_map.begin()->second.z;
	double x = b.x + b.ax * (z - b.z);
	double y = b.y + b.ay * (z - b.z);


	//true でArea内　falseでarea外

	//点(x,y)からx軸性の方向に直線を引き、その直線と多角形の辺が何回交わるか。
	//下から上に交わったときwn+1
	//上から下に交わったときwn-1
	int wn = 0;
	double vt;
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		// 上向きの辺、下向きの辺によって処理が分かれる。
	// 上向きの辺。点Pがy軸方向について、始点と終点の間にある。ただし、終点は含まない。(ルール1)
		if (itr->p[0].y <= y && itr->p[1].y > y) {
			// 辺は点pよりも右側にある。ただし、重ならない。(ルール4)
			// 辺が点pと同じ高さになる位置を特定し、その時のxの値と点pのxの値を比較する。
			vt = (y - itr->p[0].y) / (itr->p[1].y - itr->p[0].y);
			if (x < itr->p[0].x + vt * (itr->p[1].x - itr->p[0].x)) {
				++wn;  //ここが重要。上向きの辺と交差した場合は+1
			}
		}
		// 下向きの辺。点Pがy軸方向について、始点と終点の間にある。ただし、始点は含まない。(ルール2)
		else if (itr->p[0].y > y && itr->p[1].y <= y) {
			// 辺は点pよりも右側にある。ただし、重ならない。(ルール4)
			// 辺が点pと同じ高さになる位置を特定し、その時のxの値と点pのxの値を比較する。
			vt = (y - itr->p[0].y) / (itr->p[1].y - itr->p[0].y);
			if (x < itr->p[0].x + vt * (itr->p[1].x - itr->p[0].x)) {
				--wn;  //ここが重要。下向きの辺と交差した場合は-1
			}
		}
	}
	if (wn >= 1)return true;
	return false;
	//*/
}

int divide_edge_out_fwd(Momentum_recon::Mom_basetrack& btrk, std::map<int, std::vector<Fiducial_Area_Local>>& area, double edge_cut, int ex_pl_max) {


	double xmin, xmax, ymin, ymax;


	int up_pl;
	double up_z, up_x, up_y, up_ax, up_ay, ex_z, ex_x, ex_y;
	int flg = 2;

	up_pl = btrk.pl;// most downsttream pl number
	//up_z = z_map.at(up_pl);
	//up_x = itr->basetracks.rbegin()->x;
	//up_y = itr->basetracks.rbegin()->y;
	//up_ax = itr->basetracks.rbegin()->ax;
	//up_ay = itr->basetracks.rbegin()->ay;
	for (int ex_pl = 0; ex_pl <= ex_pl_max; ex_pl++) {
		if (area.count(up_pl - ex_pl) == 0)continue;
		//ex_z = z_map.at(up_pl + ex_pl);
		//ex_x = up_x + up_ax * (ex_z - up_z);
		//ex_y = up_y + up_ay * (ex_z - up_z);

		if (!judge_fiducial_area_fwd(area.at(up_pl - ex_pl), btrk)) {
			//true でarea外?
			flg = 0;//sideout
		}
	}
	return flg;
}
bool judge_fiducial_area_fwd(std::vector<Fiducial_Area_Local>& area, Momentum_recon::Mom_basetrack& b) {

	std::map<double, Point> point_map;
	double ex_x, ex_y, dist;
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		// extraporate upstram direction
		ex_x = b.x - b.ax * (-itr->p[0].z + b.z);
		ex_y = b.y - b.ay * (-itr->p[0].z + b.z);
		dist = pow(ex_x - itr->p[0].x, 2) + pow(ex_y - itr->p[0].y, 2);
		point_map.insert(std::make_pair(dist, itr->p[0]));
	}
	//if (ex_x < 5000 || ex_x>250000 - 5000) return false;
	//if (ex_y < 5000 || ex_y>250000 - 5000) return false;
	//return true;

	/**/
	//外挿先から距離の一番近い点のz座標を使用
	double z = point_map.begin()->second.z;
	double x = b.x - b.ax * (-z + b.z);
	double y = b.y - b.ay * (-z + b.z);


	//true でArea内　falseでarea外

	//点(x,y)からx軸性の方向に直線を引き、その直線と多角形の辺が何回交わるか。
	//下から上に交わったときwn+1
	//上から下に交わったときwn-1
	int wn = 0;
	double vt;
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		// 上向きの辺、下向きの辺によって処理が分かれる。
	// 上向きの辺。点Pがy軸方向について、始点と終点の間にある。ただし、終点は含まない。(ルール1)
		if (itr->p[0].y <= y && itr->p[1].y > y) {
			// 辺は点pよりも右側にある。ただし、重ならない。(ルール4)
			// 辺が点pと同じ高さになる位置を特定し、その時のxの値と点pのxの値を比較する。
			vt = (y - itr->p[0].y) / (itr->p[1].y - itr->p[0].y);
			if (x < itr->p[0].x + vt * (itr->p[1].x - itr->p[0].x)) {
				++wn;  //ここが重要。上向きの辺と交差した場合は+1
			}
		}
		// 下向きの辺。点Pがy軸方向について、始点と終点の間にある。ただし、始点は含まない。(ルール2)
		else if (itr->p[0].y > y && itr->p[1].y <= y) {
			// 辺は点pよりも右側にある。ただし、重ならない。(ルール4)
			// 辺が点pと同じ高さになる位置を特定し、その時のxの値と点pのxの値を比較する。
			vt = (y - itr->p[0].y) / (itr->p[1].y - itr->p[0].y);
			if (x < itr->p[0].x + vt * (itr->p[1].x - itr->p[0].x)) {
				--wn;  //ここが重要。下向きの辺と交差した場合は-1
			}
		}
	}
	if (wn >= 1)return true;
	//std::cout << "(" << ex_x << ", " << ex_y << "), (" << b.x << ", " << b.y << ", " << b.ax << ", " << b.ay << ")" << std::endl;
	return false;
	//*/
}






//mfile0::M_Base
//basetrack-alignment mapの対応
double select_triangle_vale(corrmap_3d::align_param2* param, Point* p) {
	double x, y;
	double dist = 0;
	x = (param->corr_p[0]->x + param->corr_p[1]->x + param->corr_p[2]->x) / 3;
	y = (param->corr_p[0]->y + param->corr_p[1]->y + param->corr_p[2]->y) / 3;
	dist = (p->x - x) * (p->x - x) + (p->y - y) * (p->y - y);
	return dist;
}
corrmap_3d::align_param2* search_param(std::vector<corrmap_3d::align_param*>& param, Point* p, std::multimap<int, corrmap_3d::align_param2*>& triangles) {
	//三角形内部
	//最近接三角形
	double dist = 0;
	std::map<double, corrmap_3d::align_param* > dist_map;
	//align_paramを近い順にsort
	for (auto itr = param.begin(); itr != param.end(); itr++) {
		dist = ((*itr)->x - p->x) * ((*itr)->x - p->x) + ((*itr)->y - p->y) * ((*itr)->y - p->y);
		dist_map.insert(std::make_pair(dist, (*itr)));
	}

	double sign[3];
	bool flg = false;
	int id;

	corrmap_3d::align_param2* ret = triangles.begin()->second;
	for (auto itr = dist_map.begin(); itr != dist_map.end(); itr++) {
		if (itr != dist_map.begin())continue;


		//corrmapのID
		id = itr->second->id;
		if (triangles.count(id) == 0) {
			fprintf(stderr, "alignment triangle ID=%d not found\n", id);
			exit(1);
		}
		//idの属する三角形を探索
		auto range = triangles.equal_range(id);
		for (auto itr2 = range.first; itr2 != range.second; itr2++) {
			sign[0] = (itr2->second->corr_p[1]->x - itr2->second->corr_p[0]->x) * (p->y - itr2->second->corr_p[1]->y) - (itr2->second->corr_p[1]->y - itr2->second->corr_p[0]->y) * (p->x - itr2->second->corr_p[1]->x);
			sign[1] = (itr2->second->corr_p[2]->x - itr2->second->corr_p[1]->x) * (p->y - itr2->second->corr_p[2]->y) - (itr2->second->corr_p[2]->y - itr2->second->corr_p[1]->y) * (p->x - itr2->second->corr_p[2]->x);
			sign[2] = (itr2->second->corr_p[0]->x - itr2->second->corr_p[2]->x) * (p->y - itr2->second->corr_p[0]->y) - (itr2->second->corr_p[0]->y - itr2->second->corr_p[2]->y) * (p->x - itr2->second->corr_p[0]->x);
			//printf("point %.lf,%.1lf\n", base.x, base.y);
			//printf("triangle %.1lf %.1lf\n", itr2->second->corr_p[0]->x, itr2->second->corr_p[0]->y);
			//printf("triangle %.1lf %.1lf\n", itr2->second->corr_p[1]->x, itr2->second->corr_p[1]->y);
			//printf("triangle %.1lf %.1lf\n", itr2->second->corr_p[2]->x, itr2->second->corr_p[2]->y);
			//printf("sign %.1lf %1.lf %.1lf\n", sign[0], sign[1], sign[2]);
			//printf("  signbit %d %d %d\n", std::signbit(sign[0]), std::signbit(sign[1]), std::signbit(sign[2]));
			//printf("n signbit %d %d %d\n", !std::signbit(sign[0]), !std::signbit(sign[1]), !std::signbit(sign[2]));
			//printf("judge %d\n", (std::signbit(sign[0]) && std::signbit(sign[1]) && std::signbit(sign[2])) || (!std::signbit(sign[0]) && !std::signbit(sign[1]) && !std::signbit(sign[2])));
			//printf("\n");

			//符号が3つとも一致でtrue
			if ((std::signbit(sign[0]) && std::signbit(sign[1]) && std::signbit(sign[2])) || (!std::signbit(sign[0]) && !std::signbit(sign[1]) && !std::signbit(sign[2]))) {
				ret = itr2->second;
				flg = true;
				break;
			}
		}
		if (flg)break;
	}
	if (flg) {
		//printf("point in trianlge\n");
		return ret;
	}

	//distが最小になるcorrmapをとってくる
	dist = -1;
	for (auto itr = dist_map.begin(); itr != dist_map.end(); itr++) {
		//corrmapのID
		id = itr->second->id;
		if (triangles.count(id) == 0) {
			fprintf(stderr, "alignment triangle ID=%d not found\n", id);
			exit(1);
		}
		//idの属する三角形を探索
		auto range = triangles.equal_range(id);
		for (auto itr2 = range.first; itr2 != range.second; itr2++) {
			if (dist<0 || dist>select_triangle_vale(itr2->second, p)) {
				dist = select_triangle_vale(itr2->second, p);
				ret = itr2->second;
			}
		}
	}
	//printf("point not in trianlge\n");
	return ret;
}
std::vector <std::pair<Point*, corrmap_3d::align_param2*>>track_affineparam_correspondence(std::vector<Point*>& p, std::vector <corrmap_3d::align_param2>& param) {

	//local alignの視野中心を取り出して、位置でhash
	//local alignの視野中心の作るdelaunay三角形をmapで対応

	std::map<int, corrmap_3d::align_param*> view_center;
	std::multimap<int, corrmap_3d::align_param2*>triangles;
	double xmin = 999999, ymin = 999999, hash = 2000;
	for (auto itr = param.begin(); itr != param.end(); itr++) {
		for (int i = 0; i < 3; i++) {
			view_center.insert(std::make_pair(itr->corr_p[i]->id, (itr->corr_p[i])));
			triangles.insert(std::make_pair(itr->corr_p[i]->id, &(*itr)));
			xmin = std::min(itr->corr_p[i]->x, xmin);
			ymin = std::min(itr->corr_p[i]->y, ymin);
		}
	}
	std::multimap<std::pair<int, int>, corrmap_3d::align_param*> view_center_hash;
	std::pair<int, int>id;
	for (auto itr = view_center.begin(); itr != view_center.end(); itr++) {
		id.first = int((itr->second->x - xmin) / hash);
		id.second = int((itr->second->y - ymin) / hash);
		view_center_hash.insert(std::make_pair(id, itr->second));
	}

	std::vector < std::pair<Point*, corrmap_3d::align_param2*>> ret;
	std::vector<corrmap_3d::align_param*> param_cand;
	int loop = 0, ix, iy, count = 0;
	for (auto itr = p.begin(); itr != p.end(); itr++) {
		if (count % 100000 == 0) {
			printf("\r search correspond triangles %d/%d(%4.1lf%%)", count, p.size(), count * 100. / p.size());
		}
		count++;
		ix = ((*itr)->x - xmin) / hash;
		iy = ((*itr)->y - ymin) / hash;
		loop = 1;
		while (true) {
			param_cand.clear();
			for (int iix = ix - loop; iix <= ix + loop; iix++) {
				for (int iiy = iy - loop; iiy <= iy + loop; iiy++) {
					id.first = iix;
					id.second = iiy;
					if (view_center_hash.count(id) != 0) {
						auto range = view_center_hash.equal_range(id);
						for (auto res = range.first; res != range.second; res++) {
							param_cand.push_back(res->second);
						}
					}
				}
			}
			if (param_cand.size() > 2)break;
			loop++;
		}
		corrmap_3d::align_param2* param2 = search_param(param_cand, *itr, triangles);
		ret.push_back(std::make_pair((*itr), param2));
	}
	printf("\r search correspond triangles %d/%d(%4.1lf%%)\n", count, p.size(), count * 100. / p.size());

	return ret;
}
//変換 zshrink補正-->9para変換
void trans_base_all(std::vector < std::pair<Point*, corrmap_3d::align_param2*>>& track_pair) {
	std::map<std::tuple<int, int, int>, corrmap_3d::align_param2*> param_map;
	std::multimap<std::tuple<int, int, int>, Point*>base_map;
	std::tuple<int, int, int>id;
	//三角形ごとにbasetrackをまとめる
	for (auto itr = track_pair.begin(); itr != track_pair.end(); itr++) {
		std::get<0>(id) = itr->second->corr_p[0]->id;
		std::get<1>(id) = itr->second->corr_p[1]->id;
		std::get<2>(id) = itr->second->corr_p[2]->id;
		param_map.insert(std::make_pair(id, itr->second));
		base_map.insert(std::make_pair(id, itr->first));
	}


	//ここで三角形ごとに変換
	int count = 0;
	std::vector<Point*> t_base;
	for (auto itr = param_map.begin(); itr != param_map.end(); itr++) {
		if (count % 1000 == 0) {
			printf("\r basetrack trans num of triangles %d/%d(%4.1lf%%)", count, param_map.size(), count * 100. / param_map.size());
		}
		count++;

		t_base.clear();

		if (base_map.count(itr->first) == 0)continue;
		auto range = base_map.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			t_base.push_back(res->second);
		}
		trans_base(t_base, itr->second);

	}
	printf("\r basetrack trans num of triangles %d/%d(%4.1lf%%)\n", count, param_map.size(), count * 100. / param_map.size());

}
void trans_base(std::vector<Point*>& p, corrmap_3d::align_param2* param) {

	matrix_3D::matrix_33 x_rot_mat(0, param->x_rot), y_rot_mat(1, param->y_rot), z_rot_mat(2, param->z_rot), all_trans(0, 0), shear_mat(0, 0), shrink_mat(0, 0);

	shrink_mat.val[0][0] *= param->x_shrink;
	shrink_mat.val[1][1] *= param->y_shrink;
	//shrink_mat.val[2][2] *= param->z_shrink;
	shear_mat.val[0][1] = param->yx_shear;
	//shear_mat.val[0][2] = param->zx_shear;
	//shear_mat.val[1][2] = param->zy_shear;

	matrix_3D::vector_3D shift, center;
	center.x = param->x;
	center.y = param->y;
	center.z = param->z;
	shift.x = param->dx;
	shift.y = param->dy;
	shift.z = param->dz;

	all_trans.matrix_multiplication(shear_mat);
	all_trans.matrix_multiplication(shrink_mat);
	all_trans.matrix_multiplication(z_rot_mat);
	all_trans.matrix_multiplication(y_rot_mat);
	all_trans.matrix_multiplication(x_rot_mat);

	//all_trans.Print();
	matrix_3D::vector_3D base_p0;
	double base_thick = 210;
	for (auto itr = p.begin(); itr != p.end(); itr++) {
		base_p0.x = (*itr)->x;
		base_p0.y = (*itr)->y;
		base_p0.z = param->z;

		//base_p1.x = (*itr)->x + (*itr)->ax*base_thick;
		//base_p1.y = (*itr)->y + (*itr)->ay*base_thick;
		////角度shrink分はここでかける
		//base_p1.z = param->z + base_thick / param->z_shrink;

		//視野中心を原点に移動
		//base_p0 = matrix_3D::addition(base_p0, matrix_3D::const_multiple(center, -1));
		//base_p1 = matrix_3D::addition(base_p1, matrix_3D::const_multiple(center, -1));

		//変換の実行
		base_p0.matrix_multiplication(all_trans);
		base_p0 = matrix_3D::addition(base_p0, shift);
		//base_p1.matrix_multiplication(all_trans);
		//base_p1 = matrix_3D::addition(base_p1, shift);

		//原点をもとに戻す
		//base_p0 = matrix_3D::addition(base_p0, center);
		//base_p1 = matrix_3D::addition(base_p1, center);

		(*itr)->x = base_p0.x;
		(*itr)->y = base_p0.y;
		(*itr)->z = base_p0.z;

		//printf("ax:%.4lf --> %.4lf\n", (*itr)->ax, (base_p0.x - base_p1.x) / (base_p0.z - base_p1.z));
		//printf("ay:%.4lf --> %.4lf\n", (*itr)->ay, (base_p0.y - base_p1.y) / (base_p0.z - base_p1.z));

		//(*itr)->ax = (base_p0.x - base_p1.x) / (base_p0.z - base_p1.z) + param->zx_shear;
		//(*itr)->ay = (base_p0.y - base_p1.y) / (base_p0.z - base_p1.z) + param->zy_shear;

	}
}
