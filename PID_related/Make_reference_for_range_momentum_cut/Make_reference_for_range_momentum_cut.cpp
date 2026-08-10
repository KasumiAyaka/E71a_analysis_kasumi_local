#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
# include <iostream>
# include <random>

#include <list>

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


struct tkey {
	int eid;
	double dz;
	double w;
};
bool operator<(const tkey& lhs, const tkey& rhs) {
	return std::tie(lhs.eid, lhs.w, lhs.dz) < std::tie(rhs.eid, rhs.w, rhs.dz);
}

class stop_track {
public:
	int64_t chainid, groupid;
	int  nseg, npl, pl0, pl1, vph, rawid;
	double ax, ay, x, y, z;
	int stoppl;
	int unixtime;
	double ip;
	int pid;
	double mom, mulikelihood, pliklihoood, weight;
	double vx, vy, vz;//for drbag
	int ntrk;
	int stopflg;
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
	std::vector< std::pair<int, stop_track>> trk;
	std::vector<track_pair>pair;
	int unixtime;
	double dz;
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

void add_stop_flg(std::multimap<tkey, stop_track>& stop, std::map<int, std::vector<Fiducial_Area_Local>>& area, double edge_cut, int ex_pl_max);
int divide_edge_out(stop_track btrk, std::map<int, std::vector<Fiducial_Area_Local>>& area, double edge_cut, int ex_pl_max);
int divide_edge_out_fwd(stop_track btrk, std::map<int, std::vector<Fiducial_Area_Local>>& area, double edge_cut, int ex_pl_max);
bool judge_fiducial_area(std::vector<Fiducial_Area_Local>& area, stop_track b);
bool judge_fiducial_area_fwd(std::vector<Fiducial_Area_Local>& area, stop_track b);


std::map<int, std::vector<Fiducial_Area_Local>> read_fiducial_Area(std::string filename);
void trans_mfile_cordinate(std::vector<corrmap_3d::align_param2>& param, std::vector<Fiducial_Area_Local>& area);
void trans_base(std::vector<Point*>& p, corrmap_3d::align_param2* param);
void trans_base_all(std::vector < std::pair<Point*, corrmap_3d::align_param2*>>& track_pair);
std::vector <std::pair<Point*, corrmap_3d::align_param2*>>track_affineparam_correspondence(std::vector<Point*>& p, std::vector <corrmap_3d::align_param2>& param);
bool sort_M_base(const mfile0::M_Base& left, const mfile0::M_Base& right) {
	return left.pos < right.pos;
}
void read_stop_txt_for_true(std::string in_momch, std::multimap<tkey, stop_track>& tracks);


void read_stop_txt(std::string in_momch, std::multimap<tkey, stop_track>& tracks);
void clustering_2trk_vtx2_ver2(std::multimap<int, stop_track>& tracks, int pl, std::ofstream& ofs, int ecc, int mode);

int main(int argc, char** argv) {
	if (argc < 5) {
		fprintf(stderr, "usage:prg in_mc_momch.txt out.txt #ECC [1:true/0:recon]   [d(um)] [plcut]:If vpl > plcut then erase that event.\n");
		fprintf(stderr, "usage:prg in_mc_momch.txt out.txt #ECC [1:true/0:recon] 1 [d(um)] [plcut]\n-->does not display logs...\n");
		fprintf(stderr, "d,plcutは任意\n");
		exit(1);
	}
	std::string in_momch = argv[1];
	std::string out_txt = argv[2];
	int ecc = std::stoi(argv[3]);
	int mode = -1;
	double d = 4.000;
	int plcut = 130;
	int trktype = std::stoi(argv[4]);
	if (argc ==6) {
		mode = std::stoi(argv[5]);
	}
	else if (argc == 7) {
		d = std::stod(argv[5]);
		plcut = std::stoi(argv[6]);
	}
	else {
		mode = std::stoi(argv[5]);
		d = std::stod(argv[6]);
		plcut = std::stoi(argv[7]);
	}



	std::string file_in_ECC, file_in_area;
	if (ecc < 4) {
		file_in_ECC = "T:\\NINJA\\E71a\\ECC" + std::to_string(ecc);
		file_in_area = file_in_ECC + "\\Area0";
	}
	else if (ecc == 4) {
		file_in_ECC = "K:\\NINJA\\E71a\\ECC" + std::to_string(ecc);
		file_in_area = file_in_ECC + "\\Area0";
	}
	else if (ecc < 7) {
		file_in_ECC = "I:\\NINJA\\E71a\\ECC" + std::to_string(ecc);
		file_in_area = file_in_ECC + "\\Area0";
	}
	else {
		file_in_ECC = "K:\\NINJA\\E71a\\ECC" + std::to_string(ecc);
		file_in_area = file_in_ECC + "\\Area0";
	}
	//chk
	std::cout << " Path" << std::endl;
	std::cout << "  " << file_in_ECC << std::endl;
	std::cout << "  " << file_in_area << std::endl;

	//std::cin;




	auto start = std::chrono::system_clock::now();//for measure working time
	
												  //read momch & set track
	std::cout << "\t* Reading momch.txt" << std::endl;
	std::multimap<tkey, stop_track> stop;
	if (trktype == 0) {
		read_stop_txt(in_momch, stop);
	}
	else {
		read_stop_txt_for_true(in_momch, stop);
	}
	std::cout << "\tfin reading stoptrack." << std::endl;



	// Add partner stop flg
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
	add_stop_flg(stop, area, 0, 4);






	//file消去
	std::ofstream ofs(out_txt);
	ofs.close();

	std::set<tkey> set;
	for (auto itr = stop.begin(); itr != stop.end(); itr++) {
		set.insert(itr->first);
	}
	std::cout << "\t  #event = " << set.size() << std::endl;

	std::cout << "\t* Calc Vertex point" << std::endl;
	ofs.open(out_txt);
	int remain = 0, allnum = 0;
	for (auto ev = set.begin(); ev != set.end(); ev++) {
		auto tks = stop.equal_range(*ev);
		if (mode == 1)	printf("\t\t* event %5d, #trk = %2d\n", *ev, stop.count(*ev));
		std::multimap<int, stop_track> rid;
		for (auto itr0 = tks.first; itr0 != tks.second; itr0++) {
			rid.insert(std::make_pair(itr0->second.chainid, itr0->second));
		}
		allnum++;
		std::cout << "\t" << rid.begin()->second.stoppl%1000 << " " << rid.begin()->second.x << "," << rid.begin()->second.y << std::endl;
		if (rid.begin()->second.stoppl%1000 > plcut)continue;
		if (rid.begin()->second.x < d || rid.begin()->second.x>250000 - d || rid.begin()->second.y < d || rid.begin()->second.y>250000 - d)continue;
		clustering_2trk_vtx2_ver2(rid, rid.begin()->second.stoppl, ofs, ecc, mode);
		remain++;
		rid.clear();
	}
	std::cout << remain << "/" << allnum << std::endl;
	auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
	MeasureProcessingTime(start, end);

}
void read_stop_txt(std::string in_momch, std::multimap<tkey, stop_track>& tracks) {

	std::ifstream ifs(in_momch);
	if (!ifs) {
		std::cerr << "file open error!" << std::endl;
		exit(0);
	}
	tkey tk;
	stop_track stop_tmp;

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	int recon, tr, k, bpl, dflg, nbtk;
	int flg = 0;
	while (std::getline(ifs, str)) {
		str_v = StringSplit_with_tab(str);
		//std::cout << str << std::endl;

		if (str_v.size() == 19 && std::stoi(str_v[1]) < 0) {
			//event header
			stop_tmp.groupid = std::stoi(str_v[0]);
			stop_tmp.unixtime = std::stoi(str_v[1]);
			stop_tmp.stoppl = std::stoi(str_v[4]);
			stop_tmp.weight = std::stod(str_v[13]);
			// to check the calc of vtx
			stop_tmp.vx = std::stod(str_v[7]);
			stop_tmp.vy = std::stod(str_v[8]);
			stop_tmp.vz = std::stod(str_v[9]);
			k = 0;
			stop_tmp.ntrk = std::stoi(str_v[17]);// #trk
			tr = std::stoi(str_v[18]);
			tk.eid = stop_tmp.groupid;
			tk.dz = stop_tmp.vz;
			tk.w = stop_tmp.weight;
			//std::cout << "header " << stop_tmp.groupid << " "<<stop_tmp.ntrk<< std::endl;
		}
		if (str_v.size() == 27) {
			flg = 0;
			// chain header
			if (int(std::stoi(str_v[2]) / 10000) > 0) {//recon
			//if (int(std::stoi(str_v[2]) / 10000) <= 0) {//MC true
				flg++;
				//if (k < stop_tmp.ntrk) {
					//recon
					//std::cout << "  chain header " << std::stoi(str_v[2]) << std::endl;
				stop_tmp.chainid = std::stoi(str_v[0]);
				stop_tmp.stopflg = std::stoi(str_v[1]);
				stop_tmp.pid = std::stoi(str_v[2]);
				dflg = std::stoi(str_v[3]);
				stop_tmp.mulikelihood = std::stod(str_v[23]);
				stop_tmp.pliklihoood = std::stod(str_v[24]);
				stop_tmp.nseg = std::stoi(str_v[25]);
				stop_tmp.ip = 0;
				//recon
				/**/
				stop_tmp.mom = std::stod(str_v[7]); //ecc_mcs assume muon
				if (stop_tmp.pid % 10000 == 2212) {
					stop_tmp.mom = std::stod(str_v[8]);	//assume p
					//c.ecc_mcs_mom_error[1][0]//+ [1][1] //-
				}
				// true (momはBM rng momしか持ってないので)
				//stop_tmp.mom = std::stod(str_v[9]);//BM rng mom
				//std::cout << stop_tmp.mom << std::endl;
				//std::cout << c.particle_flg << " " << c.muon_likelihood << " " << c.proton_likelihood << std::endl;
			}
			k++;
			nbtk = 0;
		}
		if (str_v.size() == 19 && flg > 0) {
			// basetracks	
			bpl = std::stod(str_v[0]);
			//std::cout << dflg << std::endl;
			//if (bpl > stop_tmp.stoppl) {
			if (dflg == -1) {
				// bwd
				if (nbtk == 0) {
					stop_tmp.pl0 = bpl;
					stop_tmp.rawid = std::stoi(str_v[1]);
					stop_tmp.ax = std::stod(str_v[2]);
					stop_tmp.ay = std::stod(str_v[3]);
					stop_tmp.x = std::stod(str_v[4]);
					stop_tmp.y = std::stod(str_v[5]);
					stop_tmp.z = std::stod(str_v[6]);
					stop_tmp.vph = std::stoi(str_v[10]) % 10000 + std::stoi(str_v[16]) % 10000;
					//if (sqrt(stop_tmp.ax * stop_tmp.ax + stop_tmp.ay * stop_tmp.ay) > 4.0) {
					//	std::cout << stop_tmp.pid << " " << stop_tmp.stoppl << " " << stop_tmp.pl0 << " " << stop_tmp.ax << " " << stop_tmp.ay << " " << sqrt(stop_tmp.ax * stop_tmp.ax + stop_tmp.ay * stop_tmp.ay) << std::endl;
					//}
				}
				if (nbtk + 1 == stop_tmp.nseg) {
					stop_tmp.pl1 = bpl;
					stop_tmp.npl = stop_tmp.pl1 - stop_tmp.pl0 + 1;
					if (stop_tmp.pid != 0)tracks.insert(std::make_pair(tk, stop_tmp));
					//std::cout << "\tbtrk header" << std::endl;
					//std::cout << stop_tmp.mom << std::endl;
				}
			}
			else if (dflg == 1) {
				// fwd
				if (nbtk == 0) {
					stop_tmp.pl0 = bpl;
					//if (sqrt(stop_tmp.ax * stop_tmp.ax + stop_tmp.ay * stop_tmp.ay) > 4.0) {
					//	std::cout << stop_tmp.pid << " " << stop_tmp.stoppl << " " << stop_tmp.pl0 << " " << stop_tmp.ax << " " << stop_tmp.ay << " " << sqrt(stop_tmp.ax * stop_tmp.ax + stop_tmp.ay * stop_tmp.ay) << std::endl;
					//}
				}
				if (nbtk + 1 == stop_tmp.nseg) {
					stop_tmp.pl1 = bpl;
					stop_tmp.rawid = std::stoi(str_v[1]);
					stop_tmp.ax = std::stod(str_v[2]);
					stop_tmp.ay = std::stod(str_v[3]);
					stop_tmp.x = std::stod(str_v[4]);
					stop_tmp.y = std::stod(str_v[5]);
					stop_tmp.z = std::stod(str_v[6]);
					stop_tmp.vph = std::stoi(str_v[10]) % 10000 + std::stoi(str_v[16]) % 10000;
					stop_tmp.npl = stop_tmp.pl1 - stop_tmp.pl0 + 1;
					if (stop_tmp.pid != 0)tracks.insert(std::make_pair(tk, stop_tmp));
					//std::cout << "\tbtrk header" << std::endl;

				}
			}
			else {
				std::cout << "Unexpected error....?" << std::endl;
				return;
			}
			nbtk++;
			//if (sqrt(stop_tmp.ax * stop_tmp.ax + stop_tmp.ay * stop_tmp.ay) > 4.0) {
			//	std::cout << stop_tmp.pid << " " << stop_tmp.stoppl << " " << stop_tmp.pl0 << " " << stop_tmp.ax << " " << stop_tmp.ay << " " << sqrt(stop_tmp.ax * stop_tmp.ax + stop_tmp.ay * stop_tmp.ay) << std::endl;
			//}

		}
		//if (str_v.size() == 19) {//denag
		//	if (sqrt(std::stod(str_v[2]) * std::stod(str_v[2]) + std::stod(str_v[3]) * std::stod(str_v[3])) > 4)std::cout << flg<<" "<<std::stod(str_v[0]) << " " << std::stod(str_v[2]) << " " << std::stod(str_v[3]) << " " << sqrt(std::stod(str_v[2]) * std::stod(str_v[2]) + std::stod(str_v[3]) * std::stod(str_v[3])) << std::endl;

		//}	
		//if (str_v.size() == 14) {
		//	// linklet
		//}
	}
	std::cout << "\t  tracks = " << tracks.size() << std::endl;
}

void read_stop_txt_for_true(std::string in_momch, std::multimap<tkey, stop_track>& tracks) {

	std::ifstream ifs(in_momch);
	if (!ifs) {
		std::cerr << "file open error!" << std::endl;
		exit(0);
	}
	tkey tk;
	stop_track stop_tmp;

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	int recon, tr, k, bpl, dflg, nbtk;
	int flg = 0;
	while (std::getline(ifs, str)) {
		str_v = StringSplit_with_tab(str);
		//std::cout << str << std::endl;

		if (str_v.size() == 19 && std::stoi(str_v[1]) < 0) {
			//event header
			stop_tmp.groupid = std::stoi(str_v[0]);
			stop_tmp.unixtime = std::stoi(str_v[1]);
			stop_tmp.stoppl = std::stoi(str_v[4]);
			stop_tmp.weight = std::stod(str_v[13]);
			// to check the calc of vtx
			stop_tmp.vx = std::stod(str_v[7]);
			stop_tmp.vy = std::stod(str_v[8]);
			stop_tmp.vz = std::stod(str_v[9]);
			k = 0;
			stop_tmp.ntrk = std::stoi(str_v[17]);// #trk
			tr = std::stoi(str_v[18]);
			tk.eid = stop_tmp.groupid;
			tk.dz = stop_tmp.vz;
			tk.w = stop_tmp.weight;
			//std::cout << "header " << stop_tmp.groupid << " "<<stop_tmp.ntrk<< std::endl;
		}
		if (str_v.size() == 27) {
			flg = 0;
			// chain header
			if (int(std::stoi(str_v[2]) / 10000) <= 0) {//MC true
				flg++;
				stop_tmp.chainid = std::stoi(str_v[0]);
				stop_tmp.stopflg = std::stoi(str_v[1]);
				stop_tmp.pid = std::stoi(str_v[2]);
				dflg = std::stoi(str_v[3]);
				stop_tmp.mulikelihood = std::stod(str_v[23]);
				stop_tmp.pliklihoood = std::stod(str_v[24]);
				stop_tmp.nseg = std::stoi(str_v[25]);
				stop_tmp.ip = 0;

				// true (momはBM rng momしか持ってないので)
				stop_tmp.mom = std::stod(str_v[9]);//BM rng mom
			}
			k++;
			nbtk = 0;
		}
		if (str_v.size() == 19 && flg > 0) {
			// basetracks	
			bpl = std::stod(str_v[0]);
			//std::cout << dflg << std::endl;
			//if (bpl > stop_tmp.stoppl) {
			if (dflg == -1) {
				// bwd
				if (nbtk == 0) {
					stop_tmp.pl0 = bpl;
					stop_tmp.rawid = std::stoi(str_v[1]);
					stop_tmp.ax = std::stod(str_v[2]);
					stop_tmp.ay = std::stod(str_v[3]);
					stop_tmp.x = std::stod(str_v[4]);
					stop_tmp.y = std::stod(str_v[5]);
					stop_tmp.z = std::stod(str_v[6]);
					stop_tmp.vph = std::stoi(str_v[10]) % 10000 + std::stoi(str_v[16]) % 10000;
					//if (sqrt(stop_tmp.ax * stop_tmp.ax + stop_tmp.ay * stop_tmp.ay) > 4.0) {
					//	std::cout << stop_tmp.pid << " " << stop_tmp.stoppl << " " << stop_tmp.pl0 << " " << stop_tmp.ax << " " << stop_tmp.ay << " " << sqrt(stop_tmp.ax * stop_tmp.ax + stop_tmp.ay * stop_tmp.ay) << std::endl;
					//}
				}
				if (nbtk + 1 == stop_tmp.nseg) {
					stop_tmp.pl1 = bpl;
					stop_tmp.npl = stop_tmp.pl1 - stop_tmp.pl0 + 1;
					if (stop_tmp.pid != 0)tracks.insert(std::make_pair(tk, stop_tmp));
					//std::cout << "\tbtrk header" << std::endl;
					//std::cout << stop_tmp.mom << std::endl;
				}
			}
			else if (dflg == 1) {
				// fwd
				if (nbtk == 0) {
					stop_tmp.pl0 = bpl;
					//if (sqrt(stop_tmp.ax * stop_tmp.ax + stop_tmp.ay * stop_tmp.ay) > 4.0) {
					//	std::cout << stop_tmp.pid << " " << stop_tmp.stoppl << " " << stop_tmp.pl0 << " " << stop_tmp.ax << " " << stop_tmp.ay << " " << sqrt(stop_tmp.ax * stop_tmp.ax + stop_tmp.ay * stop_tmp.ay) << std::endl;
					//}
				}
				if (nbtk + 1 == stop_tmp.nseg) {
					stop_tmp.pl1 = bpl;
					stop_tmp.rawid = std::stoi(str_v[1]);
					stop_tmp.ax = std::stod(str_v[2]);
					stop_tmp.ay = std::stod(str_v[3]);
					stop_tmp.x = std::stod(str_v[4]);
					stop_tmp.y = std::stod(str_v[5]);
					stop_tmp.z = std::stod(str_v[6]);
					stop_tmp.vph = std::stoi(str_v[10]) % 10000 + std::stoi(str_v[16]) % 10000;
					stop_tmp.npl = stop_tmp.pl1 - stop_tmp.pl0 + 1;
					if (stop_tmp.pid != 0)tracks.insert(std::make_pair(tk, stop_tmp));
					//std::cout << "\tbtrk header" << std::endl;

				}
			}
			else {
				std::cout << "Unexpected error....?" << std::endl;
				return;
			}
			nbtk++;

		}
	}
	std::cout << "\t  tracks = " << tracks.size() << std::endl;
}

void clustering_2trk_vtx2_ver2(std::multimap<int, stop_track>& tracks, int pl, std::ofstream& ofs, int ecc, int mode) {
	double refz = 0; int utime;
	for (auto itr1 = tracks.begin(); itr1 != tracks.end(); itr1++) {
		if (itr1->second.pid % 10000 == 13) {
			refz = itr1->second.z;
			utime = itr1->second.unixtime;
		}
	}

	//rawid,stop
	std::vector<track_multi> ret;
	double zrange[2] = { 0,0 };
	if (pl <= 15 || (pl >= 16 && pl % 2 == 0)) {
		zrange[0] = -1000;
	}
	else if (pl % 2 == 1) {
		zrange[0] = -3200; //3500->3200
	}

	double extra[2];
	track_multi multi;
	multi.pl = pl;
	double vz;
	//全2trkのmd計算
	if (tracks.size() != 1) {
		for (auto itr1 = tracks.begin(); itr1 != tracks.end(); itr1++) {
			multi.eventid = itr1->second.groupid;
			multi.unixtime = utime;
			vz = itr1->second.vz;
			for (auto itr2 = std::next(itr1, 1); itr2 != tracks.end(); itr2++) {
				matrix_3D::vector_3D pos0, pos1, dir0, dir1;
				pos0.x = itr1->second.x;
				pos0.y = itr1->second.y;
				pos0.z = itr1->second.z;
				pos1.x = itr2->second.x;
				pos1.y = itr2->second.y;
				pos1.z = itr2->second.z;
				dir0.x = itr1->second.ax;
				dir0.y = itr1->second.ay;
				dir0.z = 1;
				dir1.x = itr2->second.ax;
				dir1.y = itr2->second.ay;
				dir1.z = 1;

				// pos0を基準にzrangeの範囲内で最近接距離をとる位置(extra)を探索
				double md = minimum_distance_fixed(pos0, pos1, dir0, dir1, zrange, extra, refz);
				track_pair pair_tmp;
				matrix_3D::vector_3D extra0 = addition(pos0, const_multiple(dir0, extra[0]));
				matrix_3D::vector_3D extra1 = addition(pos1, const_multiple(dir1, extra[1]));

				pair_tmp.x = (extra0.x + extra1.x) / 2;
				pair_tmp.y = (extra0.y + extra1.y) / 2;
				pair_tmp.z = (extra0.z + extra1.z) / 2;
				pair_tmp.dz = pair_tmp.z - refz;
				pair_tmp.eventid = multi.eventid;
				pair_tmp.md = md;
				pair_tmp.oa = matrix_3D::opening_angle(dir0, dir1);
				pair_tmp.t[0] = itr1->second;
				pair_tmp.t[1] = itr2->second;

				multi.pair.push_back(pair_tmp);
			}
		}


		matrix_3D::vector_3D p_vtx, pos, dir;
		tkey k;
		//加重平均でvtx pointの決定
		multi.x = 0;
		multi.y = 0;
		multi.z = 0;
		for (auto itr = multi.pair.begin(); itr != multi.pair.end(); itr++) {
			multi.x += itr->x;
			multi.y += itr->y;
			multi.z += itr->z;
		}
		multi.x = multi.x / multi.pair.size();
		multi.y = multi.y / multi.pair.size();
		multi.z = multi.z / multi.pair.size();
		//multi.dz = multi.z - refz;
		multi.dz = vz - refz;
		//各trkに対してIPの計算
		for (auto itr = tracks.begin(); itr != tracks.end(); itr++) {
			matrix_3D::vector_3D pos0, pos1, dir0, dir1;
			pos0.x = itr->second.x;
			pos0.y = itr->second.y;
			pos0.z = itr->second.z;
			pos1.x = multi.x;
			pos1.y = multi.y;
			pos1.z = multi.z;
			dir0.x = itr->second.ax;
			dir0.y = itr->second.ay;
			dir0.z = 1;
			itr->second.ip = matrix_3D::inpact_parameter(pos0, dir0, pos1);
			//k.eid = itr->second.groupid;
			multi.trk.push_back(std::make_pair(itr->first, itr->second));
		}
		ret.push_back(multi);
	}
	else {
		for (auto itr1 = tracks.begin(); itr1 != tracks.end(); itr1++) {
			multi.eventid = itr1->second.groupid;
			multi.unixtime = utime;
			itr1->second.ip = 0;
			multi.trk.push_back(std::make_pair(itr1->first, itr1->second));
			track_pair pair_tmp = { 0 };

			multi.pair.push_back(pair_tmp);
		}
		multi.dz = -20000;

		ret.push_back(multi);

	}

	for (auto itr0 = ret.begin(); itr0 != ret.end(); itr0++) {
		for (auto itr1 = itr0->trk.begin(); itr1 != itr0->trk.end(); itr1++) {
			if (mode == 1)std::cout << "\tCalc:( " << itr0->x << ", " << itr0->y << ", " << itr0->z << " ), MC_recon:( " << itr1->second.vx << ", " << itr1->second.vy << ", " << itr1->second.vz << " )" << std::endl;
			ofs << std::right << std::fixed
				// information of vertex
				<< std::setw(2) << std::setprecision(0) << ecc << " "
				<< std::setw(6) << std::setprecision(0) << itr0->eventid << " "
				<< std::setw(3) << std::setprecision(0) << itr0->unixtime << " "
				<< std::setw(6) << std::setprecision(0) << itr0->pl << " "
				<< std::setw(4) << std::setprecision(0) << itr0->trk.size() << " "
				//<< std::setw(10) << std::setprecision(1) << itr1->second.vx << " "
				//<< std::setw(10) << std::setprecision(1) << itr1->second.vy << " "
				//<< std::setw(10) << std::setprecision(1) << itr1->second.vz << " "
				<< std::setw(10) << std::setprecision(1) << itr0->x << " "
				<< std::setw(10) << std::setprecision(1) << itr0->y << " "
				<< std::setw(10) << std::setprecision(1) << itr0->z << " "
				<< std::setw(8) << std::setprecision(1) << itr0->dz << " "
				// information of track
				<< std::setw(3) << std::setprecision(0) << itr1->second.stopflg << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.chainid << " "
				<< std::setw(4) << std::setprecision(0) << itr1->second.pl0 << " " // downstream
				<< std::setw(4) << std::setprecision(0) << itr1->second.pl1 << " " // upstream
				<< std::setw(10) << std::setprecision(0) << itr1->second.rawid << " "//vertex side
				//<< std::setw(7) << std::setprecision(0) << itr1->second.groupid << " "
				<< std::setw(4) << std::setprecision(0) << itr1->second.nseg << " "
				<< std::setw(4) << std::setprecision(0) << itr1->second.npl << " "
				<< std::setw(6) << std::setprecision(0) << itr1->second.vph << " "
				<< std::setw(10) << std::setprecision(0) << itr1->second.pid << " "
				<< std::setw(7) << std::setprecision(4) << itr1->second.ax << " "
				<< std::setw(7) << std::setprecision(4) << itr1->second.ay << " "
				<< std::setw(10) << std::setprecision(1) << itr1->second.x << " "
				<< std::setw(10) << std::setprecision(1) << itr1->second.y << " "
				<< std::setw(10) << std::setprecision(1) << itr1->second.z << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.ip << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.mom << " "
				<< std::setw(9) << std::setprecision(4) << itr1->second.mulikelihood << " "
				<< std::setw(9) << std::setprecision(4) << itr1->second.pliklihoood << " "
				<< std::setw(9) << std::setprecision(4) << itr1->second.weight
				<< std::endl;
		}
	}


}


void add_stop_flg(std::multimap<tkey, stop_track> &stop, std::map<int, std::vector<Fiducial_Area_Local>>& area, double edge_cut, int ex_pl_max) {

	int vpl;
	int cnt = 0; int stopnum = 0;
	for (auto& ev : stop) {
		vpl = ev.second.stoppl;
		if (ev.second.pid == 13)continue;
		ev.second.stopflg = 2;// stop
		cnt++;
		if(ev.second.pl0>ev.second.stoppl)	{// bwd
			if (ev.second.pl1 > 130) {
				ev.second.stopflg = 0;// penetrate
			}
			else {
				//最上流から4PL外挿,edgeから5mm以内に入ったらedge out
				ev.second.stopflg = divide_edge_out(ev.second, area, edge_cut, ex_pl_max);
			}
		}
		else {
			if (ev.second.pl0 < 5) {
				ev.second.stopflg = 0;//penetrate
			}
			else {
				//最上流から4PL外挿,edgeから5mm以内に入ったらedge out
				ev.second.stopflg = divide_edge_out_fwd(ev.second, area, edge_cut, ex_pl_max);

			}
		}
		if (ev.second.stopflg == 0)stopnum++;

	}



	std::cout << " # of Stop partner track = " << stopnum << " / " << cnt << std::endl;

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

int divide_edge_out(stop_track btrk, std::map<int, std::vector<Fiducial_Area_Local>>& area, double edge_cut, int ex_pl_max) {


	double xmin, xmax, ymin, ymax;


	int up_pl;
	double up_z, up_x, up_y, up_ax, up_ay, ex_z, ex_x, ex_y;
	int flg = 2;

	up_pl = btrk.pl0;// most upstream pl number
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
bool judge_fiducial_area(std::vector<Fiducial_Area_Local>& area,stop_track b) {

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

int divide_edge_out_fwd(stop_track btrk, std::map<int, std::vector<Fiducial_Area_Local>>& area, double edge_cut, int ex_pl_max) {


	double xmin, xmax, ymin, ymax;


	int up_pl;
	double up_z, up_x, up_y, up_ax, up_ay, ex_z, ex_x, ex_y;
	int flg = 2;

	up_pl = btrk.stoppl;// most downsttream pl number
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
bool judge_fiducial_area_fwd(std::vector<Fiducial_Area_Local>& area,stop_track b) {

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
