#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
# include <iostream>
# include <random>

#include <list>

struct tkey {
	int eid;
	double ip;
};
class stop_track {
public:
	// event
	int groupid;
	int64_t unixtime;
	int vpl, vertex_material;
	double w;//int bunch;
	int ntrk;
	//chain
	int chainid, stopflg, pid, charge;
	double MCS, MCSerr[2];
	int BMmomflg;// 
	double BMrng, BMrngerr[2];//rng
	double BMcurv, BMcurverr[2];//curverture
	double  mulikelihood, pliklihoood;
	int nseg;
	double pb;
	//basetrack
	int npl, pl0, pl1, vph, rawid;
	double ax, ay, x, y, z;
	double ip;
	double md, oa;// muonとの
	
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
	int64_t unixtime;
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


void read_stop_txt(std::vector<Momentum_recon::Event_information>& momch, std::multimap<int, stop_track>& tracks, double tan_thr);
void clustering_2trk_vtx2_ver2(std::multimap<int, stop_track>& tracks, std::multimap<int, stop_track>& notuse, int pl, std::ofstream& ofs, int ecc);
void clustering_2trk_vtx2(std::multimap<int, stop_track>& tracks, int pl, std::ofstream& ofs, int ecc);
void read_mcmomch_txt(std::string in_momch, std::multimap<int, stop_track>& tracks, double tan_thr);
int main(int argc, char** argv) {
	if (argc < 3) {
		fprintf(stderr, "usage:prg in-momch.momch out-vtx.txt #ECC [mode = 0] [tan_thr = 15.0] [edgecut] [option]\n");
		fprintf(stderr, "---------------------------------------------------------------------\n");
		fprintf(stderr, "  * mode (default = 0) \n");
		fprintf(stderr, "     mode =  0 : use all track to calcurate vertex point\n     mode = -5 : after manchk\n     mode = -2 : after manchk & rng-mom cut\n");
		fprintf(stderr, "  * tan_thr (default = 15.0) \n");
		fprintf(stderr, "     cut chains if it's sqrt(ax**2 + ay**2) > tan_thr.\n");
		fprintf(stderr, "  * edgecut (default = 0.0 um --> no cut)  \n");
		fprintf(stderr, "     rejected area.(muon xpos,ypos)\n");
		fprintf(stderr, "  * option \n");
		fprintf(stderr, "     if input 1 --> use momch.txt\n");
		fprintf(stderr, "---------------------------------------------------------------------\n");
		exit(1);
	}
	std::string in_momch = argv[1];
	std::string out_txt = argv[2];
	int ecc = std::stoi(argv[3]);
	int mode = 0;
	if (argc >= 5) {
		mode = std::stoi(argv[4]);
	}
	double tan_thr = 15.0;
	if (argc >= 6) {
		tan_thr = std::stod(argv[5]);
	}
	double edgecut = 0.0;
	if (argc >= 7) {
		edgecut = std::stod(argv[6]);
	}
	int readtxt = -1;
	if (argc >= 8) {
		readtxt = 1;
	}
//	int a;
//	std::cout << "(mode, tan_thr, edgecut) = ( " << mode << ", " << tan_thr << ", " << edgecut << ")" << std::endl;
//	std::cin >> a;


	auto start = std::chrono::system_clock::now();//for measure working time
	std::multimap<int, stop_track> stop;

	if (readtxt < 0) {
		std::cout << "Reading binary momch file." << std::endl;
		//read momch
		std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);
		// reading stop.txt
		read_stop_txt(momch, stop, tan_thr);
	}
	else {
		std::cout << "Reading text momch file." << std::endl;
		read_mcmomch_txt(in_momch, stop, tan_thr);
	}
	std::cout << "fin reading stoptrack." << std::endl;

	// 古いfile消去
	std::ofstream ofs(out_txt);
	ofs.close();

	std::set<int> set;
	for (auto itr = stop.begin(); itr != stop.end(); itr++) {
		set.insert(itr->first);// groupid
	}
	std::cout << set.size() << std::endl;

	ofs.open(out_txt);
	double pos[2];
	if (mode != 0) {// apply partner cut
		for (auto ev = set.begin(); ev != set.end(); ev++) {
			auto tks = stop.equal_range(*ev);
			printf("event %5d, vtx search\n", *ev);
			std::multimap<int, stop_track> rid, all;
			for (auto itr0 = tks.first; itr0 != tks.second; itr0++) {
				std::cout << std::right << std::fixed
					<< std::setw(5) << std::setprecision(0) << *ev << " "
					<< std::setw(5) << std::setprecision(0) << itr0->second.chainid << " "
					<< std::setw(3) << std::setprecision(0) << itr0->second.charge;
				if (itr0->second.pid % 10000 == 13) {
					rid.insert(std::make_pair(itr0->second.chainid, itr0->second));
					std::cout << std::endl;
					pos[0] = itr0->second.x;
					pos[1] = itr0->second.y;
				}
				else {
					if (itr0->second.charge < mode) {
						//muon以外でchargeはmanual checkの結果
						// -10 -> manchk
						// - 5 -> rng-mom
						all.insert(std::make_pair(itr0->second.chainid, itr0->second));
						std::cout << " -> This chain doesn't use to calcurate vertex point." << std::endl;
					}
					else {
						rid.insert(std::make_pair(itr0->second.chainid, itr0->second));
						std::cout << std::endl;
					}
				}
			}
			if (edgecut > 0) {
				if (pos[0] < edgecut || pos[0] > 250000.0 - edgecut || pos[1] < edgecut || pos[1] > 250000.0 - edgecut) {
					std::cout << "\t( x, y ) = ( " << pos[0] << ", " << pos[1] << " )\t-->Edge Area(rejected)." << std::endl;
				}
				else {
					clustering_2trk_vtx2_ver2(rid, all, rid.begin()->second.vpl, ofs, ecc);
				}
			}
			else {
				clustering_2trk_vtx2_ver2(rid, all, rid.begin()->second.vpl, ofs, ecc);
			}
			rid.clear();
		}
	}
	else {// use all trk to calculate vtx
		for (auto ev = set.begin(); ev != set.end(); ev++) {
			auto tks = stop.equal_range(*ev);
			printf("event %5d, vtx search\n", *ev);
			std::multimap<int, stop_track> rid, all;
			for (auto itr0 = tks.first; itr0 != tks.second; itr0++) {
				std::cout << std::right << std::fixed
					<< std::setw(5) << std::setprecision(0) << *ev << " "
					<< std::setw(5) << std::setprecision(0) << itr0->second.chainid << " "
					<< std::setw(3) << std::setprecision(0) << itr0->second.charge;
				rid.insert(std::make_pair(itr0->second.rawid, itr0->second));
				if (itr0->second.pid%10000 == 13) {
					pos[0] = itr0->second.x;
					pos[1] = itr0->second.y;
					std::cout << "\t( x, y ) = ( " << itr0->second.x << ", " << itr0->second.y << " )" << std::endl;
				}
				else {
					std::cout << std::endl;
				}
			}
			if (edgecut > 0.0) {
				if (pos[0] < edgecut || pos[0] > 250000.0 - edgecut || pos[1] < edgecut || pos[1] > 250000.0 - edgecut) {
					std::cout << "\t--> Edge Area(rejected)." << std::endl;
				}
				else {
					std::cout << "use" << std::endl;
					clustering_2trk_vtx2(rid, rid.begin()->second.vpl, ofs, ecc);
					
				}
			}
			else {
				clustering_2trk_vtx2(rid, rid.begin()->second.vpl, ofs, ecc);
				std::cout << std::endl;
			}
			rid.clear();
		}

	}
	auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
	MeasureProcessingTime(start, end);

}

void read_stop_txt(std::vector<Momentum_recon::Event_information>& momch, std::multimap<int, stop_track>& tracks, double tan_thr) {
	stop_track stop_tmp;

	int pos = 0;
	int cnt = 0;
	int cutnum = 0;
	for (auto& ev : momch) {
		// event header
		//stop_tmp.vpl = ev.vertex_pl;
		stop_tmp.vpl = ev.vertex_pl%1000;//formc
		stop_tmp.groupid = ev.groupid;
		stop_tmp.unixtime = ev.unix_time;
		stop_tmp.vertex_material = ev.vertex_material;
		stop_tmp.w = ev.weight;
		stop_tmp.ntrk = ev.chains.size();
		std::cout << ev.weight << std::endl;
		for (auto& c : ev.chains) {
			// chain header
			stop_tmp.chainid = c.chainid;
			stop_tmp.stopflg = c.stop_flg;
			stop_tmp.pid = c.particle_flg;
			stop_tmp.charge = c.charge_sign;// muon以外にはmanual checkの結果
			stop_tmp.MCS = c.ecc_mcs_mom[0];//assume mu 
			stop_tmp.MCSerr[0] = c.ecc_mcs_mom_error[0][0];
			stop_tmp.MCSerr[1] = c.ecc_mcs_mom_error[0][1];
			if (int(c.particle_flg %10000)== 2212) {
				stop_tmp.MCS = c.ecc_mcs_mom[1];//assume p
				stop_tmp.MCSerr[0] = c.ecc_mcs_mom_error[1][0];
				stop_tmp.MCSerr[1] = c.ecc_mcs_mom_error[1][1];
				//c.ecc_mcs_mom_error[1][0]//+ [1][1] //-
			}
			stop_tmp.pb = c.Get_muon_mcs_pb();
			stop_tmp.BMmomflg = c.bm_range_mom_error[0];
			stop_tmp.BMrng = c.bm_range_mom;
			stop_tmp.BMrngerr[0] = c.bm_range_mom_error[0];
			stop_tmp.BMrngerr[1] = c.bm_range_mom_error[1];
			stop_tmp.BMcurv = c.bm_curvature_mom;
			stop_tmp.BMcurverr[0] = c.bm_curvature_mom_error[0];
			stop_tmp.BMcurverr[1] = c.bm_curvature_mom_error[1];
			stop_tmp.mulikelihood = c.muon_likelihood;
			stop_tmp.pliklihoood = c.proton_likelihood;

			stop_tmp.nseg = c.base.size();
			//std::cout << c.base.size();
			stop_tmp.pl0 = c.base.begin()->pl;//dounstream?
			stop_tmp.pl1 = c.base.rbegin()->pl;
			stop_tmp.npl = stop_tmp.pl1 - stop_tmp.pl0 + 1;
			//std::cout << c.particle_flg << " " << c.muon_likelihood << " " << c.proton_likelihood << std::endl;

			if (stop_tmp.pl1 >= stop_tmp.vpl) {
				//stop
				//std::cout << pos << " " << c.base.rbegin()->ax << std::endl;
				pos = stop_tmp.pl1;
				stop_tmp.rawid = c.base.rbegin()->rawid;
				stop_tmp.ax = c.base.rbegin()->ax;
				stop_tmp.ay = c.base.rbegin()->ay;
				stop_tmp.x = c.base.rbegin()->x;
				stop_tmp.y = c.base.rbegin()->y;
				stop_tmp.z = c.base.rbegin()->z;
				stop_tmp.vph = c.base.rbegin()->m[0].ph % 10000 + c.base.rbegin()->m[1].ph % 10000;
			}
			if (stop_tmp.pl0 > stop_tmp.vpl) {
				//start
				stop_tmp.rawid = c.base.begin()->rawid;
				stop_tmp.ax = c.base.begin()->ax;
				stop_tmp.ay = c.base.begin()->ay;
				stop_tmp.x = c.base.begin()->x;
				stop_tmp.y = c.base.begin()->y;
				stop_tmp.z = c.base.begin()->z;
				stop_tmp.vph = c.base.begin()->m[0].ph % 10000 + c.base.begin()->m[1].ph % 10000;
			}
			stop_tmp.ip = 0;
			stop_tmp.md = 0;
			stop_tmp.oa = 0;
			std::cout << stop_tmp.ax<<","<<stop_tmp.ay << std::endl;

			if (sqrt(stop_tmp.ax * stop_tmp.ax + stop_tmp.ay * stop_tmp.ay) <= tan_thr) {//angle cut
				tracks.insert(std::make_pair(stop_tmp.groupid, stop_tmp));
				cutnum++;
			}

		}
		cnt++;
	}

	//printf("input fin %d track\n", cnt);
	std::cout << " * Angle Cut : " << cnt << "-->" << cutnum << std::endl;
}
void read_mcmomch_txt(std::string in_momch, std::multimap<int, stop_track>& tracks, double tan_thr) {

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
	int all = 0;
	int cutnum = 0;
	int addgid = 0;
	while (std::getline(ifs, str)) {
		str_v = StringSplit_with_tab(str);
		//std::cout << str << std::endl;

		if (str_v.size() == 19 && std::stoi(str_v[1]) < 0) {//event header
			//stop_tmp.groupid = std::stoi(str_v[0]);
			stop_tmp.groupid = addgid;
			addgid++;
			stop_tmp.unixtime = std::stoi(str_v[1]);
			stop_tmp.vpl = std::stoi(str_v[4]) % 1000;
			stop_tmp.vertex_material = std::stoi(str_v[6]);

			// vertex pos(recon)
			//stop_tmp.vx = std::stod(str_v[7]);
			//stop_tmp.vy = std::stod(str_v[8]);
			//stop_tmp.vz = std::stod(str_v[9]);

			stop_tmp.w = std::stod(str_v[13]);
			k = 0;
			stop_tmp.ntrk = std::stoi(str_v[17]);// # of recon chain
			tr = std::stoi(str_v[18]);// # of true chain
			//tk.eid = stop_tmp.groupid;
			//tk.dz = stop_tmp.vz;
			//tk.w = stop_tmp.weight;

		}
		if (str_v.size() == 27) {// chain header
			all++;
			flg = 0;
			stop_tmp.chainid = std::stoi(str_v[0]);
			stop_tmp.stopflg = std::stoi(str_v[1]);
			stop_tmp.pid = std::stoi(str_v[2]);
			dflg = std::stoi(str_v[3]);
			stop_tmp.charge = std::stoi(str_v[4]);// muon以外にはmanual checkの結果
			stop_tmp.MCS = std::stod(str_v[7]);//assume mu 
			stop_tmp.MCSerr[0] = std::stod(str_v[15]);
			stop_tmp.MCSerr[1] = std::stod(str_v[16]);
			if (stop_tmp.pid % 10000 == 2212) {
				stop_tmp.MCS = std::stod(str_v[8]);	//assume p
				stop_tmp.MCSerr[0] = std::stod(str_v[17]);
				stop_tmp.MCSerr[1] = std::stod(str_v[18]);
			}
			stop_tmp.pb = std::stod(str_v[5]);// ecc range
			stop_tmp.BMmomflg = 0;
			stop_tmp.BMrng = std::stod(str_v[9]);
			stop_tmp.BMrngerr[0] = std::stod(str_v[19]);
			stop_tmp.BMrngerr[1] = std::stod(str_v[20]);
			stop_tmp.BMcurv = std::stod(str_v[11]);
			stop_tmp.BMcurverr[0] = std::stod(str_v[21]);
			stop_tmp.BMcurverr[1] = std::stod(str_v[22]);
			stop_tmp.mulikelihood = std::stod(str_v[23]);
			stop_tmp.pliklihoood = std::stod(str_v[24]);
			stop_tmp.nseg = std::stoi(str_v[25]);

			stop_tmp.ip = 0;
			stop_tmp.md = 0;
			stop_tmp.oa = 0;
			if (int(std::stoi(str_v[2]) / 10000) > 0) {//recon
				flg++;
			}
			k++;
			nbtk = 0;
		}
		if (str_v.size() == 19 && flg > 0) {// && k<=stop_tmp.ntrk) {// basetracks	
			bpl = std::stod(str_v[0]);

			if (dflg == -1) {// bwd
				if (nbtk == 0) {// downstream
					stop_tmp.pl0 = bpl;
					stop_tmp.rawid = std::stoi(str_v[1]);
					stop_tmp.ax = std::stod(str_v[2]);
					stop_tmp.ay = std::stod(str_v[3]);
					stop_tmp.x = std::stod(str_v[4]);
					stop_tmp.y = std::stod(str_v[5]);
					stop_tmp.z = std::stod(str_v[6]);
					stop_tmp.vph = std::stoi(str_v[10]) % 10000 + std::stoi(str_v[16]) % 10000;
				}
				if (nbtk + 1 == stop_tmp.nseg) {// upstream
					stop_tmp.pl1 = bpl;
					stop_tmp.npl = stop_tmp.pl1 - stop_tmp.pl0 + 1;
					if (sqrt(stop_tmp.ax * stop_tmp.ax + stop_tmp.ay * stop_tmp.ay) <= tan_thr) {
						tracks.insert(std::make_pair(stop_tmp.groupid, stop_tmp));
					}
					else {
						cutnum++;
					}
				}
			}
			else if (dflg == 1) {// fwd
				if (nbtk == 0) {
					stop_tmp.pl0 = bpl;
					stop_tmp.vph = std::stoi(str_v[10]) % 10000 + std::stoi(str_v[16]) % 10000;
				}
				if (nbtk + 1 == stop_tmp.nseg) {
					stop_tmp.pl1 = bpl;
					stop_tmp.rawid = std::stoi(str_v[1]);
					stop_tmp.ax = std::stod(str_v[2]);
					stop_tmp.ay = std::stod(str_v[3]);
					stop_tmp.x = std::stod(str_v[4]);
					stop_tmp.y = std::stod(str_v[5]);
					stop_tmp.z = std::stod(str_v[6]);
					stop_tmp.npl = stop_tmp.pl1 - stop_tmp.pl0 + 1;
					if (std::stoi(str_v[10]) % 10000 + std::stoi(str_v[16]) % 10000 > 0) {
						stop_tmp.vph = std::stoi(str_v[10]) % 10000 + std::stoi(str_v[16]) % 10000;

					}
					if (sqrt(stop_tmp.ax * stop_tmp.ax + stop_tmp.ay * stop_tmp.ay) <= tan_thr) {//angle cut
						tracks.insert(std::make_pair(stop_tmp.groupid, stop_tmp));
					}
					else {
						cutnum++;
					}
				}
			}
			else {
				std::cout << "Unexpected error....?" << std::endl;
				return;
			}
			nbtk++;
		}
		//if (str_v.size() == 19) {//denag
		//	if (sqrt(std::stod(str_v[2]) * std::stod(str_v[2]) + std::stod(str_v[3]) * std::stod(str_v[3])) > 4)std::cout << flg<<" "<<std::stod(str_v[0]) << " " << std::stod(str_v[2]) << " " << std::stod(str_v[3]) << " " << sqrt(std::stod(str_v[2]) * std::stod(str_v[2]) + std::stod(str_v[3]) * std::stod(str_v[3])) << std::endl;

		//}	
		//if (str_v.size() == 14) {
		//	// linklet
		//}
	}
		std::cout << "\t  tracks : " << all << " --> " << tracks.size() << " ( rejected : " << cutnum << ")" << std::endl;
}

void clustering_2trk_vtx2_ver2(std::multimap<int, stop_track>& tracks, std::multimap<int, stop_track>& notuse, int pl, std::ofstream& ofs, int ecc) {
	double refz = 0; int utime;
	for (auto itr1 = tracks.begin(); itr1 != tracks.end(); itr1++) {
		if (itr1->second.chainid == 0) {
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
	//全2trkのmd計算
	std::cout << tracks.size() << std::endl;
	if (tracks.size() != 1) {
		for (auto itr1 = tracks.begin(); itr1 != tracks.end(); itr1++) {
			multi.eventid = itr1->second.groupid;
			multi.unixtime = utime;
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
				if (itr1->second.pid != 13 && itr2->second.pid == 13) {
					itr1->second.md = md;
					itr1->second.oa = pair_tmp.oa;
				}
				if (itr2->second.pid != 13 && itr1->second.pid == 13) {
					itr2->second.md = md;
					itr2->second.oa = pair_tmp.oa;
				}
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
		multi.dz = multi.z - refz;
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
			ofs << std::right << std::fixed
				// information of vertex
				<< std::setw(1) << std::setprecision(0) << ecc << " "
				<< std::setw(5) << std::setprecision(0) << itr0->eventid << " "
				<< std::setw(10) << std::setprecision(0) << itr0->unixtime << " "
				<< std::setw(10) << std::setprecision(4) << itr1->second.w << " "
				<< std::setw(3) << std::setprecision(0) << itr0->pl << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.vertex_material << " "
				<< std::setw(3) << std::setprecision(0) << itr0->trk.size() << " "//finalchk後のvtxを構成するtrk数
				<< std::setw(3) << std::setprecision(0) << itr1->second.ntrk << " "//もともとの
				<< std::setw(10) << std::setprecision(1) << itr0->x << " "
				<< std::setw(10) << std::setprecision(1) << itr0->y << " "
				<< std::setw(10) << std::setprecision(1) << itr0->z << " "
				<< std::setw(8) << std::setprecision(1) << itr0->dz << " "
				// information of chain
				<< std::setw(3) << std::setprecision(0) << itr1->second.chainid << " "
				<< std::setw(8) << std::setprecision(0) << itr1->second.pid << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.stopflg << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.charge << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.MCS << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.MCSerr[0] << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.MCSerr[1] << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.BMmomflg << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.BMrng << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.BMrngerr[0] << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.BMrngerr[1] << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.BMcurv << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.BMcurverr[0] << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.BMcurverr[1] << " "
				<< std::setw(9) << std::setprecision(4) << itr1->second.mulikelihood << " "
				<< std::setw(9) << std::setprecision(4) << itr1->second.pliklihoood << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.pb << " "
				// information of trk
				<< std::setw(3) << std::setprecision(0) << itr1->second.nseg << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.npl << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.pl0 << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.pl1 << " "
				<< std::setw(10) << std::setprecision(0) << itr1->second.rawid << " "
				<< std::setw(6) << std::setprecision(0) << itr1->second.vph << " "
				<< std::setw(7) << std::setprecision(4) << itr1->second.ax << " "
				<< std::setw(7) << std::setprecision(4) << itr1->second.ay << " "
				<< std::setw(10) << std::setprecision(1) << itr1->second.x << " "
				<< std::setw(10) << std::setprecision(1) << itr1->second.y << " "
				<< std::setw(10) << std::setprecision(1) << itr1->second.z << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.ip << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.md << " "
				<< std::setw(8) << std::setprecision(4) << itr1->second.oa << " "

				<< std::endl;

		}
		for (auto itr2 = notuse.begin(); itr2 != notuse.end(); itr2++) {
			// information of vertex
			ofs << std::setw(1) << std::setprecision(0) << ecc << " "
				<< std::setw(5) << std::setprecision(0) << itr0->eventid << " "
				<< std::setw(10) << std::setprecision(0) << itr0->unixtime << " "
				<< std::setw(10) << std::setprecision(4) << itr2->second.w << " "
				<< std::setw(3) << std::setprecision(0) << itr0->pl << " "
				<< std::setw(3) << std::setprecision(0) << itr2->second.vertex_material << " "
				<< std::setw(3) << std::setprecision(0) << itr0->trk.size() << " "//finalchk後のvtxを構成するtrk数
				<< std::setw(3) << std::setprecision(0) << itr2->second.ntrk << " "//もともとの
				<< std::setw(10) << std::setprecision(1) << itr0->x << " "
				<< std::setw(10) << std::setprecision(1) << itr0->y << " "
				<< std::setw(10) << std::setprecision(1) << itr0->z << " "
				<< std::setw(8) << std::setprecision(1) << itr0->dz << " "
				// information of chain
				<< std::setw(3) << std::setprecision(0) << itr2->second.chainid << " "
				<< std::setw(8) << std::setprecision(0) << itr2->second.pid << " "
				<< std::setw(3) << std::setprecision(0) << itr2->second.stopflg << " "
				<< std::setw(3) << std::setprecision(0) << itr2->second.charge << " "
				<< std::setw(8) << std::setprecision(1) << itr2->second.MCS << " "
				<< std::setw(8) << std::setprecision(1) << itr2->second.MCSerr[0] << " "
				<< std::setw(8) << std::setprecision(1) << itr2->second.MCSerr[1] << " "
				<< std::setw(3) << std::setprecision(0) << itr2->second.BMmomflg << " "
				<< std::setw(8) << std::setprecision(1) << itr2->second.BMrng << " "
				<< std::setw(8) << std::setprecision(1) << itr2->second.BMrngerr[0] << " "
				<< std::setw(8) << std::setprecision(1) << itr2->second.BMrngerr[1] << " "
				<< std::setw(8) << std::setprecision(1) << itr2->second.BMcurv << " "
				<< std::setw(8) << std::setprecision(1) << itr2->second.BMcurverr[0] << " "
				<< std::setw(8) << std::setprecision(1) << itr2->second.BMcurverr[1] << " "
				<< std::setw(9) << std::setprecision(4) << itr2->second.mulikelihood << " "
				<< std::setw(9) << std::setprecision(4) << itr2->second.pliklihoood << " "
				<< std::setw(8) << std::setprecision(1) << itr2->second.pb << " "
				// information of trk
				<< std::setw(3) << std::setprecision(0) << itr2->second.nseg << " "
				<< std::setw(3) << std::setprecision(0) << itr2->second.npl << " "
				<< std::setw(3) << std::setprecision(0) << itr2->second.pl0 << " "
				<< std::setw(3) << std::setprecision(0) << itr2->second.pl1 << " "
				<< std::setw(10) << std::setprecision(0) << itr2->second.rawid << " "
				<< std::setw(6) << std::setprecision(0) << itr2->second.vph << " "
				<< std::setw(7) << std::setprecision(4) << itr2->second.ax << " "
				<< std::setw(7) << std::setprecision(4) << itr2->second.ay << " "
				<< std::setw(10) << std::setprecision(1) << itr2->second.x << " "
				<< std::setw(10) << std::setprecision(1) << itr2->second.y << " "
				<< std::setw(10) << std::setprecision(1) << itr2->second.z << " "
				<< std::setw(8) << std::setprecision(1) << itr2->second.ip << " "
				<< std::setw(8) << std::setprecision(1) << itr2->second.md << " "
				<< std::setw(8) << std::setprecision(4) << itr2->second.oa << " "

				<< std::endl;
		}
	}


}


void clustering_2trk_vtx2(std::multimap<int, stop_track>& tracks, int pl, std::ofstream& ofs, int ecc) {
	double refz = 0; int utime;
	for (auto itr1 = tracks.begin(); itr1 != tracks.end(); itr1++) {
		if (itr1->second.chainid == 0) {
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
	//全2trkのmd計算
	std::cout << tracks.size() << std::endl;
	if (tracks.size() != 1) {
		for (auto itr1 = tracks.begin(); itr1 != tracks.end(); itr1++) {
			multi.eventid = itr1->second.groupid;
			multi.unixtime = utime;
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
				if (itr1->second.pid != 13 && itr2->second.pid == 13) {
					itr1->second.md = pair_tmp.md;
					itr1->second.oa = pair_tmp.oa;
				}
				if (itr2->second.pid != 13 && itr1->second.pid == 13) {
					itr2->second.md = pair_tmp.md;
					itr2->second.oa = pair_tmp.oa;
				}
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
		multi.dz = multi.z - refz;
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
			ofs << std::right << std::fixed
				// information of vertex
				<< std::setw(1) << std::setprecision(0) << ecc << " "
				<< std::setw(5) << std::setprecision(0) << itr0->eventid << " "
				<< std::setw(10) << std::setprecision(4) << itr0->unixtime << " "
				<< std::setw(10) << std::setprecision(4) << itr1->second.w << " "
				<< std::setw(3) << std::setprecision(0) << itr0->pl << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.vertex_material << " "
				<< std::setw(3) << std::setprecision(0) << itr0->trk.size() << " "				//finalchk後のvtxを構成するtrk数
				<< std::setw(3) << std::setprecision(0) << itr1->second.ntrk << " "				//もともとの
				<< std::setw(10) << std::setprecision(1) << itr0->x << " "
				<< std::setw(10) << std::setprecision(1) << itr0->y << " "
				<< std::setw(10) << std::setprecision(1) << itr0->z << " "
				<< std::setw(8) << std::setprecision(1) << itr0->dz << " "
				// information of chain
				<< std::setw(3) << std::setprecision(0) << itr1->second.chainid << " "
				<< std::setw(8) << std::setprecision(0) << itr1->second.pid << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.stopflg << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.charge << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.MCS << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.MCSerr[0] << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.MCSerr[1] << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.BMmomflg << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.BMrng << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.BMrngerr[0] << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.BMrngerr[1] << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.BMcurv << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.BMcurverr[0] << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.BMcurverr[1] << " "
				<< std::setw(9) << std::setprecision(4) << itr1->second.mulikelihood << " "
				<< std::setw(9) << std::setprecision(4) << itr1->second.pliklihoood << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.pb << " "
				// information of trk
				<< std::setw(3) << std::setprecision(0) << itr1->second.nseg << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.npl << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.pl0 << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.pl1 << " "
				<< std::setw(10) << std::setprecision(0) << itr1->second.rawid << " "
				<< std::setw(6) << std::setprecision(0) << itr1->second.vph << " "
				<< std::setw(7) << std::setprecision(4) << itr1->second.ax << " "
				<< std::setw(7) << std::setprecision(4) << itr1->second.ay << " "
				<< std::setw(10) << std::setprecision(1) << itr1->second.x << " "
				<< std::setw(10) << std::setprecision(1) << itr1->second.y << " "
				<< std::setw(10) << std::setprecision(1) << itr1->second.z << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.ip << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.md << " "
				<< std::setw(8) << std::setprecision(4) << itr1->second.oa << " "

				<< std::endl;

		}
	}


}