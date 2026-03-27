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
	int ecc;
	// event
	int groupid;
	int64_t unixtime;
	int bunch;
	int vpl, vertex_material;
	int ntrk,ntrk_all;
	double vx, vy, vz, dz;
	//chain
	int chainid, pid, stopflg, charge;
	double MCS, MCSerr[2];
	int BMmomflg;// 
	double BMrng, BMrngerr[2];//rng
	double BMcurv, BMcurverr[2];//curverture
	double  mulikelihood, pliklihoood;
	double pb;
	int nseg;
	//basetrack
	int npl, pl0, pl1, rawid, vph;
	double ax, ay, x, y, z;
	double ip;
	double md, oa;// muonとの
	//all
	double al_md, al_dz;
	int selflg;
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


void set_chain(std::string in, std::multimap<int, stop_track>& stop);
void clustering_2trk_vtx2_ver2(std::multimap<int, stop_track>& tracks, std::multimap<int, stop_track>& notuse, int pl, std::ofstream& ofs);

int main(int argc, char** argv) {
	if (argc < 9) {
		fprintf(stderr, "usage:prg in.txt out.txt tan_thr md_thr npl_thr bwdtrkcutflg fillfactor_thr VPH_thr [mode] [p_fillfactor_thr] [p_md_thr]\n");
		fprintf(stderr, "---------------------------------------------------------------------\n");
		fprintf(stderr, "For mip, if \n");
		fprintf(stderr, "  tan > tan_thr => Doesn't use to calcurate vtx point\n");
		fprintf(stderr, "  md  > md_thr  => Doesn't use to calcurate vtx point\n");
		fprintf(stderr, "  npl < npl_thr => Doesn't use to calcurate vtx point\n");
		fprintf(stderr, "  bwdtrkcut = 1 => Bwd trk isn't used to calcurate vtx point(-1-->off)\n");
		fprintf(stderr, "  ffactor < thr => Doesn't use to calcurate vtx point\n");
		fprintf(stderr, "  VPH < VPH_thr => Doesn't use to calcurate vtx point\n");
		fprintf(stderr, "\nIf you input mode == 1, this cut-off criterion will also be applied to protons.\n");
		fprintf(stderr, "For protons, if \n");
		fprintf(stderr, "  ffactor < p_fillfactor_thr => Doesn't use to calcurate vtx point(Default = 0.68)\n");
		fprintf(stderr, "  md      > p_md_thr         => Doesn't use to calcurate vtx point(Default = 80)\n");
		fprintf(stderr, "this cut-off criterion will always be applied to protons\n");
		fprintf(stderr, "---------------------------------------------------------------------\n");
		std::cout << argc << std::endl;
		exit(1);
	}
	std::string in_txt = argv[1];
	std::string out_txt = argv[2];
	double tan_thr = std::stoi(argv[3]);
	double md_thr = std::stoi(argv[4]);
	int	npl_thr = std::stod(argv[5]);
	int	mode = std::stoi(argv[6]);
	double fillfactor_thr = std::stod(argv[7]);
	int vph_thr = std::stoi(argv[8]);
	int	mode_p = -1;
	if (argc == 10) {
		mode_p = std::stoi(argv[9]);
	}
	double p_fillfactor_thr = 0.68;
	double p_md_thr = 80;
	if (argc == 11) {
		p_fillfactor_thr = std::stoi(argv[10]);
	}
	if (argc == 12) {
		p_md_thr = std::stoi(argv[11]);
	}

	// reading stop.txt
	std::multimap<int, stop_track> stop;
	set_chain(in_txt, stop);
	std::cout << "fin reading stoptrack." << std::endl;

	// 古いfile消去
	std::ofstream ofs(out_txt);
	ofs.close();

	std::set<int> set;
	for (auto itr = stop.begin(); itr != stop.end(); itr++) {
		set.insert(itr->first);
	}
	std::cout << set.size() << std::endl;

	ofs.open(out_txt);
	double tan, fillfactor;
	int flg;
	int p_flg;
	for (auto ev = set.begin(); ev != set.end(); ev++) {
		auto tks = stop.equal_range(*ev);
		printf(" * event %5d, vtx search\n", *ev);
		std::multimap<int, stop_track> sig, noi;
		for (auto itr0 = tks.first; itr0 != tks.second; itr0++) {
			flg = 0;
			p_flg = 1;
			std::cout << std::right << std::fixed
				<< std::setw(10) << std::setprecision(0) << *ev << " "
				<< std::setw(5) << std::setprecision(0) << itr0->second.chainid << " "
				<< std::setw(5) << std::setprecision(0) << itr0->second.charge << " "
				<< std::setw(5) << std::setprecision(0) << itr0->second.pid << " ";
			fillfactor = double(itr0->second.nseg) / double(itr0->second.npl);
			std::cout << std::setw(5) << std::setprecision(4) << fillfactor;

			if (itr0->second.pid == 13) {
				sig.insert(std::make_pair(itr0->second.rawid, itr0->second));
				std::cout << std::endl;
			}
			else {
				tan = sqrt(itr0->second.ax * itr0->second.ax + itr0->second.ay * itr0->second.ay);
				if (itr0->second.al_md > md_thr) {// cut
					flg++;
				}
				if (itr0->second.npl < npl_thr) {// cut
					flg++;
				}
				if (tan > tan_thr) {//cut 
					flg++;
				}
				if (mode > 0 && itr0->second.pl0 > itr0->second.vpl) {// cut bwd trk
					flg++;
				}
				if (fillfactor < fillfactor_thr) {// fillfactor cut
					flg++;
				}
				if (itr0->second.vph < vph_thr) {// vph cut
					flg++;
				}

				if (itr0->second.pid == 2212 && mode_p < 0 && fillfactor >= p_fillfactor_thr && itr0->second.al_md <= p_md_thr) {// && itr0->second.al_md <= md_thr) {// カットしないproton
					p_flg = -1;
				}
				//if (itr0->second.pid == 2212 && mode_p < 0 ) {// カットしないproton
				//	p_flg = -1;
				//}

				if (p_flg > 0 && flg > 0) {
					// vertexの計算に使わない
					noi.insert(std::make_pair(itr0->second.chainid, itr0->second));
					//<< std::setw(5) << std::setprecision(0) << itr0->second.nseg << " "
					//<< std::setw(5) << std::setprecision(0) << itr0->second.npl << " "
					std::cout << " : flg = " << flg << " > 0 ";
					std::cout << " -> This chain doesn't use to calcurate vertex point." << std::endl;
				}
				else {
					// vertexの計算に使う
					sig.insert(std::make_pair(itr0->second.rawid, itr0->second));
					std::cout << std::endl;
				}
			}
		}
		clustering_2trk_vtx2_ver2(sig, noi, sig.begin()->second.vpl, ofs);
		sig.clear();
		noi.clear();
	}

}
void set_chain(std::string in, std::multimap<int, stop_track> &stop) {
		std::ifstream ifs(in);
	if (!ifs) {
		std::cerr << " File open error ! " << std::endl;
		return;
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;
	stop_track s;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);
		// information of vertex
		s.ecc = std::stoi(str_v[0]);
		s.groupid = std::stoi(str_v[1]);
		s.unixtime = std::stoll(str_v[2]);
		s.bunch = std::stoi(str_v[3]);
		s.vpl = std::stoi(str_v[4]);
		s.vertex_material = std::stoi(str_v[5]);
		s.ntrk = std::stoi(str_v[6]);
		s.ntrk_all = std::stoi(str_v[7]);
		s.vx = std::stod(str_v[8]);
		s.vy = std::stod(str_v[9]);
		s.vz = std::stod(str_v[10]);
		s.dz = 0;
		// information of chain
		s.chainid = std::stoi(str_v[12]);
		s.pid = std::stoi(str_v[13]);
		s.stopflg = std::stoi(str_v[14]);
		s.charge = std::stoi(str_v[15]);
		s.MCS = std::stod(str_v[16]);
		s.MCSerr[0] = std::stod(str_v[17]);
		s.MCSerr[1] = std::stod(str_v[18]);
		s.BMmomflg = std::stoi(str_v[19]);
		s.BMrng = std::stod(str_v[20]);
		s.BMrngerr[0] = std::stod(str_v[21]);
		s.BMrngerr[1] = std::stod(str_v[22]);
		s.BMcurv = std::stod(str_v[23]);
		s.BMcurverr[0] = std::stod(str_v[24]);
		s.BMcurverr[1] = std::stod(str_v[25]);
		s.mulikelihood = std::stod(str_v[26]);
		s.pliklihoood = std::stod(str_v[27]);
		s.pb = std::stod(str_v[28]);
		// information of trk
		s.nseg = std::stoi(str_v[29]);
		s.npl = std::stoi(str_v[30]);
		s.pl0 = std::stoi(str_v[31]);
		s.pl1 = std::stoi(str_v[32]);
		s.rawid = std::stoi(str_v[33]);
		s.vph = std::stoi(str_v[34]);
		s.ax = std::stod(str_v[35]);
		s.ay = std::stod(str_v[36]);
		s.x = std::stod(str_v[37]);
		s.y = std::stod(str_v[38]);
		s.z = std::stod(str_v[39]);
		s.ip = 0;
		s.md = 0;
		s.oa = 0;
		s.al_md = std::stod(str_v[41]);
		s.al_dz = std::stod(str_v[11]);
		s.selflg = 0;
		stop.insert(std::make_pair(s.groupid, s));
	}


}
void clustering_2trk_vtx2_ver2(std::multimap<int, stop_track>& tracks, std::multimap<int, stop_track>& notuse, int pl, std::ofstream& ofs) {
	double refz = 0;
	for (auto itr1 = tracks.begin(); itr1 != tracks.end(); itr1++) {
		if (itr1->second.chainid == 0) {
			refz = itr1->second.z;
		}
	}

	//rawid,stop
	std::vector<track_multi> ret;
	double zrange[2] = { 0,0 };
	if (pl <= 15 || (pl >= 16 && pl % 2 == 0)) {
		zrange[0] = -1000;
	}
	else if (pl % 2 == 1) {
		zrange[0] = -3200; 
	}

	double extra[2];
	track_multi multi;
	multi.pl = pl;
	//全2trkのmd計算
	std::cout << "\tntrk = " << tracks.size() << std::endl;
	if (tracks.size() != 1) {
		for (auto itr1 = tracks.begin(); itr1 != tracks.end(); itr1++) {
			multi.eventid = itr1->second.groupid;
			multi.unixtime = itr1->second.unixtime;
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
					itr1->second.dz = pair_tmp.dz;
				}
				if (itr2->second.pid != 13 && itr1->second.pid == 13) {
					itr2->second.md = md;
					itr2->second.oa = pair_tmp.oa;
					itr2->second.dz = pair_tmp.dz;
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
			itr->second.selflg = 1;
			//k.eid = itr->second.groupid;
			multi.trk.push_back(std::make_pair(itr->first, itr->second));
		}
		for (auto itr = notuse.begin(); itr != notuse.end(); itr++) {
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
			itr->second.selflg = -1;
			multi.trk.push_back(std::make_pair(itr->first, itr->second));
		}

		ret.push_back(multi);
	}
	else {
		for (auto itr1 = tracks.begin(); itr1 != tracks.end(); itr1++) {
			multi.eventid = itr1->second.groupid;
			multi.unixtime = itr1->second.unixtime;;
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
				<< std::setw(1) << std::setprecision(0) << itr1->second.ecc << " "
				<< std::setw(5) << std::setprecision(0) << itr0->eventid << " "
				<< std::setw(10) << std::setprecision(0) << itr0->unixtime << " "
				<< std::setw(1) << std::setprecision(0) << itr1->second.bunch << " "
				<< std::setw(3) << std::setprecision(0) << itr0->pl << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.vertex_material << " "
				// information of vertex new!
				<< std::setw(3) << std::setprecision(0) << tracks.size() << " "
				<< std::setw(10) << std::setprecision(1) << itr0->x << " "
				<< std::setw(10) << std::setprecision(1) << itr0->y << " "
				<< std::setw(10) << std::setprecision(1) << itr0->z << " "
				<< std::setw(8) << std::setprecision(1) << itr0->dz << " "
				// information of vertex org
				<< std::setw(3) << std::setprecision(0) << itr1->second.ntrk_all << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.al_md << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.al_dz << " "
				// information of chain
				<< std::setw(3) << std::setprecision(0) << itr1->second.chainid << " "
				<< std::setw(4) << std::setprecision(0) << itr1->second.pid << " "
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
				<< std::setw(3) << std::setprecision(0) << itr1->second.nseg << " "
				<< std::setw(3) << std::setprecision(0) << itr1->second.npl << " "
				// information of btrk
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
				<< std::setw(8) << std::setprecision(4) << itr1->second.dz << " "
				<< std::setw(2) << std::setprecision(0) << itr1->second.selflg 
				<< std::endl;

		}
	}


}

