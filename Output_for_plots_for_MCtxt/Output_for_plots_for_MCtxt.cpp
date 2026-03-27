#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
# include <iostream>
# include <random>

#include <list>

struct tkey {
	int eid;
	double dz;
	double w;
};
bool operator<(const tkey& lhs, const tkey& rhs) {
	return std::tie(lhs.eid, lhs.w, lhs.dz) < std::tie(rhs.eid,rhs.w, rhs.dz);
}

class stop_track {
public:
	int64_t chainid, groupid;
	int  nseg, npl, pl0, pl1, vph, rawid;
	double ax, ay, x, y, z;
	// ph-->pid
	int stoppl;
	int unixtime;
	double ip;
	int pid;
	double mom, mulikelihood, pliklihoood,weight;
	double vx, vy, vz;//for drbag
	int ntrk;
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


void read_stop_txt(std::string in_momch, std::multimap<tkey, stop_track>& tracks);
void clustering_2trk_vtx2_ver2(std::multimap<int, stop_track>& tracks, int pl, std::ofstream& ofs, int ecc, int mode);

int main(int argc, char** argv) {
	if (argc < 3) {
		fprintf(stderr, "usage:prg in-mfile.all out-vtx.txt #ECC\n");
		fprintf(stderr, "usage:prg in-mfile.all out-vtx.txt #ECC 1\n-->does not display logs...\n");
		exit(1);
	}
	std::string in_momch = argv[1];
	std::string out_txt = argv[2];
	int ecc = std::stoi(argv[3]);
	int mode = -1;
	if (argc == 5) {
		mode = std::stoi(argv[4]);
	}
	//read momch

	// reading stop.txt
	std::multimap<tkey, stop_track> stop;
	read_stop_txt(in_momch, stop);
	std::cout << "fin reading stoptrack." << std::endl;

	//file消去
	std::ofstream ofs(out_txt);
	ofs.close();

	std::set<tkey> set;
	for (auto itr = stop.begin(); itr != stop.end(); itr++) {
		set.insert(itr->first);
	}
	std::cout << set.size() << std::endl;

	std::cout << "Calc Vertex point" << std::endl;
	ofs.open(out_txt);
	for (auto ev = set.begin(); ev != set.end(); ev++) {
		auto tks = stop.equal_range(*ev);
		if (mode == 1)	printf("  * event %5d, #trk = %2d\n", *ev, stop.count(*ev));
		std::multimap<int, stop_track> rid;
		for (auto itr0 = tks.first; itr0 != tks.second; itr0++) {
			rid.insert(std::make_pair(itr0->second.rawid, itr0->second));
		}
		clustering_2trk_vtx2_ver2(rid, rid.begin()->second.stoppl, ofs, ecc, mode);
		rid.clear();
	}

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

		if (str_v.size() == 19&& std::stoi(str_v[1]) <0) {
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
			//if (int(std::stoi(str_v[2]) / 10000) <= 0) {//recon true
					flg++;
			//if (k < stop_tmp.ntrk) {
				//recon
				//std::cout << "  chain header " << std::stoi(str_v[2]) << std::endl;
				stop_tmp.chainid = std::stoi(str_v[0]);
				stop_tmp.pid = std::stoi(str_v[2]);
				dflg = std::stoi(str_v[3]);
				stop_tmp.mulikelihood = std::stod(str_v[23]);
				stop_tmp.pliklihoood = std::stod(str_v[24]);
				stop_tmp.nseg = std::stoi(str_v[25]);
				stop_tmp.ip = 0;
				
				//recon
				/**/
				stop_tmp.mom = std::stod(str_v[7]); //ecc_mcs assume muon
				if (stop_tmp.pid%10000 == 2212) {
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
		if (str_v.size() == 19 && flg>0) {
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
	std::cout << "tracks = " << tracks.size() << std::endl;
}
void clustering_2trk_vtx2_ver2(std::multimap<int, stop_track>& tracks, int pl, std::ofstream& ofs, int ecc,int mode) {
	double refz = 0; int utime;
	for (auto itr1 = tracks.begin(); itr1 != tracks.end(); itr1++) {
		if (itr1->second.pid%10000 == 13) {
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
			if (mode == 1)std::cout << "\tData:( " << itr0->x << ", " << itr0->y << ", " << itr0->z << " ), MC:( " << itr1->second.vx << ", " << itr1->second.vy << ", " << itr1->second.vz << " )" << std::endl;
			ofs << std::right << std::fixed
				// information of vertex
				<< std::setw(2) << std::setprecision(0) << ecc << " "
				<< std::setw(6) << std::setprecision(0) << itr0->eventid << " "
				<< std::setw(3) << std::setprecision(0) << itr0->unixtime << " "
				<< std::setw(6) << std::setprecision(0) << itr0->pl << " "
				<< std::setw(4) << std::setprecision(0) << itr0->trk.size() << " "
				<< std::setw(10) << std::setprecision(1) << itr1->second.vx << " "
				<< std::setw(10) << std::setprecision(1) << itr1->second.vy << " "
				<< std::setw(10) << std::setprecision(1) << itr1->second.vz << " "
				<< std::setw(8) << std::setprecision(1) << itr0->dz << " "
				// information of track
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
				<< std::setw(9) << std::setprecision(4) << itr1->second.pliklihoood <<" "
				<<std::setw(9)<<std::setprecision(4)<<itr1->second.weight
				<< std::endl;
		}
	}


}