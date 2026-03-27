// 026/3/2
#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
# include <iostream>
# include <random>

#include <list>

struct tkey {
	int utime;
	double bunch;
};
bool operator<(const tkey& lhs, const tkey& rhs) {
	return std::tie(lhs.utime,lhs.bunch) < std::tie(rhs.utime, rhs.bunch);
}

class stop_track {
public:
	int64_t chainid, groupid;
	int  nseg, npl, pl0, pl1, vph, rawid;
	double ax, ay, x, y, z;
	// ph-->pid
	int stoppl;
	int64_t unixtime;
	int bunch;
	double ip;
	int pid;
	double  mulikelihood, pliklihoood, weight;
	double vx, vy, vz;//for drbag
	int ntrk;
	int charge;
	int BMmomflg;// 
	double BMrng, BMrngerr[2];//rng
	double BMcurv, BMcurverr[2];//curverture
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


struct Shartingfile_version2 {
	int timestamp, bunch;
	double charge;
	int BMmomflg;// 
	double BMrng, BMrngerr[2];//rng
	double BMcurv, BMcurverr[2];//curverture
	double eccoss, ossfix, fixtss, tssst;//chi2
	int eventid;
};

double minimum_distance_fixed(matrix_3D::vector_3D pos0, matrix_3D::vector_3D pos1, matrix_3D::vector_3D dir0, matrix_3D::vector_3D dir1, double z_range[2], double extra[2], double refz) {
	double extra0_distance, extra1_distance, delta;
	matrix_3D::vector_3D pos;
	pos.x = pos1.x - pos0.x;
	pos.y = pos1.y - pos0.y;
	pos.z = pos1.z - pos0.z;
	//ÇŸÇ⁄ïΩçsÇ»èÍçá
	if (opening_angle(dir0, dir1) < 0.0001) {
		extra0_distance = (pos1.z + pos0.z) / 2 - pos0.z;
		extra1_distance = (pos1.z + pos0.z) / 2 - pos1.z;
	}
	else {
		delta = dot(dir0, dir0) * dot(dir1, dir1) - pow(dot(dir0, dir1), 2.);
		extra0_distance = (+1 * dot(pos, dir0) * dot(dir1, dir1) - dot(dir0, dir1) * dot(pos, dir1)) / delta;
		extra1_distance = (-1 * dot(pos, dir1) * dot(dir0, dir0) + dot(dir0, dir1) * dot(pos, dir0)) / delta;
	}
	//range[0]:è¨,range[1]:ëÂ
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

void read_sharingfile(std::string in_sf, std::map<int, Shartingfile_version2>& tracks);
void read_stop(std::vector<Momentum_recon::Event_information>& momch, std::map<int, stop_track>& tracks, double tan_thr);
void Matching(std::map<int, stop_track>& tracks, std::map<int, Shartingfile_version2>& sfv2, std::ofstream& ofs, int ecc, int mode);

int main(int argc, char** argv) {
	if (argc < 5) {
		fprintf(stderr, "usage:prg in.momch in_sf_ver2.txt out.txt #ECC [output log=1]\n");
		exit(1);
	}
	std::string in_momch = argv[1];
	std::string in_sf = argv[2];
	std::string out_txt = argv[3];
	int ecc = std::stoi(argv[4]);
	int mode = -1;
	if (argc == 6) {
		mode = std::stoi(argv[5]);
	}


	// read momch
	//std::cout << "\tRead\n" << in_momch << std::endl;
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);
	std::map<int, stop_track>momchtrk;
	read_stop(momch, momchtrk,4);

	// read sharing file
	//std::cout << "\tRead\n" << in_sf << std::endl;
	std::map<int, Shartingfile_version2>sftrk;
	read_sharingfile(in_sf, sftrk);

	std::ofstream ofs(out_txt);

	Matching(momchtrk, sftrk, ofs, ecc, mode);


}
void read_sharingfile(std::string in_sf, std::map<int, Shartingfile_version2>& tracks) {

	std::ifstream ifs(in_sf);
	if (!ifs) {
		std::cerr << "\tfile open error!" << std::endl;
		exit(0);
	}

	tkey tk;
	Shartingfile_version2 sfv2;

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	int recon, tr, k, bpl, dflg, nbtk;
	int flg = 0;
	while (std::getline(ifs, str)) {
		str_v = StringSplit_with_tab(str);
		//std::cout << str << std::endl;
		if (str_v.size() != 30) {
			std::cerr << "\tUse sharing file version2 by matsuo." << std::endl;
			std::cerr << str_v.size() << std::endl;
			break;
		}


		sfv2.timestamp = std::stoll(str_v[10]);
		sfv2.bunch = std::stoi(str_v[12]);
		sfv2.charge = std::stoi(str_v[15]);
		sfv2.BMmomflg = std::stoi(str_v[16]);
		sfv2.BMrng = std::stod(str_v[17]);
		sfv2.BMrngerr[0] = std::stod(str_v[18]);
		sfv2.BMrngerr[1] = std::stod(str_v[19]);
		sfv2.BMcurv = std::stod(str_v[20]);
		sfv2.BMcurverr[0] = std::stod(str_v[21]);
		sfv2.BMcurverr[1] = std::stod(str_v[22]);
		sfv2.eccoss = std::stod(str_v[23]);
		sfv2.ossfix = std::stod(str_v[24]);
		sfv2.fixtss = std::stod(str_v[25]);
		sfv2.tssst = std::stod(str_v[26]);
		sfv2.eventid = std::stoi(str_v[27]);

		tk.utime = sfv2.timestamp;
		tk.bunch = sfv2.bunch;

		//tracks.insert(std::make_pair(tk, sfv2));
		tracks.insert(std::make_pair(sfv2.eventid, sfv2));
	}

	std::cout << " * # of Muon prediction(sf ver2) = " << tracks.size() << std::endl;
}

void read_stop(std::vector<Momentum_recon::Event_information>& momch, std::map<int, stop_track>& tracks, double tan_thr) {
	stop_track stop_tmp;

	int pos = 0;
	int cnt = 0;
	int cutnum = 0;
	std::cout << " * # of momch track = " << momch.size() << std::endl;
	for (auto& ev : momch) {
		// event header
		stop_tmp.stoppl = ev.vertex_pl;
		stop_tmp.groupid = ev.groupid;
		stop_tmp.unixtime = ev.unix_time;
		stop_tmp.bunch = ev.weight;

		for (auto& c : ev.chains) {
			if (c.chainid == 0) {
				// chain header
				stop_tmp.chainid = c.chainid;
				stop_tmp.nseg = c.base.size();
				stop_tmp.pl0 = c.base.begin()->pl;//downstream
				stop_tmp.pl1 = c.base.rbegin()->pl;// upstream
				stop_tmp.npl = stop_tmp.pl1 - stop_tmp.pl0 + 1;
				stop_tmp.pid = c.particle_flg;
				stop_tmp.BMrng = c.ecc_mcs_mom[0];//assume mu 
				stop_tmp.BMrngerr[0] = c.ecc_mcs_mom_error[0][0];//assume mu 
				stop_tmp.BMrngerr[1] = c.ecc_mcs_mom_error[0][1];//assume mu 
				if (c.particle_flg == 2212) {
					stop_tmp.BMrng = c.ecc_mcs_mom[1];//assume p
					//c.ecc_mcs_mom_error[1][0]//+ [1][1] //-
					stop_tmp.BMrngerr[0] = c.ecc_mcs_mom_error[1][0];//assume mu 
					stop_tmp.BMrngerr[1] = c.ecc_mcs_mom_error[1][1];//assume mu 
				}
				stop_tmp.mulikelihood = c.muon_likelihood;
				stop_tmp.pliklihoood = c.proton_likelihood;
				//std::cout << c.particle_flg << " " << c.muon_likelihood << " " << c.proton_likelihood << std::endl;

				if (stop_tmp.pl1 >= stop_tmp.stoppl) {
					//stop
					pos = stop_tmp.pl1;
					stop_tmp.rawid = c.base.rbegin()->rawid;
					stop_tmp.ax = c.base.rbegin()->ax;
					stop_tmp.ay = c.base.rbegin()->ay;
					stop_tmp.x = c.base.rbegin()->x;
					stop_tmp.y = c.base.rbegin()->y;
					stop_tmp.z = c.base.rbegin()->z;
					stop_tmp.vph = c.base.rbegin()->m[0].ph % 10000 + c.base.rbegin()->m[1].ph % 10000;
				}
				if (stop_tmp.pl0 > stop_tmp.stoppl) {
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
				if (sqrt(stop_tmp.ax * stop_tmp.ax + stop_tmp.ay * stop_tmp.ay) <= tan_thr) {
					tracks.insert(std::make_pair(stop_tmp.groupid, stop_tmp));
					cutnum++;
				}
			}
		}
		cnt++;
	}

	//printf("input fin %d track\n", cnt);
	std::cout << " \t* Angle Cut (tan<=4): " << cnt << "-->" << cutnum << std::endl;
}

void Matching(std::map<int, stop_track>& tracks, std::map<int, Shartingfile_version2>& sfv2, std::ofstream& ofs, int ecc, int mode) {
	double refz = 0; int utime;

	tkey k;
	// for muon selection
	int cnt = 0;
	for (auto itr0 = tracks.begin(); itr0 != tracks.end(); itr0++) {

		//if (itr1->second.pid % 10000 != 13)continue;
		//std::cout << "\tPID : " << itr1->second.pid % 10000 << std::endl;

		//if (itr0->second.x < 40000 || itr0->second.x > 210000 || itr0->second.y < 40000 || itr0->second.y > 210000) {
		//	//std::cout << " \t* Area Cut : " << itr0->second.groupid << std::endl;
		//	continue;
		//}
		cnt++;

		k.bunch = itr0->second.bunch;
		k.utime = itr0->second.unixtime;

		auto itr_sf = sfv2.find((itr0->first));
		if (itr_sf == sfv2.end())continue;

		//if (itr_sf->second.eventid != itr0->eventid) {
		//	std::cout << "MOMCH : " << itr0->eventid << ", ninjasf : " << itr_sf->second.eventid << std::endl;
		//}
		//std::cout << " * Matched EVENT : " << itr0->second.groupid << std::endl;
		if (itr_sf->second.timestamp != itr0->second.unixtime && itr_sf->second.bunch != itr0->second.bunch) {
			std::cout << "\tEVENT "<< itr0->second.groupid<< "MOMCH   : " << itr0->second.unixtime << ", " << itr0->second.bunch << ";ninjasf : " << itr_sf->second.timestamp << ", " << itr_sf->second.bunch << std::endl;
		}

		ofs << std::right << std::fixed
			// information of vertex
			<< std::setw(10) << std::setprecision(0) << itr_sf->second.timestamp << " "
			<< std::setw(2) << std::setprecision(4) << itr_sf->second.bunch << " "// bunch
			<< std::setw(5) << std::setprecision(0) << itr0->second.groupid << " "
			<< std::setw(3) << std::setprecision(0) << itr0->second.ntrk << " "

			<< std::setw(3) << std::setprecision(0) << itr0->second.pl1 << " "//vpl
			// information of track
			<< std::setw(3) << std::setprecision(0) << itr0->second.nseg << " "
			<< std::setw(3) << std::setprecision(0) << itr0->second.npl << " "

			<< std::setw(5) << std::setprecision(0) << itr0->second.pid << " "
			<< std::setw(8) << std::setprecision(4) << itr0->second.ax << " "
			<< std::setw(8) << std::setprecision(4) << itr0->second.ay << " "
			<< std::setw(10) << std::setprecision(1) << itr0->second.x << " "
			<< std::setw(10) << std::setprecision(1) << itr0->second.y << " "
			//mcs
			<< std::setw(10) << std::setprecision(1) << itr0->second.BMrng << " "
			<< std::setw(10) << std::setprecision(1) << itr0->second.BMrngerr[0] << " "
			<< std::setw(10) << std::setprecision(1) << itr0->second.BMrngerr[1] << " "
			// sf information
			<< std::setw(2) << std::setprecision(0) << itr_sf->second.charge << " "
			<< std::setw(2) << std::setprecision(0) << itr_sf->second.BMmomflg << " "
			<< std::setw(10) << std::setprecision(1) << itr_sf->second.BMrng << " "
			<< std::setw(10) << std::setprecision(1) << itr_sf->second.BMrngerr[0] << " "
			<< std::setw(10) << std::setprecision(1) << itr_sf->second.BMrngerr[1] << " "
			<< std::setw(10) << std::setprecision(1) << itr_sf->second.BMcurv << " "
			<< std::setw(10) << std::setprecision(1) << itr_sf->second.BMcurverr[0] << " "
			<< std::setw(10) << std::setprecision(1) << itr_sf->second.BMcurverr[1] << " "

			<< std::setw(6) << std::setprecision(3) << itr_sf->second.eccoss << " "
			<< std::setw(6) << std::setprecision(3) << itr_sf->second.ossfix << " "
			<< std::setw(6) << std::setprecision(3) << itr_sf->second.fixtss << " "
			<< std::setw(6) << std::setprecision(3) << itr_sf->second.tssst << " "

			<< std::setw(10) << std::setprecision(0) << itr_sf->second.timestamp%1500000000 << 0<<itr_sf->second.bunch

			<< std::endl;



	}
	//std::cout << " * Area Cut : " << tracks.size()<< "-->" << cnt << std::endl;
	std::cout << " * output : " <<  cnt << std::endl;

}