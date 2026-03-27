// 2026/02/09
// kasumi
// Add_tumestamp_to_momch
// Input SharingFile and Momch, then add timestamp (and momentum,charge) to Momch.

#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
# include <iostream>
# include <random>

#include <list>
#include <cassert>
#include <filesystem>

class stop_track {
public:
	int chainid, groupid;
	int  nseg, npl, pl0, pl1, ph, rawid, pid;
	double ax, ay, x, y, z;
	// ph-->pid
	int stoppl;
	int64_t unixtime;



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
	int timestamp, bunch;
	double BMmom, charge;
	//2026/3/2 sdd
	int BMmomflg;// 
	double BMerrp, BMerrm;//rng
	double BMcurv, BMcurverrp, BMcurverrm;

};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.timestamp, lhs.bunch, lhs.BMmom, lhs.charge, lhs.BMmomflg, lhs.BMcurv, lhs.BMerrp, lhs.BMerrm, lhs.BMcurverrp, lhs.BMcurverrm)
		< std::tie(rhs.timestamp, rhs.bunch, rhs.BMmom, rhs.charge, rhs.BMmomflg, rhs.BMcurv, rhs.BMerrp, rhs.BMerrm, rhs.BMcurverrp, rhs.BMcurverrm);
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

void read_sf(std::string in, std::map<int, Lst>& map);
void Add_information_of_sharingfile(std::vector<Momentum_recon::Event_information>& momch, std::map<int, Lst>& sfinf);



int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "usage:prg in.momch in_sf.txt  output.momch\n");
		exit(1);
	}
	std::string in_momch = argv[1];// input momch
	std::string in_sf = argv[2];
	std::string out_momch = argv[3];// output momch

	//bool result = std::filesystem::create_directories(out_momch);
	//assert(result);
	//assert(std::filesystem::exists(out_momch));//directryの存在確認
	//assert(std::filesystem::is_directory(out_momch));//指定されたパスがディレクトリを指しているかを確認する

	auto start = std::chrono::system_clock::now();//for measure working time


	// read momch
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);

	// read sharingfile
	std::map<int, Lst> sfinf;	
	read_sf(in_sf, sfinf);

	// add timestamp, muon_charge and BM_momentum
	Add_information_of_sharingfile(momch, sfinf);

	// write out
	Momentum_recon::Write_Event_information_extension(out_momch, momch);

	auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
	MeasureProcessingTime(start, end);
}

void read_sf(std::string in,std::map<int,Lst>&map) {

	std::ifstream ifs(in);
	if (!ifs) {
		std::cerr << " File open error ! " << std::endl;
		return;
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;
	int eid;
	double mom;
	Lst l;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);

		l.timestamp = std::stoll(str_v[10]);// unixtime
		l.bunch = std::stoi(str_v[12]);//bunch
		l.charge = std::stod(str_v[15]);// muon charge

		if (str_v.size() == 24) {
			l.BMmomflg = -1;
			l.BMmom = std::stod(str_v[16]);// BM momentum
			l.BMerrm = -2;
			l.BMerrp = -2;
			l.BMcurv = -2;
			l.BMcurverrm - 2;
			l.BMcurverrp = -2;
			eid = std::stoi(str_v[21]);
		}
		else if (str_v.size() == 30) {
			l.BMmomflg = std::stoi(str_v[16]);
			l.BMmom = std::stod(str_v[17]);// BM momentum
			l.BMerrm = std::stod(str_v[18]);
			l.BMerrp = std::stod(str_v[19]);
			l.BMcurv = std::stod(str_v[20]);
			l.BMcurverrm= std::stod(str_v[21]);
			l.BMcurverrp = std::stod(str_v[22]);
			eid = std::stoi(str_v[27]);
		}
		//std::cout << l.timestamp << " " << l.charge << " " << l.BMmom << std::endl;
		map.insert(std::make_pair(eid, l));

	}

	std::cout << "\t* The Number of muon prediction is " << map.size() << std::endl;

}


void Add_information_of_sharingfile(std::vector<Momentum_recon::Event_information>& momch, std::map<int, Lst>& sfinf) {
	std::cout << "\t* The Number of momch events " << momch.size() << std::endl;

	int cnt = 0;
	for (auto& ev : momch) {

		auto itr = sfinf.find(ev.groupid);
		
		ev.unix_time = itr->second.timestamp;
		//ev.true_vertex_position[0] = itr->second.bunch;
		ev.weight = itr->second.bunch;

		for (auto& c : ev.chains) {
			//if (c.particle_flg == 13) {
			if (c.chainid == 0) {
				// muon charge
				c.charge_sign = itr->second.charge;
				// BM momentum flg
				c.bm_range_mom_error[0] = itr->second.BMmomflg;
				// bm運動量のエラーは今sfは持っていないので，rangeかcurvertureどっちを使ったかのflgを代わりに入れる。
				// 0->rng;1->cuv

				// BM range momentum
				c.bm_range_mom = itr->second.BMmom;
				//c.bm_range_mom_error[0] = itr->second.BMerrm;
				c.bm_range_mom_error[1] = itr->second.BMerrp;
				// BM curverture momentum
				c.bm_curvature_mom = itr->second.BMcurv;
				c.bm_curvature_mom_error[0] = itr->second.BMcurverrm;
				c.bm_curvature_mom_error[1] = itr->second.BMcurverrp;

				cnt++;
			}
		}
	}
	std::cout << "\t* The Number of add info events " << cnt << std::endl;

}