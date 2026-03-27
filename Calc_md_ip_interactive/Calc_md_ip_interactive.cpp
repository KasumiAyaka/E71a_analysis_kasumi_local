#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
# include <iostream>
# include <random>

#include <list>

struct Key {
	int pl, cid, gid;
};
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.gid, lhs.cid, lhs.pl) < std::tie(rhs.gid, rhs.cid, rhs.pl);
}

struct Btrk {
	int vph, rid;
	double x, y, z, ax, ay;
};

struct VTXInfo {
	int gid,ntrk,vpl;
	double vx, vy, vz, dz;
};
struct ChainInfo {
	int cid, nseg, npl, stop_flg;
};

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
	double mom, rng;
	double mu_md;
	int stop_flg;
	double dl, dal, dr, dar;
	int d_pl;
	int vph2;//íºãﬂÇÃéüÇÃbasetrack
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

void read_stop_txt(std::vector<Momentum_recon::Event_information>& momch, std::multimap<Key, Btrk>& trk);
void Culc_MD(std::multimap<Key, Btrk>& trk, Key k1, Key k2);
void Calc_MinimumDistance(std::multimap<Key, Btrk>& trk, Key k[2]);
void Calc_MinimumDistance_zrange(std::multimap<Key, Btrk>& trk, Key k[2]);
void Calc_AngDiff(std::multimap<Key, Btrk>& trk, Key k[2]);
void Calc_PosDiff(std::multimap<Key, Btrk>& trk, Key k[2]);


int main(int argc, char** argv) {
	if (argc !=2) {
		fprintf(stderr, "usage:prg in-mfile.momch mode\n");
		exit(1);
	}
	std::string in_momch = argv[1];
	int mode = 1;
	//read momch
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);

	// reading stop.txt
	std::multimap<Key, Btrk> trk;
	read_stop_txt(momch, trk);
	std::cout << "fin set tracks." << std::endl;


	int a = 0;
	Key k[2];

	std::cout << std::setfill(' ');
	std::cout << "\n Input number : \n"
		<< "   finish    --> -1\n"
		<< "   continue  -->  1"
		<< std::endl;
	std::cin >> a;
	std::cout << std::endl;

	do {

		if (a > 0) {
			std::cout << std::endl;
			std::cout << " Select mode : \n"
				<< "  finish                         --> -1\n"
				<< " ---------------------------------------\n"
				<< "  Calculation\n"
				<< "    * MD,IP of pair and track    -->  1\n"
				<< "    * Angle Distane              -->  2\n"
				<< "    * Position Distane           -->  3\n"
				<< "    * Minimum Distane            -->  4\n"
				<< "    * MD with z-range striction  -->  5\n"
				<< std::endl;

			std::cin >> mode;
			a = mode;
			std::cout << std::endl;

			if (mode > 0 && mode < 5) {
				std::cout << " Input basetrack1 :  [group-id] [chain-id] [pl]" << std::endl;
				std::cin >> k[0].gid >> k[0].cid >> k[0].pl;
				std::cout << "\n Input basetrack2 : [group-id] [chain-id] [pl]" << std::endl;
				std::cin >> k[1].gid >> k[1].cid >> k[1].pl;
				std::cout << std::endl;

				if (trk.find(k[0]) == trk.end() || trk.find(k[1]) == trk.end()) {
					std::cout << "\n   not exist that basetrack!\n" << std::endl;
					mode = 100;
				}

			}

			if (mode == 1) {
				std::cout << " * IP of pair and track " << std::endl;
				Culc_MD(trk, k[0], k[1]);
			}
			if (mode == 2) {
				std::cout << " * Calc md " << std::endl;
				Calc_MinimumDistance(trk, k);
			}
			else if (mode == 3) {
				std::cout << " * Calc andle diff " << std::endl;
				Calc_AngDiff(trk, k);
			}
			else if (mode == 4) {
				std::cout << " * Calc position diff " << std::endl;
				Calc_PosDiff(trk, k);
			}
			else if (mode == 5) {
				std::cout << " * Calc md (with Z-Range) " << std::endl;
				Culc_MD(trk, k[0], k[1]);
			}
			else if (mode < 0) {
				a == -1;
			}
			else {
				mode = 100;
			}

		}

	} while (a > 0);

	std::cout << " Finish process...!" << std::endl;

}
void read_stop_txt(std::vector<Momentum_recon::Event_information>& momch, std::multimap<Key, Btrk>& trk) {
	Btrk btrk;
	Key k;

	int pos = 0;
	int cnt = 0;
	double ax, ay, angle, dax, day, dal, dar, dl, dr, dx, dy;
	for (auto& ev : momch) {
		k.gid = ev.groupid;
		for (auto& c : ev.chains) {
			k.cid = c.chainid;
			for (auto& b : c.base) {
				k.pl = b.pl;
				btrk.ax = b.ax;
				btrk.ay = b.ay;
				btrk.rid = b.rawid;
				btrk.vph = b.m[0].ph % 10000 + b.m[1].ph % 10000;
				btrk.x = b.x;
				btrk.y = b.y;
				btrk.z = b.z;
				trk.insert(std::make_pair(k, btrk));
			}
		}
	}
	//printf("input fin %d track\n", cnt);
}

void Culc_MD(std::multimap<Key, Btrk>& trk, Key k1, Key k2) {
	Key k3;
	int pl;
	double refz;
	if (k1.pl > k2.pl) {
		pl = k2.pl;
		k3 = k2;
		k3.cid = 0;
	}
	else {
		pl = k1.pl;
		k3 = k1;
		k3.cid = 0;
	}
	refz = trk.find(k3)->second.z;

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

	auto itr1 = trk.find(k1);
	auto itr2 = trk.find(k2);
	// 2trkÇÃmdåvéZ
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

	// pos0ÇäÓèÄÇ…zrangeÇÃîÕàÕì‡Ç≈ç≈ãﬂê⁄ãóó£ÇÇ∆ÇÈà íu(extra)ÇíTçı
	track_pair pair_tmp;
	pair_tmp.md = minimum_distance_fixed(pos0, pos1, dir0, dir1, zrange, extra, refz);
	matrix_3D::vector_3D extra0 = addition(pos0, const_multiple(dir0, extra[0]));
	matrix_3D::vector_3D extra1 = addition(pos1, const_multiple(dir1, extra[1]));

	pair_tmp.x = (extra0.x + extra1.x) / 2;
	pair_tmp.y = (extra0.y + extra1.y) / 2;
	pair_tmp.z = (extra0.z + extra1.z) / 2;
	pair_tmp.dz = pair_tmp.z - refz;
	pair_tmp.oa = matrix_3D::opening_angle(dir0, dir1);

	std::cout << std::fixed << std::right
		<< "   trk1   : GroupID"
		<< std::setw(5) << std::setprecision(0) << k1.gid << ", ChainID"
		<< std::setw(5) << std::setprecision(0) << k1.cid << ", PL"
		<< std::setw(5) << std::setprecision(1) << k1.pl << ", RawID"
		<< std::setw(5) << std::setprecision(0) << itr1->second.rid << "\n"
		<< "   trk2   : GroupID" << std::setw(10) << std::setprecision(5) << pair_tmp.oa << "\n"
		<< std::setw(5) << std::setprecision(0) << k2.gid << ", ChainID"
		<< std::setw(5) << std::setprecision(0) << k2.cid << ", PL"
		<< std::setw(5) << std::setprecision(1) << k2.pl << ", RawID"
		<< std::setw(5) << std::setprecision(0) << itr2->second.rid << "\n"
		<< std::endl;

	std::cout << std::fixed << std::right
		<< "   md    = " << std::setw(10) << std::setprecision(1) << pair_tmp.md << "\n"
		<< "   oa    = " << std::setw(10) << std::setprecision(5) << pair_tmp.oa << "\n"
		<< "   dz    = " << std::setw(10) << std::setprecision(1) << pair_tmp.dz << "\n"
		<< "   point = ("
		<< std::setw(10) << std::setprecision(1) << pair_tmp.x << ","
		<< std::setw(10) << std::setprecision(1) << pair_tmp.y << ","
		<< std::setw(10) << std::setprecision(1) << pair_tmp.z << ")\n"
		<< std::endl;

	int a = 1;

	do {

		std::cout << " Calculate ip to this vtx? \n"
			<< " ---------------------------------------\n"
			<< "  Yes  -->   1\n"
			<< "  No   -->  -1\n"
			<< std::endl;

		std::cin >> a;
		if (a > 0) {
			std::cout << std::endl;
			std::cout << " input groupid chainid pl " << std::endl;
			std::cin >> k3.gid >> k3.cid >> k3.pl;
			auto trk3 = trk.find(k3);
			std::cout << std::endl;

			matrix_3D::vector_3D p_vtx, pos, dir;
			//äetrkÇ…ëŒÇµÇƒIPÇÃåvéZ
			matrix_3D::vector_3D vtx;
			double ip;
			pos0.x = trk3->second.x;
			pos0.y = trk3->second.y;
			pos0.z = trk3->second.z;
			vtx.x = pair_tmp.x;
			vtx.y = pair_tmp.y;
			vtx.z = pair_tmp.z;
			dir0.x = trk3->second.ax;
			dir0.y = trk3->second.ay;
			dir0.z = 1;
			ip = matrix_3D::inpact_parameter(pos0, dir0, vtx);
			std::cout << "   ip = " << ip << std::endl;
			std::cout << std::endl;
		}

	} while (a > 0);

}

void Calc_MinimumDistance(std::multimap<Key, Btrk>& trk, Key k[2]) {

	if (trk.count(k[0]) != 1) {
		std::cerr << "duplicated!" << std::endl;
		return;
	}

	auto itr1 = trk.find(k[0]);
	auto itr2 = trk.find(k[1]);

	std::cout << " ( PL, rawid )" << std::endl;
	std::cout << std::fixed << std::right
		<< "trk1 :  (" << std::setw(3) << std::setprecision(0) << itr1->first.pl << ", "
		<< std::setw(12) << std::setprecision(0) << itr1->second.rid << ")\n"
		<< "trk2 :  (" << std::setw(3) << std::setprecision(0) << itr2->first.pl << ", "
		<< std::setw(12) << std::setprecision(0) << itr2->second.rid << ")\n"
		<< std::endl;

	matrix_3D::vector_3D pos0, pos1, dir0, dir1;

	pos0.x = itr1->second.x;
	pos0.y = itr1->second.y;
	pos0.z = itr1->second.z;
	dir0.x = itr1->second.ax;
	dir0.y = itr1->second.ay;
	dir0.z = 1;

	pos1.x = itr2->second.x;
	pos1.y = itr2->second.y;
	pos1.z = itr2->second.z;
	dir1.x = itr2->second.ax;
	dir1.y = itr2->second.ay;
	dir1.z = 1;


	double oa, md;
	double point[3];
	md = matrix_3D::minimum_distance(pos0, pos1, dir0, dir1, point);
	oa = matrix_3D::opening_angle(dir0, dir1);

	std::cout << std::fixed << std::right
		<< "   md    =  " << std::setw(10) << std::setprecision(1) << md << "\n"
		<< "   oa    =      " << std::setw(10) << std::setprecision(5) << oa << "\n"
		<< "   point = ("
		<< std::setw(10) << std::setprecision(1) << point[0] << ", "
		<< std::setw(10) << std::setprecision(1) << point[1] << ", "
		<< std::setw(10) << std::setprecision(1) << point[2] << ")\n"
		<< std::endl;

}

void Calc_MinimumDistance_zrange(std::multimap<Key, Btrk>& trk, Key k[2]) {

	if (trk.count(k[0]) != 1) {
		std::cerr << "duplicated!" << std::endl;
		return;
	}

	auto itr1 = trk.find(k[0]);
	auto itr2 = trk.find(k[1]);
	std::cout << "PL :  " << itr1->first.pl << ", " << itr2->first.pl << std::endl;

	matrix_3D::vector_3D pos0, pos1, dir0, dir1;

	pos0.x = itr1->second.x;
	pos0.y = itr1->second.y;
	pos0.z = itr1->second.z;
	dir0.x = itr1->second.ax;
	dir0.y = itr1->second.ay;
	dir0.z = 1;

	pos1.x = itr2->second.x;
	pos1.y = itr2->second.y;
	pos1.z = itr2->second.z;
	dir1.x = itr2->second.ax;
	dir1.y = itr2->second.ay;
	dir1.z = 1;


	double oa, md, dz;
	double extra[2], z_range[2];
	z_range[0] = pos0.z;
	z_range[1] = pos1.z;
	oa = matrix_3D::opening_angle(dir0, dir1);
	md = matrix_3D::minimum_distance(pos0, pos1, dir0, dir1, z_range, extra);
	dz = (fabs(extra[0]) + fabs(extra[1])) / 2;

	std::cout << std::fixed << std::right
		<< "   md    = " << std::setw(10) << std::setprecision(1) << md << "\n"
		<< "   oa    =     " << std::setw(10) << std::setprecision(5) << oa << "\n"
		<< "   dz    = " << std::setw(10) << std::setprecision(1) << dz
		<< std::endl;

}

void Calc_AngDiff(std::multimap<Key, Btrk>& trk, Key k[2]) {

	if (trk.count(k[0]) != 1) {
		std::cerr << "duplicated!" << std::endl;
		return;
	}

	auto itr1 = trk.find(k[0]);
	auto itr2 = trk.find(k[1]);
	std::cout << "PL :  " << itr1->first.pl << ", " << itr2->first.pl << std::endl;

	double ax, ay, dax, day, dlat, drad, angle;
	// dax,day
	dax = itr2->second.ax - itr1->second.ax;
	day = itr2->second.ay - itr1->second.ay;
	std::cout << std::fixed << std::right
		<< "   dax   = " << std::setw(10) << std::setprecision(5) << dax << "\n"
		<< "   day   = " << std::setw(10) << std::setprecision(5) << day
		<< std::endl;


	// d(lat)
	/*
	ax = itr1->second.ax;
	ay = itr1->second.ay;
	angle = sqrt(ax * ax + ay * ay);

	dax = itr2->second.ax - ax;
	day = itr2->second.ay - ay;
	if (angle < 0.01) {
		dlat = dax;
	}
	else {
		dlat = (dax * ay - day * ax) / angle;
	}

	// d(rad)
	if (angle < 0.01) {
		drad = day;
	}
	else {
		drad = (dax * ax + day * ay) / angle;
	}
	*/
	ax = itr2->second.ax + itr1->second.ax;
	ay = itr2->second.ay + itr1->second.ay;
	angle = sqrt(ax * ax + ay * ay);
	dlat = (-ay * dax + ax * day) / angle;
	drad = (ax * dax + ay * day) / angle;
	if (angle < 0.01) {
		dlat = dax;
		drad = day;
	}
	std::cout << std::fixed << std::right
		<< "   dal   = " << std::setw(10) << std::setprecision(5) << dlat << "\n"
		<< "   dar   = " << std::setw(10) << std::setprecision(5) << drad
		<< std::endl;

	// d(tan)
	double tan1, tan2, dtan;
	tan1 = sqrt(itr1->second.ax * itr1->second.ax + itr1->second.ay * itr1->second.ay);
	tan2 = sqrt(itr2->second.ax * itr2->second.ax + itr2->second.ay * itr2->second.ay);

	tan1 = atan(tan1);// theta [rad]
	tan2 = atan(tan2);// theta [rad]
	dtan = tan2 - tan1;// d(theta) [rad]



	std::cout << std::fixed << std::right
		<< "   dtan  = " << std::setw(10) << std::setprecision(5) << dtan << " [rad], ( tan = " << tan(dtan) << " )\n"
		<< std::endl;

}

void Calc_PosDiff(std::multimap<Key, Btrk>& trk, Key k[2]) {

	if (trk.count(k[0]) != 1) {
		std::cerr << "duplicated!" << std::endl;
		return;
	}

	auto itr1 = trk.find(k[0]);
	auto itr2 = trk.find(k[1]);
	std::cout << "PL :  " << itr1->first.pl << ", " << itr2->first.pl << std::endl;

	double extra[2], z_range[2];
	matrix_3D::vector_3D pos0, pos1, point;

	pos0.x = -itr1->second.x;
	pos0.y = -itr1->second.y;
	pos0.z = -itr1->second.z;

	pos1.x = itr2->second.x;
	pos1.y = itr2->second.y;
	pos1.z = itr2->second.z;

	point = matrix_3D::addition(pos0, pos1);
	double d = matrix_3D::distance(pos0, pos1);

	std::cout << std::fixed << std::right
		<< "   distance   =  " << std::setw(10) << std::setprecision(1) << d << "\n"
		<< "   (dx,dy,dz) = ("
		<< std::setw(10) << std::setprecision(1) << point.x << ","
		<< std::setw(10) << std::setprecision(1) << point.y << ","
		<< std::setw(10) << std::setprecision(1) << point.z << ")\n"
		<< std::endl;

}
