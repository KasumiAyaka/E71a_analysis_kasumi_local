// 2025/12/12
// Calc_vertex_information_interactive
// kasumi 


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
bool operator<(const tkey& lhs, const tkey& rhs) {
	return std::tie(lhs.eid, lhs.ip) < std::tie(rhs.eid, rhs.ip);
};
struct Track{
	int cid, pid, stopflg, veflg1, veflg2, veflg3, nseg, npl, pl0, pl1;
	double tan, mcs, rng;
	int vph, dpl, vph2, dal, dl, mu_md, mu_dz;
	int rawid;
	double ax, ay, x, y, z;
	double dar;
	double pb;
};
bool operator<(const Track& lhs, const Track& rhs) {
	return std::tie(lhs.pid, lhs.cid, lhs.stopflg, lhs.veflg1 , lhs.veflg2 , lhs.veflg3, lhs.nseg, lhs.npl, lhs.pl0, lhs.pl1, lhs.tan, lhs.mcs, lhs.rng,
		lhs.vph, lhs.dpl, lhs.vph2, lhs.dal, lhs.dar, lhs.dl, lhs.mu_md, lhs.mu_dz, lhs.rawid, lhs.ax, lhs.ay, lhs.x, lhs.y, lhs.z, lhs.pb)
		< std::tie(rhs.pid, rhs.cid, rhs.stopflg, rhs.veflg1, rhs.veflg2, rhs.veflg3, rhs.nseg, rhs.npl, rhs.pl0, rhs.pl1, rhs.tan, rhs.mcs, rhs.rng,
			rhs.vph, rhs.dpl, rhs.vph2, rhs.dal, rhs.dar, rhs.dl, rhs.mu_md, rhs.mu_dz, rhs.rawid, rhs.ax, rhs.ay, rhs.x, rhs.y, rhs.z, rhs.pb);
};

class track_pair {
public:
	double x, y, z, md, oa;
	double dz;
	Track t[2];
};

class track_multi {
public:
	int pl;
	double x, y, z;
	std::vector< std::pair<int, Track>> trk;
	std::vector<track_pair>pair;
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


void set_basetrack(std::set<Track>& map);
void clustering_2trk_vtx2_ver3(std::set<Track>& tracks);
void set_basetrack(std::vector<std::string>& str_v, std::set<Track>& map);

int main(int argc, char** argv) {
	//if (argc < 3) {
	//	fprintf(stderr, "usage:prg in-mfile.momch out-vtx.txt VElist.txt\n");
	//	exit(1);
	//}

	std::set<Track>set;
	std::cout << " Input Track information " << std::endl;
	set_basetrack(set);

	int a = 0;
	do {

		std::cout << "Quit? --> Yes :    -1\n"
			<< "      --> No  :     Input Track Information\n" << std::endl;

		std::string str;					//1strein into
		std::vector<std::string> str_v;		//input 1 ward
		std::getline(std::cin, str);
		str_v = StringSplit_with_tab(str);
		if (str_v.size() > 1) {
			a = 1;
		}
		else {
			a = -1;
		}
		std::cout << std::endl;

		if (a > 0) {
			//std::cout << " Input Track information " << std::endl;

			set_basetrack(str_v, set);
			std::cout << std::endl;

		}
	} while (a > 0);
	std::cout << " Input track : " << set.size() << std::endl;
	clustering_2trk_vtx2_ver3(set);

}

void set_basetrack(std::vector<std::string>& str_v, std::set<Track>& map) {
	Track t;
	//std::cout << "Input trk" << std::endl;

	
	std::cout << str_v.size() << std::endl;
	if (str_v.size() != 27) {
		std::cout << "argc = " << map.size() << "; 27 is correct." << std::endl;
		std::cout << "Please input Correct track infomation." << std::endl;
		return;
	}


	t.cid = std::stoi(str_v[0]);
	t.pid = std::stoi(str_v[1]);
	t.stopflg= std::stoi(str_v[2]);
	t.veflg1 = std::stoi(str_v[3]);
	t.veflg2 = std::stoi(str_v[4]);
	t.veflg3 = std::stoi(str_v[5]);
	t.nseg = std::stoi(str_v[6]);
	t.npl = std::stoi(str_v[7]);
	t.pl0 = std::stoi(str_v[8]);
	t.pl1 = std::stoi(str_v[9]);

	t.tan = std::stod(str_v[10]);
	t.mcs = std::stod(str_v[11]);
	t.rng = std::stod(str_v[12]);
	t.vph = std::stoi(str_v[13]);
	t.dpl = std::stoi(str_v[14]);
	t.vph2 = std::stoi(str_v[15]);
	t.dal = std::stod(str_v[16]);
	t.dar = std::stod(str_v[17]);
	t.dl = std::stod(str_v[18]);
	t.mu_md = std::stod(str_v[19]);
	t.mu_dz = std::stod(str_v[20]);

	t.rawid = std::stoi(str_v[21]);
	t.ax = std::stod(str_v[22]);
	t.ay = std::stod(str_v[23]);
	t.x = std::stod(str_v[24]);
	t.y = std::stod(str_v[25]);
	t.z = std::stod(str_v[26]);
	
	//std::cin >> t.cid >> t.pid >> t.stopflg >> t.veflg >> t.nseg >> t.npl >> t.pl0 >> t.pl1
	//	>> t.tan >> t.mcs >> t.rng >> t.vph >> t.dpl >> t.vph2 >> t.dal >> t.dl >> t.mu_md >> t.mu_dz
	//	>> t.rawid >> t.ax >> t.ay >> t.x >> t.y >> t.z;
	map.insert(t);
	std::cout << "ok." << std::endl;
	//std::cout << t.ax << " " << t.ay << " " << t.x << " "<<t.y << " " << t.z << std::endl;
}
void set_basetrack(std::set<Track>& map) {
	Track t;
	//std::cout << "Input trk" << std::endl;

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;
	std::getline(std::cin, str);
	str_v = StringSplit_with_tab(str);

	std::cout << str_v.size() << std::endl;
	if (str_v.size() != 27) {
		std::cout << "argc = " << map.size() << "; 27 is correct." << std::endl;
		std::cout << "Please input Correct track infomation." << std::endl;
		return;
	}


	t.cid = std::stoi(str_v[0]);
	t.pid = std::stoi(str_v[1]);
	t.stopflg = std::stoi(str_v[2]);
	t.veflg1 = std::stoi(str_v[3]);
	t.veflg2 = std::stoi(str_v[4]);
	t.veflg3 = std::stoi(str_v[5]);
	t.nseg = std::stoi(str_v[6]);
	t.npl = std::stoi(str_v[7]);
	t.pl0 = std::stoi(str_v[8]);
	t.pl1 = std::stoi(str_v[9]);

	t.tan = std::stod(str_v[10]);
	t.mcs = std::stod(str_v[11]);
	t.rng = std::stod(str_v[12]);
	t.vph = std::stoi(str_v[13]);
	t.dpl = std::stoi(str_v[14]);
	t.vph2 = std::stoi(str_v[15]);
	t.dal = std::stod(str_v[16]);
	t.dar = std::stod(str_v[17]);
	t.dl = std::stod(str_v[18]);
	t.mu_md = std::stod(str_v[19]);
	t.mu_dz = std::stod(str_v[20]);

	t.rawid = std::stoi(str_v[21]);
	t.ax = std::stod(str_v[22]);
	t.ay = std::stod(str_v[23]);
	t.x = std::stod(str_v[24]);
	t.y = std::stod(str_v[25]);
	t.z = std::stod(str_v[26]);

	//std::cin >> t.cid >> t.pid >> t.stopflg >> t.veflg >> t.nseg >> t.npl >> t.pl0 >> t.pl1
	//	>> t.tan >> t.mcs >> t.rng >> t.vph >> t.dpl >> t.vph2 >> t.dal >> t.dl >> t.mu_md >> t.mu_dz
	//	>> t.rawid >> t.ax >> t.ay >> t.x >> t.y >> t.z;
	map.insert(t);
	std::cout << "ok." << std::endl;
	//std::cout << t.ax << " " << t.ay << " " << t.x << " "<<t.y << " " << t.z << std::endl;
}

void clustering_2trk_vtx2_ver3(std::set<Track>& tracks) {

	//rawid,stop
	Track t;
	int mucnt = 0;
	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	for (auto itr = tracks.begin(); itr != tracks.end(); itr++) {
		if (itr->pid == 13) {
			t.pl1 = itr->pl1;
			t.z = itr->z;
			mucnt++;
		}
	}
	if (mucnt < 1) {
		std::cout << " Input Muon Track information " << std::endl;
		std::getline(std::cin, str);
		str_v = StringSplit_with_tab(str);
		t.pl1 = std::stoi(str_v[9]);
		t.z = std::stod(str_v[26]);
		std::cout << std::endl;
	}

	std::vector<track_multi> ret;
	double zrange[2] = { -3200,100 };
	if (t.pl1 <= 15 || (t.pl1 >= 16 && t.pl1 % 2 == 0)) {
		zrange[0] = -1000;
	}
	else if (t.pl1 % 2 == 1) {
		zrange[0] = -3200; //3500->3200
	}


	double extra[2];
	track_multi multi;
	multi.pl = t.pl1;
	//全2trkのmd計算
	std::cout << tracks.size() << std::endl;
	matrix_3D::vector_3D pos0, pos1, dir0, dir1;
	double ip;

	if (tracks.size() != 1) {
		for (auto itr1 = tracks.begin(); itr1 != tracks.end(); itr1++) {
			for (auto itr2 = std::next(itr1, 1); itr2 != tracks.end(); itr2++) {
				pos0.x = itr1->x;
				pos0.y = itr1->y;
				pos0.z = itr1->z;
				pos1.x = itr2->x;
				pos1.y = itr2->y;
				pos1.z = itr2->z;
				dir0.x = itr1->ax;
				dir0.y = itr1->ay;
				dir0.z = 1;
				dir1.x = itr2->ax;
				dir1.y = itr2->ay;
				dir1.z = 1;

				// pos0を基準にzrangeの範囲内で最近接距離をとる位置(extra)を探索
				double md = minimum_distance_fixed(pos0, pos1, dir0, dir1, zrange, extra, t.z);
				track_pair pair_tmp;
				matrix_3D::vector_3D extra0 = addition(pos0, const_multiple(dir0, extra[0]));
				matrix_3D::vector_3D extra1 = addition(pos1, const_multiple(dir1, extra[1]));

				pair_tmp.x = (extra0.x + extra1.x) / 2;
				pair_tmp.y = (extra0.y + extra1.y) / 2;
				pair_tmp.z = (extra0.z + extra1.z) / 2;
				pair_tmp.dz = pair_tmp.z - t.z;
				pair_tmp.md = md;
				pair_tmp.oa = matrix_3D::opening_angle(dir0, dir1);
				pair_tmp.t[0] = *itr1;
				pair_tmp.t[1] = *itr2;

				std::cout << std::fixed << std::right
					<< std::setw(3) << std::setprecision(0) << itr1->cid << " - " << std::setw(3) << std::setprecision(0) << itr2->cid << "\n"
					<< "   md    = " << std::setw(10) << std::setprecision(1) << pair_tmp.md << "\n"
					<< "   oa    = " << std::setw(10) << std::setprecision(5) << pair_tmp.oa << "\n"
					<< "   dz    = " << std::setw(10) << std::setprecision(1) << pair_tmp.dz << ", "
					<< "   point = ("
					<< std::setw(10) << std::setprecision(1) << pair_tmp.x << ","
					<< std::setw(10) << std::setprecision(1) << pair_tmp.y << ","
					<< std::setw(10) << std::setprecision(1) << pair_tmp.z << ")\n"
					<< std::endl;

				multi.pair.push_back(pair_tmp);
			}
		}


		matrix_3D::vector_3D p_vtx, pos, dir;
		tkey k;
		double ip;
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
		multi.dz = multi.z - t.z;
		//各trkに対してIPの計算

		for (auto itr = tracks.begin(); itr != tracks.end(); itr++) {
			pos0.x = itr->x;
			pos0.y = itr->y;
			pos0.z = itr->z;
			pos1.x = multi.x;
			pos1.y = multi.y;
			pos1.z = multi.z;
			dir0.x = itr->ax;
			dir0.y = itr->ay;
			dir0.z = 1;
			ip = matrix_3D::inpact_parameter(pos0, dir0, pos1);

			std::cout << std::fixed << std::right
				<< "   ip    = " << std::setw(10) << std::setprecision(1) << ip << "   (" << std::setw(3) << std::setprecision(0) << itr->cid << "," << std::setw(4) << std::setprecision(0) << itr->pid << ")" << std::endl;
		}
	}
	std::cout << std::fixed << std::right
		<< "   point = ("
		<< std::setw(10) << std::setprecision(1) << multi.x << ","
		<< std::setw(10) << std::setprecision(1) << multi.y << ","
		<< std::setw(10) << std::setprecision(1) << multi.z << ")\n"
		<< "   dz    = "
		<< std::setw(10) << std::setprecision(1) << multi.z - t.z 
		<< std::endl;

	//std::cout << " Calculate ip to this vtx point? \n"
	//	<< " ---------------------------------------\n"
	//	<< "  Yes  -->   1\n"
	//	<< "  No   -->  -1\n"
	//	<< std::endl;
	int a = 1;
	//std::getline(std::cin, str);
	//str_v = StringSplit_with_tab(str);
	////std::cout << str_v.size() << std::endl;
	//a = std::stoi(str_v[0]);
	//std::cout << std::endl;


	if (a == 1) {
		pos1.x = multi.x;
		pos1.y = multi.y;
		pos1.z = multi.z;
		std::cout << " Input  Track information " << std::endl;

		std::getline(std::cin, str);
		str_v = StringSplit_with_tab(str);
		t.cid = std::stoi(str_v[0]);
		t.pid = std::stoi(str_v[1]);
		t.stopflg = std::stoi(str_v[2]);
		t.veflg1 = std::stoi(str_v[3]);
		t.veflg2 = std::stoi(str_v[4]);
		t.veflg3 = std::stoi(str_v[5]);
		t.nseg = std::stoi(str_v[6]);
		t.npl = std::stoi(str_v[7]);
		t.pl0 = std::stoi(str_v[8]);
		t.pl1 = std::stoi(str_v[9]);

		t.tan = std::stod(str_v[10]);
		t.mcs = std::stod(str_v[11]);
		t.rng = std::stod(str_v[12]);
		t.vph = std::stoi(str_v[13]);
		t.dpl = std::stoi(str_v[14]);
		t.vph2 = std::stoi(str_v[15]);
		t.dal = std::stod(str_v[16]);
		t.dar = std::stod(str_v[17]);
		t.dl = std::stod(str_v[18]);
		t.mu_md = std::stod(str_v[19]);
		t.mu_dz = std::stod(str_v[20]);

		t.rawid = std::stoi(str_v[21]);
		t.ax = std::stod(str_v[22]);
		t.ay = std::stod(str_v[23]);
		t.x = std::stod(str_v[24]);
		t.y = std::stod(str_v[25]);
		t.z = std::stod(str_v[26]);
		std::cout << std::endl;

		pos0.x = t.x;
		pos0.y = t.y;
		pos0.z = t.z;
		dir0.x = t.ax;
		dir0.y = t.ay;
		dir0.z = 1;
		ip = matrix_3D::inpact_parameter(pos0, dir0, pos1);
		//std::cout << std::fixed << std::right
		//	<< "   point  = ("
		//	<< std::setw(10) << std::setprecision(1) << multi.x << ","
		//	<< std::setw(10) << std::setprecision(1) << multi.y << ","
		//	<< std::setw(10) << std::setprecision(1) << multi.z << ")"
		//	<< std::endl;
		std::cout << std::fixed << std::right
			<< "   ip = " << std::setw(10) << std::setprecision(1) << ip << ", \n"
			<< std::endl;

	}
}