#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
# include <iostream>
# include <random>

#include <list>

struct tkey {
	int pl;
	double vx,vy;
};
bool operator<(const tkey& lhs, const tkey& rhs) {
	return std::tie(lhs.pl, lhs.vx, lhs.vy) < std::tie(rhs.pl, rhs.vx, rhs.vy);
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
	double mom, rng;
	double mu_md, mu_dz;
	int stop_flg;
	double dl, dal, dr, dar;
	int d_pl;
	int vph2;//直近の次のbasetrack
	double pb;
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

std::set<tkey>mc_event_list(std::string input);
std::set<tkey>mc_event_list(std::vector<Momentum_recon::Event_information>& momch);

void Output_mc_true_brock(std::vector<Momentum_recon::Event_information>& momch, std::set<tkey>& key, std::string output);
void read_muon_mc_txt(std::vector<Momentum_recon::Event_information>& momch, std::set<tkey>& key, std::multimap<int, stop_track>& tracks);
void clustering_2trk_vtx2_ver4(std::multimap<int, stop_track>& tracks, int pl, std::ofstream& ofs);
void Check_true_and_recon(std::vector<Momentum_recon::Event_information>& momch, std::set<tkey>& key, std::string output);

int main(int argc, char** argv) {
	if (argc < 3) {
		fprintf(stderr, "usage:prg in_mc_event.momch in-mc_org.momch out-vtx.txt out-true-partner.txt\n");
		exit(1);
	}
	std::string in_momch = argv[1];//momch
	std::string in_event = argv[2];// event list
	std::string output_vertex = argv[3];//vtx
	std::string output_true_partner = argv[4];//trule 
	std::cout << "\t argc = " << argc << "\n\n" << std::endl;

	//read momch
	std::cout << " * Read momch(mc event)" << std::endl;
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);
	std::set<tkey> list = mc_event_list(momch);
	momch.clear();
	std::cout << "Finish reading momch. ( # of Event(list) = " << list.size() << ")\n" << std::endl;

	// output true track information
	std::cout << " * Read momch(mc original)" << std::endl;
	//std::set<tkey> list = mc_event_list(in_event);
	//Output_mc_true_brock(momch, list, output_true_partner);	
	//std::cout << " Finish output true track information." << std::endl;
	// From momch
	std::vector<Momentum_recon::Event_information> momch_org = Momentum_recon::Read_Event_information_extension(in_event);
	std::cout << "Finish reading momch. ( # of Event(mc_org)  = " << momch_org.size() << ")\n" << std::endl;


	// reading stop.txt
	std::cout << " * Set tracks" << std::endl;
	std::multimap<int, stop_track> stop;
	read_muon_mc_txt(momch_org, list, stop);// use all chain(no cut)

	// Matching
	Output_mc_true_brock(momch_org, list, output_true_partner);
	//Check_true_and_recon(momch_org, list, output_true_partner);


	// make vertex
	std::cout << " * Make vertex" << std::endl;
	//file消去
	std::ofstream ofs(output_vertex);
	ofs.close();
	ofs.open(output_vertex);
	//ofs << "stopflg 0:penetrate/sideout, 2:ecc stop, VEflg 0:remain,  >0:erase(#of matched basetrk)" << std::endl;



	std::set<int> set;
	for (auto itr = stop.begin(); itr != stop.end(); itr++) {
		set.insert(itr->first);
		//std::cout << "\t" << itr->first << std::endl;
	}
	std::cout << "\t #of event : " << set.size() << std::endl;

	for (auto ev = set.begin(); ev != set.end(); ev++) {
		auto tks = stop.equal_range(*ev);
		printf("event %5d, TrackNum = ", *ev);
		std::multimap<int, stop_track> rid;
		for (auto itr0 = tks.first; itr0 != tks.second; itr0++) {
			rid.insert(std::make_pair(itr0->second.chainid, itr0->second));
		}
		//std::cout << rid.size() << std::endl;
		//clustering_2trk_vtx2_ver3(rid, rid.begin()->second.stoppl, ofs, ve);
		clustering_2trk_vtx2_ver4(rid, rid.begin()->second.stoppl, ofs);
		rid.clear();
	}
	std::cout << " Fin." << std::endl;

}

std::set<tkey>mc_event_list(std::string input) 
{
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error...\n" << input << std::endl;
		exit(0);
	}

	tkey k;
	std::set<tkey>ret;

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	while (std::getline(ifs, str)) {
		str_v = StringSplit_with_tab(str);
		
		if (str_v.size() == 3) {
			k.pl = std::stoi(str_v[0]);
			k.vx = std::stod(str_v[1]);
			k.vy - std::stod(str_v[2]);
			ret.insert(k);// groupid chainid
		}
	}

	std::cout << "\t* The Number of listed chain is " << ret.size() << std::endl;
	return ret;
	
}

std::set<tkey>mc_event_list(std::vector<Momentum_recon::Event_information>& momch)
{
	tkey k;
	std::set<tkey> ret;
	for (auto& ev : momch) {

		k.pl = ev.vertex_pl;
		k.vx = ev.true_vertex_position[0];
		k.vy = ev.true_vertex_position[1];
		//std::cout << std::right << std::fixed
		//	<< std::setw(10) << std::setprecision(1) << k.pl << " "
		//	<< std::setw(10) << std::setprecision(1) << k.vx << " "
		//	<< std::setw(10) << std::setprecision(1) << k.vy << " "
		//	<< std::endl;

		ret.insert(k);
	}
	return ret;
}

void Output_mc_true_brock(std::vector<Momentum_recon::Event_information>& momch, std::set<tkey>& key, std::string output)
{
	std::ofstream ofs(output);
	if (!ofs) {
		std::cerr << "\tFile Open error...\n" << output << std::endl;
		return;
	}

	tkey k;
	for (auto& ev : momch) {

		k.pl = ev.vertex_pl;
		k.vx = ev.true_vertex_position[0];
		k.vy = ev.true_vertex_position[1];
		if (key.count(k) > 0) {
			// find

			//ofs << std::right << std::fixed
			//	<< std::setw(10) << std::setprecision(1) << k.pl << " "
			//	<< std::setw(10) << std::setprecision(1) << k.vx << " "
			//	<< std::setw(10) << std::setprecision(1) << k.vy << " "
			//	<< std::setw(8) << std::setprecision(1) << ev.weight << std::endl;

			ofs << std::right << std::fixed
				<< std::setw(12) << "groupid" << " "
				<< std::setw(12) << "unix_time" << " "
				<< std::setw(4) << "vpl" << " "
				<< std::setw(10) << "v_material" << " "
				<< std::setw(10) << "vx " << " "
				<< std::setw(10) << "vy" << " "
				<< std::setw(10) << "vz" << " "
				<< std::setw(5) << "recon" << " "
				<< std::setw(5) << "true" << " "
				<< std::setw(8) << "weight" << std::endl;

			ofs << std::right << std::fixed
				<< std::setw(12) << std::setprecision(0) << ev.groupid << " "
				<< std::setw(12) << std::setprecision(0) << ev.unix_time << " "
				<< std::setw(4) << std::setprecision(0) << ev.vertex_pl%1000 << " "
				<< std::setw(10) << std::setprecision(0) << ev.vertex_material << " "
				<< std::setw(10) << std::setprecision(1) << ev.vertex_position[0] << " "
				<< std::setw(10) << std::setprecision(1) << ev.vertex_position[1] << " "
				<< std::setw(10) << std::setprecision(1) << ev.vertex_position[2] << " "
				<< std::setw(5) << std::setprecision(0) << ev.chains.size() << " "
				<< std::setw(5) << std::setprecision(0) << ev.true_chains.size() << " "
				<< std::setw(8) << std::setprecision(1) << ev.weight << std::endl;

			ofs << std::right << std::fixed
				<< std::setw(12) << "chainid" << " "
				<< std::setw(12) << "particle_flg" << " "
				//<< std::setw(12) << "mu rng" << " "
				//<< std::setw(12) << "p rng" << " "
				//<< std::setw(12) << "mu mcs" << " "
				//<< std::setw(12) << "p mcs" << " "
				<< std::setw(12) << "mom" << " "//				<< std::setw(12) << "bm rng" << " "
				//<< std::setw(12) << "bm pi" << " "
				//<< std::setw(12) << "muon pb" << " "
				<< std::endl;
			for (auto& c : ev.true_chains) {

				ofs << std::right << std::fixed
					<< std::setw(12) << std::setprecision(0) << c.chainid << " "
					<< std::setw(12) << std::setprecision(0) << c.particle_flg << " "
					//<< std::setw(12) << std::setprecision(1) << c.ecc_range_mom[0] << " "
					//<< std::setw(12) << std::setprecision(1) << c.ecc_range_mom[1] << " "
					//<< std::setw(12) << std::setprecision(1) << c.ecc_mcs_mom[0] << " "
					//<< std::setw(12) << std::setprecision(1) << c.ecc_mcs_mom[1] << " "
					<< std::setw(12) << std::setprecision(1) << c.bm_range_mom << " "
					//<< std::setw(12) << std::setprecision(1) << c.bm_curvature_mom << " "
					//<< std::setw(12) << std::setprecision(1) << c.Get_muon_mcs_pb() << " "
					<< std::endl;
				//std::cout << "\tVPH " << c.base.begin()->m[0].ph << ", " << c.base.begin()->m[1].ph << std::endl;

			}
			ofs << std::endl;

		}
	}



	//printf("\t* input fin.\n\t #of track : %d --> %d\n", cnt, use);

}

void read_muon_mc_txt(std::vector<Momentum_recon::Event_information>& momch, std::set<tkey>& key, std::multimap<int, stop_track>& tracks)
{
	stop_track stop_tmp;

	int pos = 0;
	int cnt = 0;
	int use = 0;
	double ax, ay, angle, dax, day, dal, dar, dl, dr, dx, dy;
	tkey k;
	for (auto& ev : momch) {

		k.pl = ev.vertex_pl;
		k.vx = ev.true_vertex_position[0];
		k.vy = ev.true_vertex_position[1];
		cnt++;
		//std::cout << std::right << std::fixed
			//<< std::setw(10) << std::setprecision(1) << k.pl << " "
			//<< std::setw(10) << std::setprecision(1) << k.vx << " "
			//<< std::setw(10) << std::setprecision(1) << k.vy << " "
			//<< std::endl;


		if (key.count(k) > 0) {
			// find
			stop_tmp.stoppl = ev.vertex_pl % 1000;
			stop_tmp.groupid = ev.groupid;
			stop_tmp.unixtime = 10;
			use++;

			//for (auto& c : ev.true_chains) {}
			for (auto& c : ev.chains) {//recon


				//if (int(stop_tmp.pid%10000) != 13)continue;

				//std::cout << stop_tmp.groupid << " " << c.chainid << std::endl;
				if (c.base.size() == 1)continue;
				stop_tmp.chainid = c.chainid;
				stop_tmp.nseg = c.base.size();
				stop_tmp.pl0 = c.base.begin()->pl;//dounstream
				stop_tmp.pl1 = c.base.rbegin()->pl;//upstream
				stop_tmp.npl = stop_tmp.pl1 - stop_tmp.pl0 + 1;
				stop_tmp.pid = int(c.particle_flg % 10000);
				stop_tmp.mom = c.ecc_mcs_mom[0];
				stop_tmp.rng = c.ecc_range_mom[0];
				stop_tmp.stop_flg = c.stop_flg;
				if (stop_tmp.pid != 13) {
					stop_tmp.stop_flg = 2;// ECC stop
				}

				if (stop_tmp.pid == 2212) {
					stop_tmp.mom = c.ecc_mcs_mom[1];
					stop_tmp.rng = c.ecc_range_mom[1];
				}
				stop_tmp.pb = c.Get_muon_mcs_pb();
				stop_tmp.vph = c.base.begin()->m[0].imager % 10000 + c.base.begin()->m[1].imager % 10000;
				//for (auto& b : c.base) {
				//	std::cout << "\tVPH " << b.m[0].ph<<","<<b.m[0].view << "," << b.m[0].zone << "," << b.m[0].pixelnum << "," << b.m[0].imager << "," << b.m[0].hitnum<< ", " << b.m[1].ph << std::endl;
				//}
				if (stop_tmp.pl1 <= stop_tmp.stoppl) {//fwd
					// stoptrack
					if (stop_tmp.pid !=13&&c.base.begin()->pl < 4) {//penetrate
						stop_tmp.stop_flg = 0;
					}
					if (stop_tmp.pid != 13 && c.base.begin()->x < 10000 || c.base.begin()->x > 240000) {//side out
						stop_tmp.stop_flg = 0;
					}
					if (stop_tmp.pid != 13 && c.base.begin()->y < 10000 || c.base.begin()->y > 240000) {//side out
						stop_tmp.stop_flg = 0;
					}
					pos = stop_tmp.pl1;
					stop_tmp.rawid = c.base.rbegin()->rawid;
					stop_tmp.ax = c.base.rbegin()->ax;
					stop_tmp.ay = c.base.rbegin()->ay;
					stop_tmp.x = c.base.rbegin()->x;
					stop_tmp.y = c.base.rbegin()->y;
					stop_tmp.z = c.base.rbegin()->z;
					//stop_tmp.vph = c.base.rbegin()->m[0].ph % 10000 + c.base.rbegin()->m[1].ph % 10000;
					stop_tmp.vph2 = std::next(c.base.rbegin(), 1)->m[0].ph % 10000 + std::next(c.base.rbegin(), 1)->m[1].ph % 10000;

					//std::cout << c.base_pair.size() << std::endl;
					if (c.base_pair.size() != 0) {
						auto itr = c.base_pair.rbegin();
						ax = itr->first.ax;
						ay = itr->first.ay;
						angle = sqrt(ax * ax + ay * ay);
						dax = itr->second.ax - ax;
						day = itr->second.ay - ay;
						dx = itr->second.x - itr->first.x;
						dy = itr->second.y - itr->second.y;

						if (angle < 0.01) {
							dal = dax;
						}
						else {
							dal = (dax * ay - day * ax) / angle;
						}
						if (angle < 0.01) {
							dar = day;
						}
						else {
							dar = (dax * ax + day * ay) / angle;
						}
						stop_tmp.dal = dal;
						stop_tmp.dl = (dx * ay - dy * ax) / angle;
						stop_tmp.dar = dar;
						stop_tmp.dr = (dx * ax + dy * ay) / angle;
						stop_tmp.d_pl = itr->second.pl - itr->first.pl;
						//std::cout << " Angle diff : " << dal << ", " << dar << std::endl;
					}

				}
				if (stop_tmp.pl0 > stop_tmp.stoppl) {//bwd
					//start
									// bwd
					if (stop_tmp.pid != 13 && c.base.rbegin()->pl > 130) {//penetrate
						stop_tmp.stop_flg = 0;
					}
					if (stop_tmp.pid != 13 && c.base.rbegin()->x < 10000 || c.base.rbegin()->x  >240000) {//side out
						stop_tmp.stop_flg = 0;
					}
					if (stop_tmp.pid != 13 && c.base.rbegin()->y < 10000 || c.base.rbegin()->y > 240000) {//side out
						stop_tmp.stop_flg = 0;
					}
					stop_tmp.rawid = c.base.begin()->rawid;
					stop_tmp.ax = c.base.begin()->ax;
					stop_tmp.ay = c.base.begin()->ay;
					stop_tmp.x = c.base.begin()->x;
					stop_tmp.y = c.base.begin()->y;
					stop_tmp.z = c.base.begin()->z;
					//stop_tmp.vph = c.base.begin()->m[0].ph % 10000 + c.base.begin()->m[1].ph % 10000;
					stop_tmp.vph2 = std::next(c.base.begin(), 1)->m[0].ph % 10000 + std::next(c.base.begin(), 1)->m[1].ph % 10000;
					if (c.base_pair.size() != 0) {
						auto itr = c.base_pair.begin();
						ax = itr->first.ax;
						ay = itr->first.ay;
						angle = sqrt(ax * ax + ay * ay);
						dax = itr->second.ax - ax;
						day = itr->second.ay - ay;

						if (angle < 0.01) {
							dal = dax;
						}
						else {
							dal = (dax * ay - day * ax) / angle;
						}
						if (angle < 0.01) {
							dar = day;
						}
						else {
							dar = (dax * ax + day * ay) / angle;
						}
						stop_tmp.dal = dal;
						stop_tmp.dl = (dx * ay - dy * ax) / angle;
						stop_tmp.dar = dar;
						stop_tmp.dr = (dx * ax + dy * ay) / angle;
						stop_tmp.d_pl = itr->second.pl - itr->first.pl;
					}
				}

				stop_tmp.ip = 0;
				tracks.insert(std::make_pair(stop_tmp.groupid, stop_tmp));
			}
			//std::cout << use << " / " << cnt << std::endl;

		}
	}


	printf("\t* input fin.\n\t #of event : %d --> %d\n", cnt, use);
	if (use == 0) {
		std::cout << "This momch doesn't have indicated events." << std::endl;
		exit(0);
	}
}
void clustering_2trk_vtx2_ver4(std::multimap<int, stop_track>& tracks, int pl, std::ofstream& ofs) 
{
	double refz = 0; int utime;
	for (auto itr1 = tracks.begin(); itr1 != tracks.end(); itr1++) {
		if (itr1->second.pid == 13) {
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

	tkey k;
	int veflg = 0;
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
				pair_tmp.t[0] = itr1->second;
				pair_tmp.t[1] = itr2->second;

				if (itr1->second.pid == 13) {
					itr2->second.mu_md = md;
					itr2->second.mu_dz = pair_tmp.dz;
				}
				else if (itr2->second.pid == 13) {
					itr1->second.mu_md = md;
					itr1->second.mu_dz = pair_tmp.dz;
				}
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
		ofs << std::right << std::fixed
			<< std::setw(12) << std::setprecision(0) << itr0->eventid << " "
			<< std::setw(12) << std::setprecision(0) << itr0->unixtime << " "
			<< std::setw(4) << std::setprecision(0) << itr0->pl << " "
			<< std::setw(4) << std::setprecision(0) << itr0->trk.size() << " "
			<< std::setw(10) << std::setprecision(1) << itr0->x << " "
			<< std::setw(10) << std::setprecision(1) << itr0->y << " "
			<< std::setw(10) << std::setprecision(1) << itr0->z << " "
			<< std::setw(8) << std::setprecision(1) << itr0->dz << std::endl;
		for (auto itr1 = itr0->pair.begin(); itr1 != itr0->pair.end(); itr1++) {
			ofs << std::right << std::fixed
				//<< std::setw(4) << std::setprecision(0) << itr1->t[0].pl1 << " "
				<< std::setw(4) << std::setprecision(0) << itr1->t[0].chainid << " "
				//<< std::setw(4) << std::setprecision(0) << itr1->t[1].pl1 << " "
				<< std::setw(4) << std::setprecision(0) << itr1->t[1].chainid << " "
				<< std::setw(10) << std::setprecision(1) << itr1->x << " "
				<< std::setw(10) << std::setprecision(1) << itr1->y << " "
				<< std::setw(8) << std::setprecision(1) << itr1->dz << " "
				<< std::setw(8) << std::setprecision(4) << itr1->oa << " "
				<< std::setw(6) << std::setprecision(1) << itr1->md << std::endl;
		}
		ofs << std::right << std::fixed
			<< std::setw(4) << "cid" << " "
			<< std::setw(4) << "pid" << " "
			<< std::setw(4) << "Stop" << " "
			<< std::setw(4) << "VE1" << " "
			<< std::setw(4) << "VE2" << " "
			<< std::setw(4) << "VE3" << " "

			<< std::setw(4) << "nseg" << " "
			<< std::setw(4) << "npl" << " "
			<< std::setw(4) << "pl0" << " "
			<< std::setw(4) << "pl1" << " "
			<< std::setw(8) << "tan" << " "
			//<< std::setw(8) << "pb" << " "
			<< std::setw(8) << "mcs" << " "
			<< std::setw(8) << "pb" << " "//<< std::setw(8) << "rng" << " "
			<< std::setw(4) << "vph" << " "
			<< std::setw(4) << "dpl" << " "
			<< std::setw(4) << "vph2" << " "

			<< std::setw(10) << "dal" << " "
			<< std::setw(10) << "dar" << " "
			<< std::setw(10) << "dl" << " "
			//<< std::setw(8) << "dr" << " "
			<< std::setw(6) << "mu_md" << " "
			<< std::setw(8) << "mu_dz" << " "
			//<< std::setw(6) << "ip" << " "
			<< std::setw(10) << "rawid" << " "
			<< std::setw(7) << "ax" << " "
			<< std::setw(7) << "ay" << " "
			<< std::setw(8) << "x" << " "
			<< std::setw(8) << "y" << " "
			<< std::setw(8) << "z" << " "
			<< std::endl;//itr->first=ip

		for (auto itr1 = itr0->trk.begin(); itr1 != itr0->trk.end(); itr1++) {
			veflg = 0;
			/*auto vetrk=ve.find(k);
			if (vetrk != ve.end()) {
				veflg = 1;
			}
			if (itr1->second.pid == 13) {
				veflg = 0;
			}*/

			ofs << std::right << std::fixed
				<< std::setw(4) << std::setprecision(0) << itr1->second.chainid << " "
				<< std::setw(4) << std::setprecision(0) << itr1->second.pid << " "
				<< std::setw(4) << std::setprecision(0) << itr1->second.stop_flg << " ";
				ofs << std::setw(4) << std::setprecision(0) << -1 << " "
					<< std::setw(4) << std::setprecision(0) << -1 << " "
					<< std::setw(4) << std::setprecision(0) << -1 << " ";

			ofs
				<< std::setw(4) << std::setprecision(0) << itr1->second.nseg << " "
				<< std::setw(4) << std::setprecision(0) << itr1->second.npl << " "
				<< std::setw(4) << std::setprecision(0) << itr1->second.pl0 << " "
				<< std::setw(4) << std::setprecision(0) << itr1->second.pl1 << " "
				<< std::setw(8) << std::setprecision(4) << sqrt(itr1->second.ax * itr1->second.ax + itr1->second.ay * itr1->second.ay) << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.mom << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.pb << " "//tmp
				<< std::setw(4) << std::setprecision(0) << itr1->second.vph << " "
				<< std::setw(4) << std::setprecision(0) << itr1->second.d_pl << " "
				<< std::setw(4) << std::setprecision(0) << itr1->second.vph2 << " "
				<< std::setw(10) << std::setprecision(4) << itr1->second.dal << " "
				<< std::setw(10) << std::setprecision(4) << itr1->second.dar << " "
				<< std::setw(10) << std::setprecision(1) << itr1->second.dl << " "
				//<< std::setw(8) << std::setprecision(1) << itr1->second.dr << " "
				<< std::setw(6) << std::setprecision(1) << itr1->second.mu_md << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.mu_dz << " "
				//<< std::setw(6) << std::setprecision(1) << itr1->second.ip << " "
				<< std::setw(10) << std::setprecision(0) << itr1->second.rawid << " "
				<< std::setw(7) << std::setprecision(4) << itr1->second.ax << " "
				<< std::setw(7) << std::setprecision(4) << itr1->second.ay << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.x << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.y << " "
				<< std::setw(8) << std::setprecision(1) << itr1->second.z << " "
				<< std::endl;//itr->first=ip
		}
		ofs << std::endl;
	}

}


void Check_true_and_recon(std::vector<Momentum_recon::Event_information>& momch, std::set<tkey>& key, std::string output)
{
	std::ofstream ofs(output);
	if (!ofs) {
		std::cerr << "\tFile Open error...\n" << output << std::endl;
		return;
	}

	tkey k;
	for (auto& ev : momch) {

		k.pl = ev.vertex_pl;
		k.vx = ev.true_vertex_position[0];
		k.vy = ev.true_vertex_position[1];
		if (key.count(k) > 0) {
			// find

			ofs << std::right << std::fixed
				<< std::setw(10) << std::setprecision(1) << k.pl << " "
				<< std::setw(10) << std::setprecision(1) << k.vx << " "
				<< std::setw(10) << std::setprecision(1) << k.vy << " "
				<< std::setw(8) << std::setprecision(1) << ev.weight << std::endl;

			ofs
				<< std::right << std::fixed
				<< std::setw(12) << "groupid" << " "
				<< std::setw(12) << "unix_time" << " "
				<< std::setw(4) << "vpl" << " "
				<< std::setw(10) << "v_material" << " "
				<< std::setw(10) << "vx " << " "
				<< std::setw(10) << "vy" << " "
				<< std::setw(10) << "vz" << " "
				<< std::setw(5) << "recon" << " "
				<< std::setw(5) << "true" << " "
				<< std::setw(8) << "weight" << std::endl;

			ofs << std::right << std::fixed
				<< std::setw(12) << std::setprecision(0) << ev.groupid << " "
				<< std::setw(12) << std::setprecision(0) << ev.unix_time << " "
				<< std::setw(4) << std::setprecision(0) << ev.vertex_pl << " "
				<< std::setw(10) << std::setprecision(0) << ev.vertex_material << " "
				<< std::setw(10) << std::setprecision(1) << ev.vertex_position[0] << " "
				<< std::setw(10) << std::setprecision(1) << ev.vertex_position[1] << " "
				<< std::setw(10) << std::setprecision(1) << ev.vertex_position[2] << " "
				<< std::setw(5) << std::setprecision(0) << ev.chains.size() << " "
				<< std::setw(5) << std::setprecision(0) << ev.true_chains.size() << " "
				<< std::setw(8) << std::setprecision(1) << ev.weight << std::endl;

			ofs << std::right << std::fixed
				<< std::setw(12) << "chainid" << " "
				<< std::setw(12) << "particle_flg" << " "
				<< std::setw(12) << "mu rng" << " "
				<< std::setw(12) << "p rng" << " "
				<< std::setw(12) << "mu mcs" << " "
				<< std::setw(12) << "p mcs" << " "
				<< std::setw(12) << "bm rng" << " "
				<< std::setw(12) << "bm pi" << " "
				<< std::endl;
			std::cout << "true" << std::endl;
			for (auto& c : ev.true_chains) {
				
				ofs << std::right << std::fixed
					<< std::setw(12) << std::setprecision(0) << c.chainid << " "
					<< std::setw(12) << std::setprecision(0) << c.particle_flg << " "
					<< std::setw(12) << std::setprecision(1) << c.ecc_range_mom[0] << " "
					<< std::setw(12) << std::setprecision(1) << c.ecc_range_mom[1] << " "
					<< std::setw(12) << std::setprecision(1) << c.ecc_mcs_mom[0] << " "
					<< std::setw(12) << std::setprecision(1) << c.ecc_mcs_mom[1] << " "
					<< std::setw(12) << std::setprecision(1) << c.bm_range_mom << " "
					<< std::setw(12) << std::setprecision(1) << c.bm_curvature_mom << " "
					<< std::endl;

				for (auto& b : c.base) {
					ofs << std::right << std::fixed
						<< std::setw(12) << std::setprecision(0) << b.pl << " "
						<< std::setw(12) << std::setprecision(0) << b.rawid << " "
						<< std::setw(12) << std::setprecision(4) << b.ax << " "
						<< std::setw(12) << std::setprecision(4) << b.ay << " "
						<< std::setw(12) << std::setprecision(1) << b.x << " "
						<< std::setw(12) << std::setprecision(1) << b.y << " "
						<< std::setw(12) << std::setprecision(1) << b.z << " "
						<< std::endl;

				}
				for (auto& b : c.base_pair) {
					ofs << std::right << std::fixed
						<< std::setw(12) << std::setprecision(0) << b.first.pl << " "
						<< std::setw(12) << std::setprecision(0) << b.first.rawid << " "
						<< std::setw(12) << std::setprecision(4) << b.first.ax << " "
						<< std::setw(12) << std::setprecision(4) << b.first.ay << " "
						<< std::setw(12) << std::setprecision(1) << b.first.x << " "
						<< std::setw(12) << std::setprecision(1) << b.first.y << " "
						<< std::setw(12) << std::setprecision(1) << b.first.z << " "
						<< std::setw(12) << std::setprecision(0) << b.second.pl << " "
						<< std::setw(12) << std::setprecision(0) << b.second.rawid << " "
						<< std::setw(12) << std::setprecision(4) << b.second.ax << " "
						<< std::setw(12) << std::setprecision(4) << b.second.ay << " "
						<< std::setw(12) << std::setprecision(1) << b.second.x << " "
						<< std::setw(12) << std::setprecision(1) << b.second.y << " "
						<< std::setw(12) << std::setprecision(1) << b.second.z << " "
						<< std::endl;

				}

			}
			std::cout << "recon" << std::endl;

			for (auto& c : ev.chains) {

				ofs << std::right << std::fixed
					<< std::setw(12) << std::setprecision(0) << c.chainid << " "
					<< std::setw(12) << std::setprecision(0) << c.particle_flg << " "
					<< std::setw(12) << std::setprecision(1) << c.ecc_range_mom[0] << " "
					<< std::setw(12) << std::setprecision(1) << c.ecc_range_mom[1] << " "
					<< std::setw(12) << std::setprecision(1) << c.ecc_mcs_mom[0] << " "
					<< std::setw(12) << std::setprecision(1) << c.ecc_mcs_mom[1] << " "
					<< std::setw(12) << std::setprecision(1) << c.bm_range_mom << " "
					<< std::setw(12) << std::setprecision(1) << c.bm_curvature_mom << " "
					<< std::endl;

				for (auto& b : c.base) {
					ofs << std::right << std::fixed
						<< std::setw(12) << std::setprecision(0) << b.pl << " "
						<< std::setw(12) << std::setprecision(0) << b.rawid << " "
						<< std::setw(12) << std::setprecision(4) << b.ax << " "
						<< std::setw(12) << std::setprecision(4) << b.ay << " "
						<< std::setw(12) << std::setprecision(1) << b.x << " "
						<< std::setw(12) << std::setprecision(1) << b.y << " "
						<< std::setw(12) << std::setprecision(1) << b.z << " "
						<< std::endl;

				}
			}


		}
	}



	//printf("\t* input fin.\n\t #of track : %d --> %d\n", cnt, use);

}

