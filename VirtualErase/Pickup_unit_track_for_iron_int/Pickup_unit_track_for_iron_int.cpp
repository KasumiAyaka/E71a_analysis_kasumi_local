// 2025/09/10
// Pickup_unit_track.cpp
// kasumi

// 2025/12/02
// 大幅に変更。元のsrcは下記に。
// T:\NINJA\E71a\work\kasumi\VirtualErase\prg\src


#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>


struct Key {
	int pl, rid;
};
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.pl, lhs.rid) < std::tie(rhs.pl, rhs.rid);
}

struct UnitList {
	int ve;
	Key t[2];
};
bool operator<(const UnitList& lhs, const UnitList& rhs) {
	return std::tie(lhs.ve, lhs.t[0], lhs.t[1]) < std::tie(rhs.ve, rhs.t[0], rhs.t[1]);
}
struct Event {
	int gid, cid, pid;
	UnitList u1, u2, u3;
};
bool operator<(const Event& lhs, const Event& rhs) {
	return std::tie(lhs.gid, lhs.cid, lhs.u1, lhs.u2, lhs.u3) < std::tie(rhs.gid, rhs.cid, rhs.u1, rhs.u2, rhs.u3);
}


void set_list(std::string input, std::multimap<int, int>& map);
void select_ve_tracks(std::vector<Momentum_recon::Event_information>& momch0, int npl_thr, std::string output);

int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "usage:input.momch output.txt npl_thr(<=npl_thr)\n");
		exit(1);
	}
	std::string file_in_momch = argv[1];//input momch
	std::string file_out_txt = argv[2];
	int npl_thr = std::stoi(argv[3]);
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch);

	select_ve_tracks(momch, npl_thr, file_out_txt);
	//Momentum_recon::Write_Event_information_extension(file_out_momch, new_momch);

}
void set_list(std::string input, std::multimap<int, int>& map) {
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error!" << std::endl;
		exit(0);
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	while (std::getline(ifs, str)) {
		// gid cid pid nseg npl upl dpl
		str_v = StringSplit_with_tab(str);
		//std::cout << str << std::endl;
		map.insert(std::make_pair(std::stoi(str_v[0]), std::stoi(str_v[1])));
	}
	std::cout << "   *  Total Track number : " << map.size() << std::endl;
}
void select_ve_tracks(std::vector<Momentum_recon::Event_information>& momch0, int npl_thr, std::string output) {
	Momentum_recon::Mom_chain chain;

	//std::multimap<int, int> map;//gid,cid
	std::ofstream ofs(output);

	std::set<Event> set;
	int cnt = 0;
	int npl, nseg, vpl;
	Event e;
	for (auto& ev : momch0) {
		vpl = ev.vertex_pl;
		for (auto& c : ev.chains) {
			e = { -1 };
			e.gid = ev.groupid;
			e.cid = c.chainid;
			e.pid = c.particle_flg;
			nseg = c.base.size();//nseg
			npl = c.base.rbegin()->pl - c.base.begin()->pl + 1;//npl
			if (c.particle_flg == 13)continue;
			if (npl > npl_thr)continue;
			e.u1.ve = 1;
			e.u2.ve = 2;
			e.u3.ve = 3;

			e.u1.t[0].pl = -1;
			e.u1.t[0].rid = -1;
			e.u1.t[1].pl = -1;
			e.u1.t[1].rid = -1;
			e.u2.t[0].pl = -1;
			e.u2.t[0].rid = -1;
			e.u2.t[1].pl = -1;
			e.u2.t[1].rid = -1;
			e.u3.t[0].pl = -1;
			e.u3.t[0].rid = -1;
			e.u3.t[1].pl = -1;
			e.u3.t[1].rid = -1;
			if (c.base.rbegin()->pl == vpl) {
				// fwd
				for (auto itr = c.base.rbegin(); itr != c.base.rend(); itr++) {
					
					if (itr->pl == vpl) {
						e.u1.t[1].pl = itr->pl;
						e.u1.t[1].rid = itr->rawid;
					}
					if (itr->pl == vpl - 1) {
						// vtx 直下の次のunit
						e.u2.t[0].pl = itr->pl;
						e.u2.t[0].rid = itr->rawid;
					}
					if (itr->pl == vpl - 2) {
						e.u2.t[1].pl = itr->pl;
						e.u2.t[1].rid = itr->rawid;
					}
					if (itr->pl == vpl - 3) {
						// vtx直下の2つ次のunit
						e.u3.t[0].pl = itr->pl;
						e.u3.t[0].rid = itr->rawid;
					}
					if (itr->pl == vpl - 4) {
						e.u3.t[1].pl = itr->pl;
						e.u3.t[1].rid = itr->rawid;
					}
				}
			}
			else {
				// bwd
				for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
					if (itr->pl == vpl + 1) {
						e.u1.t[0].pl = itr->pl;
						e.u1.t[0].rid = itr->rawid;
					}
					if (itr->pl == vpl + 3) {
						e.u2.t[0].pl = itr->pl;
						e.u2.t[0].rid = itr->rawid;
					}
					if (itr->pl == vpl + 2) {
						e.u2.t[1].pl = itr->pl;
						e.u2.t[1].rid = itr->rawid;
					}
					if (itr->pl == vpl + 5) {
						e.u3.t[0].pl = itr->pl;
						e.u3.t[0].rid = itr->rawid;
					}
					if (itr->pl == vpl + 4) {
						e.u3.t[1].pl = itr->pl;
						e.u3.t[1].rid = itr->rawid;
					}
				}

			}
			set.insert(e);
			cnt++;
		}
	}
	std::cout << " npl <= " << npl_thr << " : " << cnt << " trks" << std::endl;

	for (auto itr = set.begin(); itr != set.end(); itr++) {
		ofs << std::right << std::fixed //<< std::setfill(' ')
			<< std::setw(5) << std::setprecision(0) << itr->gid << " "
			<< std::setw(4) << std::setprecision(0) << itr->cid << " "
			<< std::setw(4) << std::setprecision(0) << itr->pid << " "
			<< std::setw(4) << std::setprecision(0) << itr->u1.ve << " "
			<< std::setw(4) << std::setprecision(0) << itr->u1.t[0].pl << " "
			<< std::setw(10) << std::setprecision(0) << itr->u1.t[0].rid << " "
			<< std::setw(4) << std::setprecision(0) << itr->u1.t[1].pl << " "
			<< std::setw(10) << std::setprecision(0) << itr->u1.t[1].rid
			<< std::endl;
		ofs << std::right << std::fixed //<< std::setfill(' ')
			<< std::setw(5) << std::setprecision(0) << itr->gid << " "
			<< std::setw(4) << std::setprecision(0) << itr->cid << " "
			<< std::setw(4) << std::setprecision(0) << itr->pid << " "
			<< std::setw(4) << std::setprecision(0) << itr->u2.ve << " "
			<< std::setw(4) << std::setprecision(0) << itr->u2.t[0].pl << " "
			<< std::setw(10) << std::setprecision(0) << itr->u2.t[0].rid << " "
			<< std::setw(4) << std::setprecision(0) << itr->u2.t[1].pl << " "
			<< std::setw(10) << std::setprecision(0) << itr->u2.t[1].rid
			<< std::endl;
		ofs << std::right << std::fixed //<< std::setfill(' ')
			<< std::setw(5) << std::setprecision(0) << itr->gid << " "
			<< std::setw(4) << std::setprecision(0) << itr->cid << " "
			<< std::setw(4) << std::setprecision(0) << itr->pid << " "
			<< std::setw(4) << std::setprecision(0) << itr->u3.ve << " "
			<< std::setw(4) << std::setprecision(0) << itr->u3.t[0].pl << " "
			<< std::setw(10) << std::setprecision(0) << itr->u3.t[0].rid << " "
			<< std::setw(4) << std::setprecision(0) << itr->u3.t[1].pl << " "
			<< std::setw(10) << std::setprecision(0) << itr->u3.t[1].rid
			<< std::endl;
	}
}



