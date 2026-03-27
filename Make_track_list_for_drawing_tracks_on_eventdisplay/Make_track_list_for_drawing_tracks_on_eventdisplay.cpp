// Make_track_list_for_drawing_tracks_on_eventdisplay

// 2025/12/11
// kasumi

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}

void Output_chain_information(std::vector<Momentum_recon::Event_information>& momch, std::string output);
void Output_chain_information(std::vector<Momentum_recon::Event_information>& momch, std::string output, std::string list);

int main(int argc, char** argv) {
	if (argc < 3) {
		fprintf(stderr, "usage: in.momch out.txt [event list]");
		exit(1);
	}
	std::string file_in_momch_ori = argv[1];//input momch
	std::string file_out_momch = argv[2];

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch_ori);
	if (argc == 3) {
		Output_chain_information(momch, file_out_momch);
		//Momentum_recon::Write_Event_information_extension(file_out_momch, momch);
	}
	else if(argc==4){
		std::string event_list = argv[3];
		Output_chain_information(momch, file_out_momch, event_list);

	}

}
void Output_chain_information(std::vector<Momentum_recon::Event_information>& momch, std::string output) {
	Momentum_recon::Mom_chain chain;
	int cnt;
	int vpl, upl, dpl;
	std::ofstream ofs(output);

	int ut, eid, cid;
	double x0, y0, z0, x1, y1, z1;
	int pid, dir;

	for (auto& ev : momch) {
		ut = ev.unix_time;
		eid = ev.groupid;
		for (auto& c : ev.chains) {
			cid = c.chainid;;
			// most upstream pl
			x0 = c.base.rbegin()->x;
			y0 = c.base.rbegin()->y;
			z0 = c.base.rbegin()->z;
			// most downstream pl
			x1 = c.base.begin()->x;
			y1 = c.base.begin()->y;
			z1 = c.base.begin()->z;

			pid = c.particle_flg;
			dir = c.direction;
			ofs << std::right << std::fixed
				<< std::setw(12) << std::setprecision(0) << ut << " "
				<< std::setw(5) << std::setprecision(0) << eid << " "
				<< std::setw(3) << std::setprecision(0) << cid << " "
				<< std::setw(10) << std::setprecision(1) << x0 << " "
				<< std::setw(10) << std::setprecision(1) << y0 << " "
				<< std::setw(10) << std::setprecision(1) << z0 << " "
				<< std::setw(10) << std::setprecision(1) << x1 << " "
				<< std::setw(10) << std::setprecision(1) << y1 << " "
				<< std::setw(10) << std::setprecision(1) << z1 << " "
				<< std::setw(4) << std::setprecision(0) << pid << " "
				<< std::setw(2) << std::setprecision(0) << dir << " "
				<< std::setw(3) << std::setprecision(0) << ev.vertex_pl << " "
				<< std::endl;

		}
	}

}

void Output_chain_information(std::vector<Momentum_recon::Event_information>& momch, std::string output, std::string list) {
	std::ifstream ifs(list);
	int l;
	std::set<int> set;
	while (ifs >> l) {
		set.insert(l);
	}
	std::cout << " # of Event " << set.size() << std::endl;
	Momentum_recon::Mom_chain chain;
	int cnt = 0; int cnt2 = 0;;
	int vpl, upl, dpl;
	std::ofstream ofs(output);

	int ut, eid, cid;
	double x0, y0, z0, x1, y1, z1;
	int pid, dir;

	for (auto& ev : momch) {
		ut = ev.unix_time;
		eid = ev.groupid;
		for (auto& c : ev.chains) {
			cid = c.chainid;;
			// most upstream pl
			x1 = c.base.rbegin()->x;
			y1 = c.base.rbegin()->y;
			z1 = c.base.rbegin()->z;
			// most downstream pl
			x0 = c.base.begin()->x;
			y0 = c.base.begin()->y;
			z0 = c.base.begin()->z;

			pid = c.particle_flg;
			dir = c.direction;

			if (set.find(eid) != set.end()) {
				//std::cout << " found." << std::endl;
				ofs << std::right << std::fixed
					<< std::setw(12) << std::setprecision(0) << ut << " "
					<< std::setw(5) << std::setprecision(0) << eid << " "
					<< std::setw(3) << std::setprecision(0) << cid << " "
					<< std::setw(10) << std::setprecision(1) << x0 << " "
					<< std::setw(10) << std::setprecision(1) << y0 << " "
					<< std::setw(10) << std::setprecision(1) << z0 << " "
					<< std::setw(10) << std::setprecision(1) << x1 << " "
					<< std::setw(10) << std::setprecision(1) << y1 << " "
					<< std::setw(10) << std::setprecision(1) << z1 << " "
					<< std::setw(4) << std::setprecision(0) << pid << " "
					<< std::setw(2) << std::setprecision(0) << dir << " "
					//<< std::setw(3) << std::setprecision(0) << ev.vertex_pl << " "
					<< std::endl;

				cnt2++;
			}
		}
			cnt++;

	}
	std::cout << "( # of Event, #chains) = (" << cnt <<", "<<cnt2 <<")"<< std::endl;

}