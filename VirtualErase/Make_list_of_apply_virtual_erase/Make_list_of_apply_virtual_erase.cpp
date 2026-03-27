// 2025/09/10
// Pickup_unit_track.cpp
// kasumi

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}

struct VE_LIST {
	int gid, cid, pl0, pl1, pl2, rid0, rid1, rid2,vpl;
};
bool operator<(const VE_LIST& lhs, const VE_LIST& rhs) {
	return std::tie(lhs.gid, lhs.cid, lhs.vpl, lhs.pl2, lhs.pl1, lhs.pl0, lhs.rid2, lhs.rid1, lhs.rid0)
		< std::tie(rhs.gid, rhs.cid, rhs.vpl, rhs.pl2, rhs.pl1, rhs.pl0, rhs.rid2, rhs.rid1, rhs.rid0);
}

void select_momch(std::vector<Momentum_recon::Event_information>& momch0, std::string output, std::string output1, int thr_npl);

int main(int argc, char** argv) {
	if (argc < 5 || argc>5) {
		fprintf(stderr, "usage:input.momch output.txt rootlist.txt thr_npl\n");
		exit(1);
	}
	std::string file_in_momch_ori = argv[1];//input momch
	std::string file_out_txt = argv[2];//unit_btk
	std::string file_list_txt = argv[3];//unit_btk
	int thr_npl = std::stoi(argv[4]);

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch_ori);
	std::vector<Momentum_recon::Event_information> new_momch;

	select_momch(momch, file_out_txt, file_list_txt, thr_npl);
	//Momentum_recon::Write_Event_information_extension(file_out_momch, new_momch);

}
void select_momch(std::vector<Momentum_recon::Event_information>& momch0, std::string output, std::string output1, int thr_npl) {
	Momentum_recon::Mom_chain chain;

	// Select 3 basetracks which are most close  plate of vpl.
	// If partner track is "fwd" track, pickup vpl, vpl-1, vpl-2.
	// If partner track is "bwd" track, pickup vpl+3, vpl+2, vpl1.

	//std::multimap<int, int> map;//gid,cid
	std::ofstream ofs(output);
	std::ofstream ofs1(output1);
	VE_LIST vl;

	int cnt = 0;
	int flg = 0;
	int vph = 0;
	int vpl = 0;
	int all_partner_track = 0;
	int all_event_num = 0;
	int rcnt;
	for (auto& ev : momch0) {
		vpl = ev.vertex_pl;
		all_event_num++;


		for (auto& c : ev.chains) {
			all_partner_track++;
			// partner-track, nseg<50
			if (c.particle_flg != 13 && c.base.rbegin()->pl - c.base.begin()->pl < thr_npl) {
				//std::cout << c.base.rbegin()->pl - c.base.begin()->pl << " " << c.base.begin()->pl << "," << c.base.rbegin()->pl << std::endl;
				cnt++;
				rcnt = 0;
				vl.gid = ev.groupid;
				vl.cid = c.chainid;
				vl.vpl = ev.vertex_pl;
				vl.pl0 = vl.pl1 = vl.pl2 = -1;
				vl.rid0 = vl.rid1 = vl.rid2 = -1;
				for (auto& b : c.base) {
					// fwd track
					if (vpl - b.pl > -1 && vpl - b.pl < 3) {
						ofs << std::right << std::fixed
							<< std::setw(5) << std::setprecision(0) << ev.groupid << " "
							<< std::setw(5) << std::setprecision(0) << c.chainid << " "
							<< std::setw(3) << std::setprecision(0) << b.pl << " "
							<< std::setw(12) << std::setprecision(0) << b.rawid << " "
							<< std::setw(8) << std::setprecision(4) << b.ax << " "
							<< std::setw(8) << std::setprecision(4) << b.ay << " "
							<< std::setw(10) << std::setprecision(1) << b.x << " "
							<< std::setw(10) << std::setprecision(1) << b.y << " "
							<< std::setw(8) << std::setprecision(0) << b.m[0].ph << " "
							<< std::setw(8) << std::setprecision(0) << b.m[1].ph << " "
							<< std::setw(3) << std::setprecision(0) << ev.vertex_pl << " "
							<< std::setw(3) << std::setprecision(0) << fabs(ev.vertex_pl - b.pl) << " "
							<< std::setw(6) << std::setprecision(1) << c.ecc_mcs_mom[0] << " "
							<< std::setw(6) << std::setprecision(4) << c.ecc_range_mom_error[0][0] << " "
							<< std::setw(6) << std::setprecision(4) << c.ecc_range_mom_error[0][1] << " "
							//<< std::setw(3) << std::setprecision(0) << c.ecc_mcs_mom << " "
							<< std::endl;
						if (vpl - b.pl == 0) {
							vl.rid0 = b.rawid;
							vl.pl0 = b.pl;
						}
						if (vpl - b.pl == 1) {
							vl.rid1 = b.rawid;
							vl.pl1 = b.pl;
						}
						if (vpl - b.pl == 2) {
							vl.rid2 = b.rawid;
							vl.pl2 = b.pl;
						}
					}

					// bwd trk
					if (b.pl - vpl > 0 && b.pl - vpl < 4) {
						ofs << std::right << std::fixed
							<< std::setw(5) << std::setprecision(0) << ev.groupid << " "
							<< std::setw(5) << std::setprecision(0) << c.chainid << " "
							<< std::setw(3) << std::setprecision(0) << b.pl << " "
							<< std::setw(12) << std::setprecision(0) << b.rawid << " "
							<< std::setw(8) << std::setprecision(4) << b.ax << " "
							<< std::setw(8) << std::setprecision(4) << b.ay << " "
							<< std::setw(10) << std::setprecision(1) << b.x << " "
							<< std::setw(10) << std::setprecision(1) << b.y << " "
							<< std::setw(8) << std::setprecision(0) << b.m[0].ph << " "
							<< std::setw(8) << std::setprecision(0) << b.m[1].ph << " "
							<< std::setw(3) << std::setprecision(0) << ev.vertex_pl << " "
							<< std::setw(3) << std::setprecision(0) << fabs(ev.vertex_pl - b.pl) << " "
							<< std::setw(6) << std::setprecision(1) << c.ecc_mcs_mom[0] << " "
							<< std::setw(6) << std::setprecision(4) << c.ecc_range_mom_error[0][0] << " "
							<< std::setw(6) << std::setprecision(4) << c.ecc_range_mom_error[0][1] << " "
							//<< std::setw(3) << std::setprecision(0) << c.ecc_mcs_mom << " "
							<< std::endl;

						if (b.pl - vpl == 1) {
							vl.pl2 = b.pl;
							vl.rid2 = b.rawid;
						}
						if (b.pl - vpl == 2) {
							vl.pl1 = b.pl;
							vl.rid1 = b.rawid;
						}
						if (b.pl - vpl == 3) {
							vl.pl0 = b.pl;
							vl.rid0 = b.rawid;
						}
					}
				}
				ofs1 << std::right << std::fixed
					<< std::setw(5) << std::setprecision(0) << ev.groupid << " "
					<< std::setw(5) << std::setprecision(0) << c.chainid << " "
					<< std::setw(3) << std::setprecision(0) << vl.vpl << " "
					<< std::setw(3) << std::setprecision(0) << vl.pl0 << " "
					<< std::setw(12) << std::setprecision(0) << vl.rid0 << " "
					<< std::setw(3) << std::setprecision(0) << vl.pl1 << " "
					<< std::setw(12) << std::setprecision(0) << vl.rid1 << " "
					<< std::setw(3) << std::setprecision(0) << vl.pl2 << " "
					<< std::setw(12) << std::setprecision(0) << vl.rid2 << " "
					<< std::endl;
			}
		}
	}
	std::cout << " # of events =  " << all_event_num << std::endl;
	std::cout << " # of all partner track = " << all_partner_track << std::endl;
	std::cout << " # of partner track (npl < " << thr_npl << ") " << cnt << std::endl;
}


