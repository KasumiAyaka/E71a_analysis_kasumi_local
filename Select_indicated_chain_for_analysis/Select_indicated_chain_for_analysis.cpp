// 2025/11/11
// kasumi
// chain analysis: output information of linklet


#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}

struct BaseTrack {
	int rid, pl, ph;
	double ax, ay, x, y, z;
};

struct Linklet {
	BaseTrack b0, b1;
};

void Output_chain_information(std::vector<Momentum_recon::Event_information>& momch, std::string output, int nseg_thr, int npl_thr);

int main(int argc, char** argv) {
	if (argc < 3) {
		fprintf(stderr, "usage: in.momch out.txt nseg_thr npl_thr");
		exit(1);
	}
	std::string file_in_momch_ori = argv[1];//input momch
	std::string file_out_momch = argv[2];
	int nseg_thr = std::stoi(argv[3]);
	int npl_thr = std::stoi(argv[4]);

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch_ori);
	Output_chain_information(momch, file_out_momch, nseg_thr, npl_thr);
	//Momentum_recon::Write_Event_information_extension(file_out_momch, momch);

}
void Output_chain_information(std::vector<Momentum_recon::Event_information>& momch, std::string output,int nseg_thr,int npl_thr) {
	Momentum_recon::Mom_chain chain;
	int cnt,flg;
	int vpl, upl, dpl;
	std::ofstream ofs(output);

	Linklet l;
	for (auto& ev : momch) {

		for (auto& c : ev.chains) {
			if (c.base.size() > nseg_thr) {
				// nseg cut
				flg++;
			}

			if (c.base.rbegin()->pl - c.base.begin()->pl + 1 > npl_thr) {
				// npl cut
				flg++;
			}

			if (flg < 2)continue;
			flg = 0;
			cnt = 0;

			for (auto& b : c.base) {
				// request continuity
				if (b.pl > ev.vertex_pl) {
					// bwd
					if (b.pl == ev.vertex_pl + 1) {
						l.b0.ax = b.ax;
						l.b0.ay = b.ay;
						l.b0.x = b.x;
						l.b0.y = b.y;
						l.b0.z = b.z;
						l.b0.ph = b.m[0].ph + b.m[1].ph;
						l.b0.pl = b.pl;
						l.b0.rid = b.rawid;
						cnt++;
					}
					if (b.pl == ev.vertex_pl + 2) {
						cnt++;
					}
					if (b.pl == ev.vertex_pl + 3) {
						l.b1.ax = b.ax;
						l.b1.ay = b.ay;
						l.b1.x = b.x;
						l.b1.y = b.y;
						l.b1.z = b.z;
						l.b1.ph = b.m[0].ph + b.m[1].ph;
						l.b1.pl = b.pl;
						l.b1.rid = b.rawid;
						cnt++;
					}
				}
				else {
					//fwd
					if (b.pl == ev.vertex_pl) {
						l.b0.ax = b.ax;
						l.b0.ay = b.ay;
						l.b0.x = b.x;
						l.b0.y = b.y;
						l.b0.z = b.z;
						l.b0.ph = b.m[0].ph + b.m[1].ph;
						l.b0.pl = b.pl;
						l.b0.rid = b.rawid;
						cnt++;
					}
					if (b.pl == ev.vertex_pl - 1) {
						cnt++;
					}
					if (b.pl == ev.vertex_pl - 2) {
						l.b1.ax = b.ax;
						l.b1.ay = b.ay;
						l.b1.x = b.x;
						l.b1.y = b.y;
						l.b1.z = b.z;
						l.b1.ph = b.m[0].ph + b.m[1].ph;
						l.b1.pl = b.pl;
						l.b1.rid = b.rawid;
						cnt++;
					}
				}

				if (cnt == 3) {

					ofs << ev.groupid << " "
						<< ev.vertex_pl << " "
						<< c.chainid << " "
						<< c.particle_flg << " "
						<< c.base.size() << " "
						<< c.base.rbegin()->pl - c.base.begin()->pl + 1 << " "
						<< c.ecc_mcs_mom[0] << " "
						<< l.b0.pl << " "
						<< l.b0.rid << " "
						<< l.b0.ph << " "
						<< l.b0.ax << " "
						<< l.b0.ay << " "
						<< l.b0.x << " "
						<< l.b0.y << " "
						<< l.b0.z << " "

						<< l.b1.pl << " "
						<< l.b1.rid << " "
						<< l.b1.ph << " "
						<< l.b1.ax << " "
						<< l.b1.ay << " "
						<< l.b1.x << " "
						<< l.b1.y << " "
						<< l.b1.z << " "
						<< std::endl;
					cnt = 0;

				}

			}
		}
		
	}

}
