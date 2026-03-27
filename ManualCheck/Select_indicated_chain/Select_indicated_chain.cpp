#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}


void select_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, int thr_nseg, int thr_npl, int thr_vph);

int main(int argc, char** argv) {
	if (argc < 6) {
		fprintf(stderr, "usage:input.momch output.momch thr_nseg thr_npl thr_vph \n");
		exit(1);
	}
	std::string file_in_momch_ori = argv[1];//input momch
	std::string file_out_momch = argv[2];//output
	int thr_nseg = std::stoi(argv[3]);
	int thr_npl = std::stoi(argv[4]);
	int thr_vph = std::stoi(argv[5]);

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch_ori);
	std::vector<Momentum_recon::Event_information> new_momch;


	select_momch(momch, new_momch, thr_nseg, thr_npl, thr_vph);

	Momentum_recon::Write_Event_information_extension(file_out_momch, new_momch);
}
void select_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, 
	int thr_nseg, int thr_npl, int thr_vph) {
	Momentum_recon::Mom_chain chain;

	int flg, cnt, vph, nseg, npl;
	for (auto& ev : momch0) {
		//std::cout << ev.chains.size() << std::endl;
		for (auto& c : ev.chains) {
			npl = c.base.rbegin()->pl - c.base.begin()->pl + 1;
			nseg = c.base.size();
			// debag
			//std::cout << "\tnseg / npl = " << c.base.size() << " / " << c.base.rbegin()->pl - c.base.begin()->pl + 1 << std::endl;
			if (npl == thr_npl)continue;
			if (nseg >= thr_nseg) continue;
			flg = 0; cnt = 0;
			vph = 0;
			for (auto& b : c.base) {
				vph = b.m[0].ph % 10000 + b.m[1].ph % 10000 + vph;
				cnt++;
			}

			vph = vph / cnt;
			if (vph > thr_vph) {
				flg++;
			}
		}
		if (flg > 0) {
			
			momch.push_back(ev);
			//std::cout << "\tnseg / npl = " << nseg << " / " << npl << std::endl;
		}
	}

}



