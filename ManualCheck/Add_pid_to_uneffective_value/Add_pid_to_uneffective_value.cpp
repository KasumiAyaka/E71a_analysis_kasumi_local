#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}


void select_momch(std::vector<Momentum_recon::Event_information>& momch0);

int main(int argc, char** argv) {
	if (argc < 3 || argc>3) {
		fprintf(stderr, "usage:input.momch output.momch\n");
		exit(1);
	}
	std::string file_in_momch_ori = argv[1];//input momch
	std::string file_out_momch = argv[2];//output

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch_ori);

	select_momch(momch);

	Momentum_recon::Write_Event_information_extension(file_out_momch, momch);
}
void select_momch(std::vector<Momentum_recon::Event_information>& momch0){

	int flg = 0;
	int vph = 0;
	for (auto& ev : momch0) {

		for (auto& c : ev.chains) {
			if (c.particle_flg < 0) {
				c.particle_flg = 0;
			}
		}
	}
}



