// 2024/06/14
//kasumi
//I:\NINJA\E71a\work\kasumi\ECC\MuonAnalysis\x64\Release\erase_cahin_from_momch.exe  event_water_fin3.momch ..\temp_vtx_info\elaselist.txt event_water_fin3_ECC6.momch
//2024/07/12  list Çgid cidÇ…ïœçX
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}

void chainge_pid(std::vector<Momentum_recon::Event_information>& momch);

int main(int argc, char** argv) {
	if (argc < 2) {
		fprintf(stderr, "usage: in.momch  out.momch");
		exit(1);
	}
	std::string file_in_momch_ori = argv[1];//input momch
	std::string file_out_momch = file_in_momch_ori;
	if (argc == 3) {
		file_out_momch = argv[2];
	}

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch_ori);
	chainge_pid(momch);
	Momentum_recon::Write_Event_information_extension(file_out_momch, momch);
}
void chainge_pid(std::vector<Momentum_recon::Event_information>& momch){
	Momentum_recon::Mom_chain chain;
	for (auto& ev : momch) {
		for (auto& c : ev.chains) {
			if (c.particle_flg == -1) {
				c.particle_flg = 0;
			}
		}
	}

}
