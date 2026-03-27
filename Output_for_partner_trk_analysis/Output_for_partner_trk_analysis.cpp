// 2025/11/05
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

int main(int argc, char** argv) {
	if (argc < 3) {
		fprintf(stderr, "usage: in.momch out.txt");
		exit(1);
	}
	std::string file_in_momch_ori = argv[1];//input momch
	std::string file_out_momch = argv[2];

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch_ori);
	Output_chain_information(momch, file_out_momch);
	//Momentum_recon::Write_Event_information_extension(file_out_momch, momch);

}
void Output_chain_information(std::vector<Momentum_recon::Event_information>& momch, std::string output) {
	Momentum_recon::Mom_chain chain;
	int cnt;
	int vpl, upl, dpl;
	std::ofstream ofs(output);

	for (auto& ev : momch) {
		for (auto itr = ev.chains.begin(); itr != ev.chains.end();itr++) {
			vpl = ev.vertex_pl;
			upl = itr->base.rbegin()->pl;
			dpl = itr->base.begin()->pl;

			ofs << ev.groupid << " "
				<< ev.vertex_pl << " "
				<< itr->chainid << " "
				<< itr->particle_flg << " "
				<< itr->base.size() << " "
				<< itr->base.rbegin()->pl - itr->base.begin()->pl + 1 << " ";
			std::cout << vpl << " " << upl << " " << dpl << std::endl;
			if (dpl - vpl > 0) {
				//bwd
				ofs << dpl << " "
					<< itr->base.begin()->m[0].ph + itr->base.begin()->m[1].ph << " "
					<< itr->base.begin()->ax << " "
					<< itr->base.begin()->ay << " "
					<< itr->base.begin()->x << " "
					<< itr->base.begin()->y << " "
					<< itr->base.begin()->z << " ";
				std::cout << "bwd" << std::endl;
			}else
			{
				//fwd
				ofs << upl << " "
					<< itr->base.rbegin()->m[0].ph + itr->base.rbegin()->m[1].ph << " "
					<< itr->base.rbegin()->ax << " "
					<< itr->base.rbegin()->ay << " "
					<< itr->base.rbegin()->x << " "
					<< itr->base.rbegin()->y << " "
					<< itr->base.rbegin()->z << " ";
				std::cout << "fwd" << std::endl;
			}
			ofs << itr->Get_muon_mcs_pb() << " "
				<< itr->ecc_mcs_mom[0] << " "
				<< std::endl;
		}
	}

}
