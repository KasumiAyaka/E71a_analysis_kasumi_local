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
	//Momentum_recon::Write_Event_information_extension(file_out_momch, momch);
}
void chainge_pid(std::vector<Momentum_recon::Event_information>& momch) {
	Momentum_recon::Mom_chain chain;
	double vph;
	double vph_thr = 60.0;
	int cnt;
	for (auto& ev : momch) {
		for (auto& c : ev.chains) {
			if (c.particle_flg < 1) {
				// â^ìÆó ÇÃÇ¬Ç¢ÇƒÇ¢Ç»Ç¢shortÇ»proton
				vph = 0;
				cnt = 0;
				//for (auto& b : c.base) {
				//	vph += b.m[0].ph % 10000;
				//	vph += b.m[1].ph % 10000;
				//	cnt++;
				//}
				//vph = vph / cnt;
				if (abs(c.base.rbegin()->pl - ev.vertex_pl) > abs(c.base.begin()->pl - ev.vertex_pl)) {
					vph += c.base.begin()->m[0].ph % 10000;
					vph += c.base.begin()->m[1].ph % 10000;
				}
				else {
					vph += c.base.rbegin()->m[0].ph % 10000;
					vph += c.base.rbegin()->m[1].ph % 10000;
				}

				std::cout << std::right << std::fixed
					<< std::setw(5) << std::setprecision(0) << ev.groupid << " "
					<< std::setw(3) << std::setprecision(0) << c.chainid << " "
					<< std::setw(4) << std::setprecision(0) << c.particle_flg << " "
					<< std::setw(5) << std::setprecision(0) << vph << " "
					<< std::setw(5) << std::setprecision(1) << sqrt(c.base.begin()->ax * c.base.begin()->ax + c.base.begin()->ay * c.base.begin()->ay) << " "
					<< std::setw(3) << std::setprecision(0) << c.base.size() << " /"
					<< std::setw(3) << std::setprecision(0) << c.base.rbegin()->pl-c.base.begin()->pl+1 << " "
					<< std::setw(8) << std::setprecision(1) << c.ecc_mcs_mom[1]
					<< std::endl;
				if (vph > vph_thr) {
					//c.particle_flg = 2212;

				}
				
			}
		}
	}

}
