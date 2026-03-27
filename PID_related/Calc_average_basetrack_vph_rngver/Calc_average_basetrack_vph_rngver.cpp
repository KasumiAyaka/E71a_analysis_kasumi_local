// 2025/11/28
// kasumi
// event‚Ìpartner‚Ì•]‰¿‚Ìreference‚Ì‚½‚ß
// C:\Users\kasumi\source\repos\PID_related\x64\Release\Calc_average_basetrack_vph.exe
// vph_under1500.momch average.txt average.bin
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <set>
#include <map>
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
class output_format {
public:
	int groupid, chainid, count_vph, sensorid, nseg, npl,pid;
	float average_vph, sd_vph, mcs_p,mcs_mip, ax, ay, angle;
};

std::vector<output_format> Calc_average_momch(std::vector<Momentum_recon::Event_information>& momch);
bool Calc_average(Momentum_recon::Mom_chain& c, int& count_vph, float& average_vph, float& sd_vph);
void output_bin(std::string filename, std::vector<output_format>& ret);


int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "usage:filename\n");
		fprintf(stderr, "Calc_average_basetrack_vph.exe vph_under1500.momch average.bin");
		exit(1);
	}
	std::string file_in_momch = argv[1];
	std::string file_out_bin  = argv[2];

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch);
	std::vector<output_format>ave = Calc_average_momch(momch);
	output_bin(file_out_bin, ave);

	//for (int i_mom = 0; i_mom < 10; i_mom++) {
	//	for (int i_ang = 0; i_ang < 30; i_ang++) {
	//		output_bin(file_out_bin, ave, i_ang, i_mom);
	//	}
	//}
	//output(file_out, ave);

	exit(0);
}
std::vector<output_format> Calc_average_momch(std::vector<Momentum_recon::Event_information>& momch) {
	std::vector<output_format> ret;

	for (auto& ev : momch) {
		for (auto& c : ev.chains) {
			output_format out;
			out.groupid = ev.groupid;
			out.chainid = c.chainid;
			out.mcs_mip = c.ecc_mcs_mom[0];//mip
			out.mcs_p = c.ecc_mcs_mom[1];//p
			out.nseg = c.base.size();
			out.npl = c.base.rbegin()->pl - c.base.begin()->pl + 1;
			out.ax = 0;
			out.ay = 0;
			out.sensorid = 0;
			out.pid = c.particle_flg;

			int count = 0;
			for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
				out.ax += itr->ax;
				out.ay += itr->ay;
				count++;
			}
			out.ax /= count;
			out.ay /= count;
			out.angle = sqrt(out.ax * out.ax + out.ay * out.ay);
			if (Calc_average(c, out.count_vph, out.average_vph, out.sd_vph)) {
				ret.push_back(out);
			}
		}
	}
	return ret;

}
bool Calc_average(Momentum_recon::Mom_chain& c, int& count_vph, float& average_vph, float& sd_vph) {
	count_vph = 0;
	average_vph = 0;
	double average_vph2 = 0;
	double average_pixel2 = 0;
	int vph = 0;
	for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
		vph = itr->m[0].ph % 10000 + itr->m[1].ph % 10000;
		average_vph = average_vph + vph;
		count_vph += 1;
		average_vph2 += pow(vph, 2);
		//std::cout << count_vph<< " " << itr->m[0].ph % 10000 << " " << itr->m[1].ph % 10000 << " " << vph << std::endl;
	}
	//std::cout << average_vph << " " << count_vph << std::endl;
	if (count_vph <= 2)return false;

	average_vph /= count_vph;
	sd_vph = average_vph2 / count_vph - pow(average_vph, 2);

	if (sd_vph <= 0.1)return false;

	sd_vph = sqrt(sd_vph) * sqrt(count_vph) / sqrt(count_vph - 1);
	return true;

}

void output_bin(std::string filename, std::vector<output_format>& ret) {
	std::ofstream ofs(filename, std::ios::binary);
	if (!ofs) {
		//file open s
		fprintf(stderr, "File[%s] is not exist!!\n", filename.c_str());
		exit(1);
	}
	if (ret.size() == 0) {
		fprintf(stderr, "target linklet ... null\n");
		fprintf(stderr, "File[%s] has no text\n", filename.c_str());
	}
	int64_t count = 0;
	int64_t max = ret.size();
	for (int i = 0; i < ret.size(); i++) {
		if (count % 10000 == 0) {
			std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%";
		}
		count++;
		ofs.write((char*)&ret[i], sizeof(output_format));
	}
	std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%" << std::endl;
	ofs.close();

}
