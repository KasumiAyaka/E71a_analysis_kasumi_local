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
	int groupid, chainid, count_vph, sensorid, nseg, npl;
	float average_vph, sd_vph, ecc_mcs_mom, ax, ay, angle;
};

std::vector<output_format> Calc_average_momch(std::vector<Momentum_recon::Event_information>& momch);
void output(std::string file_out, std::vector<output_format>& ret);
bool Calc_average(Momentum_recon::Mom_chain& c, int& count_vph, float& average_vph, float& sd_vph);
void output_bin(std::string filename, std::vector<output_format>& ret);
void output_bin(std::string filename, std::vector<output_format>& ret, int i_ang, int i_mom);


int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "usage:filename\n");
		fprintf(stderr, "Calc_average_basetrack_vph.exe vph_under1500.momch average.txt average.bin");
		exit(1);
	}
	std::string file_in_momch = argv[1];
	std::string file_out = argv[2];
	std::string file_out_bin = argv[3];
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
			out.ecc_mcs_mom = c.Get_muon_mcs_pb();//pbeta
			out.nseg = c.base.size();
			out.npl = c.base.rbegin()->pl - c.base.begin()->pl + 1;
			out.ax = 0;
			out.ay = 0;
			out.sensorid = 0;
			
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

void output(std::string file_out, std::vector<output_format>& ret) {
	std::ofstream ofs(file_out);
	int count = 0, all = ret.size();
	for (auto itr = ret.begin(); itr != ret.end(); itr++) {
		if (count % 10000 == 0) {
			fprintf(stderr, "\r write file %10d/%10d(%4.1lf%%)", count, all, count * 100. / all);
		}
		count++;
		ofs << std::right << std::fixed
			<< std::setw(10) << std::setprecision(0) << itr->groupid << " "
			<< std::setw(10) << std::setprecision(0) << itr->chainid << " "
			<< std::setw(3) << std::setprecision(0) << itr->sensorid << " "
			<< std::setw(3) << std::setprecision(0) << itr->nseg << " "
			<< std::setw(3) << std::setprecision(0) << itr->npl << " "
			<< std::setw(7) << std::setprecision(4) << itr->ax << " "
			<< std::setw(7) << std::setprecision(4) << itr->ay << " "
			<< std::setw(7) << std::setprecision(4) << itr->angle << " "
			<< std::setw(10) << std::setprecision(3) << itr->ecc_mcs_mom << " "
			<< std::setw(5) << std::setprecision(0) << itr->count_vph << " "// == nseg‚Å‚ÍH
			<< std::setw(8) << std::setprecision(3) << itr->average_vph << " "
			<< std::setw(8) << std::setprecision(4) << itr->sd_vph << std::endl;
		
	}
	fprintf(stderr, "\r write file %10d/%10d(%4.1lf%%)\n", count, all, count * 100. / all);

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
void output_bin(std::string filename, std::vector<output_format>& ret, int i_ang, int i_mom) {
	double angle_min, angle_max;
	double mom_min, mom_max;
	if (i_ang < 10) {
		angle_min = i_ang * 0.1;
		angle_max = (i_ang + 1) * 0.1;
	}
	else if (i_ang < 15) {
		angle_min = (i_ang - 10) * 0.2 + 1.0;
		angle_max = (i_ang - 10 + 1) * 0.2 + 1.0;
	}
	else {
		angle_min = (i_ang - 15) * 0.3 + 2.0;
		angle_max = (i_ang - 15 + 1) * 0.3 + 2.0;
	}
	mom_min = i_mom * 100;
	mom_max = (i_mom + 1) * 100;
	std::stringstream file_out;
	file_out << filename
		<< "_" << std::setw(2) << std::setfill('0') << i_ang
		<< "_" << std::setw(2) << std::setfill('0') << i_mom;
	fprintf(stderr, "write [%s]\n", file_out.str().c_str());
	std::ofstream ofs(file_out.str(), std::ios::binary);
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
		if (ret[i].angle < angle_min)continue;
		if (angle_max <= ret[i].angle)continue;
		if (ret[i].ecc_mcs_mom < mom_min)continue;
		if (mom_max <= ret[i].ecc_mcs_mom)continue;
		ofs.write((char*)&ret[i], sizeof(output_format));
	}
	std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%" << std::endl;
	ofs.close();

}