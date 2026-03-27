// 2025/12/04
// kasumi
// input 
// K:\NINJA\E71a\work\suzuki\PID\VPHcorrection\gel\vph\corr_plot\vph_pl_correction_mip.txt
//ofs << i_pl << " " << face << " " << i_ang << " "
//<< fit_func->GetParameter(1) * fit_func->GetParameter(2) << " "
//<< sqrt(pow(fit_func->GetParameter(1) * fit_func->GetParError(2), 2) + pow(fit_func->GetParameter(2) * fit_func->GetParError(1), 2)) << " "
//<< gaus->GetParameter(1) << " " << gaus->GetParError(1) << std::endl;

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <map>
#include <fstream>
#include <sstream>

struct Key {
	int pl, i_ang;
};
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.pl, lhs.i_ang) < std::tie(rhs.pl, rhs.i_ang);
}

struct VPH_MIP {
	Key k;
	int  face;
	double vphcenter, sigma;
	double vphcentre_ave, sigma_ave;
};

namespace {
	std::vector<std::string> StringSplit(std::string str) {
		std::stringstream ss{ str };
		std::vector<std::string> v;
		std::string buf;
		while (std::getline(ss, buf, ' ')) {
			if (buf != "") {
				v.push_back(buf);
			}
		}
		return v;
	}
	std::vector<std::string> StringSplit_with_tab(std::string str) {
		std::vector<std::string> v;

		std::vector<std::string> str_v = StringSplit(str);
		for (int i = 0; i < str_v.size(); i++) {
			std::stringstream ss{ str_v[i] };
			std::string buf;
			while (std::getline(ss, buf, '\t')) {
				if (buf != "") {
					v.push_back(buf);
				}
			}
		}
		return v;
	}
}
int angle_to_iang(double& angle) {
	int i_ang = 0;
	if (angle < 1.0) {
		i_ang = angle / 0.1;
	}
	else if (angle < 2.0) {
		i_ang = (angle - 1) / 0.2 + 10;
	}
	else {
		i_ang = (angle - 2) / 0.3 + 15;
	}
	return i_ang;
}
void i_ang_to_angle(int& i_ang, double& angle_min, double& angle_max) {
	if (i_ang < 7) {
		angle_min = i_ang * 0.1;
		angle_max = (i_ang + 1) * 0.1;
	}
	else if (i_ang < 11) {
		angle_min = (i_ang - 7) * 0.2 + 0.7;
		angle_max = (i_ang - 7 + 1) * 0.2 + 0.7;
	}
	else if (i_ang < 15) {
		angle_min = (i_ang - 11) * 0.4 + 1.5;
		angle_max = (i_ang - 11 + 1) * 0.4 + 1.5;
	}
	else if (i_ang < 16) {
		angle_min = (i_ang - 15) * 0.6 + 3.1;
		angle_max = (i_ang - 15 + 1) * 0.6 + 3.1;
	}
	else {
		angle_min = (i_ang - 15) * 0.6 + 3.1;
		angle_max = 20;
	}

}

int angle_to_iang(float& angle) {
	int i_ang = 0;
	if (angle < 1.0) {
		i_ang = angle / 0.1;
	}
	else if (angle < 2.0) {
		i_ang = (angle - 1) / 0.2 + 10;
	}
	else {
		i_ang = (angle - 2) / 0.3 + 15;
	}
	return i_ang;
}
void i_ang_to_angle(int& i_ang, float& angle_min, float& angle_max) {
	if (i_ang < 7) {// <=0.7
		angle_min = i_ang * 0.1;
		angle_max = (i_ang + 1) * 0.1;
	}
	else if (i_ang < 11) {// <= 1.5
		angle_min = (i_ang - 7) * 0.2 + 0.7;
		angle_max = (i_ang - 7 + 1) * 0.2 + 0.7;
	}
	else if (i_ang < 15) {// <= 3.1
		angle_min = (i_ang - 11) * 0.4 + 1.5;
		angle_max = (i_ang - 11 + 1) * 0.4 + 1.5;
	}
	else if (i_ang < 16) {// <= 3.7
		angle_min = (i_ang - 15) * 0.6 + 3.1;
		angle_max = (i_ang - 15 + 1) * 0.6 + 3.1;
	}
	else {// <=20.0
		angle_min = (i_ang - 15) * 0.6 + 3.1;
		angle_max = 20;
	}

}
void set_microtrack_vph(std::string input, std::multimap<Key, VPH_MIP>& map);
void Output_VPH_thr(std::map<Key, VPH_MIP>& btrk, std::string output);
void calcurate_basetrack_vph(std::multimap<Key, VPH_MIP>& map, std::map<Key, VPH_MIP>& btrk);

int main(int argc, char** argv) {

	if (argc < 3) {
		fprintf(stderr, "usage : prg vph_pl_correction_mip.txt output.txt\n");
		exit(1);
	}

	std::string input = argv[1];// VPH of microtrack
	std::string output = argv[2];// VPH thr of basetrack

	std::multimap<Key, VPH_MIP> map;
	set_microtrack_vph(input, map);
	std::map<Key, VPH_MIP>btrk;
	calcurate_basetrack_vph(map, btrk);
	Output_VPH_thr(btrk, output);


}
void set_microtrack_vph(std::string input, std::multimap<Key, VPH_MIP>& map) {

	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error!" << std::endl;
		std::cout << input << std::endl;
		return;
	}

	VPH_MIP m;
	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	int cnt = 0;
	m.k.pl = 3;
	m.face = 1;
	while (std::getline(ifs, str)) {
		// gid cid pid nseg npl upl dpl
		str_v = StringSplit_with_tab(str);
		//std::cout << str << std::endl;
		m.k.pl = std::stoi(str_v[0]);
		m.face = std::stoi(str_v[1]);
		m.k.i_ang = std::stoi(str_v[2]);
		m.vphcenter = std::stod(str_v[3]);
		m.sigma = std::stod(str_v[4]);
		m.vphcentre_ave = std::stod(str_v[5]);
		m.sigma_ave = std::stod(str_v[6]);
		map.insert(std::make_pair(m.k, m));

		if (m.k.i_ang == 14) {
			m.k.i_ang = m.k.i_ang + 1;
			map.insert(std::make_pair(m.k, m));
		}

	}
	std::cout << map.size() << std::endl;
	Key k13, k14, k15, k16;
	for (int i = 3; i < 134; i++) {
		k16.i_ang = 16; k15.pl = i;
		k15.i_ang = 15; k14.pl = i;
		k14.i_ang = 14; k13.pl = i;
		k13.i_ang = 13; k16.pl = i;

		auto i13 = map.equal_range(k13);
		auto i14 = map.equal_range(k14);
		auto i15 = map.equal_range(k15);
		auto i16 = map.equal_range(k16);

		m = map.find(k15)->second;
		map.insert(std::make_pair(k16, map.find(k15)->second));// i:16<--15

		map.find(k15)->second = map.find(k14)->second;// i:15<--14
		m = map.find(k14)->second;
		map.find(k14)->second.vphcenter = 0.5 * (m.vphcenter + map.find(k13)->second.vphcenter);
		map.find(k14)->second.vphcentre_ave = 0.5 * (m.vphcentre_ave + map.find(k13)->second.vphcentre_ave);
	}
}

void calcurate_basetrack_vph(std::multimap<Key, VPH_MIP>& map, std::map<Key, VPH_MIP>& btrk) {

	VPH_MIP m;
	for (int i = 3; i < 134; i++) {
		for (int j = 0; j < 17; j++) {
			m.k.pl = i;
			m.k.i_ang = j;

			m.face = 0;
			m.vphcenter = 0;
			m.sigma = 0;
			m.vphcentre_ave = 0;
			m.sigma_ave = 0;
			//std::cout << map.count(m.k) << std::endl;

			auto p = map.equal_range(m.k);
			for (auto itr = p.first; itr != p.second; itr++) {
				m.face = 0;
				m.vphcenter = m.vphcenter + itr->second.vphcenter;
				m.sigma = m.sigma;
				m.vphcentre_ave = m.vphcentre_ave + itr->second.vphcentre_ave;
				m.sigma_ave = m.sigma_ave;
			}
			btrk.insert(std::make_pair(m.k, m));
		}
	}
}

void Output_VPH_thr(std::map<Key, VPH_MIP>& btrk,std::string output) {
	std::ofstream ofs(output);
	if (!ofs) {
		std::cerr << "File open error!" << std::endl;
		std::cout << output << std::endl;
		return;
	}

	Key k;
	double vph_thr,vph_dist,ang_min,ang_max;
	for (int i = 3; i < 134; i++) {
		for (int j = 0; j < 17; j++) {
			k.i_ang = j;
			k.pl = i;
			vph_dist = btrk.find(k)->second.vphcenter;
			vph_thr = 6 * sqrt(vph_dist) + vph_dist;
			i_ang_to_angle(j, ang_min, ang_max);
			ofs << std::fixed << std::right
				<< std::setw(4) << std::setprecision(0) << k.pl << " "
				<< std::setw(4) << std::setprecision(1) << ang_min << " "
				<< std::setw(4) << std::setprecision(1) << ang_max << " "
				<< std::setw(10) << std::setprecision(4) << btrk.find(k)->second.vphcenter<<" "
				<< std::setw(10) << std::setprecision(4) << sqrt(vph_dist) << " "
				<< std::setw(8) << std::setprecision(1) << vph_thr << " "
				<< std::endl;
		}
	}
}
