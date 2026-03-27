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

struct AngleRange {
	double ang_min, ang_max;
};
bool operator<(const AngleRange& lhs, const AngleRange& rhs) {
	return std::tie(lhs.ang_min, lhs.ang_max) < std::tie(rhs.ang_min, rhs.ang_max);
}

struct VPH_thr {
	int pl;
	AngleRange a;
	double vph, vph_sig, vph_thr;
};

struct Results {
	int ecc;
	int gid, cid, pid, venum;
	int pl0, rid0;
	double tan0, dal0, dl0;
	int vph0, vph00;
	int pl1, rid1;
	double tan1, dal1, dl1;
	int vph1, vph11;
	double dar0, dar1;
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

void set_vph_thr(std::string input, std::multimap<int, VPH_thr>& map);
void judge_black(std::string input, std::vector<Results>& res);
void add_black_flg(std::vector<Results>& res, std::multimap<int, VPH_thr>& thr);
void write_results(std::vector<Results>& res, std::string output);

int main(int argc, char** argv) {

	if (argc < 4) {
		fprintf(stderr, "usage : prg ve.txt mip_thr.txt output.txt\n");
		exit(1);
	}

	std::string in_ve = argv[1];// ve
	std::string in_vph_thr = argv[2];// thr
	std::string output = argv[3];// VPH thr of basetrack

	std::multimap<int, VPH_thr> map;
	std::vector<Results> res;
	set_vph_thr(in_vph_thr, map);
	judge_black(in_ve, res);
	add_black_flg(res, map);
	write_results(res, output);
	std::cout << "fin." << std::endl;
}

void set_vph_thr(std::string input, std::multimap<int, VPH_thr>& map) {

	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error!" << std::endl;
		std::cout << input << std::endl;
		return;
	}

	VPH_thr m;
	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	while (std::getline(ifs, str)) {
		// gid cid pid nseg npl upl dpl
		str_v = StringSplit_with_tab(str);
		//std::cout << str << std::endl;

		m.pl = std::stoi(str_v[0]);
		m.a.ang_min = std::stod(str_v[1]);
		m.a.ang_max = std::stod(str_v[2]);
		m.vph = std::stod(str_v[3]);
		m.vph_sig= std::stod(str_v[4]);
		m.vph_thr = std::stod(str_v[5]);
		map.insert(std::make_pair(m.pl, m));
	}

	std::cout << map.size() << std::endl;
}
void judge_black(std::string input, std::vector<Results>& res) {
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error!" << std::endl;
		std::cout << input << std::endl;
		return;
	}

	Results r;
	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	while (std::getline(ifs, str)) {
		// gid cid pid nseg npl upl dpl
		str_v = StringSplit_with_tab(str);
		//std::cout << str << std::endl;

		r.ecc = std::stoi(str_v[0]);
		r.gid = std::stoi(str_v[1]);
		r.cid= std::stoi(str_v[2]);
		r.pid = std::stoi(str_v[3]);
		r.venum = std::stoi(str_v[4]);

		r.pl0 = std::stoi(str_v[5]);
		r.rid0 = std::stoi(str_v[6]); 
		r.tan0 = std::stod(str_v[7]);
		r.dal0 = std::stod(str_v[8]);
		r.dl0 = std::stod(str_v[9]);
		r.vph0 = std::stoi(str_v[10]);
		r.vph00 = std::stoi(str_v[11]);

		r.pl1 = std::stoi(str_v[12]);
		r.rid1 = std::stoi(str_v[13]);
		r.tan1 = std::stod(str_v[14]);
		r.dal1 = std::stod(str_v[15]);
		r.dl1 = std::stod(str_v[16]);
		r.vph1 = std::stoi(str_v[17]);
		r.vph11 = std::stoi(str_v[18]);
		
		r.dar0 = std::stod(str_v[19]);
		r.dar1 = std::stod(str_v[20]);

		res.push_back(r);
	}

}

void add_black_flg(std::vector<Results>& res, std::multimap<int, VPH_thr>& thr) {

	int blkflg;
	int pl00, pl11;
	for (auto itr = res.begin(); itr != res.end(); itr++) {
		blkflg = 0;

		// upstream
		if (itr->pl0 > 0 && itr->tan0 != -10) {
			if (itr->pl0 % 2 == 1) {//even
				pl00 = itr->pl0 - 3;
			}
			else {
				pl00 = itr->pl0 + 3;
			}

			auto p0 = thr.equal_range(pl00);
			for (auto p = p0.first; p != p0.second; p++) {
				if (itr->tan0 > p->second.a.ang_min && itr->tan0 < p->second.a.ang_max) {
					if (itr->vph00 > p->second.vph_thr) {
						// blkflg
						blkflg = 1;
					}
					else {
						blkflg = -1;
					}
					std::cout << itr->cid << " " << itr->venum << " " << itr->vph00 << " " << p->second.vph_thr << " " << blkflg << "  " << itr->tan0 << " " << p->second.a.ang_min << std::endl;
					itr->vph00 = blkflg;
				}
			}
		}
		blkflg = 0;

		// downstream
		if (itr->pl1 > 0 && itr->tan1 != -10) {
			if (itr->pl1 % 2 == 1) {//even
				pl11 = itr->pl1 - 3;
			}
			else {
				pl11 = itr->pl1 + 3;
			}
			auto p1 = thr.equal_range(pl11);
			for (auto p = p1.first; p != p1.second; p++) {
				if (itr->tan1 > p->second.a.ang_min && itr->tan1 < p->second.a.ang_max) {
					if (itr->vph11 > p->second.vph_thr) {
						// blkflg
						blkflg = 1;
					}
					else {
						blkflg = -1;
					}
					itr->vph11 = blkflg;
				}
			}
		}

	}
}
void write_results(std::vector<Results>& res,std::string output) {
	std::ofstream ofs(output);
	for (auto itr = res.begin(); itr != res.end(); itr++) {
		ofs << std::fixed << std::right//<< std::setfill(' ')
			<< std::setw(2) << std::setprecision(0) << itr->ecc << " "
			<< std::setw(5) << std::setprecision(0) << itr->gid << " "
			<< std::setw(3) << std::setprecision(0) << itr->cid << " "
			<< std::setw(4) << std::setprecision(0) << itr->pid << " "
			<< std::setw(2) << std::setprecision(0) << itr->venum << " "

			<< std::setw(3) << std::setprecision(0) << itr->pl0 << " "
			<< std::setw(10) << std::setprecision(0) << itr->rid0 << " "
			<< std::setw(8) << std::setprecision(4) << itr->tan0 << " "
			<< std::setw(8) << std::setprecision(4) << itr->dal0 << " "
			<< std::setw(8) << std::setprecision(1) << itr->dl0 << " "
			<< std::setw(4) << std::setprecision(0) << itr->vph0 << " "
			<< std::setw(4) << std::setprecision(0) << itr->vph00 << " "

			<< std::setw(3) << std::setprecision(0) << itr->pl1 << " "
			<< std::setw(10) << std::setprecision(0) << itr->rid1 << " "
			<< std::setw(8) << std::setprecision(4) << itr->tan1 << " "
			<< std::setw(8) << std::setprecision(4) << itr->dal1 << " "
			<< std::setw(8) << std::setprecision(1) << itr->dl1 << " "
			<< std::setw(4) << std::setprecision(0) << itr->vph1 << " "
			<< std::setw(4) << std::setprecision(0) << itr->vph11 

			<< std::setw(10) << std::setprecision(4) << itr->dar0 << " "
			<< std::setw(10) << std::setprecision(4) << itr->dar1 << " "
			<< std::endl;
	}
}
