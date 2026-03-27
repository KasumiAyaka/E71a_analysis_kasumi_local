// 2025/12/05
// kasumi
// Make_VE_flg_list


#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <map>
#include <fstream>
#include <sstream>
#include <set>

struct AngleRange {
	double min, max;
};
bool operator<(const AngleRange& lhs, const AngleRange& rhs) {
	return std::tie(lhs.min, lhs.max) < std::tie(rhs.min, rhs.max);
}

struct LateralThr {
	AngleRange tan;
	double mean, sig;
};

struct Results {
	int ecc;
	int gid, cid, pid, venum;
	int pl0, rid0;
	double tan0, dal0, dl0;
	int vph0, blkflg0;
	int pl1, rid1;
	double tan1, dal1, dl1;
	int vph1, blkflg1;
	double dar0, dar1;
};
struct Label {
	int ecc;
	int gid, cid;
};
bool operator<(const Label& lhs, const Label& rhs) {
	return std::tie(lhs.ecc, lhs.gid, lhs.cid) < std::tie(rhs.ecc, rhs.gid, rhs.cid);
}

struct Criterias {
	int dl, dal, dar, pid;
};
struct VE_flg {
	int venum;
	int veflg;
	//int blkflg;
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

void set_ve_results(std::string input, std::vector<Results>& res);
void set_lateral_thr(std::string input, std::vector<LateralThr>& ang);
void decide_ve_flg(std::vector<Results>& res, std::vector<LateralThr>& lat_thr, std::vector<LateralThr>& rad_thr, double d_dl, double d_dal, double d_dar, std::multimap<Label, VE_flg>& output);
int judge_lateral_pos(double dl, double d);
int judge_lateral_ang(double dal, LateralThr dal_thr, double d);
int judge_radial_ang(double dar, LateralThr dar_thr, double d);
int judge_black(int pid, int blkflg);
LateralThr decide_angle(double tan, std::vector<LateralThr>& ang);
int judge_ve_flg(Criterias c);
	void write_results(std::multimap<Label, VE_flg>& res, std::string output);

	int main(int argc, char** argv) {

		if (argc < 4) {
			fprintf(stderr, "usage : VElist.txt lat_thr.txt output.txt [dl_thr=10.] [dalsig=5.] [dar_sig=6.]\n");
			exit(1);
		}
		std::cout << " argc = " << argc << std::endl;

		std::string in_ve = argv[1];// VPH of microtrack
		std::string in_lat_thr = argv[2];// VPH thr of basetrack
		std::string in_rad_thr = argv[3];// VPH thr of basetrack
		std::string output = argv[4];// VPH thr of basetrack
		double d_dl = 10;
		double d_dal = 5.0;
		double d_dar = 6.0;
		if (argc == 6) {
			d_dl = std::stod(argv[5]);
			d_dal = std::stod(argv[6]);
			d_dar = std::stod(argv[7]);
		}

		std::vector<Results> res;
		set_ve_results(in_ve, res);
		std::vector<LateralThr> lat_thr, rad_thr;
		set_lateral_thr(in_lat_thr, lat_thr); std::cout << "\t* Lat_thr = " << lat_thr.size() << std::endl;
		set_lateral_thr(in_rad_thr, rad_thr); std::cout << "\t* Rad_thr = " << rad_thr.size() << std::endl;
		std::multimap<Label, VE_flg> out;
		decide_ve_flg(res, lat_thr, rad_thr, d_dl, d_dal,d_dar, out);
		write_results(out, output);

		std::cout << "fin." << std::endl;
	}

void set_ve_results(std::string input, std::vector<Results>& res) {
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
		r.cid = std::stoi(str_v[2]);
		r.pid = std::stoi(str_v[3]);
		r.venum = std::stoi(str_v[4]);

		r.pl0 = std::stoi(str_v[5]);
		r.rid0 = std::stoi(str_v[6]);
		r.tan0 = std::stod(str_v[7]);
		r.dal0 = std::stod(str_v[8]);
		r.dl0 = std::stod(str_v[9]);
		r.vph0 = std::stoi(str_v[10]);
		r.blkflg0 = std::stoi(str_v[11]);

		r.pl1 = std::stoi(str_v[12]);
		r.rid1 = std::stoi(str_v[13]);
		r.tan1 = std::stod(str_v[14]);
		r.dal1 = std::stod(str_v[15]);
		r.dl1 = std::stod(str_v[16]);
		r.vph1 = std::stoi(str_v[17]);
		r.blkflg1 = std::stoi(str_v[18]);

		r.dar0 = std::stod(str_v[19]);
		r.dar1 = std::stod(str_v[20]);

		res.push_back(r);
	}

}

void set_lateral_thr(std::string input, std::vector<LateralThr>& ang) {
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error!" << std::endl;
		std::cout << input << std::endl;
		return;
	}

	LateralThr l;
	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	while (std::getline(ifs, str)) {
		// gid cid pid nseg npl upl dpl
		str_v = StringSplit_with_tab(str);
		
	if (str_v.size() == 4) {
			l.mean = std::stod(str_v[2]);
			l.sig = std::stod(str_v[3]);

			l.tan.min = std::stod(str_v[0]);
			l.tan.max = std::stod(str_v[1]);
		}
		else {
			std::cout << input << " has unexpected file format?" << std::endl;
		}
		ang.push_back(l);
	}
}

void decide_ve_flg(std::vector<Results>& res, std::vector<LateralThr>& lat_thr, std::vector<LateralThr>& rad_thr, double d_dl, double d_dal, double d_dar, std::multimap<Label, VE_flg>& output) {
	std::ofstream ofs("Make_VE_flg.log.txt");

	LateralThr dal_thr,dar_thr;
	VE_flg out; Label lb;
	Criterias ct0,ct1;
	int tmp0, tmp1;
	for (auto itr = res.begin(); itr != res.end(); itr++) {
		ct0 = { 0 };
		ct1 = { 0 };
		tmp0 = tmp1 = 0;
		if (itr->pl0 > 0 && itr->tan0 > 0) {
			dal_thr = decide_angle(itr->tan0, lat_thr);
			dar_thr = decide_angle(itr->tan0, rad_thr);
			ct0.dl= judge_lateral_pos(itr->dl0, d_dl);//1/-1
			ct0.dal = judge_lateral_ang(itr->dal0, dal_thr, d_dal);//1/-1
			ct0.dar = judge_radial_ang(itr->dar0, dar_thr, d_dar);//1/-1
			ct0.pid = judge_black(itr->pid, itr->blkflg0);// 1/0/-1

			tmp0 = judge_ve_flg(ct0);//1 or 0
		}
		else {
			ct0.dl = ct0.dal = ct0.dar = ct0.pid;
		}

		if (itr->pl1 > 0 && itr->tan1 > 0) {
			dal_thr = decide_angle(itr->tan1, lat_thr);
			dar_thr = decide_angle(itr->tan1, rad_thr);
			ct1.dl = judge_lateral_pos(itr->dl1, d_dl);//1/-1
			ct1.dal = judge_lateral_ang(itr->dal1, dal_thr, d_dal);//1/-1
			ct1.dar = judge_radial_ang(itr->dar1, dar_thr, d_dar);//1/-1
			ct1.pid = judge_black(itr->pid, itr->blkflg1);// 1/0/-1

			tmp1 = judge_ve_flg(ct1);//1 or 0
		}
		else {
			ct1.dl = ct1.dal = ct1.dar = ct1.pid = 0;
		}

		ofs << std::fixed << std::right//<< std::setfill(' ')
			<< std::setw(2) << std::setprecision(0) << itr->ecc << " "
			<< std::setw(5) << std::setprecision(0) << itr->gid << " "
			<< std::setw(3) << std::setprecision(0) << itr->cid << " "
			<< std::setw(4) << std::setprecision(0) << itr->pid << " "
			<< std::setw(4) << std::setprecision(0) << itr->blkflg0 << " "
			<< std::setw(4) << std::setprecision(0) << itr->blkflg1 << " "
			<< std::setw(2) << std::setprecision(0) << itr->venum << "  upside:"
			<< std::setw(2) << std::setprecision(0) << ct0.dl << " "
			<< std::setw(2) << std::setprecision(0) << ct0.dal << " "
			<< std::setw(2) << std::setprecision(0) << ct0.dar << " "
			<< std::setw(2) << std::setprecision(0) << ct0.pid << " "
			<< std::setw(2) << std::setprecision(0) << ct1.dl << " "
			<< std::setw(2) << std::setprecision(0) << ct1.dal << " "
			<< std::setw(2) << std::setprecision(0) << ct1.dar << " "
			<< std::setw(2) << std::setprecision(0) << ct1.pid <<  "  ve result:"
			<< std::setw(2) << std::setprecision(0) << tmp0 << " "
			<< std::setw(2) << std::setprecision(0) << tmp1 << "  FLG:"
			<< std::setw(2) << std::setprecision(0) << tmp0+tmp1 << " "
			<< std::endl;


		lb.ecc = itr->ecc; lb.gid = itr->gid; lb.cid = itr->cid;
		out.venum = itr->venum;
		out.veflg = tmp0 + tmp1;
		output.insert(std::make_pair(lb, out));
	}
}

int judge_lateral_pos(double dl, double d) {
	int flg = 0;
	if (dl*dl>d*d) {
		// remain
		flg = -1;
	}
	else {
		// ve
		flg = 1;
	}

	return flg;
}
int judge_lateral_ang(double dal, LateralThr dal_thr, double d) {
	int flg = 0;
	//if (dal < dl_thr.mean - d * dl_thr.sig || dal > dl_thr.mean + d * dl_thr.sig) {
	if (fabs( dal) > d * dal_thr.sig) {
		// remain
		flg = -1;
	}
	else {
		// ve
		flg = 1;
	}

	return flg;
}
int judge_black(int pid, int blkflg) {
	int flg = 0;
	if (pid == 2212) {//blk
		if (blkflg == 1) {//blk-->ok
			flg = 1;
		}
		else if (blkflg == -1) {//mip-->x
			flg = -1;
		}
		else {
			std::cout << "judge of vph if not going well...?" << std::endl;
		}
	}

	if (pid == 211) {
		if (blkflg == -1) {
			flg = 1;
		}
		else  if (blkflg == 1) {
			flg = -1;
		}
		else {
			std::cout << "judge of vph if not going well...?" << std::endl;
		}

	}

	return flg;
}
LateralThr decide_angle(double tan, std::vector<LateralThr>& ang) {
	int flg = 0;
	LateralThr ret;
	ret = ang[0];
	for (auto i = 0; i < ang.size(); i++) {
		if (tan > ang[i].tan.min && tan < ang[i].tan.max) {
			flg++;
			ret = ang[i];
		}
	}
	if (flg != 1) {
		std::cout << "unexpected tan " << tan << std::endl;
	}
	return ret;
}
int judge_radial_ang(double dar, LateralThr dar_thr, double d) {
	int flg = 0;
	if (fabs(dar) >  d * dar_thr.sig) {
		// remain
		flg = -1;
	}
	else {
		// ve
		flg = 1;
	}

	return flg;
}



int judge_ve_flg(Criterias c) {
	int ret = 1;
	if (c.dl == -1 || c.dal == -1 || c.dar == -1||c.pid==-1) {
		ret = 0;
	}
	//if (flg[0] + flg[1]==2) {
	//	ret = 1;
	//}
	//if (ret == 1 && flg[2] == -1) {
	//	ret = 0;
	//}
	//if (ret == 1 && flg[2] == 0) {
	//	ret = 10;
	//}

	return ret;
}

void write_results(std::multimap<Label, VE_flg>& res, std::string output) {
	std::ofstream ofs(output);
	std::set<Label>l;
	for (auto itr = res.begin(); itr != res.end(); itr++) {
		l.insert(itr->first);
	}
	int c1, c2, c3;
	for (auto itr = l.begin(); itr != l.end(); itr++) {
		ofs << std::fixed << std::right//<< std::setfill(' ')
			<< std::setw(2) << std::setprecision(0) << itr->ecc << " "
			<< std::setw(5) << std::setprecision(0) << itr->gid << " "
			<< std::setw(3) << std::setprecision(0) << itr->cid << " ";

		c1 = 0;
		auto p = res.equal_range(*itr);
		for (auto q = p.first; q != p.second; q++) {
			if (q->second.venum == 1) {
				ofs << std::fixed << std::right//<< std::setfill(' ')
					<< std::setw(2) << std::setprecision(0) << q->second.veflg << " ";
				c1++;
			}
		}
		if (c1 == 0) ofs << std::fixed << std::right << std::setw(2) << std::setprecision(0) << 0 << " ";
		
		c2 = 0;
		auto p2 = res.equal_range(*itr);
		for (auto q = p2.first; q != p2.second; q++) {
			if (q->second.venum == 2) {
				ofs << std::fixed << std::right//<< std::setfill(' ')
					<< std::setw(2) << std::setprecision(0) << q->second.veflg << " ";
				c2++;
			}
		}
		if (c2 == 0) ofs << std::fixed << std::right << std::setw(2) << std::setprecision(0) << 0 << " ";
		
		c3 = 0;
		auto p3 = res.equal_range(*itr);
		for (auto q = p3.first; q != p3.second; q++) {
			if (q->second.venum == 3) {
				ofs << std::fixed << std::right//<< std::setfill(' ')
					<< std::setw(2) << std::setprecision(0) << q->second.veflg
					<< std::endl;
				c3++;
			}
		}
		if (c3 == 0) ofs << std::fixed << std::right << std::setw(2) << std::setprecision(0) << 0 << "\n";
	}
}
