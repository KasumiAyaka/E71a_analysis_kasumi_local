#pragma comment(lib, "FILE_structure.lib")
#include <FILE_structure.hpp>


#include <vector>
#include <map>
#include <fstream>
#include <ios>     // std::left, std::right
#include <iomanip> 
#include <math.h>
#include <algorithm>

class output {
public:
	int groupid, chainid, pid;
	double pb, vph, angle;
};

class PID_point {
public:
	int i_ang;
	double angle_min, angle_max, pb, vph_thr;

};
class Chain_inf {
public:
	int groupid, chainid, pid;
	double pb, vph, angle, proton_likelihood, pion_likelihood, likelihood_ratio;
};

std::map<int, std::map<double, PID_point >> Read_PID_point(std::string filename);
void Add_PID(std::vector<Momentum_recon::Event_information>& momch, std::map<int, std::map<double, PID_point >>& pid_point_map);
	void outputfile(std::string filename, std::vector<Chain_inf>& out);
std::map<int, std::vector<Chain_inf>> Read_PID_Lattice(std::string filename);
//void Add_Likelihood(std::vector<Chain_inf>& cid, std::map<int, std::vector<Chain_inf>>& likelihood_point_map);
void outputfile2(std::string filename, std::vector<Chain_inf>& out);
int angle_to_iangle(double angle);
bool judege_iron_penetrate(int pl0, int pl1);

int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "uasge:output-file output-file pion_thr pion_vph_file likelihood-point out-likelihood\n");
		exit(1);
	}
	std::string file_in_momch = argv[1];
	std::string file_in_pid = argv[2];
	std::string file_out_momch = argv[3];
	//std::string file_in_likelihood_point = argv[4];

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch);
	std::map<int, std::map<double, PID_point >> pid_point_map = Read_PID_point(file_in_pid);
	//std::map<int, std::vector<Chain_inf>> likelihood_point = Read_PID_Lattice(file_in_likelihood_point);

	Add_PID(momch, pid_point_map);
	Momentum_recon::Write_Event_information_extension(file_out_momch, momch);

	//outputfile(file_out_pid, c_inf);

	//std::string file_out_likelihood = argv[5];
	//Add_Likelihood(c_inf, likelihood_point);
	//outputfile2(file_out_likelihood, c_inf);

}
std::vector<output> Calc_inf(std::vector<Momentum_recon::Event_information>& momch) {
	std::vector<output> ret;
	for (auto& ev : momch) {
		for (auto& c : ev.chains) {
			output out;
			out.chainid = c.chainid;
			out.groupid = ev.groupid;
			out.pb = c.Get_muon_mcs_pb();
			//out.pb = c.ecc_mcs_mom[0];
			if (!isfinite(out.pb))out.pb = -1;
			out.pid = 0;
			out.vph = 0;
			//if (c.chainid == 0) {
			//	out.pid = c.particle_flg;
			//}
			int count = 0;
			for (auto& b : c.base) {
				out.vph += b.m[0].ph % 10000;
				out.vph += b.m[1].ph % 10000;
				count += 2;
			}
			out.vph /= count;
			if (c.base.rbegin()->pl <= ev.vertex_pl) {
				out.angle = sqrt(pow(c.base.rbegin()->ax, 2) + pow(c.base.rbegin()->ay, 2));
			}
			else {
				out.angle = sqrt(pow(c.base.begin()->ax, 2) + pow(c.base.begin()->ay, 2));
			}
			if (!judege_iron_penetrate(c.base.begin()->pl, c.base.rbegin()->pl))continue;
			ret.push_back(out);
		}
	}
	return ret;
}
bool judege_iron_penetrate(int pl0, int pl1) {
	int num_film = 0;
	int num_iron = 0;
	int num_water = 0;
	for (int pl = pl0; pl < pl1; pl++) {
		num_film++;
		if (pl == 3)continue;
		else if (pl <= 14)num_iron++;
		else if (pl == 15)continue;
		else if (pl % 2 == 0)num_iron++;
		else num_water++;
	}
	return num_iron > 0;

}


std::map<int, std::map<double, PID_point >> Read_PID_point(std::string filename) {
	std::ifstream ifs(filename);
	std::map<int, std::map<double, PID_point >> ret;
	PID_point data;
	while (ifs >> data.i_ang >> data.angle_min >> data.angle_max >> data.pb >> data.vph_thr) {
		int i_ang = int(data.angle_max * 10);
		auto res = ret.find(i_ang);
		if (res == ret.end()) {
			std::map<double, PID_point > map_tmp;
			map_tmp.insert(std::make_pair(data.pb, data));
			ret.insert(std::make_pair(i_ang, map_tmp));
		}
		else {
			res->second.insert(std::make_pair(data.pb, data));
		}
	}
	return ret;
}
std::map<int, std::vector<Chain_inf>> Read_PID_Lattice(std::string filename) {
	Chain_inf p;
	int i_ang;
	double angle;
	std::map<int, std::vector<Chain_inf>> ret;
	std::ifstream ifs(filename);
	int count = 0;
	while (ifs >> p.groupid >> p.chainid >> p.angle >> p.pb >> p.vph >> p.proton_likelihood >> p.pion_likelihood >> p.likelihood_ratio >> p.pid) {
		if (count % 10000 == 0) {
			fprintf(stderr, "\r Read point ... %d", count);
		}
		count++;

		i_ang = angle_to_iangle(p.angle);
		auto res = ret.find(i_ang);
		if (res == ret.end()) {
			std::vector<Chain_inf> vec;
			vec.push_back(p);
			ret.insert(std::make_pair(i_ang, vec));
		}
		else {
			res->second.push_back(p);
		}
	}
	fprintf(stderr, "\r Read point ... %d\n", count);

	return ret;
}

void Add_PID(std::vector<Momentum_recon::Event_information>& momch,std::map<int, std::map<double, PID_point >>& pid_point_map) {
	double a;
	int i_ang, iangcnt = 0;


	for (auto& ev : momch) {
		for (auto& c : ev.chains) {
			output out;
			out.chainid = c.chainid;
			out.groupid = ev.groupid;
			out.pb = c.Get_muon_mcs_pb();
			//out.pb = c.ecc_mcs_mom[0];
			if (!isfinite(out.pb))out.pb = -1;
			out.pid = 0;
			out.vph = 0;
			//if (c.chainid == 0) {
			//	out.pid = c.particle_flg;
			//}
			int count = 0;
			for (auto& b : c.base) {
				out.vph += b.m[0].ph % 10000;
				out.vph += b.m[1].ph % 10000;
				count += 2;
			}
			out.vph /= count;
			if (c.base.rbegin()->pl <= ev.vertex_pl) {
				out.angle = sqrt(pow(c.base.rbegin()->ax, 2) + pow(c.base.rbegin()->ay, 2));
			}
			else {
				out.angle = sqrt(pow(c.base.begin()->ax, 2) + pow(c.base.begin()->ay, 2));
			}
			if (!judege_iron_penetrate(c.base.begin()->pl, c.base.rbegin()->pl))continue;

			// add pid
			//if(typeid(itr->pb).name()!="double"&& typeid(itr->pb).name() != "int"){
			//if (std::isnan(itr->pb)) {//add kasumi 2025/03/16
			//	std::cout << "eventid = " << itr->groupid << " :pb-value is nan!" << std::endl;
			//	itr->pid = -1;
			//	continue;
			//}
			//if (itr->pb > 700) {//add kasumi 2024/12/24
			//	itr->pid = 0;
			//	continue;
			//}

			i_ang = int(out.angle * 10);
			std::map<double, PID_point > pid_point;
			auto res = pid_point_map.lower_bound(i_ang);
			if (res == pid_point_map.end()) {
				res = std::next(res, -1);
			}
			pid_point = res->second;

			PID_point point;
			auto res2 = pid_point.upper_bound(out.pb);
			//’ápb‚Ìedge
			if (res2 == pid_point.begin()) {
				point = res2->second;
			}
			//‚pb‚Ìedge
			else if (res2 == pid_point.end()) {
				point = std::next(res2, -1)->second;
			}
			else {
				PID_point point2[2];
				point2[0] = std::next(res2, -1)->second;
				point2[1] = res2->second;

				if (fabs(point2[0].pb - out.pb) < fabs(point2[1].pb - out.pb)) {
					point = point2[0];
				}
				else {
					point = point2[1];
				}
			}

			if (point.vph_thr > out.vph) {
				out.pid = 211;
			}
			else {
				out.pid = 2212;
			}

			a = out.pb;
			if (a == 10) {
				// "pb" is "nun" ==> pb=10 : pid-> -1
				std::cout << " " << out.groupid << " " << out.chainid << " " << out.pb << " " << "MeV " << out.vph << " " << out.pid << " " << out.angle << std::endl;
				out.pid = 0;
			}
			/*
			i_ang = angle_to_iangle(out.angle);
			std::vector<Chain_inf> likelihood_point;
			auto res2 = likelihood_point_map.find(i_ang);
			
			if (res2 == likelihood_point_map.end()) {
				res2 = std::next(res, -1);
			}
			likelihood_point = res->second;
			double pb_width = -1, vph_width = -1;
			for (int i = 1; i < likelihood_point.size(); i++) {
				if (pb_width < 0) {
					if (likelihood_point[i].pb - likelihood_point[i - 1].pb > 0) {
						pb_width = likelihood_point[i].pb - likelihood_point[i - 1].pb;
					}
				}
				if (vph_width < 0) {
					if (likelihood_point[i].vph - likelihood_point[i - 1].vph > 0) {
						vph_width = likelihood_point[i].vph - likelihood_point[i - 1].vph;
					}
				}
				if (vph_width > 0 && pb_width > 0)break;
			}
			//printf("\nvph width:%lf, pb widht:%lf\n", vph_width, pb_width);

			Chain_inf likelihood_val = *likelihood_point.begin();
			double dist = DBL_MAX;
			for (auto& p : likelihood_point) {
				double dist_now = pow((p.pb - itr->pb) / pb_width, 2) + pow((p.vph - itr->vph) / vph_width, 2);
				if (dist_now < dist) {
					dist = dist_now;
					likelihood_val = p;
				}
			}
			itr->pion_likelihood = likelihood_val.pion_likelihood;
			itr->proton_likelihood = likelihood_val.proton_likelihood;
			itr->likelihood_ratio = likelihood_val.pion_likelihood / (likelihood_val.pion_likelihood + likelihood_val.proton_likelihood);
			*/
			c.particle_flg = out.pid;


		}
	}





}
void Add_Likelihood(std::vector<Momentum_recon::Event_information>& momch, std::vector<Chain_inf>& cid, std::map<int, std::vector<Chain_inf>>& likelihood_point_map) {
	int i_ang, cnt = 0;
	for (auto& ev : momch) {
		for (auto& c : ev.chains) {

		}
	}

}
void iang_to_angle_range(int i_ang, double& angle_min, double& angle_max) {
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
	else {
		angle_min = (i_ang - 15) * 0.6 + 3.1;
		angle_max = (i_ang - 15 + 1) * 0.6 + 3.1;
	}
}
int angle_to_iangle(double angle) {
	if (angle < 0.7) {
		return int(angle / 0.1 + 0.1);
	}
	else if (angle < 1.5) {
		return int((angle - 0.7) / 0.2 + 0.1) + 7;
	}
	else if (angle < 3.1) {
		return int((angle - 1.5) / 0.4 + 0.1) + 11;
	}
	else {
		return int((angle - 3.1) / 0.6 + 0.1) + 15;
	}
	return -1;
}

void outputfile(std::string filename, std::vector<Chain_inf>& out) {
	int all = out.size(), cnt = 0;
	std::ofstream ofs(filename);
	for (auto itr = out.begin(); itr != out.end(); itr++) {
		if (cnt % 10000 == 0) {
			fprintf(stderr, "\r Write file %d/%d(%4.1lf%%)", cnt, all, cnt * 100. / all);
		}
		cnt++;
		ofs << std::right << std::fixed
			<< std::setw(5) << std::setprecision(0) << itr->groupid << " "
			<< std::setw(5) << std::setprecision(0) << itr->chainid << " "
			<< std::setw(5) << std::setprecision(0) << itr->pid << " "
			<< std::setw(7) << std::setprecision(4) << itr->angle << " "
			<< std::setw(8) << std::setprecision(1) << itr->pb << " "
			<< std::setw(6) << std::setprecision(1) << itr->vph << std::endl;
	}
	fprintf(stderr, "\r Write file %d/%d(%4.1lf%%)\n", cnt, all, cnt * 100. / all);


}
void outputfile2(std::string filename, std::vector<Chain_inf>& out) {
	int all = out.size(), cnt = 0;
	std::ofstream ofs(filename);
	for (auto itr = out.begin(); itr != out.end(); itr++) {
		if (cnt % 10000 == 0) {
			fprintf(stderr, "\r Write file %d/%d(%4.1lf%%)", cnt, all, cnt * 100. / all);
		}
		cnt++;
		ofs << std::right << std::fixed
			<< std::setw(5) << std::setprecision(0) << itr->groupid << " "
			<< std::setw(5) << std::setprecision(0) << itr->chainid << " "
			<< std::setw(6) << std::setprecision(4) << itr->proton_likelihood << " "
			<< std::setw(6) << std::setprecision(4) << itr->pion_likelihood << " "
			<< std::setw(6) << std::setprecision(4) << itr->likelihood_ratio << std::endl;
	}
	fprintf(stderr, "\r Write file %d/%d(%4.1lf%%)\n", cnt, all, cnt * 100. / all);


}