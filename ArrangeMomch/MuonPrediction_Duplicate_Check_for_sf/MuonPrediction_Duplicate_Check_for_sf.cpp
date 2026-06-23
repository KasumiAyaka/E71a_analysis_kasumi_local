// 2026/02/09
// kasumi
// Add_tumestamp_to_momch
// Input SharingFile and Momch, then add timestamp (and momentum,charge) to Momch.

# include <iostream>
# include <random>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>
#include <chrono>
#include <list>
#include <cassert>
#include <filesystem>

#include <set>
#include <map>
#include <tuple>

struct SF {
	int PL;
	int64_t ISS, OSS, FIX, NT;//rawid
	int Spot;//stoptA*100+SpotB
	int zone1, raw1, zone2, raw2;//micro
	int64_t utime;
	int t_trackid, bunch, entry_in_daily_file, nplane, charge;
	double mom;//BM
	double chi2[4];// shifter
	int eid;//event id
	int PMtype;//0:ecclike,1:sandmu,2:upstreamint
	int ECCtype;//0 stop, 1 pene, 2 edgeout, -2 shower, -1 notuse
};


struct StartEnd {
	int PL;
	int64_t ISS;//rawid
	int64_t utime;
	int t_trackid, bunch, entry_in_daily_file;
};
bool operator<(const StartEnd& lhs, const StartEnd& rhs) {
	return std::tie(lhs.utime, lhs.bunch, lhs.PL, lhs.ISS, lhs.t_trackid, lhs.entry_in_daily_file)
		< std::tie(rhs.utime, rhs.bunch, rhs.PL, rhs.ISS, rhs.t_trackid, rhs.entry_in_daily_file);
}

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


void MeasureProcessingTime(std::chrono::system_clock::time_point& start, std::chrono::system_clock::time_point& end) {
	auto dur = end - start;        // 要した時間を計算
	auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
	// 要した時間をミリ秒（1/1000秒）に変換して表示
	std::cout << msec << " milli sec \n";
	if (msec / 1000 < 60) {
		std::cout << msec / 1000 << "sec\n";
	}
	else if (msec / 1000 / 60 < 60) {
		std::cout << msec / 1000 / 60 << "min\n";
	}
	else if (msec / 1000 / 3600 < 24) {
		std::cout << msec / 1000 / 3600 << "h\n";
	}
	else {
		std::cout << (msec / 1000 / 3600) / 24 << "day" << (msec / 1000 / 3600) % 24 << "h\n";
	}
};


void read_sf(std::string in);



int main(int argc, char** argv) {
	if (argc != 2) {
		fprintf(stderr, "usage:prg in_sf.txt\n");
		exit(1);
	}

	std::string in_sf = argv[1];

	auto start = std::chrono::system_clock::now();//for measure working time

	read_sf(in_sf);

	auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
	MeasureProcessingTime(start, end);
}



void read_sf(std::string in) {

	std::ifstream ifs(in);
	if (!ifs) {
		std::cerr << " File open error ! " << std::endl;
		return;
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;
	SF s;
	StartEnd tmp;
	std::set<uint16_t> timestamp;
	std::set<std::pair<uint16_t, int>> tsbn;
	std::set<StartEnd> se;
	int cnt = 0;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);
		if (str_v.size() == 24) {
			s.PL = std::stoi(str_v[0]);
			s.ISS = std::stoll(str_v[1]);
			s.OSS = std::stoll(str_v[2]);
			s.FIX = std::stoll(str_v[3]);
			s.NT = std::stoll(str_v[4]);
			s.Spot = std::stoi(str_v[5]);
			s.zone1 = std::stoi(str_v[6]);
			s.raw1 = std::stoi(str_v[7]);
			s.zone2 = std::stoi(str_v[8]);
			s.raw2 = std::stoi(str_v[9]);
			s.utime = std::stoll (str_v[10]);
			s.t_trackid = std::stoi(str_v[11]);
			s.bunch = std::stoi(str_v[12]);
			s.entry_in_daily_file = std::stoi(str_v[13]);
			s.nplane = std::stoi(str_v[14]);
			s.charge = std::stoi(str_v[15]);
			s.mom = std::stod(str_v[16]);
			s.chi2[0] = std::stod(str_v[17]);
			s.chi2[1] = std::stod(str_v[18]);
			s.chi2[2] = std::stod(str_v[19]);
			s.chi2[3] = std::stod(str_v[20]);
			s.eid = std::stoi(str_v[21]);
			s.PMtype = std::stoi(str_v[22]);
			s.ECCtype = std::stoi(str_v[23]);
			tmp.bunch = s.bunch;
			tmp.entry_in_daily_file = s.entry_in_daily_file;
			tmp.ISS = s.ISS;
			tmp.PL = s.PL;
			tmp.t_trackid = s.t_trackid;
			tmp.utime = s.utime;
			timestamp.insert(s.utime);
			tsbn.insert(std::make_pair(s.utime, s.bunch));
			se.insert(tmp);
			cnt++;

		}
		else if (str_v.size() == 30) {
		}

	}

	std::cout << "\t* The Number of muon prediction is " <<cnt << std::endl;
	std::cout << "timestamp uniqe       : " << timestamp.size() << std::endl;
	std::cout << "timestamp-bunch uniqe : " << tsbn.size() << std::endl;
	std::cout << "StartEnd              : " << se.size() << std::endl;



}
