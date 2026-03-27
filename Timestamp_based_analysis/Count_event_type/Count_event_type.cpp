// 2026/2/16
// kasumi
// Count_event_type
// タイムスタンプ・bunchベースでmuon predictionがどういう結果になったのかを調べる。
// *** output ***
// unixtime bunch ecc muonid vpl Areaflg materialflg #stop #penetrate #edgeout #others
// Areaflg : Central = 1;Edge = -1;
// materialflg : water = 0;iron = 2; FeECC = 1; out of ecc = -1;

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <chrono>
#include <filesystem>
#include <iomanip>


struct Key {
	int utime, bunch;
};
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.utime, lhs.bunch) < std::tie(rhs.utime, rhs.bunch);
}

struct Rfile {
	int eccnum, eid, vpl, areaflg, muonflg;
	// areaflg : Central = 1 / Edge = -1
	// muonflg : Stop = 1 / located = 11,12,3 / Edgeout = -2 / pene = -1 / not connected = -10
};
bool operator<(const Rfile& lhs, const Rfile& rhs) {
	return std::tie(lhs.eccnum, lhs.eid, lhs.vpl, lhs.areaflg, lhs.muonflg)
		< std::tie(rhs.eccnum, rhs.eid, rhs.vpl, rhs.areaflg, rhs.muonflg);
}
struct OutputFormat {
	int utime, bunch, eccnum, eid, vpl, areaflg, materialflg, located, stopnum, penenum, edgenum, othnum;
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

void Read_result_file(std::string input, std::multimap<Key, Rfile>& map);
void Count_events(std::multimap<Key, Rfile>& map, std::map<Key, OutputFormat>& out);
void Output(std::string output, std::map<Key, OutputFormat>& out);

int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "usage:prg in.momch output.momch\n");
		exit(1);
	}
	std::string in_txt = argv[1];// input momch
	std::string out_txt = argv[2];// output momch

	//bool result = std::filesystem::create_directories(out_momch);
	//assert(result);
	//assert(std::filesystem::exists(out_momch));//directryの存在確認
	//assert(std::filesystem::is_directory(out_momch));//指定されたパスがディレクトリを指しているかを確認する

	auto start = std::chrono::system_clock::now();//for measure working time
	std::multimap<Key, Rfile> mp;
	Read_result_file(in_txt, mp);
	std::map<Key, OutputFormat> out;
	Count_events(mp,out);
	Output(out_txt, out);


	auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
	MeasureProcessingTime(start, end);

}

void Read_result_file(std::string input, std::multimap<Key, Rfile>& map) {

	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << " File open error ! " << std::endl;
		return;
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	Key k;
	Rfile r;
	while (std::getline(ifs, str)) {
		str_v = StringSplit_with_tab(str);
		if (str_v.size() != 7) continue;
		k.utime = std::stoi(str_v[0]);
		k.bunch = std::stoi(str_v[1]);
		r.eccnum = std::stoi(str_v[2]);
		r.eid = std::stoi(str_v[3]);
		r.vpl = std::stoi(str_v[4]);
		r.muonflg = std::stoi(str_v[5]);
		r.areaflg = std::stoi(str_v[6]);
		map.insert(std::make_pair(k,r));
	}

	std::cout << "\t* The Number of Events(muon prediction base) =" << map.size() << std::endl;

}

void Count_events(std::multimap<Key, Rfile>& map, std::map<Key, OutputFormat>& out) {
	std::set<Key>tmp;
	for (auto itr = map.begin(); itr != map.end(); itr++) {
		tmp.insert(itr->first);
	}

	OutputFormat o;
	int flg = 0;
	for (auto itr = tmp.begin(); itr != tmp.end(); itr++) {
		auto p = map.equal_range(*itr);
		o.stopnum = o.edgenum = o.penenum = o.located = o.othnum = 0;
		flg = 0;
		for (auto itr0 = p.first; itr0 != p.second; itr0++) {
			if (itr0->second.muonflg > 1) {
				o.stopnum++;
				o.utime = itr0->first.utime;
				o.bunch = itr0->first.bunch;
				o.eccnum = itr0->second.eccnum;
				o.eid = itr0->second.eid;
				o.vpl = itr0->second.vpl;
				o.areaflg = itr0->second.areaflg;
				flg++;
			}
			else if (itr0->second.muonflg == 1) {
				o.located++;
			}
			else if (itr0->second.muonflg == -2) {
				o.edgenum++;
			}
			else if (itr0->second.muonflg == -1) {
				o.penenum++;
			}
			else if (itr0->second.muonflg == -10) {
				o.othnum++;
			}
			else {
				std::cout << "\t unextoected format: \n" << itr0->first.utime << " " << itr0->first.bunch << " "
					<< itr0->second.eid << " " << itr0->second.muonflg << std::endl;

			}
			o.materialflg = itr0->second.muonflg;

			if (flg == 0) {
				o.utime = itr0->first.utime;
				o.bunch = itr0->first.bunch;
				o.eccnum = itr0->second.eccnum;
				o.eid = itr0->second.eid;
				o.vpl = itr0->second.vpl;
				o.areaflg = itr0->second.areaflg;
			}

		}
		out.insert(std::make_pair(*itr, o));
	}

	std::cout << "\t* The number of Events(timestamp&bunch base) = " << out.size() << std::endl;
}
void Output(std::string output,std::map<Key, OutputFormat>& out) {
	std::ofstream ofs(output);
	if (!ofs) {
		std::cerr << "Faild to open " << output << std::endl;
		exit(0);
	}
	int nu = 0;//located neutrino int
	int stop = 0;// stop but doesn't exist in the final output
	int pene = 0;
	int ot = 0;// missed events
	for (auto itr = out.begin(); itr != out.end(); itr++) {
		ofs << std::right << std::fixed
			<< std::setw(12) << std::setprecision(0) << itr->first.utime << " "
			<< std::setw(3) << std::setprecision(0) << itr->first.bunch << " "
			<< std::setw(6) << std::setprecision(0) << itr->second.eccnum << " "
			<< std::setw(4) << std::setprecision(0) << itr->second.eid << " "
			<< std::setw(4) << std::setprecision(0) << itr->second.vpl << " "
			<< std::setw(4) << std::setprecision(0) << itr->second.areaflg << " "
			<< std::setw(4) << std::setprecision(0) << itr->second.materialflg << " "
			<< std::setw(4) << std::setprecision(0) << itr->second.located << " "
			<< std::setw(4) << std::setprecision(0) << itr->second.stopnum << " "
			<< std::setw(4) << std::setprecision(0) << itr->second.penenum << " "
			<< std::setw(4) << std::setprecision(0) << itr->second.edgenum << " "
			<< std::setw(4) << std::setprecision(0) << itr->second.othnum
			<< std::endl;
		
		// unixtime bunch ecc muonid vpl Areaflg materialflg #stop #penetrate #edgeout #others
		if (itr->second.located > 0) {// exist in final output file
			nu++;
		}
		else{
			if (itr->second.stopnum > 0) {// ecc stop
				stop++;
			}
			else {
				if (itr->second.penenum + itr->second.edgenum > 0) {// pene
					pene++;
				}
				else {
					if (itr->second.othnum > 0) {
						ot++;
					}
				}
			}
		}
	}


	std::cout << "\t* Counts" << std::endl;
	std::cout << "\t    Located neutrino interaction = " << nu << std::endl;
	std::cout << "\t    ECC stop but missed event    = " << stop << std::endl;
	std::cout << "\t    Penetrate / Side out event   = " << pene << std::endl;
	std::cout << "\t    Missed events                = " << ot << std::endl;
	std::cout << "\t    ------------------------------> In Total = " << nu + stop + pene + ot << std::endl;

}