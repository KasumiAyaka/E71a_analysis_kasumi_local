//“r’†
// 2026/06/26


#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <iomanip>
#include <set>

struct Timestamp {
	int utime, bunch, ecc, eid;
};

struct EventInfo {
	int ecc, eid, pl;
	double x, y;
};

struct EventSummary {
	Timestamp t;
	int pl;
	int areaflg, muonflg;
	//area : central==1,edge==-1
	//muon : stop==1,edgeout==-2,penetrate==-1
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

void SetTimeStamp(std::string input, std::map<int, EventSummary>& eve);
void SetEventInfo(std::string input, std::multimap<int, EventInfo>& eve);
void SetEventSummary(std::string input, std::map<int, EventSummary>& eve);
void AddEventInfo(std::map<int, EventSummary>& sf, std::multimap<int, EventInfo>& eve, int mflg);
void Write(std::string output, std::map<int, EventSummary>& map);
void SetEventSummary_utime(std::string input, std::multimap<int, EventSummary>& eve);
void DiscardSandMuon(std::string input, std::multimap<int, EventSummary>& map, std::map<int, EventSummary>& out);
void Count(std::multimap<int, EventSummary>& eve, std::string output);

int main(int argc, char** argv) {

	if (argc < 1 || argc>6) {
		fprintf(stderr, " *** usage *** \n");
		fprintf(stderr, " * timestamp-info-->eventSummary\nusage:input-timestamp.txt output.txt\n\n");
		fprintf(stderr, " * add event-info to timestamp-info\nusage:input-timestamp.txt in-event.txt output.txt muonflg\nmuonflg=1:stop, =-2:edgeout, =-1:penetrate\n\n");
		fprintf(stderr, " * add event-info to EventSummary\nusage:input-eventSummary.txt in-event.txt output.txt muonflg 1\n\n");
		fprintf(stderr, " * sandmuon erase\nusage:input-eventSummary.txt in-sandmuon_timestamplist.txt output.txt \n\n");
		exit(1);
	}

	std::string input1, input2, output;
	int mode, flg;

	if (argc == 3) {
		std::cout << " input:sf-->output:EventSummary " << std::endl;
		input1 = argv[1];
		output = argv[2];

		std::map<int, EventSummary> sf;
		SetTimeStamp(input1, sf);
		Write(output, sf);
	}
	else if (argc == 5) {
		std::cout << " add event-info to sf" << std::endl;
		//add event info
		input1 = argv[1];
		input2 = argv[2];
		output = argv[3];
		flg = std::stoi(argv[4]);
		std::map<int, EventSummary> sf;
		std::multimap<int, EventInfo> eve;
		SetTimeStamp(input1, sf);
		SetEventInfo(input2, eve);

		AddEventInfo(sf, eve, flg);
		Write(output, sf);
	}
	else if (argc == 6) {
		std::cout << " add event-info to EventSumaary" << std::endl;
		//add event info
		input1 = argv[1];
		input2 = argv[2];
		output = argv[3];
		flg = std::stoi(argv[4]);

		std::map<int, EventSummary> sf;
		std::multimap<int, EventInfo> eve;
		SetEventSummary(input1, sf);
		SetEventInfo(input2, eve);
		AddEventInfo(sf, eve, flg);
		Write(output, sf);
	}
	else if (argc == 4) {
		input1 = argv[1];
		input2 = argv[2];
		output = argv[3];
		std::multimap<int, EventSummary> es;
		std::map<int, EventSummary> sf;
		SetEventSummary_utime(input1, es);
		DiscardSandMuon(input2, es, sf);
		Write(output, sf);
	}
	else if (argc == 2) {
		input1 = argv[1];
		std::multimap<int, EventSummary> es;
		std::cout << "input output file name :";
		std::cin >> output;
		SetEventSummary_utime(input1, es);
		Count(es, output);
	}
	else {
		std::cout << "Use Correct mode." << std::endl;
	}
}

void SetTimeStamp(std::string input, std::map<int, EventSummary>& eve) {

	std::ifstream ifs(input);
	EventSummary e;
	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;
	if (!ifs) {
		std::cerr << "File Open Error!  " << input << std::endl;
	}
	while (std::getline(ifs, str)) {
		str_v = StringSplit_with_tab(str);
		//std::cout << "ok?" << std::endl;

		e.t.utime = std::stoi(str_v[0]);
		e.t.bunch = std::stoi(str_v[1]);
		e.t.ecc = std::stoi(str_v[2]);
		e.t.eid = std::stoi(str_v[3]);

		e.areaflg = -10;
		e.muonflg = -10;
		e.pl = -10;
		eve.insert(std::make_pair(e.t.eid, e));
	}
	std::cout << " * Finish input Timestamp, size : " << eve.size() << std::endl;
}

void SetEventInfo(std::string input, std::multimap<int, EventInfo>& eve) {

	std::ifstream ifs(input);
	EventInfo e;
	if (!ifs) {
		std::cerr << "File Open Error!  " << input << std::endl;
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;
	while (std::getline(ifs, str)) {
		str_v = StringSplit_with_tab(str);
		//e.ecc= std::stoi(str_v[0]);
		e.eid = std::stoi(str_v[0]);
		e.pl = std::stoi(str_v[1]);
		e.x = std::stod(str_v[2]);
		e.y = std::stod(str_v[3]);
		eve.insert(std::make_pair(e.eid, e));
		//std::cout << "ok?" << std::endl;

	}
	std::cout << " * Finish input EventInfo, size : " << eve.size() << std::endl;

}

void SetEventSummary(std::string input, std::map<int, EventSummary>& eve) {

	std::ifstream ifs(input);
	EventSummary e;
	if (!ifs) {
		std::cerr << "File Open Error!  " << input << std::endl;
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;
	while (std::getline(ifs, str)) {
		str_v = StringSplit_with_tab(str);

		e.t.utime = std::stoi(str_v[0]);
		e.t.bunch = std::stoi(str_v[1]);
		e.t.ecc = std::stoi(str_v[2]);
		e.t.eid = std::stoi(str_v[3]);

		e.areaflg = std::stoi(str_v[6]);
		e.muonflg = std::stoi(str_v[5]);
		e.pl = std::stoi(str_v[4]);
		eve.insert(std::make_pair(e.t.eid, e));
	}
	std::cout << " * Finish input EventSummary, size : " << eve.size() << std::endl;

}
void SetEventSummary_utime(std::string input, std::multimap<int, EventSummary>& eve) {

	std::ifstream ifs(input);
	EventSummary e;
	if (!ifs) {
		std::cerr << "File Open Error!  " << input << std::endl;
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;
	while (std::getline(ifs, str)) {
		str_v = StringSplit_with_tab(str);

		e.t.utime = std::stoi(str_v[0]);
		e.t.bunch = std::stoi(str_v[1]);
		e.t.ecc = std::stoi(str_v[2]);
		e.t.eid = std::stoi(str_v[3]);

		e.areaflg = std::stoi(str_v[6]);
		e.muonflg = std::stoi(str_v[5]);
		e.pl = std::stoi(str_v[4]);
		eve.insert(std::make_pair(e.t.utime, e));
	}
	std::cout << " * Finish input EventSummary, size : " << eve.size() << std::endl;

}

void AddEventInfo(std::map<int, EventSummary>& sf, std::multimap<int, EventInfo>& eve, int mflg) {

	for (auto itr = sf.begin(); itr != sf.end(); itr++) {

		auto p = eve.equal_range(itr->first);
		if (eve.count(itr->first) == 0) {
			//std::cout << "This event is already judged." << std::endl;
			continue;
		}
		double x = p.first->second.x, y = p.first->second.y;
		int pl = p.first->second.pl;
		for (auto itr0 = p.first; itr0 != p.second; itr0++) {
			if (itr0->second.pl > pl) {
				x = itr0->second.x;
				y = itr0->second.y;
				pl = itr0->second.pl;
			}
			if (eve.count(itr->first) > 1) {
				std::cout << itr0->second.eid << " " << itr0->second.x << " " << itr0->second.y << " " << itr0->second.pl << std::endl;
			}

		}

		//if (itr->second.pl > pl)continue;
		if (x < 40000.0001 || x>210000.0001 || y < 40000.0001 || y>210000.0001) {
			itr->second.areaflg = -1;
		}
		else {
			itr->second.areaflg = 1;
		}

		itr->second.muonflg = mflg;
		itr->second.pl = pl;
		if (pl > 130) {
			itr->second.muonflg = -1;
		}
	}

	std::cout << " * Finish Merge timestamp and event information" << std::endl;
}
void DiscardSandMuon(std::string input, std::multimap<int, EventSummary>& map, std::map<int, EventSummary>& out) {

	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File Open Error!  " << input << std::endl;
	}
	int ut;
	int cnt = 0;
	while (ifs >> ut) {
		if (map.count(ut) > 0) {
			auto p = map.equal_range(ut);
			for (auto itr = p.first; itr != p.second; itr++) {
				itr->second.muonflg = -1;
				cnt++;
			}
		}
	}

	std::cout << " * Finish ovewriting muon-classification-flg : # = " << cnt << std::endl;
	for (auto itr = map.begin(); itr != map.end(); itr++) {
		out.insert(std::make_pair(itr->second.t.eid, itr->second));
	}
	map.clear();
}
void Write(std::string output, std::map<int, EventSummary>& map) {

	std::ofstream ofs(output);
	if (!ofs) {
		std::cerr << "file open error!" << std::endl;
	}
	for (auto itr = map.begin(); itr != map.end(); itr++) {
		ofs << std::right << std::fixed
			<< std::setw(10) << std::setprecision(0) << itr->second.t.utime << " "
			<< std::setw(3) << std::setprecision(0) << itr->second.t.bunch << " "
			<< std::setw(3) << std::setprecision(0) << itr->second.t.ecc << " "
			<< std::setw(6) << std::setprecision(0) << itr->second.t.eid << " "
			<< std::setw(4) << std::setprecision(0) << itr->second.pl << " "
			<< std::setw(4) << std::setprecision(0) << itr->second.muonflg << " "
			<< std::setw(4) << std::setprecision(0) << itr->second.areaflg
			<< std::endl;
	}
	std::cout << " * Finish writting." << std::endl;
}

void Count(std::multimap<int, EventSummary>& eve, std::string output) {
	std::set<int> ut;
	for (auto itr = eve.begin(); itr != eve.end(); itr++) {
		ut.insert(itr->first);
	}

	std::ofstream ofs(output);
	int stop_c, stop_e, edge, pene, non, iron, fe;
	int st = 0, oth = 0, a = 0, b = 0, dis = 0;
	for (auto itr = ut.begin(); itr != ut.end(); itr++) {
		auto p = eve.equal_range(*itr);
		stop_c = 0; stop_e = 0; pene = 0; edge = 0; non = 0;
		iron = 0; fe = 0;
		for (auto itr0 = p.first; itr0 != p.second; itr0++) {
			if (itr0->second.muonflg > 0) {
				if (itr0->second.pl % 2 == 0 && itr0->second.pl > 15) {
					iron++;
				}
				else if (itr0->second.pl % 2 == 1 && itr0->second.pl > 15) {
					if (itr0->second.areaflg > 0) {
						stop_c++;
					}
					else {
						stop_e++;
					}
				}
				else if (itr0->second.pl < 16) {
					fe++;
				}
			}
			else if (itr0->second.muonflg == -1) {
				pene++;
			}
			else if (itr0->second.muonflg == -2) {
				edge++;
			}
			else {
				non++;
			}
		}

		ofs << std::right << std::fixed
			<< std::setw(10) << std::setprecision(0) << *itr << " "
			<< std::setw(5) << std::setprecision(0) << stop_c << " "
			<< std::setw(5) << std::setprecision(0) << stop_e << " "
			<< std::setw(5) << std::setprecision(0) << iron << " "
			<< std::setw(5) << std::setprecision(0) << fe << " "
			<< std::setw(5) << std::setprecision(0) << pene << " "
			<< std::setw(5) << std::setprecision(0) << edge << " "
			<< std::setw(5) << std::setprecision(0) << non
			<< std::endl;

		if (stop_c > 0 || stop_e > 0) {
			st++;
		}
		else if (iron > 0) {
			a++;
		}
		else if (fe > 0) {
			b++;
		}
		else if (pene > 0 || edge > 0) {
			oth++;
		}
		else {
			dis++;
		}
	}
	std::cout << "* timestamp based count\nStop Total: " << st + a + b
		<< "\n  Water = " << st << ", iron = " << a << ", FeECC = " << b
		<< "\npenetrate or edgeout : " << oth
		<< "\nevaporate : " << dis << std::endl;
}