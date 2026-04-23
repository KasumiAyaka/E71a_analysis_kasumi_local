// 2026/4/23
// kasumi
// Check_eventnum_and_partnernum
// To display the number of interactions for each substance and the number for each type of partner.

#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
# include <iostream>
# include <random>

#include <list>
#include <cassert>
#include <filesystem>

void MeasureProcessingTime(std::chrono::system_clock::time_point& start, std::chrono::system_clock::time_point& end) {
	auto dur = end - start;        // 要した時間を計算
	auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
	// 要した時間をミリ秒（1/1000秒）に変換して表示
	if (msec < 1000) {
		std::cout << msec << " milli sec \n";
	}
	else if (msec / 1000 < 60) {
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

struct Materials {
	int p, pi, mu, oth;
	int pid;
};

void count_events(std::vector<Momentum_recon::Event_information>& momch, std::multimap<int, Materials>& mat);
void output(std::string output, std::multimap<int, Materials>& mat);


int main(int argc, char** argv) {
	if (argc < 2 || argc>4) {
		fprintf(stderr, "===============================================================================\n");
		fprintf(stderr, " usage : prg in.momch [output.txt]\n\n");
		fprintf(stderr, " If you want to output the number of particles for each event, please enter the output file name as an argument.\n The output file format is :\n group ID, vertex material, number of muons, number of pions, number of protons, and number of others.\n");
		fprintf(stderr, "===============================================================================\n");
		exit(1);
	}
	std::string in_momch = argv[1];// input momch
	std::string out_txt;

	//bool result = std::filesystem::create_directories(out_momch);
	//assert(result);
	//assert(std::filesystem::exists(out_momch));//directryの存在確認
	//assert(std::filesystem::is_directory(out_momch));//指定されたパスがディレクトリを指しているかを確認する

	auto start = std::chrono::system_clock::now();//for measure working time


	// read momch
	std::cout << " Reading momch" << std::endl;
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);
	std::multimap<int, Materials>mat;
	count_events(momch, mat);

	if (argc > 2) {
		output(argv[2], mat);
	}

	// write out
	//Momentum_recon::Write_Event_information_extension(out_momch, momch);

	auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
	MeasureProcessingTime(start, end);
}


void count_events(std::vector<Momentum_recon::Event_information>& momch, std::multimap<int, Materials>& mat) {

	int p, pi, mu, oth;
	int wt, fe, bs, emu, pene, ed, dup;
	Materials m;
	int cnt = 0;
	mu = 0;
	wt = bs = fe = emu = pene = ed = dup = 0;
	std::cout << "\n Counting the number of events, partner tracks." << std::endl;
	for (auto& ev : momch) {
		if (ev.vertex_material == 0) {//water
			wt++;
		}
		else if (ev.vertex_material == 1) {//base
			bs++;
		}
		else if (ev.vertex_material == 2) {//iron
			fe++;
		}
		else if (ev.vertex_material == 5) {//emulsion
			emu++;
		}
		else if (ev.vertex_material == -2) {//pene
			pene++;
		}
		else if (ev.vertex_material == -3) {//sideout
			ed++;
		}
		else if (ev.vertex_material == -5) {// timestamp duplicated
			dup++;
		}
		else {
			std::cout << "unecpected material..." << std::endl;
		}
		mu = 0;
		pi = 0;
		p = 0;
		oth = 0;

		for (auto& c : ev.chains) {
			if (c.particle_flg % 10000 == 13) {//muon
				mu++;
			}
			else if (c.particle_flg % 10000 == 211) {//pion
				pi++;
			}
			else if (c.particle_flg % 10000 == 2212) {//proton
				p++;
			}
			else {// cannot distinguish
				oth++;
			}
		}

		m.mu = mu;
		m.pi = pi;
		m.p = p;
		m.oth = oth;
		m.pid = ev.groupid;
		mat.insert(std::make_pair(ev.vertex_material, m));
		cnt++;
	}

	std::cout << "\n *** Event Summary ***" << std::endl;
	std::cout << "  # of events = " << wt + bs + fe + emu + pene + ed + dup << " /" << cnt << std::endl;
	std::cout << "  Material : ( water, base, Iron, emulsion, penetrate, sideout, duplicated ) \n           = ( "
		<< std::setw(5) << wt << ", " << std::setw(4) << bs << ", " << std::setw(4) << fe << ", " << std::setw(8) << emu << ", "
		<< std::setw(9) << pene << ", " << std::setw(7) << ed << ", " << std::setw(10) << dup << " )\n" << std::endl;

	mu = 0;
	pi = 0;
	p = 0;
	oth = 0;
	for (auto itr = mat.equal_range(0).first; itr != mat.equal_range(0).second; itr++) {
		mu = itr->second.mu + mu;
		p = itr->second.p + p;
		pi = itr->second.pi + pi;
		oth = itr->second.oth + oth;
	}
	std::cout << "  water    : ( mu, pi, p, oth ) = ( " << std::setw(3) << mu << ", " << std::setw(3) << pi << ", " << std::setw(3) << p << ", " << std::setw(3) << oth << " )" << std::endl;
	mu = 0;
	pi = 0;
	p = 0;
	oth = 0;
	for (auto itr = mat.equal_range(1).first; itr != mat.equal_range(1).second; itr++) {
		mu = itr->second.mu + mu;
		p = itr->second.p + p;
		pi = itr->second.pi + pi;
		oth = itr->second.oth + oth;
	}
	std::cout << "  base     : ( mu, pi, p, oth ) = ( " << std::setw(3) << mu << ", " << std::setw(3) << pi << ", " << std::setw(3) << p << ", " << std::setw(3) << oth << " )" << std::endl;
	mu = 0;
	pi = 0;
	p = 0;
	oth = 0;
	for (auto itr = mat.equal_range(2).first; itr != mat.equal_range(2).second; itr++) {
		mu = itr->second.mu + mu;
		p = itr->second.p + p;
		pi = itr->second.pi + pi;
		oth = itr->second.oth + oth;
	}
	std::cout << "  Iron     : ( mu, pi, p, oth ) = ( " << std::setw(3) << mu << ", " << std::setw(3) << pi << ", " << std::setw(3) << p << ", " << std::setw(3) << oth << " )" << std::endl;
	mu = 0;
	pi = 0;
	p = 0;
	oth = 0;
	for (auto itr = mat.equal_range(5).first; itr != mat.equal_range(5).second; itr++) {
		mu = itr->second.mu + mu;
		p = itr->second.p + p;
		pi = itr->second.pi + pi;
		oth = itr->second.oth + oth;
	}
	std::cout << "  emulsion : ( mu, pi, p, oth ) = ( " << std::setw(3) << mu << ", " << std::setw(3) << pi << ", " << std::setw(3) << p << ", " << std::setw(3) << oth << " )\n" << std::endl;

}
void output(std::string output, std::multimap<int, Materials>& mat) {
	std::ofstream ofs(output);
	if (!ofs) {
		std::cerr << "\nFailed to create the output file. :" <<output<< std::endl;
		return;
	}

	for (auto itr = mat.begin(); itr != mat.end(); itr++) {
		ofs << std::right << std::fixed
			<< std::setw(5) << std::setprecision(0) << itr->second.pid << " "// groupid
			<< std::setw(5) << std::setprecision(0) << itr->first << " "// material
			<< std::setw(5) << std::setprecision(0) << itr->second.mu << " "//munum
			<< std::setw(5) << std::setprecision(0) << itr->second.pi << " "//pinum
			<< std::setw(5) << std::setprecision(0) << itr->second.p << " "//pnum
			<< std::setw(5) << std::setprecision(0) << itr->second.oth << " "//oth num
			<< std::endl;
	}
	std::cout << "Writing to the output file is complete.\n" << std::endl;
}