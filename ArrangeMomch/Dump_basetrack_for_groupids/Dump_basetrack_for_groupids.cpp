// 2026/11/10
// kasumi
// Dump_basetrack_for_groupids
// Given a momch file and a list of groupids, dump the basetrack points
// (pl, rawid, x, y, z, ax, ay) of the muon chain (particle_flg==13, chainid==0)
// for each requested groupid, plus the event's vertex_material, so that side
// view / top view plots can be drawn from outside this program (e.g. matplotlib).

#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
#include <iostream>
#include <fstream>
#include <set>
#include <string>

int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "===============================================================================\n");
		fprintf(stderr, " usage : prg in.momch groupid_list.txt out.txt\n\n");
		fprintf(stderr, " groupid_list.txt : one groupid per line\n");
		fprintf(stderr, " out.txt columns  : groupid, pl, rawid, x, y, z, ax, ay, vertex_material\n");
		fprintf(stderr, "                    (muon chain : particle_flg==13, chainid==0)\n");
		fprintf(stderr, "===============================================================================\n");
		exit(1);
	}
	std::string in_momch = argv[1];
	std::string groupid_list_file = argv[2];
	std::string out_txt = argv[3];

	std::set<int> target_groupids;
	{
		std::ifstream ifs(groupid_list_file);
		if (!ifs) {
			std::cerr << "Failed to open groupid list file : " << groupid_list_file << std::endl;
			exit(1);
		}
		int g;
		while (ifs >> g) {
			target_groupids.insert(g);
		}
	}
	std::cout << " # of target groupids : " << target_groupids.size() << std::endl;

	std::cout << " Reading momch" << std::endl;
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);

	std::ofstream ofs(out_txt);
	if (!ofs) {
		std::cerr << "Failed to create the output file : " << out_txt << std::endl;
		exit(1);
	}

	int found_groupids = 0;
	int basetrack_num = 0;
	for (auto& ev : momch) {
		if (target_groupids.find(ev.groupid) == target_groupids.end()) continue;

		for (auto& c : ev.chains) {
			if (c.particle_flg != 13) continue;
			if (c.chainid != 0) continue;

			found_groupids++;
			for (auto& b : c.base) {
				ofs << ev.groupid << "\t" << b.pl << "\t" << b.rawid << "\t"
					<< b.x << "\t" << b.y << "\t" << b.z << "\t" << b.ax << "\t" << b.ay << "\t"
					<< ev.vertex_material << std::endl;
				basetrack_num++;
			}
		}
	}

	std::cout << " # of groupids found in momch (muon chain) : " << found_groupids << std::endl;
	std::cout << " # of basetracks written                   : " << basetrack_num << std::endl;
	std::cout << "Fin." << std::endl;
}
