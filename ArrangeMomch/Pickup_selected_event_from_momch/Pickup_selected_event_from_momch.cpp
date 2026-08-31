// groupid ���L�q���ꂽ txt ��ǂݍ��݁A�Y������C�x���g�� pid(particle_flg)==13 �� chain �����o�͂���
// kasumi

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#include <set>

void read_groupid_list(std::string filename, std::set<int>& list);
void select_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, std::set<int>& groupid_list);

int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "usage:prg in.momch groupid_list.txt out.momch\n");
		exit(1);
	}
	std::string file_in_momch = argv[1];//input momch
	std::string file_in_list = argv[2];//groupid list (1 groupid per line)
	std::string file_out_momch = argv[3];//output

	std::vector<Momentum_recon::Event_information> momch0 = Momentum_recon::Read_Event_information_extension(file_in_momch);
	std::vector<Momentum_recon::Event_information> momch;

	std::set<int> groupid_list;
	read_groupid_list(file_in_list, groupid_list);

	select_momch(momch0, momch, groupid_list);

	Momentum_recon::Write_Event_information_extension(file_out_momch, momch);
}

void read_groupid_list(std::string filename, std::set<int>& list) {
	std::ifstream ifs(filename);
	if (!ifs) {
		std::cerr << filename << " could not be opened." << std::endl;
		exit(1);
	}

	int gid;
	while (ifs >> gid) {
		list.insert(gid);
	}
	std::cout << "Groupid list size : " << list.size() << std::endl;
}

void select_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, std::set<int>& groupid_list) {
	int event_count = 0;
	int remain_chain = 0;
	int erase_chain = 0;

	for (auto& ev : momch0) {
		if (groupid_list.count(ev.groupid) == 0)continue;

		for (auto c = ev.chains.begin(); c != ev.chains.end();) {
			if (c->particle_flg != 13) {
				c = ev.chains.erase(c);
				erase_chain++;
			}
			else {
				c++;
				remain_chain++;
			}
		}
		momch.push_back(ev);
		event_count++;
	}

	std::cout << "Selected event num : " << event_count << std::endl;
	std::cout << "Remain chain num   : " << remain_chain << std::endl;
	std::cout << "Erase  chain num   : " << erase_chain << std::endl;
}
