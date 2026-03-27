// 2025/12/08 
// kasumi

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}


void add_momch(std::vector<Momentum_recon::Event_information>& momch, std::vector<Momentum_recon::Event_information>& momch_add, int add_event, int add_chain);
void read_list(std::string in, std::set<Lst>& list);
void erase_momch(std::vector<Momentum_recon::Event_information>& momch, std::set<Lst>& set);

int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "usage: in.momch remain_chain_list.txt out.momch\nerase_list\n\t[group-id] [chain-id]\n");
		exit(1);
	}
	std::string file_in_momch_ori = argv[1];//input momch
	std::string file_in_momch_add = argv[2];//erase list
	//int add_event = std::stoi(argv[3]);
	//int add_chain = std::stoi(argv[4]);
	std::string file_out_momch = argv[3];//output

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch_ori);
	//std::vector<Momentum_recon::Event_information> momch_add = Momentum_recon::Read_Event_information_extension(file_in_momch_add);

	//add_momch(momch, momch_add, add_event, add_chain);
	std::set<Lst>list;
	read_list(file_in_momch_add, list);
	erase_momch(momch, list);

	Momentum_recon::Write_Event_information_extension(file_out_momch, momch);
}
void add_momch(std::vector<Momentum_recon::Event_information>& momch, std::vector<Momentum_recon::Event_information>& momch_add, int add_event, int add_chain) {
	Momentum_recon::Mom_chain chain;
	for (auto& ev : momch_add) {
		if (ev.groupid != add_event)continue;
		for (auto& c : ev.chains) {
			if (c.chainid != add_chain)continue;
			chain = c;
		}
	}

	for (auto& ev : momch) {
		if (ev.groupid != add_event)continue;
		ev.chains.push_back(chain);
	}

}
void read_list(std::string in, std::set<Lst>& list) {
	std::ifstream ifs(in);
	if (!ifs) {
		std::cerr << "reading file error" << std::endl;
	}

	Lst l;
	while (ifs >> l.gid >> l.cid) {
		//	while (ifs >> l.gid >> l.cid >> l.pl >> l.rawid) {
		list.insert(l);
	}
	std::cout << "Reamin Track size:" << list.size() << std::endl;
}
void erase_momch(std::vector<Momentum_recon::Event_information>& momch, std::set<Lst>& set) {
	Momentum_recon::Mom_chain chains;

	std::ofstream ofs("erased_chain_log.txt");
	int count = 0;
	int remain = 0;
	Lst k;

	for (auto& ev : momch) {
		for (auto itr = ev.chains.begin(); itr != ev.chains.end();) {
			k.gid = ev.groupid;
			k.cid = itr->chainid;
			int num = set.count(k);
			//std::cout << k.gid << " " << k.cid << " " << ev.chains.size() << std::endl;
			// 条件一致した要素を削除する
			if (num == 0) {
				// 削除された要素の次を指すイテレータが返される。
				itr = ev.chains.erase(itr);
				ofs << k.gid << " " << k.cid << std::endl;
				count++;
			}
			// 要素削除をしない場合に、イテレータを進める
			else {
				remain++;
				++itr;
			}
		}
	}
	std::cout << "Remain chain num : " << remain << std::endl;
	std::cout << "Erase  chain num : " << count << std::endl;
}