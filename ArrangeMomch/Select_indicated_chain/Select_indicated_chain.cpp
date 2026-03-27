// 2025/11/04
// kasumi

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}

void Select_length(std::vector<Momentum_recon::Event_information>& momch, int npl_thr);
void Select_volume(std::vector<Momentum_recon::Event_information>& momch, int vph_thr);
void Select_angle(std::vector<Momentum_recon::Event_information>& momch, double tan_thr);
void Declease_event_num(std::vector<Momentum_recon::Event_information>& momch, int n);

int main(int argc, char** argv) {
	if (argc < 5) {
		fprintf(stderr, "usage: in.momch  <npl_thr <=vph >tan_thr n(divide num) out.momch");
		exit(1);
	}
	std::string file_in_momch_ori = argv[1];//input momch
	int npl_thr = std::stoi(argv[2]);
	int vph_thr = std::stoi(argv[3]);
	double tan_thr = std::stod(argv[4]);
	int n = std::stoi(argv[5]);
	if (n < 2) {
		n = 1;
	}

	std::string file_out_momch = file_in_momch_ori;
	if (argc == 7) {
		file_out_momch = argv[6];
	}

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch_ori);
	Select_length(momch, npl_thr);
	Select_volume(momch, vph_thr);
	Select_angle(momch, tan_thr);
	if (n > 1) {
		Declease_event_num(momch, n);
	}
	Momentum_recon::Write_Event_information_extension(file_out_momch, momch);
}
void Select_length(std::vector<Momentum_recon::Event_information>& momch, int npl_thr) {
	Momentum_recon::Mom_chain chain;
	//std::cout << " * cut npl >= " << npl_thr << std::endl;
	//for (auto& ev : momch) {
	//	for (auto itr = ev.chains.begin(); itr != ev.chains.end();) {
	//		if (itr->particle_flg == 13)continue;

	//		if (itr->base.rbegin()->pl - itr->base.begin()->pl + 1 < npl_thr) {
	//			//std::cout << "remain : npl = " << itr->base.rbegin()->pl - itr->base.begin()->pl + 1 << std::endl;
	//			// 要素削除をしない場合に、イテレータを進める
	//			++itr;

	//		}
	//		else {
	//			//std::cout << "erase : npl = " << itr->base.rbegin()->pl - itr->base.begin()->pl + 1 << std::endl;
	//			// 条件一致した要素を削除する
	//			// 削除された要素の次を指すイテレータが返される。
	//			itr = ev.chains.erase(itr);
	//		}
	//	}
	//}
	std::cout << " * cut npl < " << npl_thr << std::endl;
	for (auto& ev : momch) {
		for (auto itr = ev.chains.begin(); itr != ev.chains.end();) {
			if (itr->particle_flg == 13) {
				itr++;
				continue;
			}

			if (itr->base.rbegin()->pl - itr->base.begin()->pl + 1 >= npl_thr) {
				//std::cout << "remain : npl = " << itr->base.rbegin()->pl - itr->base.begin()->pl + 1 << std::endl;
				// 要素削除をしない場合に、イテレータを進める
				++itr;

			}
			else {
				//std::cout << "erase : npl = " << itr->base.rbegin()->pl - itr->base.begin()->pl + 1 << std::endl;
				// 条件一致した要素を削除する
				// 削除された要素の次を指すイテレータが返される。
				itr = ev.chains.erase(itr);

			}
		}
	}

}
void Select_volume(std::vector<Momentum_recon::Event_information>& momch, int npl_thr) {
	Momentum_recon::Mom_chain chain;
	int flg,vph;
	std::cout << " * cut VPH <= " << npl_thr << std::endl;

	for (auto& ev : momch) {
		for (auto itr = ev.chains.begin(); itr != ev.chains.end();) {
			if (itr->particle_flg == 13) {
				itr++;
				continue;
			}
			if (abs(itr->base.rbegin()->pl - ev.vertex_pl) > abs(itr->base.begin()->pl - ev.vertex_pl)) {
				vph += itr->base.begin()->m[0].ph % 10000;
				vph += itr->base.begin()->m[1].ph % 10000;
			}
			else {
				vph += itr->base.rbegin()->m[0].ph % 10000;
				vph += itr->base.rbegin()->m[1].ph % 10000;
			}

			if (vph > npl_thr) {
				//std::cout << "remain : vph = " << vph << std::endl;
				// 要素削除をしない場合に、イテレータを進める
				++itr;

			}
			else {
				// 条件一致した要素を削除する
				// 削除された要素の次を指すイテレータが返される。
				itr = ev.chains.erase(itr);

			}
		}
	}

}
void Select_angle(std::vector<Momentum_recon::Event_information>& momch, double tan_thr) {
	Momentum_recon::Mom_chain chain;
	double tan;
	std::cout << " * cut tan > " << tan_thr << std::endl;

	for (auto& ev : momch) {
		for (auto itr = ev.chains.begin(); itr != ev.chains.end();) {
			if (itr->particle_flg == 13) {
				itr++;
				continue;
			}
			if (abs(itr->base.rbegin()->pl - ev.vertex_pl) > abs(itr->base.begin()->pl - ev.vertex_pl)) {
				tan =sqrt( itr->base.begin()->ax* itr->base.begin()->ax+ itr->base.begin()->ay* itr->base.begin()->ay);
			}
			else {
				tan = sqrt(itr->base.rbegin()->ax * itr->base.rbegin()->ax + itr->base.rbegin()->ay * itr->base.rbegin()->ay);
			}

			if (tan <= tan_thr) {
				//std::cout << "remain : tan = " << tan << std::endl;
				// 要素削除をしない場合に、イテレータを進める
				++itr;

			}
			else {
				// 条件一致した要素を削除する
				// 削除された要素の次を指すイテレータが返される。
				itr = ev.chains.erase(itr);

			}
		}
	}

}
void Declease_event_num(std::vector<Momentum_recon::Event_information>& momch, int n) {
	Momentum_recon::Mom_chain chain;
	int cnt=0;
	std::cout << " * cut  " << n << std::endl;

	for (auto& ev : momch) {
		for (auto itr = ev.chains.begin(); itr != ev.chains.end();) {
			cnt++;
			if (itr->particle_flg == 13) {
				itr++;
				continue;
			}
			if (cnt % n == 0) {
				//std::cout << "remain : cnt % n = " << cnt % n << std::endl;
				// 要素削除をしない場合に、イテレータを進める
				++itr;
			}
			else {
				// 条件一致した要素を削除する
				// 削除された要素の次を指すイテレータが返される。
				itr = ev.chains.erase(itr);

			}
		}
	}

}

