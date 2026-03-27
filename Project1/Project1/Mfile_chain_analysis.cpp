// Mfile_chain_analysis
// output:
//   m_upstream.bmf
//	 m_edgeout
//	 m_short_thin.bmf
//	 m_stop.bmf
//	 m_start.bmf
#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
#include <set>

class Mfile_Area {
public:
	std::map<int, double> x_min, y_min, x_max, y_max, z;
};

uint64_t mfile_size(const mfile1::MFile_minimum& m_all);
Mfile_Area read_area(std::string filename);
std::vector<uint64_t> chain_id_list(const mfile1::MFile_minimum& m_all);

void divide_upstream_chain(const mfile1::MFile_minimum& m_all, const std::vector<uint64_t>& all, std::vector<uint64_t>& chain_up, std::vector<uint64_t>& chain_down, int PL_thr);
void divide_edge_out_chain(const mfile1::MFile_minimum& m_all, const std::vector<uint64_t>& all, std::vector<uint64_t>& chain_out, std::vector<uint64_t>& chain_in, const Mfile_Area& area);
bool judge_edgeout(const std::vector<mfile1::MFileBase>& base, const Mfile_Area& area);
void stop_track_selection(const mfile1::MFile_minimum& m_all, const std::vector<uint64_t>& all, std::vector<uint64_t>& remain, std::vector<uint64_t>& selected);
bool VPH_cut(const std::vector<mfile1::MFileBase>& base);
void mfile_wrtie(std::string filename, mfile1::MFile_minimum& m_all, const std::vector<uint64_t>& chain_id);

void divide_downstream_chain(const mfile1::MFile_minimum& m_all, const std::vector<uint64_t>& all, std::vector<uint64_t>& chain_up, std::vector<uint64_t>& chain_down, int PL_thr);
void divide_edge_out_chain_start(const mfile1::MFile_minimum& m_all, const std::vector<uint64_t>& all, std::vector<uint64_t>& chain_out, std::vector<uint64_t>& chain_in, const Mfile_Area& area);
bool judge_edgeout_start(const std::vector<mfile1::MFileBase>& base, const Mfile_Area& area);

void main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "usage:prg in-mfile area-file out-stop\n");
		exit(1);
	}
	std::string file_in_mfile = argv[1];
	std::string file_in_area = argv[2];
	Mfile_Area area = read_area(file_in_area);

	mfile1::MFile_minimum m;
	mfile1::read_mfile(file_in_mfile, m);

	printf("read mfile fin\n");
	printf("size:%.3lf[GB]\n", mfile_size(m) / 1000. / 1000. / 1000);

	//stoptrackの抽出
	int max_pl = area.z.rbegin()->first;
	printf("max PL = %d\n", max_pl);
	std::vector<uint64_t>chain_all, chain_up, chain_down;
	chain_all = chain_id_list(m);
	//max PL=133の場合 130,131,132,133をVETO ==> max_pl - 3
	//2021/10/31 PL132,133を除外-->max_pl -1
	divide_upstream_chain(m, chain_all, chain_up, chain_down, max_pl - 1);
	std::vector<uint64_t>chain_in, chain_out;
	divide_edge_out_chain(m, chain_down, chain_out, chain_in, area);
	std::vector<uint64_t>selected, remain;
	stop_track_selection(m, chain_in, remain, selected);

	//start trackの抽出
	int min_pl = area.z.begin()->first;
	printf("min PL = %d\n", min_pl);
	std::vector<uint64_t>schain_up, schain_down;
	chain_all = chain_id_list(m);
	//max PL=4の場合 4,5,6,7をVETO ==> min_pl + 3
	//2021/10/31 PL003~133でよい？->min_pl
	divide_downstream_chain(m, chain_all, schain_up, schain_down, min_pl + 0);
	std::vector<uint64_t>schain_in, schain_out;
	divide_edge_out_chain_start(m, schain_up, schain_out, schain_in, area);
	std::vector<uint64_t>sselected, sremain;
	stop_track_selection(m, schain_in, sremain, sselected);

	mfile_wrtie("m_upstream.bmf", m, chain_up);
	mfile_wrtie("m_edgeout.bmf", m, chain_out);
	mfile_wrtie("m_short_thin.bmf", m, remain);
	mfile_wrtie("m_stop.bmf", m, selected);
	mfile_wrtie("m_start.bmf", m, sselected);
}
uint64_t mfile_size(const mfile1::MFile_minimum& m_all) {
	uint64_t size = 0;
	size += sizeof(mfile1::MFileHeader);
	size += sizeof(mfile1::MFileInfoHeader);
	uint64_t chain_size = sizeof(mfile1::MFileChain);
	uint64_t base_size = sizeof(mfile1::MFileBase);
	size += chain_size * m_all.info_header.Nchain;
	size += base_size * m_all.info_header.Nbasetrack;
	return size;
}
Mfile_Area read_area(std::string filename) {
	std::ifstream ifs(filename);
	Mfile_Area ret;
	int pl;
	double x_min, y_min, x_max, y_max, z;
	while (ifs >> pl >> x_min >> x_max >> y_min >> y_max >> z) {
		ret.x_min.insert(std::make_pair(pl, x_min));
		ret.x_max.insert(std::make_pair(pl, x_max));
		ret.y_min.insert(std::make_pair(pl, y_min));
		ret.y_max.insert(std::make_pair(pl, y_max));
		ret.z.insert(std::make_pair(pl, z));
	}
	return ret;
}
std::vector<uint64_t> chain_id_list(const mfile1::MFile_minimum& m_all) {
	std::vector<uint64_t> ret;
	for (int i = 0; i < m_all.chains.size(); i++) {
		ret.push_back(i);
	}
	return ret;
}
void divide_upstream_chain(const mfile1::MFile_minimum& m_all, const std::vector<uint64_t>& all, std::vector<uint64_t>& chain_up, std::vector<uint64_t>& chain_down, int PL_thr) {
	for (int i = 0; i < all.size(); i++) {
		//PL_thr以上の飛跡をup stream VETO
		if (m_all.chains[all[i]].pos1 / 10 >= PL_thr) {
			chain_up.push_back(all[i]); // pl133,132
		}
		else {
			chain_down.push_back(all[i]);
		}
	}
	printf("all chain %lld\n", all.size());
	printf("\tup  stream track %lld\n", chain_up.size());
	printf("\tdownstream track %lld\n", chain_down.size());
	return;
}
void divide_edge_out_chain(const mfile1::MFile_minimum& m_all, const std::vector<uint64_t>& all, std::vector<uint64_t>& chain_out, std::vector<uint64_t>& chain_in, const Mfile_Area& area) {
	//3PL外挿して1つでもout of range ならedege out
	for (int i = 0; i < all.size(); i++) {
		if (judge_edgeout(m_all.all_basetracks[all[i]], area)) {
			chain_out.push_back(all[i]);
		}
		else {
			chain_in.push_back(all[i]);
		}
	}
	printf("all chain %lld\n", all.size());
	printf("\tedge out track %lld\n", chain_out.size());
	printf("\tstop     track %lld\n", chain_in.size());
	return;
}
bool judge_edgeout(const std::vector<mfile1::MFileBase>& base, const Mfile_Area& area) {
	int pl = base.rbegin()->pos / 10;
	double z = area.z.at(pl);
	double edge_cut = 5000;//minimumから5mm内側のlineでカット
	bool flg = false;
	double x_up, y_up;
	//1-3PL分外挿 上流に向かって
	for (int i_pl = 1; i_pl <= 3; i_pl++) {
		auto res = area.z.find(pl + i_pl);
		if (res == area.z.end()) {
			fprintf(stderr, "pl out of range %d\n", pl + i_pl);
			continue;
		}
		x_up = base.rbegin()->x + base.rbegin()->ax * (res->second - z);
		y_up = base.rbegin()->y + base.rbegin()->ay * (res->second - z);
		if (area.x_min.at(pl + i_pl) + edge_cut > x_up)flg = true;
		if (area.x_max.at(pl + i_pl) - edge_cut < x_up)flg = true;
		if (area.y_min.at(pl + i_pl) + edge_cut > y_up)flg = true;
		if (area.y_max.at(pl + i_pl) - edge_cut < y_up)flg = true;
	}
	return flg;
}
void stop_track_selection(const mfile1::MFile_minimum& m_all, const std::vector<uint64_t>& all, std::vector<uint64_t>& remain, std::vector<uint64_t>& selected) {

	std::vector<mfile0::M_Chain> ret;

	bool flg = false;
	int nPL = 0;
	double eff;
	for (auto& i : all) {
		flg = false;
		nPL = (m_all.chains[i].pos1 - m_all.chains[i].pos0) / 10 + 1;
		eff = m_all.chains[i].nseg * 1.0 / nPL;
		if (eff < 0.500000000000000001)continue;
		if (m_all.chains[i].nseg <= 5) {
			flg = VPH_cut(m_all.all_basetracks[i]);
		}
		else {
			flg = true;
		}
		if (flg) {
			selected.push_back(i);
		}
		else {
			//PH cutはかけといて，PartnerSerchで拾う,this is better way than doing no cut??<==because there are too many tracks,
			selected.push_back(i); //ph_cutかけない?
			remain.push_back(i); // 一応分離しておく
		}
	}
	printf("all chain %lld\n", all.size());
	printf("\tselected     chain %lld\n", selected.size());
	printf("\tnot selected chain %lld\n", remain.size());
	return;
}
bool VPH_cut(const std::vector<mfile1::MFileBase>& base) {
	double ax = 0, ay = 0, VPH = 0;
	int count = 0;
	for (const auto& b : base) {
		ax += b.ax;
		ay += b.ay;
		VPH += b.ph % 10000;
		count++;
	}
	ax = ax / count;
	ay = ay / count;
	VPH = VPH / count;
	double angle = sqrt(ax * ax + ay * ay);

	if (angle < 1.0) {
		return VPH > 190 - 70 * angle;
	}
	else if (1.0 <= angle && angle < 2.0) {
		return VPH > 150 - 30 * angle;
	}
	else if (2.0 <= angle && angle < 3.0) {
		return VPH > 110 - 10 * angle;
	}
	else if (3.0 <= angle) {
		return VPH > 95 - 5 * angle;
	}
	return false;
}

void divide_downstream_chain(const mfile1::MFile_minimum& m_all, const std::vector<uint64_t>& all, std::vector<uint64_t>& chain_up, std::vector<uint64_t>& chain_down, int PL_thr) {
	for (int i = 0; i < all.size(); i++) {
		//PL_thr以下の飛跡をdown stream VETO
		if (m_all.chains[all[i]].pos0 / 10 <= PL_thr) {
			chain_down.push_back(all[i]);
		}
		else {
			chain_up.push_back(all[i]);
		}
	}
	printf("all chain %lld\n", all.size());
	printf("\tup  stream track %lld\n", chain_up.size());
	printf("\tdownstream track %lld\n", chain_down.size());
	return;
}
void divide_edge_out_chain_start(const mfile1::MFile_minimum& m_all, const std::vector<uint64_t>& all, std::vector<uint64_t>& chain_out, std::vector<uint64_t>& chain_in, const Mfile_Area& area) {
	//3PL外挿して1つでもout of range ならedege out 下流に向かって
	for (int i = 0; i < all.size(); i++) {
		if (judge_edgeout_start(m_all.all_basetracks[all[i]], area)) {
			chain_out.push_back(all[i]);
		}
		else {
			chain_in.push_back(all[i]);
		}
	}
	printf("all chain %lld\n", all.size());
	printf("\tedge out track %lld\n", chain_out.size());
	printf("\tstart     track %lld\n", chain_in.size());
	return;
}
bool judge_edgeout_start(const std::vector<mfile1::MFileBase>& base, const Mfile_Area& area) {
	int pl = base.begin()->pos / 10;
	double z = area.z.at(pl);
	double edge_cut = 5000;// edgeよりも5mm内側
	bool flg = false;
	double x_down, y_down;
	//1-3PL分外挿(下流へ)
	for (int i_pl = 1; i_pl <= 3; i_pl++) {
		auto res = area.z.find(pl - i_pl);
		if (res == area.z.end()) {
			fprintf(stderr, "pl out of range %d\n", pl - i_pl);
			continue;
		}
		x_down = base.begin()->x + base.begin()->ax * (res->second - z);
		y_down = base.begin()->y + base.begin()->ay * (res->second - z);
		if (area.x_min.at(pl - i_pl) + edge_cut > x_down)flg = true;
		if (area.x_max.at(pl - i_pl) - edge_cut < x_down)flg = true;
		if (area.y_min.at(pl - i_pl) + edge_cut > y_down)flg = true;
		if (area.y_max.at(pl - i_pl) - edge_cut < y_down)flg = true;
	}
	return flg;
}



void mfile_wrtie(std::string filename, mfile1::MFile_minimum& m_all, const std::vector<uint64_t>& chain_id) {

	std::ofstream ofs(filename, std::ios::binary);

	m_all.info_header.Nchain = chain_id.size();

	size_t Nbasetrack = 0;
	for (auto& p : chain_id) {
		Nbasetrack += m_all.all_basetracks[p].size();
	}
	m_all.info_header.Nbasetrack = Nbasetrack;

	ofs.write((char*)&m_all.header, sizeof(mfile1::MFileHeader));
	ofs.write((char*)&m_all.info_header, sizeof(mfile1::MFileInfoHeader));

	int count = 0;
	int64_t max = m_all.info_header.Nchain;

	for (auto& p : chain_id) {
		if (count % 10000 == 0) {
			std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%";
		}
		count++;
		ofs.write((char*)&m_all.chains[p], sizeof(mfile1::MFileChain));
		assert(m_all.chains[p].nseg == m_all.all_basetracks[p].size());
		for (int b = 0; b < m_all.chains[p].nseg; b++) {
			ofs.write((char*)&m_all.all_basetracks[p][b], sizeof(mfile1::MFileBase));
		}

	}

	std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%" << std::endl;

}

