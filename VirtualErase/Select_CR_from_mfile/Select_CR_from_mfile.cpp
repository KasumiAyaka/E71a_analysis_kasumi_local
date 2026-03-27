// 2025/11/23
// Select_CR_from_mfile
// kasumi
// to estimate the amount of CR

#define _CRT_SECURE_NO_WARNINGS
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#include <set>
//for read jsonfile
//#include <boost/property_tree/ptree.hpp>
//#include <boost/property_tree/json_parser.hpp>
//#include <boost/foreach.hpp>
//#include <boost/optional.hpp>
//#include <boost/array.hpp>
//mfileÇÃç≈è„ó¨PLÇ∆ç≈â∫ó¨PLÇèoóÕÇ∑ÇÈ
struct PLRange {
	int max, min;
};
struct Margin {
	double x, y;
};
struct AngleCut {
	double x, y;
};

struct Param {
	int nseg;
	double fillfactor;
	double dlateral;
	PLRange rpl;
	Margin edgecut;
	AngleCut ang;
};


//using namespace boost::property_tree;
//void read_param(std::string param, Param& par);

struct Prediction {
	int PL, flg, rawid;
	double ax, ay, x, y;
	double vph;
};
struct Output {
	Prediction pre, hit;
};
struct Info {
	double ax, ay, x, y;
	int ph;
};
std::vector<mfile0::M_Chain> chain_nseg_selection(std::vector<mfile0::M_Chain>& chain, int nseg);
std::vector<mfile0::M_Chain> chain_dlat_selection(std::vector<mfile0::M_Chain>& chain, double threshold);
std::vector<mfile0::M_Chain> group_clustering(std::vector<mfile0::M_Chain>& chain);
std::vector<mfile0::M_Chain> chain_dlat_selection(std::vector<mfile0::M_Chain>& chain, double threshold);
std::vector<mfile0::M_Chain> chain_angle_selection(std::vector<mfile0::M_Chain>& chain, double thr_ax, double thr_ay);
//void CalcEfficiency(std::ofstream& ofs_eff, std::ofstream& ofs_pred, std::ofstream& ofs_id, std::vector<mfile0::M_Chain>  chain, int PL, std::map<int, double> zmap, Param& p);
void output_chain(std::ofstream& ofs, std::vector<mfile0::M_Chain>& chain);
void output_chain2(std::ofstream& ofs, std::vector<mfile0::M_Chain>& chain, int npl_thr);
int angle_divide(double angle);
std::vector<mfile0::M_Chain> chain_VPH_cut(std::vector<mfile0::M_Chain>& chain);
std::vector<mfile0::M_Chain> chain_fillfactor(std::vector<mfile0::M_Chain>& chain, double nseg_thr);


int main(int argc, char** argv) {
	if (argc < 3 || argc>5) {
		fprintf(stderr, "usage:prg in-mfile out-file\n");
		fprintf(stderr, "usage:prg in-mfile out-file npl_thr\n");
		//fprintf(stderr, "usage:prg in-mfile in_param out-file(eff) out-file(prediction) out_exist_base_id_list\n");
		exit(1);
	}
	//std::cout << "***Causion!***\nEdgecut is not avalable.\n*************" << std::endl;
	int npl_thr = -1;
	//input value
	std::string file_in_mfile = argv[1];
	//std::string file_in_param = argv[2];
	std::string file_out = argv[2];
	if (argc > 3) {
		npl_thr = std::stoi(argv[3]);
	}
	//mfile
	mfile0::Mfile m;

	//read param
	//Param p;
	//read_param(file_in_param, p);

	//read mfile
	mfile1::read_mfile_extension(file_in_mfile, m);

	//zmap
	//std::map<int, double> zmap;
	//for (auto itr = m.chains.begin(); itr != m.chains.end(); itr++) {
	//	for (auto itr2 = itr->basetracks.begin(); itr2 != itr->basetracks.end(); itr2++) {
	//		zmap.insert(std::make_pair(itr2->pos / 10, itr2->z));
	//	}
	//}


	//m.chains = chain_nseg_selection(m.chains, p.nseg);//(m.chains, 15)
	//m.chains = chain_angle_selection(m.chains, p.ang.x, p.ang.y);//(m.chains, 4.0, 4.0);
	//m.chains = chain_dlat_selection(m.chains, p.dlateral);//(m.chains, 0.005);
	//m.chains = chain_fillfactor(m.chains, p.fillfactor);

	//m.chains = group_clustering(m.chains);

	std::ofstream ofs_eff(file_out);

	if (npl_thr > 0) {
		output_chain2(ofs_eff, m.chains,npl_thr);
	}
	else {
		output_chain(ofs_eff, m.chains);
	}

	std::cout << "finished!!" << std::endl;
}
//void read_param(std::string param, Param& par) {
//	std::cout << "Set Param" << std::endl;
//	par = { 0 };
//
//	ptree pt;
//	read_json(param, pt);
//
//
//
//
//	// PL range
//	if (boost::optional<int> upl = pt.get_optional<int>("PLrange.Max")) {
//		par.rpl.max = upl.get();
//	}
//	else {
//		std::cout << "maxPL is nothing" << std::endl;
//	}
//
//	if (boost::optional<int> dpl = pt.get_optional<int>("PLrange.min")) {
//		par.rpl.min = dpl.get();
//	}
//	else {
//		std::cout << "minPL is nothing" << std::endl;
//	}
//
//
//	// nseg
//	if (boost::optional<int> upl = pt.get_optional<int>("Nseg")) {
//		par.nseg = upl.get();
//	}
//
//	// fill factor
//	if (boost::optional<double> upl = pt.get_optional<double>("FillFactor")) {
//		par.fillfactor = upl.get();
//	}
//
//
//	// dlateral
//	if (boost::optional<double> upl = pt.get_optional<double>("LateralAngleDiff")) {
//		par.dlateral = upl.get();
//	}
//
//
//	// margin
//	if (boost::optional<double> ax = pt.get_optional<double>("FilmEdgeCut.x")) {
//		par.edgecut.x = ax.get();
//	}
//	if (boost::optional<double> ay = pt.get_optional<double>("FilmEdgeCut.y")) {
//		par.edgecut.y = ay.get();
//	}
//
//	// angle
//	if (boost::optional<double> ax = pt.get_optional<double>("AngleRange.ax")) {
//		par.ang.x = ax.get();
//	}
//	if (boost::optional<double> ay = pt.get_optional<double>("AngleRange.ay")) {
//		par.ang.y = ay.get();
//	}
//
//
//
//
//	/*
//		basic_ptree<std::string, std::string>::const_iterator iter = pt.get_child("SkipPL").begin(), iterEnd = pt.get_child("SkipPL").end();
//		for (; iter != iterEnd; ++iter)
//		{
//			//->first;  // Key.  Array elements have no names
//			//->second; // The object at each step
//			//std::cout << "=> " << iter->second.get_value<std::string>() << std::endl;
//			spl.insert(std::stoi(iter->second.get_value<std::string>()));
//		}*/
//
//	std::cout << "param" << std::endl;
//	std::cout << "****************************************" << std::endl;
//	std::cout << "maxPL      : " << std::setw(6) << par.rpl.max << std::endl;
//	std::cout << "minPL      : " << std::setw(6) << par.rpl.min << std::endl;
//	std::cout << "FillFactor : " << std::setw(6) << std::fixed << std::setprecision(4) << par.fillfactor << std::endl;
//	std::cout << "nseg       : " << std::setw(6) << par.nseg << std::endl;
//	std::cout << "d(lateral) : " << std::setw(6) << std::fixed << std::setprecision(4) << par.dlateral << std::endl;
//	std::cout << "angle cut  : " << std::endl;
//	std::cout << "             ax <= " << std::setw(7) << std::fixed << std::setprecision(4) << par.ang.x << ", ay <= " << std::fixed << std::setw(7) << std::setprecision(4) << par.ang.y << std::endl;
//	std::cout << "edge cut   : " << std::endl;
//	std::cout << "             x side = " << std::fixed << std::setw(1) << par.edgecut.x << ", y side = " << std::fixed << std::setw(1) << par.edgecut.y << std::endl;
//	std::cout << "****************************************" << std::endl;
//	std::cout << std::endl;
//}
std::vector<mfile0::M_Chain> chain_nseg_selection(std::vector<mfile0::M_Chain>& chain, int nseg) {
	std::vector<mfile0::M_Chain> ret;
	for (auto itr = chain.begin(); itr != chain.end(); itr++) {
		if (itr->nseg < nseg)continue;
		ret.push_back(*itr);
	}
	printf("chain nseg >= %d: %d --> %d (%4.1lf%%)\n", nseg, chain.size(), ret.size(), ret.size() * 100. / chain.size());
	return ret;
}
std::vector<mfile0::M_Chain> group_clustering(std::vector<mfile0::M_Chain>& chain) {

	std::vector<mfile0::M_Chain> ret;
	std::set<int>gid;
	std::multimap<int, mfile0::M_Chain*>group;
	int nseg_max;
	for (auto itr = chain.begin(); itr != chain.end(); itr++) {
		gid.insert(itr->basetracks[0].group_id);
		group.insert(std::make_pair(itr->basetracks[0].group_id, &(*itr)));
	}
	for (auto itr = gid.begin(); itr != gid.end(); itr++) {
		if (group.count(*itr) == 1) {
			ret.push_back(*(group.find(*itr)->second));
		}
		else {
			auto range = group.equal_range(*itr);
			nseg_max = 1;
			for (auto itr2 = range.first; itr2 != range.second; itr2++) {
				nseg_max = std::max(itr2->second->nseg, nseg_max);
			}
			for (auto itr2 = range.first; itr2 != range.second; itr2++) {
				if (itr2->second->nseg == nseg_max) {
					ret.push_back(*(itr2->second));
					break;
				}
			}
		}
	}
	printf("chain group clustering: %d --> %d (%4.1lf%%)\n", chain.size(), ret.size(), ret.size() * 100. / chain.size());
	return ret;

}
std::vector<mfile0::M_Chain> chain_dlat_selection(std::vector<mfile0::M_Chain>& chain, double threshold) {
	std::vector<mfile0::M_Chain> ret;
	double ax, ay;
	for (auto itr = chain.begin(); itr != chain.end(); itr++) {
		ax = mfile0::chain_ax(*itr);
		ay = mfile0::chain_ay(*itr);
		if (mfile0::angle_diff_dev_lat(*itr, ax, ay) > threshold)continue;
		ret.push_back(*itr);
	}
	printf("chain lateral selection <= %5.4lf : %d --> %d (%4.1lf%%)\n", threshold, chain.size(), ret.size(), ret.size() * 100. / chain.size());
	return ret;

}
std::vector<mfile0::M_Chain> chain_angle_selection(std::vector<mfile0::M_Chain>& chain, double thr_ax, double thr_ay) {
	std::vector<mfile0::M_Chain> ret;
	double ax, ay;
	for (auto itr = chain.begin(); itr != chain.end(); itr++) {
		ax = mfile0::chain_ax(*itr);
		ay = mfile0::chain_ay(*itr);
		if (thr_ax - fabs(ax) < 0)continue;
		if (thr_ay - fabs(ay) < 0)continue;
		ret.push_back(*itr);
	}
	printf("chain average angle selection |ax|<=%3.1lf |ay|<=%3.1lf : %d --> %d (%4.1lf%%)\n", thr_ax, thr_ay, chain.size(), ret.size(), ret.size() * 100. / chain.size());
	return ret;
}
std::vector<mfile0::M_Chain> chain_fillfactor(std::vector<mfile0::M_Chain>& chain, double nseg_thr) {
	std::vector<mfile0::M_Chain> ret;
	double tmp;
	int count = 0, npl;
	for (auto itr = chain.begin(); itr != chain.end(); itr++) {
		npl = itr->pos1 / 10 - itr->pos0 / 10 + 1;
		for (auto itr2 = itr->basetracks.begin(); itr2 != itr->basetracks.end(); itr2++) {
			count++;
		}
		tmp = count / npl;
		//std::cout << "check " << tmp << " vs thr: " << nseg_thr << std::endl;
		if (tmp > nseg_thr) {
			ret.push_back(*itr);
			count = 0;
		}
	}
	printf("fillfactor selection : %d --> %d (%4.1lf%%)\n", chain.size(), ret.size(), ret.size() * 100. / chain.size());
	return ret;
}

void output_chain(std::ofstream& ofs, std::vector<mfile0::M_Chain>& chain) {
	int num = 0;
	int count = 0;
	for (auto itr = chain.begin(); itr != chain.end(); itr++) {
		for (auto itr2 = itr->basetracks.begin(); itr2 != itr->basetracks.end(); itr2++) {
			count++;
			if (itr2->pos / 10 == itr->pos0 / 10) {
				ofs << std::right << std::fixed << std::setfill(' ')
					<< std::setw(3) << std::setprecision(0) << itr->pos1 / 10 - itr->pos0 / 10 + 1 << " "
					<< std::setw(3) << std::setprecision(0) << itr2->pos / 10 << " "
					<< std::setw(7) << std::setprecision(0) << itr2->rawid << " "
					<< std::setw(7) << std::setprecision(4) << itr2->ax << " "
					<< std::setw(7) << std::setprecision(4) << itr2->ay << " "
					<< std::setw(8) << std::setprecision(1) << itr2->x << " "
					<< std::setw(8) << std::setprecision(1) << itr2->y << " "
					<< std::setw(5) << std::setprecision(0) << itr2->ph << " "
					<< std::setw(2) << std::setprecision(1) << itr2->z << " ";

			}
			if (itr2->pos / 10 == itr->pos1 / 10) {
				ofs << std::right << std::fixed << std::setfill(' ')
					<< std::setw(3) << std::setprecision(0) << itr2->pos / 10 << " "
					<< std::setw(7) << std::setprecision(0) << itr2->rawid << " "
					<< std::setw(7) << std::setprecision(4) << itr2->ax << " "
					<< std::setw(7) << std::setprecision(4) << itr2->ay << " "
					<< std::setw(8) << std::setprecision(1) << itr2->x << " "
					<< std::setw(8) << std::setprecision(1) << itr2->y << " "
					<< std::setw(5) << std::setprecision(0) << itr2->ph << " "
					<< std::setw(2) << std::setprecision(1) << itr2->z << " "
					<< std::setw(3) << std::setprecision(0) << count << std::endl;
				count = 0;
			}
		}
	}


}
void output_chain2(std::ofstream& ofs, std::vector<mfile0::M_Chain>& chain, int npl) {
	int num = 0;
	int count = 0;
	double tan, dal, dl;
	for (auto itr = chain.begin(); itr != chain.end(); itr++) {

		//if (itr->basetracks.size() < 6)continue;
		//if (itr->basetracks.rbegin()->pos - itr->basetracks.begin()->pos / 10 + 1 < 6)continue;
		if (itr->basetracks.size() != npl)continue;
		if (itr->pos1 / 10 - itr->pos0 / 10 + 1 != npl)continue;

		std::vector<mfile0::M_Base>::iterator tk0, tk1;
		if ((itr->basetracks.begin()->pos / 10 - std::next(itr->basetracks.begin(), 1)->pos/10)%2 == 1) {
			//oo oo oo
 		tk0 = std::next(itr->basetracks.begin(), npl / 2);//3
		tk1 = std::next(itr->basetracks.begin(), npl / 2 + 1);//4
		}



		tan = sqrt(tk0->ax * tk0->ax + tk0->ay * tk0->ay);
		if (tan < 0.00001) {
			dal = 0;
		}
		else {
			dal = ((tk1->ax - tk0->ax) * tk0->ay - (tk1->ay - tk0->ay) * tk0->ax) / tan;
		}
		dl = ((tk1->x - tk0->x) * tk0->ay - (tk1->y - tk0->y) * tk0->ax) / tan;

		ofs << std::right << std::fixed << std::setfill(' ')
			<< std::setw(12) << std::setprecision(0) << itr->chain_id << " "
			<< std::setw(3) << std::setprecision(0) << itr->basetracks.size() << " "
			<< std::setw(3) << std::setprecision(0) << itr->pos1 / 10 - itr->pos0 / 10 + 1 << " "
			<< std::setw(8) << std::setprecision(4) << tan << " "
			<< std::setw(8) << std::setprecision(4) << dal << " "
			<< std::setw(8) << std::setprecision(1) << dl << " "
			<< std::setw(8) << std::setprecision(4) << tk0->ph % 10000 << " "
			<< std::setw(8) << std::setprecision(1) << tk1->ph % 10000 << " "
			<< std::endl;
	}
}
int angle_divide(double angle) {
	if (angle < 0.1)return 0;
	if (0.1 <= angle && angle < 0.3)return 1;
	if (0.3 <= angle && angle < 0.5)return 2;
	if (0.5 <= angle && angle < 0.7)return 3;
	if (0.7 <= angle && angle < 0.9)return 4;
	if (0.9 <= angle && angle < 1.1)return 5;
	if (1.1 <= angle && angle < 1.3)return 6;
	if (1.3 <= angle && angle < 1.5)return 7;
	if (1.5 <= angle && angle < 2.0)return 8;
	if (2.0 <= angle && angle < 2.5)return 9;
	if (2.5 <= angle && angle < 3.0)return 10;
	if (3.0 <= angle && angle < 3.5)return 11;
	if (3.5 <= angle && angle < 4.0)return 12;
	if (4.0 <= angle && angle < 4.5)return 13;
	if (4.5 <= angle && angle < 5.0)return 14;
	if (5.0 <= angle)return 15;
	return -1;
}
