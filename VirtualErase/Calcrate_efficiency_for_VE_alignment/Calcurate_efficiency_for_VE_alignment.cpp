// 2025/11/23
// Calcurate_efficiency_for_VE_alignment
// kasumi
// to calcurate efficiency of unittrack and basetrack
//C:\Users\kasumi\source\repos\VirtualErase\x64\Release\Calcrate_efficiency_for_VE_alignment.exe

#define _CRT_SECURE_NO_WARNINGS
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.h>
#include <set>
#include <map>

//for read jsonfile
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <boost/array.hpp>
struct PLRange {
	int max, min;
};
struct Margin {
	double xmin, ymin, xmax, ymax;
};
struct AngleCut {
	double x, y;
};

struct Param {
	int nseg;
	double fillfactor;
	double dlateral;
	//PLRange rpl;
	Margin edgecut;
	AngleCut ang;
	std::map<int,int> pl;
};


using namespace boost::property_tree;
void read_param(std::string param, Param& par);

struct Prediction {
	int PL, flg, rawid;
	double ax, ay, x, y;
	double vph;
	int plnum;
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
void CalcEfficiency(std::ofstream& ofs_eff, std::ofstream& ofs_pred, std::vector<mfile0::M_Chain>  chain, int PL, std::map<int, double> zmap, Param& p);
int angle_divide(double angle);
std::vector<mfile0::M_Chain> chain_VPH_cut(std::vector<mfile0::M_Chain>& chain);
std::vector<mfile0::M_Chain> chain_fillfactor(std::vector<mfile0::M_Chain>& chain, double nseg_thr);

void CalcEfficiency_unit_linklet(std::ofstream& ofs_eff, std::ofstream& ofs_pred, std::vector<mfile0::M_Chain>  chain, int targetpl, std::map<int, double> zmap, Param& p);
void CalcEfficiency_unit_linklet2(std::ofstream& ofs_eff, std::ofstream& ofs_pred, std::vector<mfile0::M_Chain>  chain, int targetpl, std::map<int, double> zmap, Param& p);

int main(int argc, char** argv) {
	if (argc < 5) {
		fprintf(stderr, "usage:prg in-mfile.all in_param.json out-file(eff) out-file(prediction) [mode]\n");
		fprintf(stderr, "mode 1 --> oo o? oo o- / -o oo ?o oo\n");
		fprintf(stderr, "mode 2 --> oo o? oo -- / -- oo ?o oo\n");
		fprintf(stderr, "mode 3 --> -o o? o- -- / -- -o ?o o-\n");
		//fprintf(stderr, "usage:prg in-mfile in_param out-file(eff) out-file(prediction) out_exist_base_id_list\n");
		exit(1);
	}

	//input value
	std::string file_in_mfile = argv[1];// input.all
	std::string file_in_param = argv[2];// input_paran.json
	std::string file_out_eff = argv[3];/// efficiency
	std::string file_out_prediction = argv[4];// pred
	int mode = 1;
	if (argc == 6) {
		mode == std::stoi(argv[5]);
	}
	std::cout << " argc = " << argc << std::endl;
	//mfile
	mfile0::Mfile m;

	//read param
	Param p;
	read_param(file_in_param, p);

	//read mfile
	mfile1::read_mfile_extension(file_in_mfile, m);

	//zmap
	std::map<int, double> zmap;
	for (auto itr = m.chains.begin(); itr != m.chains.end(); itr++) {
		for (auto itr2 = itr->basetracks.begin(); itr2 != itr->basetracks.end(); itr2++) {
			zmap.insert(std::make_pair(itr2->pos / 10, itr2->z));
		}
	}


	m.chains = chain_nseg_selection(m.chains, p.nseg);//(m.chains, 15)
	m.chains = chain_angle_selection(m.chains, p.ang.x, p.ang.y);//(m.chains, 4.0, 4.0);
	m.chains = chain_dlat_selection(m.chains, p.dlateral);//(m.chains, 0.005);
	m.chains = chain_fillfactor(m.chains, p.fillfactor);

	m.chains = group_clustering(m.chains);

	std::ofstream ofs_eff(file_out_eff);
	std::ofstream ofs_pred(file_out_prediction);

	if (mode == 1) {
		for (auto itr = p.pl.begin(); itr != p.pl.end(); itr++) {//pl
			// --> -o oo ?o oo / --> oo o? oo o-
			CalcEfficiency(ofs_eff, ofs_pred, m.chains, itr->first, zmap, p);
			//CalcEfficiency(ofs_eff, ofs_pred, ofs_id, m.chains, PL, zmap, p);	}
		}
	}
	else if (mode == 2) {
		for (auto itr = p.pl.begin(); itr != p.pl.end(); itr++) {//pl
	// --> oo ?o oo / --> oo o? oo
			CalcEfficiency_unit_linklet(ofs_eff, ofs_pred, m.chains, itr->first, zmap, p);
		}
	}
	else if (mode == 3) {
		for (auto itr = p.pl.begin(); itr != p.pl.end(); itr++) {//pl
	// --> oo ?o oo / --> oo o? oo
			CalcEfficiency_unit_linklet2(ofs_eff, ofs_pred, m.chains, itr->first, zmap, p);
		}
	}
	else {
		std::cout << "input correct number" << std::endl;
	}

	std::cout << "finished!!" << std::endl;
}
void read_param(std::string param, Param& par) {
	std::cout << "Set Param" << std::endl;
	par = { 0 };

	ptree pt;
	read_json(param, pt);

	// nseg
	if (boost::optional<int> upl = pt.get_optional<int>("Npl")) {
		par.nseg = upl.get();
	}

	// fill factor
	if (boost::optional<double> upl = pt.get_optional<double>("FillFactor")) {
		par.fillfactor = upl.get();
	}


	// dlateral
	if (boost::optional<double> upl = pt.get_optional<double>("LateralAngleDiff")) {
		par.dlateral = upl.get();
	}


	// margin
	if (boost::optional<double> ax = pt.get_optional<double>("FilmEdgeCut.xmin")) {
		par.edgecut.xmin = ax.get();
	}
	if (boost::optional<double> ay = pt.get_optional<double>("FilmEdgeCut.ymin")) {
		par.edgecut.ymin = ay.get();
	}
	if (boost::optional<double> ax = pt.get_optional<double>("FilmEdgeCut.xmax")) {
		par.edgecut.xmax = ax.get();
	}
	if (boost::optional<double> ay = pt.get_optional<double>("FilmEdgeCut.ymax")) {
		par.edgecut.ymax = ay.get();
	}

	// angle
	if (boost::optional<double> ax = pt.get_optional<double>("AngleRange.ax")) {
		par.ang.x = ax.get();
	}
	if (boost::optional<double> ay = pt.get_optional<double>("AngleRange.ay")) {
		par.ang.y = ay.get();
	}

	// pl list
	basic_ptree<std::string, std::string>::const_iterator iter = pt.get_child("PL_list").begin(), iterEnd = pt.get_child("PL_list").end();
	int cnt = 0;
	for (; iter != iterEnd; ++iter) {
		//->first;  // Key.  Array elements have no names
		//->second; // The object at each step
		//std::cout << "=> " << iter->second.get_value<std::string>() << std::endl;
		par.pl.insert(std::make_pair(std::stoi(iter->second.get_value<std::string>()), cnt));// pl, í Çµî‘çÜ
		cnt++;
	}
	/*
		basic_ptree<std::string, std::string>::const_iterator iter = pt.get_child("SkipPL").begin(), iterEnd = pt.get_child("SkipPL").end();
		for (; iter != iterEnd; ++iter)
		{
			//->first;  // Key.  Array elements have no names
			//->second; // The object at each step
			//std::cout << "=> " << iter->second.get_value<std::string>() << std::endl;
			spl.insert(std::stoi(iter->second.get_value<std::string>()));
		}*/

	std::cout << " param" << std::endl;
	std::cout << "****************************************" << std::endl;
	std::cout << " * FillFactor :  nseg / npl >= " << std::setw(6) << std::fixed << std::setprecision(4) << par.fillfactor << std::endl;
	std::cout << " * # of p     :   npl >" << std::setw(6) << par.nseg << std::endl;
	std::cout << " * d(lateral) :   dl  <" << std::setw(6) << std::fixed << std::setprecision(4) << par.dlateral << std::endl;
	std::cout << " * d(ang-lat) :   dal <" << std::setw(6) << std::fixed << std::setprecision(4) << par.dlateral << std::endl;
	std::cout << " * angle cut  : " << std::endl;
	std::cout << "       ax <= " << std::setw(7) << std::fixed << std::setprecision(4) << par.ang.x << ", ay <= " << std::fixed << std::setw(7) << std::setprecision(4) << par.ang.y << std::endl;
	std::cout << " * Film area  : " << std::endl;
	std::cout << "      (xmax,ymax) = ( " << std::setw(7) << std::fixed << std::setprecision(1) << par.edgecut.xmax << " ," << std::setw(7) << std::fixed << std::setprecision(1) << par.edgecut.ymax << ")" << std::endl;
	std::cout << "      (xmin,ymin) = ( " << std::setw(7) << std::fixed << std::setprecision(1) << par.edgecut.xmin << " ," << std::setw(7) << std::fixed << std::setprecision(1) << par.edgecut.ymin << ")" << std::endl;
	std::cout << " * list of PL ( nu--> )   ( total : " << par.pl.size() << " tracks )" << std::endl;
	int n = 0;
	for (auto itr = par.pl.begin(); itr != par.pl.end(); itr++) {
		std::cout << "(" << std::setw(3) << itr->first << "," << std::setw(3) << itr->second << "),";
		if (n % 6 == 5) {
			std::cout << std::endl;
		}
		n++;
	}
	std::cout << std::endl;
	std::cout << "****************************************" << std::endl;
	std::cout << std::endl;
}

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

std::vector<mfile0::M_Chain> chain_VPH_cut(std::vector<mfile0::M_Chain>& chain) {
	std::vector<mfile0::M_Chain> ret;
	double ax, ay;
	for (auto itr = chain.begin(); itr != chain.end(); itr++) {
		for (auto itr2 = itr->basetracks.begin(); itr2 != itr->basetracks.end(); itr2++) {
			if (itr2->pos / 10 == 27 || itr2->pos / 10 == 28 || itr2->pos / 10 == 29 || itr2->pos / 10 == 39 || itr2->pos / 10 == 40 || itr2->pos / 10 == 41 || itr2->pos / 10 == 42 || itr2->pos / 10 == 42 || itr2->pos / 10 == 44 || itr2->pos / 10 == 45) {
				if (itr2->ph % 10000 < 15 - 7.5 * sqrt(itr2->ax * itr2->ax + itr2->ay * itr2->ay)) {
					continue;
				}
				continue;
			}
		}
		ret.push_back(*itr);
	}
	printf("chain VPH(PL027-029,039-045) selection : %d --> %d (%4.1lf%%)\n", chain.size(), ret.size(), ret.size() * 100. / chain.size());
	return ret;
}

void CalcEfficiency(std::ofstream& ofs_eff, std::ofstream& ofs_pred, std::vector<mfile0::M_Chain>  chain, int targetpl, std::map<int, double> zmap, Param& p){

	std::vector<mfile0::M_Chain> chain_sel;
	int flg = 0;
	double dz;
	int PL_min = 0;
	int PL_max = p.pl.size()-1;
	int PL = p.pl.find(targetpl)->second;//í Çµî‘çÜ

	for (auto itr = chain.begin(); itr != chain.end(); itr++) {
		flg = 0;
		for (auto itr2 = itr->basetracks.begin(); itr2 != itr->basetracks.end(); itr2++) {

			if (PL_min == PL) {
				if (abs(p.pl.find(itr2->pos / 10)->second - PL) != 0 && abs(p.pl.find(itr2->pos / 10)->second - PL) <= 6) {
					flg++;
				}
			}
			else if (PL_min + 1 == PL) {
				if (abs(p.pl.find(itr2->pos / 10)->second - PL) != 0 && abs(p.pl.find(itr2->pos / 10)->second - PL) <= 5) {
					flg++;
				}
			}
			else if (PL_min + 2 == PL) {
				if (abs(p.pl.find(itr2->pos / 10)->second - PL) != 0 && abs(p.pl.find(itr2->pos / 10)->second - PL) <= 4) {
					flg++;
				}
			}
			else if (PL_max - 2 == PL) {
				if (abs(p.pl.find(itr2->pos / 10)->second - PL) != 0 && abs(p.pl.find(itr2->pos / 10)->second - PL) <= 4) {
					flg++;
				}
			}

			else if (PL_max - 1 == PL) {
				if (abs(p.pl.find(itr2->pos / 10)->second - PL) != 0 && abs(p.pl.find(itr2->pos / 10)->second - PL) <= 5) {
					flg++;
				}
			}
			else if (PL_max == PL) {
				if (abs(p.pl.find(itr2->pos / 10)->second - PL) != 0 && abs(p.pl.find(itr2->pos / 10)->second - PL) <= 6) {
					flg++;
				}
			}
			else {
				if (abs(p.pl.find(itr2->pos / 10)->second - PL) != 0 && abs(p.pl.find(itr2->pos / 10)->second - PL) <= 3) {
					flg++;
				}
			}
		}
		if (flg == 6) {
			chain_sel.push_back(*itr);
		}
	}

	printf("PL=%d prediction selection: %d --> %d (%4.1lf%%)\n", targetpl, chain.size(), chain_sel.size(), chain_sel.size() * 100. / chain.size());
	if (chain_sel.size() == 0) {
		printf("PL%03d no prediction track\n", PL);
		return;
	}

	double px0, py0, px1, py1;
	int vph0, vph1;
	int fg;
	std::vector<Output> pred;
	for (auto itr = chain_sel.begin(); itr != chain_sel.end(); itr++) {
		Output out;
		out = { -1 };
		Prediction pred_tmp;
		pred_tmp.PL = PL;
		pred_tmp.ax = mfile0::chain_ax(*itr);// chain-average of ax
		pred_tmp.ay = mfile0::chain_ay(*itr);
		pred_tmp.flg = 0;
		pred_tmp.rawid = 0;
		//pred_tmp.vph = mfile0::chain_vph(*itr);
		pred_tmp.vph = 0;
		Info i0 = { 0 }, i1 = { 0 };
		fg = 0;
		for (auto itr2 = itr->basetracks.begin(); itr2 != itr->basetracks.end(); itr2++) {

			if (p.pl.find(itr2->pos / 10)->second == PL - 1) {
				// extraporate for upstream (PL0-->pl1)
				px0 = itr2->x + pred_tmp.ax * (zmap[PL] - zmap[PL - 1]);
				py0 = itr2->y + pred_tmp.ay * (zmap[PL] - zmap[PL - 1]);
				vph0 = itr2->ph;
				i0.x = px0;
				i0.y = py0;
				i0.ph = vph0;
				i0.ax = itr2->ax;
				i0.ay = itr2->ay;

				fg++;
			}
			else if (p.pl.find(itr2->pos / 10)->second == PL + 1) {
				// downstream extraporation(PL1-->PL0)
				px1 = itr2->x + pred_tmp.ax * (zmap[PL] - zmap[PL + 1]);
				py1 = itr2->y + pred_tmp.ay * (zmap[PL] - zmap[PL + 1]);
				vph1 = itr2->ph;
				i1.x = px1;
				i1.y = py1;
				i1.ph = vph1;
				i1.ax = itr2->ax;
				i1.ay = itr2->ay;
				fg++;
			}
			else if (p.pl.find(itr2->pos / 10)->second == PL) {
				// judge
				out.hit.x = itr2->x;
				out.hit.y = itr2->y;
				out.hit.ax = itr2->ax;
				out.hit.ay = itr2->ay;
				out.hit.PL = PL;
				out.hit.rawid = itr2->rawid;
				out.hit.vph = itr2->ph;
				out.hit.flg = 1;

			}

			if (fg == 2) {
				pred_tmp.x = (i0.x + i1.x) / 2;
				pred_tmp.y = (i0.y + i1.y) / 2;
				pred_tmp.ax = (i0.ax + i1.ax) / 2;
				pred_tmp.ay = (i0.ay + i1.ay) / 2;
				pred_tmp.vph = (i0.ph + i1.ph) / 2;

				fg = 0;
				//	fprintf(stderr, "PL exception PL%03d\n", PL);
			}

			if (p.pl.find(itr2->pos / 10)->second == PL_min + 1 && PL == PL_min) {
				pred_tmp.x = i1.x;
				pred_tmp.y = i1.y;
				pred_tmp.ax = i1.ax;
				pred_tmp.ay = i1.ay;
				pred_tmp.vph = i1.ph;

			}
			else if (p.pl.find(itr2->pos / 10)->second == PL_max - 1 && PL == PL_max) {
				pred_tmp.x = i0.x;
				pred_tmp.y = i0.y;
				pred_tmp.ax = i0.ax;
				pred_tmp.ay = i0.ay;
				pred_tmp.vph = i0.ph;

			}
		}
		out.pre = pred_tmp;
		pred.push_back(out);
	}


	double xmin, xmax, ymin, ymax;
	for (auto itr = pred.begin(); itr != pred.end(); itr++) {
		if (itr == pred.begin()) {
			xmin = itr->pre.x;
			xmax = itr->pre.x;
			ymin = itr->pre.y;
			ymax = itr->pre.y;
		}
		xmin = std::min(itr->pre.x, xmin);
		xmax = std::max(itr->pre.x, xmax);
		ymin = std::min(itr->pre.y, ymin);
		ymax = std::max(itr->pre.y, ymax);
	}
	//xmin = xmin + p.edgecut.xmin;
	//xmax = xmax - p.edgecut.xmax;
	//ymin = ymin + p.edgecut.ymin;
	//ymax = ymax - p.edgecut.ymax;
	std::cout << "  film size" << std::endl;
	std::cout << "      (xmax,ymax) = ( " << std::fixed << std::setw(7) << xmax << " ," << std::setw(7) << ymax << ")" << std::endl;
	std::cout << "      (xmin,ymin) = ( " << std::fixed << std::setw(7) << xmin << " ," << std::setw(7) << ymin << ")" << std::endl;

	xmin = p.edgecut.xmin;
	xmax = p.edgecut.xmax;
	ymin = p.edgecut.ymin;
	ymax = p.edgecut.ymax;

	std::cout << "  remain film size" << std::endl;
	std::cout << "      (xmax,ymax) = ( " << std::fixed << std::setw(7) << xmax << " ," << std::setw(7) << ymax << ")" << std::endl;
	std::cout << "      (xmin,ymin) = ( " << std::fixed << std::setw(7) << xmin << " ," << std::setw(7) << ymin << ")" << std::endl;

	std::vector<Output> pred_sel;
	for (auto itr = pred.begin(); itr != pred.end(); itr++) {
		if (itr->pre.x < xmin)continue;
		if (itr->pre.x > xmax)continue;
		if (itr->pre.y < ymin)continue;
		if (itr->pre.y > ymax)continue;
		pred_sel.push_back(*itr);
	}
	if (pred_sel.size() == 0) {
		printf("PL%03d no prediction track\n", targetpl);
		return;
	}


	const int NUM = 16;
	int all[NUM] = {}, count[NUM] = {};
	double eff[NUM], eff_err[NUM];
	double angle;
	for (auto itr = pred_sel.begin(); itr != pred_sel.end(); itr++) {
		angle = sqrt(itr->pre.ax * itr->pre.ax + itr->pre.ay * itr->pre.ay);
		all[angle_divide(angle)]++;
		if (itr->hit.flg == 1) {
			count[angle_divide(angle)]++;
		}
		ofs_pred << std::right << std::fixed << std::setfill(' ')
			<< std::setw(3) << std::setprecision(0) << targetpl << " "
			<< std::setw(3) << std::setprecision(0) << itr->pre.PL << " "
			<< std::setw(7) << std::setprecision(4) << itr->pre.ax << " "
			<< std::setw(7) << std::setprecision(4) << itr->pre.ay << " "
			<< std::setw(8) << std::setprecision(1) << itr->pre.x << " "
			<< std::setw(8) << std::setprecision(1) << itr->pre.y << " "
			<< std::setw(5) << std::setprecision(1) << itr->pre.vph << " "
			<< std::setw(2) << std::setprecision(0) << itr->hit.flg << " "
			<< std::setw(3) << std::setprecision(0) << itr->hit.PL << " "
			<< std::setw(7) << std::setprecision(4) << itr->hit.ax << " "
			<< std::setw(7) << std::setprecision(4) << itr->hit.ay << " "
			<< std::setw(8) << std::setprecision(1) << itr->hit.x << " "
			<< std::setw(8) << std::setprecision(1) << itr->hit.y << " "
			<< std::setw(5) << std::setprecision(1) << itr->hit.vph << std::endl;
		/*
		if (itr->flg == 1) {
			ofs_id << std::right << std::fixed << std::setfill(' ')
				<< std::setw(3) << std::setprecision(0) << itr->PL << " "
				<< std::setw(12) << std::setprecision(0) << itr->rawid << " "
				<< std::setw(7) << std::setprecision(4) << itr->ax << " "
				<< std::setw(7) << std::setprecision(4) << itr->ay << " "
				<< std::setw(8) << std::setprecision(1) << itr->x << " "
				<< std::setw(8) << std::setprecision(1) << itr->y << std::endl;

		}
		*/
	}

	for (int i = 0; i < 13; i++) {
		eff[i] = double(count[i]) / all[i];
		eff_err[i] = sqrt(count[i] * eff[i] * (1 - eff[i])) / all[i];
		ofs_eff << std::right << std::fixed << std::setfill(' ')
			<< std::setw(4) << std::setprecision(0) << targetpl << " ";

		if (i == 0)ofs_eff << "0.0 0.1 ";
		if (i == 1)ofs_eff << "0.1 0.3 ";
		if (i == 2)ofs_eff << "0.3 0.5 ";
		if (i == 3)ofs_eff << "0.5 0.7 ";
		if (i == 4)ofs_eff << "0.7 0.9 ";
		if (i == 5)ofs_eff << "0.9 1.1 ";
		if (i == 6)ofs_eff << "1.1 1.3 ";
		if (i == 7)ofs_eff << "1.3 1.5 ";
		if (i == 8)ofs_eff << "1.5 2.0 ";
		if (i == 9)ofs_eff << "2.0 2.5 ";
		if (i == 10)ofs_eff << "2.5 3.0 ";
		if (i == 11)ofs_eff << "3.0 3.5 ";
		if (i == 12)ofs_eff << "3.5 4.0 ";
		if (i == 13)ofs_eff << "4.0 4.5 ";
		if (i == 14)ofs_eff << "4.5 5.0 ";
		if (i == 15)ofs_eff << "5.0 10.0 ";

		ofs_eff << std::setw(8) << std::setprecision(0) << all[i] << " "
			<< std::setw(8) << std::setprecision(0) << count[i] << " "
			<< std::setw(5) << std::setprecision(4) << eff[i] << " "
			<< std::setw(5) << std::setprecision(4) << eff_err[i] << std::endl;

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

void CalcEfficiency_unit_linklet(std::ofstream& ofs_eff, std::ofstream& ofs_pred, std::vector<mfile0::M_Chain>  chain, int targetpl, std::map<int, double> zmap, Param& p) {

	std::vector<mfile0::M_Chain> chain_sel;
	int flg = 0;
	double dz;
	int PL_min = 0;
	int PL_max = p.pl.size() - 1;
	int PL = p.pl.find(targetpl)->second;//í Çµî‘çÜ
	int chk_pl;

	for (auto itr = chain.begin(); itr != chain.end(); itr++) {
		flg = 0;
		for (auto itr2 = itr->basetracks.begin(); itr2 != itr->basetracks.end(); itr2++) {
			chk_pl = p.pl.find(itr2->pos / 10)->second;//í Çµî‘çÜ

			if (targetpl % 2 == 0) {
				// oo o? oo
				if (PL == PL_min) {
					// --> oo oo o?
					if (chk_pl != PL) {
						if (chk_pl <= PL + 5) {
							flg++;
						}
					}
				}
				else if (PL == PL_max - 1) {
					// --> o? oo oo
					if (chk_pl != PL) {
						if (chk_pl >= PL - 4) {
							flg++;
						}
					}

				}
				else {
					// \nu --> oo o? oo
					if (chk_pl != PL) {
						if (chk_pl <= PL + 3 && chk_pl >= PL - 2) {
							flg++;
						}
					}
				}
			}
			else {
				if (PL == PL_min + 1) {
					// --> oo oo ?o
					if (chk_pl != PL) {
						if (chk_pl <= PL + 4) {
							flg++;
						}
					}
				}
				else if (PL == PL_max) {
					// --> ?o oo oo
					if (chk_pl != PL) {
						if (chk_pl >= PL - 5) {
							flg++;
						}
					}
				}
				else {
					// \nu --> oo ?o oo
					if (chk_pl != PL) {
						if (chk_pl <= PL + 2 && chk_pl >= PL - 3) {
							flg++;
						}
					}
				}
			}
		}
		if (flg == 5) {
			chain_sel.push_back(*itr);
		}
	}

	printf("PL=%d prediction selection: %d --> %d (%4.1lf%%)\n", targetpl, chain.size(), chain_sel.size(), chain_sel.size() * 100. / chain.size());
	if (chain_sel.size() == 0) {
		printf("PL%03d no prediction track\n", PL);
		return;
	}

	double px0, py0, px1, py1;
	int vph0, vph1;
	int fg;
	std::vector<Output> pred;
	for (auto itr = chain_sel.begin(); itr != chain_sel.end(); itr++) {
		Output out;
		out = { -1 };
		Prediction pred_tmp;
		pred_tmp.PL = PL;//í Çµî‘çÜ
		pred_tmp.ax = mfile0::chain_ax(*itr);// chain-average of ax
		pred_tmp.ay = mfile0::chain_ay(*itr);
		pred_tmp.flg = 0;
		pred_tmp.rawid = 0;
		//pred_tmp.vph = mfile0::chain_vph(*itr);
		pred_tmp.vph = 0;
		Info i0 = { 0 }, i1 = { 0 };
		fg = 0;
		for (auto itr2 = itr->basetracks.begin(); itr2 != itr->basetracks.end(); itr2++) {

			if (p.pl.find(itr2->pos / 10)->second == PL - 1) {
				// extraporate for upstream (PL0-->pl1)
				px0 = itr2->x + pred_tmp.ax * (zmap[PL] - zmap[PL - 1]);
				py0 = itr2->y + pred_tmp.ay * (zmap[PL] - zmap[PL - 1]);
				vph0 = itr2->ph;
				i0.x = px0;
				i0.y = py0;
				i0.ph = vph0;
				i0.ax = itr2->ax;
				i0.ay = itr2->ay;

				fg++;
			}
			else if (p.pl.find(itr2->pos / 10)->second == PL + 1) {
				// downstream extraporation(PL1-->PL0)
				px1 = itr2->x + pred_tmp.ax * (zmap[PL] - zmap[PL + 1]);
				py1 = itr2->y + pred_tmp.ay * (zmap[PL] - zmap[PL + 1]);
				vph1 = itr2->ph;
				i1.x = px1;
				i1.y = py1;
				i1.ph = vph1;
				i1.ax = itr2->ax;
				i1.ay = itr2->ay;
				fg++;
			}
			else if (p.pl.find(itr2->pos / 10)->second == PL) {
				// judge
				out.hit.x = itr2->x;
				out.hit.y = itr2->y;
				out.hit.ax = itr2->ax;
				out.hit.ay = itr2->ay;
				out.hit.PL = PL;
				out.hit.rawid = itr2->rawid;
				out.hit.vph = itr2->ph;
				out.hit.flg = 1;
			}

			if (p.pl.find(itr2->pos / 10)->second == PL_min + 1 && PL == PL_min) {
				pred_tmp.x = i1.x;
				pred_tmp.y = i1.y;
				pred_tmp.ax = i1.ax;
				pred_tmp.ay = i1.ay;
				pred_tmp.vph = i1.ph;
			}
			else if (p.pl.find(itr2->pos / 10)->second == PL_max - 1 && PL == PL_max) {
				pred_tmp.x = i0.x;
				pred_tmp.y = i0.y;
				pred_tmp.ax = i0.ax;
				pred_tmp.ay = i0.ay;
				pred_tmp.vph = i0.ph;
			}
		}
		if (fg == 2) {
			pred_tmp.x = (i0.x + i1.x) / 2;
			pred_tmp.y = (i0.y + i1.y) / 2;
			pred_tmp.ax = (i0.ax + i1.ax) / 2;
			pred_tmp.ay = (i0.ay + i1.ay) / 2;
			pred_tmp.vph = (i0.ph + i1.ph) / 2;

			fg = 0;
			//	fprintf(stderr, "PL exception PL%03d\n", PL);
		}
		out.pre = pred_tmp;
		pred.push_back(out);
	}


	double xmin, xmax, ymin, ymax;
	for (auto itr = pred.begin(); itr != pred.end(); itr++) {
		if (itr == pred.begin()) {
			xmin = itr->pre.x;
			xmax = itr->pre.x;
			ymin = itr->pre.y;
			ymax = itr->pre.y;
		}
		xmin = std::min(itr->pre.x, xmin);
		xmax = std::max(itr->pre.x, xmax);
		ymin = std::min(itr->pre.y, ymin);
		ymax = std::max(itr->pre.y, ymax);
	}
	//xmin = xmin + p.edgecut.xmin;
	//xmax = xmax - p.edgecut.xmax;
	//ymin = ymin + p.edgecut.ymin;
	//ymax = ymax - p.edgecut.ymax;
	std::cout << "  film size" << std::endl;
	std::cout << "      (xmax,ymax) = ( " << std::fixed << std::setw(7) << xmax << " ," << std::setw(7) << ymax << ")" << std::endl;
	std::cout << "      (xmin,ymin) = ( " << std::fixed << std::setw(7) << xmin << " ," << std::setw(7) << ymin << ")" << std::endl;

	xmin = p.edgecut.xmin;
	xmax = p.edgecut.xmax;
	ymin = p.edgecut.ymin;
	ymax = p.edgecut.ymax;

	std::cout << "  remain film size" << std::endl;
	std::cout << "      (xmax,ymax) = ( " << std::fixed << std::setw(7) << xmax << " ," << std::setw(7) << ymax << ")" << std::endl;
	std::cout << "      (xmin,ymin) = ( " << std::fixed << std::setw(7) << xmin << " ," << std::setw(7) << ymin << ")" << std::endl;

	std::vector<Output> pred_sel;
	for (auto itr = pred.begin(); itr != pred.end(); itr++) {
		if (itr->pre.x < xmin)continue;
		if (itr->pre.x > xmax)continue;
		if (itr->pre.y < ymin)continue;
		if (itr->pre.y > ymax)continue;
		pred_sel.push_back(*itr);
	}
	if (pred_sel.size() == 0) {
		printf("PL%03d no prediction track\n", targetpl);
		return;
	}


	const int NUM = 16;
	int all[NUM] = {}, count[NUM] = {};
	double eff[NUM], eff_err[NUM];
	double angle;
	for (auto itr = pred_sel.begin(); itr != pred_sel.end(); itr++) {
		angle = sqrt(itr->pre.ax * itr->pre.ax + itr->pre.ay * itr->pre.ay);
		all[angle_divide(angle)]++;
		if (itr->hit.flg == 1) {
			count[angle_divide(angle)]++;
		}
		ofs_pred << std::right << std::fixed << std::setfill(' ')
			<< std::setw(3) << std::setprecision(0) << targetpl << " "
			<< std::setw(3) << std::setprecision(0) << itr->pre.PL << " "
			<< std::setw(7) << std::setprecision(4) << itr->pre.ax << " "
			<< std::setw(7) << std::setprecision(4) << itr->pre.ay << " "
			<< std::setw(8) << std::setprecision(1) << itr->pre.x << " "
			<< std::setw(8) << std::setprecision(1) << itr->pre.y << " "
			<< std::setw(5) << std::setprecision(1) << itr->pre.vph << " "
			<< std::setw(2) << std::setprecision(0) << itr->hit.flg << " "
			<< std::setw(3) << std::setprecision(0) << itr->hit.PL << " "
			<< std::setw(7) << std::setprecision(4) << itr->hit.ax << " "
			<< std::setw(7) << std::setprecision(4) << itr->hit.ay << " "
			<< std::setw(8) << std::setprecision(1) << itr->hit.x << " "
			<< std::setw(8) << std::setprecision(1) << itr->hit.y << " "
			<< std::setw(5) << std::setprecision(1) << itr->hit.vph << std::endl;
		/*
		if (itr->flg == 1) {
			ofs_id << std::right << std::fixed << std::setfill(' ')
				<< std::setw(3) << std::setprecision(0) << itr->PL << " "
				<< std::setw(12) << std::setprecision(0) << itr->rawid << " "
				<< std::setw(7) << std::setprecision(4) << itr->ax << " "
				<< std::setw(7) << std::setprecision(4) << itr->ay << " "
				<< std::setw(8) << std::setprecision(1) << itr->x << " "
				<< std::setw(8) << std::setprecision(1) << itr->y << std::endl;

		}
		*/
	}

	for (int i = 0; i < 13; i++) {
		eff[i] = double(count[i]) / all[i];
		eff_err[i] = sqrt(count[i] * eff[i] * (1 - eff[i])) / all[i];
		ofs_eff << std::right << std::fixed << std::setfill(' ')
			<< std::setw(4) << std::setprecision(0) << targetpl << " ";

		if (i == 0)ofs_eff << "0.0 0.1 ";
		if (i == 1)ofs_eff << "0.1 0.3 ";
		if (i == 2)ofs_eff << "0.3 0.5 ";
		if (i == 3)ofs_eff << "0.5 0.7 ";
		if (i == 4)ofs_eff << "0.7 0.9 ";
		if (i == 5)ofs_eff << "0.9 1.1 ";
		if (i == 6)ofs_eff << "1.1 1.3 ";
		if (i == 7)ofs_eff << "1.3 1.5 ";
		if (i == 8)ofs_eff << "1.5 2.0 ";
		if (i == 9)ofs_eff << "2.0 2.5 ";
		if (i == 10)ofs_eff << "2.5 3.0 ";
		if (i == 11)ofs_eff << "3.0 3.5 ";
		if (i == 12)ofs_eff << "3.5 4.0 ";
		if (i == 13)ofs_eff << "4.0 4.5 ";
		if (i == 14)ofs_eff << "4.5 5.0 ";
		if (i == 15)ofs_eff << "5.0 10.0 ";

		ofs_eff << std::setw(8) << std::setprecision(0) << all[i] << " "
			<< std::setw(8) << std::setprecision(0) << count[i] << " "
			<< std::setw(5) << std::setprecision(4) << eff[i] << " "
			<< std::setw(5) << std::setprecision(4) << eff_err[i] << std::endl;

	}

}

void CalcEfficiency_unit_linklet2(std::ofstream& ofs_eff, std::ofstream& ofs_pred, std::vector<mfile0::M_Chain>  chain, int targetpl, std::map<int, double> zmap, Param& p) {

	std::vector<mfile0::M_Chain> chain_sel;
	int flg = 0;
	double dz;
	int PL_min = 0;
	int PL_max = p.pl.size() - 1;
	int PL = p.pl.find(targetpl)->second;//í Çµî‘çÜ
	int chk_pl;

	for (auto itr = chain.begin(); itr != chain.end(); itr++) {
		flg = 0;
		for (auto itr2 = itr->basetracks.begin(); itr2 != itr->basetracks.end(); itr2++) {
			chk_pl = p.pl.find(itr2->pos / 10)->second;//í Çµî‘çÜ

			if (targetpl % 2 == 0) {
				// -o o? o-
				if (PL == PL_min) {
					// -->  oo o?
					if (chk_pl != PL) {
						if (chk_pl <= PL + 3) {
							flg++;
						}
					}
				}
				else if (PL == PL_max - 1) {
					// --> o? oo 
					if (chk_pl != PL) {
						if (chk_pl >= PL - 2) {
							flg++;
						}
					}

				}
				else {
					// \nu --> -o o? o-
					if (chk_pl != PL) {
						if (chk_pl <= PL + 2 && chk_pl >= PL - 1) {
							flg++;
						}
					}
				}
			}
			else {
				if (PL == PL_min + 1) {
					// --> -- oo ?o
					if (chk_pl != PL) {
						if (chk_pl <= PL + 2) {
							flg++;
						}
					}
				}
				else if (PL == PL_max) {
					// --> ?o oo --
					if (chk_pl != PL) {
						if (chk_pl >= PL - 3) {
							flg++;
						}
					}
				}
				else {
					// \nu --> -o ?o o-
					if (chk_pl != PL) {
						if (chk_pl <= PL + 1 && chk_pl >= PL - 2) {
							flg++;
						}
					}
				}
			}
		}
		if (flg == 3) {
			chain_sel.push_back(*itr);
		}
	}

	printf("PL=%d prediction selection: %d --> %d (%4.1lf%%)\n", targetpl, chain.size(), chain_sel.size(), chain_sel.size() * 100. / chain.size());
	if (chain_sel.size() == 0) {
		printf("PL%03d no prediction track\n", PL);
		return;
	}

	double px0, py0, px1, py1;
	int vph0, vph1;
	int fg;
	std::vector<Output> pred;
	for (auto itr = chain_sel.begin(); itr != chain_sel.end(); itr++) {
		Output out;
		out = { -1 };
		Prediction pred_tmp;
		pred_tmp.PL = PL;//í Çµî‘çÜ
		pred_tmp.ax = mfile0::chain_ax(*itr);// chain-average of ax
		pred_tmp.ay = mfile0::chain_ay(*itr);
		pred_tmp.flg = 0;
		pred_tmp.rawid = 0;
		//pred_tmp.vph = mfile0::chain_vph(*itr);
		pred_tmp.vph = 0;
		Info i0 = { 0 }, i1 = { 0 };
		fg = 0;
		for (auto itr2 = itr->basetracks.begin(); itr2 != itr->basetracks.end(); itr2++) {

			if (p.pl.find(itr2->pos / 10)->second == PL - 1) {
				// extraporate for upstream (PL0-->pl1)
				px0 = itr2->x + pred_tmp.ax * (zmap[PL] - zmap[PL - 1]);
				py0 = itr2->y + pred_tmp.ay * (zmap[PL] - zmap[PL - 1]);
				vph0 = itr2->ph;
				i0.x = px0;
				i0.y = py0;
				i0.ph = vph0;
				i0.ax = itr2->ax;
				i0.ay = itr2->ay;

				fg++;
			}
			else if (p.pl.find(itr2->pos / 10)->second == PL + 1) {
				// downstream extraporation(PL1-->PL0)
				px1 = itr2->x + pred_tmp.ax * (zmap[PL] - zmap[PL + 1]);
				py1 = itr2->y + pred_tmp.ay * (zmap[PL] - zmap[PL + 1]);
				vph1 = itr2->ph;
				i1.x = px1;
				i1.y = py1;
				i1.ph = vph1;
				i1.ax = itr2->ax;
				i1.ay = itr2->ay;
				fg++;
			}
			else if (p.pl.find(itr2->pos / 10)->second == PL) {
				// judge
				out.hit.x = itr2->x;
				out.hit.y = itr2->y;
				out.hit.ax = itr2->ax;
				out.hit.ay = itr2->ay;
				out.hit.PL = PL;
				out.hit.rawid = itr2->rawid;
				out.hit.vph = itr2->ph;
				out.hit.flg = 1;
			}

			if (p.pl.find(itr2->pos / 10)->second == PL_min + 1 && PL == PL_min) {
				pred_tmp.x = i1.x;
				pred_tmp.y = i1.y;
				pred_tmp.ax = i1.ax;
				pred_tmp.ay = i1.ay;
				pred_tmp.vph = i1.ph;
			}
			else if (p.pl.find(itr2->pos / 10)->second == PL_max - 1 && PL == PL_max) {
				pred_tmp.x = i0.x;
				pred_tmp.y = i0.y;
				pred_tmp.ax = i0.ax;
				pred_tmp.ay = i0.ay;
				pred_tmp.vph = i0.ph;
			}
		}
		if (fg == 2) {
			pred_tmp.x = (i0.x + i1.x) / 2;
			pred_tmp.y = (i0.y + i1.y) / 2;
			pred_tmp.ax = (i0.ax + i1.ax) / 2;
			pred_tmp.ay = (i0.ay + i1.ay) / 2;
			pred_tmp.vph = (i0.ph + i1.ph) / 2;

			fg = 0;
			//	fprintf(stderr, "PL exception PL%03d\n", PL);
		}
		out.pre = pred_tmp;
		pred.push_back(out);
	}


	double xmin, xmax, ymin, ymax;
	for (auto itr = pred.begin(); itr != pred.end(); itr++) {
		if (itr == pred.begin()) {
			xmin = itr->pre.x;
			xmax = itr->pre.x;
			ymin = itr->pre.y;
			ymax = itr->pre.y;
		}
		xmin = std::min(itr->pre.x, xmin);
		xmax = std::max(itr->pre.x, xmax);
		ymin = std::min(itr->pre.y, ymin);
		ymax = std::max(itr->pre.y, ymax);
	}
	//xmin = xmin + p.edgecut.xmin;
	//xmax = xmax - p.edgecut.xmax;
	//ymin = ymin + p.edgecut.ymin;
	//ymax = ymax - p.edgecut.ymax;
	std::cout << "  film size" << std::endl;
	std::cout << "      (xmax,ymax) = ( " << std::fixed << std::setw(7) << xmax << " ," << std::setw(7) << ymax << ")" << std::endl;
	std::cout << "      (xmin,ymin) = ( " << std::fixed << std::setw(7) << xmin << " ," << std::setw(7) << ymin << ")" << std::endl;

	xmin = p.edgecut.xmin;
	xmax = p.edgecut.xmax;
	ymin = p.edgecut.ymin;
	ymax = p.edgecut.ymax;

	std::cout << "  remain film size" << std::endl;
	std::cout << "      (xmax,ymax) = ( " << std::fixed << std::setw(7) << xmax << " ," << std::setw(7) << ymax << ")" << std::endl;
	std::cout << "      (xmin,ymin) = ( " << std::fixed << std::setw(7) << xmin << " ," << std::setw(7) << ymin << ")" << std::endl;

	std::vector<Output> pred_sel;
	for (auto itr = pred.begin(); itr != pred.end(); itr++) {
		if (itr->pre.x < xmin)continue;
		if (itr->pre.x > xmax)continue;
		if (itr->pre.y < ymin)continue;
		if (itr->pre.y > ymax)continue;
		pred_sel.push_back(*itr);
	}
	if (pred_sel.size() == 0) {
		printf("PL%03d no prediction track\n", targetpl);
		return;
	}


	const int NUM = 16;
	int all[NUM] = {}, count[NUM] = {};
	double eff[NUM], eff_err[NUM];
	double angle;
	for (auto itr = pred_sel.begin(); itr != pred_sel.end(); itr++) {
		angle = sqrt(itr->pre.ax * itr->pre.ax + itr->pre.ay * itr->pre.ay);
		all[angle_divide(angle)]++;
		if (itr->hit.flg == 1) {
			count[angle_divide(angle)]++;
		}
		ofs_pred << std::right << std::fixed << std::setfill(' ')
			<< std::setw(3) << std::setprecision(0) << targetpl << " "
			<< std::setw(3) << std::setprecision(0) << itr->pre.PL << " "
			<< std::setw(7) << std::setprecision(4) << itr->pre.ax << " "
			<< std::setw(7) << std::setprecision(4) << itr->pre.ay << " "
			<< std::setw(8) << std::setprecision(1) << itr->pre.x << " "
			<< std::setw(8) << std::setprecision(1) << itr->pre.y << " "
			<< std::setw(5) << std::setprecision(0) << itr->pre.vph << " "
			<< std::setw(2) << std::setprecision(0) << itr->hit.flg << " "
			<< std::setw(3) << std::setprecision(0) << itr->hit.PL << " "
			<< std::setw(7) << std::setprecision(4) << itr->hit.ax << " "
			<< std::setw(7) << std::setprecision(4) << itr->hit.ay << " "
			<< std::setw(8) << std::setprecision(1) << itr->hit.x << " "
			<< std::setw(8) << std::setprecision(1) << itr->hit.y << " "
			<< std::setw(5) << std::setprecision(0) << itr->hit.vph << std::endl;
		/*
		if (itr->flg == 1) {
			ofs_id << std::right << std::fixed << std::setfill(' ')
				<< std::setw(3) << std::setprecision(0) << itr->PL << " "
				<< std::setw(12) << std::setprecision(0) << itr->rawid << " "
				<< std::setw(7) << std::setprecision(4) << itr->ax << " "
				<< std::setw(7) << std::setprecision(4) << itr->ay << " "
				<< std::setw(8) << std::setprecision(1) << itr->x << " "
				<< std::setw(8) << std::setprecision(1) << itr->y << std::endl;

		}
		*/
	}

	for (int i = 0; i < 13; i++) {
		eff[i] = double(count[i]) / all[i];
		eff_err[i] = sqrt(count[i] * eff[i] * (1 - eff[i])) / all[i];
		ofs_eff << std::right << std::fixed << std::setfill(' ')
			<< std::setw(4) << std::setprecision(0) << targetpl << " ";

		if (i == 0)ofs_eff << "0.0 0.1 ";
		if (i == 1)ofs_eff << "0.1 0.3 ";
		if (i == 2)ofs_eff << "0.3 0.5 ";
		if (i == 3)ofs_eff << "0.5 0.7 ";
		if (i == 4)ofs_eff << "0.7 0.9 ";
		if (i == 5)ofs_eff << "0.9 1.1 ";
		if (i == 6)ofs_eff << "1.1 1.3 ";
		if (i == 7)ofs_eff << "1.3 1.5 ";
		if (i == 8)ofs_eff << "1.5 2.0 ";
		if (i == 9)ofs_eff << "2.0 2.5 ";
		if (i == 10)ofs_eff << "2.5 3.0 ";
		if (i == 11)ofs_eff << "3.0 3.5 ";
		if (i == 12)ofs_eff << "3.5 4.0 ";
		if (i == 13)ofs_eff << "4.0 4.5 ";
		if (i == 14)ofs_eff << "4.5 5.0 ";
		if (i == 15)ofs_eff << "5.0 10.0 ";

		ofs_eff << std::setw(8) << std::setprecision(0) << all[i] << " "
			<< std::setw(8) << std::setprecision(0) << count[i] << " "
			<< std::setw(5) << std::setprecision(4) << eff[i] << " "
			<< std::setw(5) << std::setprecision(4) << eff_err[i] << std::endl;

	}

}
