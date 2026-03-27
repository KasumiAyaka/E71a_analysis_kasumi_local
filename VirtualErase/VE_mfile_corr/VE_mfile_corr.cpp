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
	double ax, ay, x, y,z;
	double vph;
};


struct Key {
	int pl, rid;
};
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.pl, lhs.rid) < std::tie(rhs.pl, rhs.rid);
}

struct UnitList {
	int ecc, gid, cid, venum;
	Key k0, k1;
};
bool operator<(const UnitList& lhs, const UnitList& rhs) {
	return std::tie(lhs.ecc, lhs.gid, lhs.cid, lhs.venum, lhs.k0, lhs.k1) < std::tie(rhs.ecc, rhs.gid, rhs.cid, rhs.venum, rhs.k0, rhs.k1);
}

struct Output {
	UnitList u;
	Prediction org, prd, hit;
};

struct Results {

	UnitList i;
	double dl0, dal0;
	int vph0, vph00;
	double dl1, dal1;
	int vph1, vph11;
	double tan0;
	double tan1;

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

void Read_VE_list(std::string input, std::map<Key, UnitList>& ret, std::vector<Results>& res);
void extraporate(std::string output, std::vector<mfile0::M_Chain>& chain, std::map<Key, UnitList>& ve, std::map<int, double>& zmap, std::map<UnitList, Output>& tmp);
void output_results(std::string output, std::map<UnitList, Output>& tmp, std::vector<Results>& mp);

int main(int argc, char** argv) {
	if (argc !=5) {
		fprintf(stderr, "usage:prg in-mfile.all veResult.txt output.txt out_ve.txt\n");
		exit(1);
	}

	std::string file_in_mfile = argv[1];
	std::string file_in_pve_txtaram = argv[2];
	std::string output0 = argv[3];
	std::string output1 = argv[4];

	//mfile
	mfile0::Mfile m;

	//read param
	//Param p;
	//read_param(file_in_param, p);

	//read mfile
	std::cout<<"read mfile"<<std::endl;
	mfile1::read_mfile_extension(file_in_mfile, m);
	std::map<Key, UnitList> ve;
	std::vector<Results> res;
	Read_VE_list(file_in_pve_txtaram, ve, res);

	//zmap
	std::map<int, double> zmap;
	for (auto itr = m.chains.begin(); itr != m.chains.end(); itr++) {
		for (auto itr2 = itr->basetracks.begin(); itr2 != itr->basetracks.end(); itr2++) {
			zmap.insert(std::make_pair(itr2->pos / 10, itr2->z));
		}
	}


	//m.chains = chain_nseg_selection(m.chains, p.nseg);//(m.chains, 15)
	//m.chains = chain_angle_selection(m.chains, p.ang.x, p.ang.y);//(m.chains, 4.0, 4.0);
	//m.chains = chain_dlat_selection(m.chains, p.dlateral);//(m.chains, 0.005);
	//m.chains = chain_fillfactor(m.chains, p.fillfactor);
	//m.chains = group_clustering(m.chains);
	std::map<UnitList, Output> tmp;
	extraporate(output0, m.chains, ve, zmap, tmp);
	output_results(output1,tmp,res);

	std::cout << "finished!!" << std::endl;
}
void Read_VE_list(std::string input, std::map<Key, UnitList>& ret, std::vector<Results>& res) {
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error!" << std::endl;
		exit(0);
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	Results r;
	while (std::getline(ifs, str)) {
		// gid cid pid nseg npl upl dpl
		str_v = StringSplit_with_tab(str);
		std::cout << str << std::endl;
		r.i.ecc == std::stoi(str_v[1]);
		r.i.gid == std::stoi(str_v[2]);
		r.i.cid == std::stoi(str_v[3]);
		r.i.venum == std::stoi(str_v[4]);
		r.i.k0.pl == std::stoi(str_v[5]);
		r.i.k0.rid == std::stoi(str_v[6]);
		r.i.k1.pl == std::stoi(str_v[12]);
		r.i.k1.rid == std::stoi(str_v[13]);

		if (std::stoi(str_v[5]) < 1) {
			ret.insert(std::make_pair(r.i.k1, r.i));
		}
		if (std::stoi(str_v[12]) < 1) {
			ret.insert(std::make_pair(r.i.k0, r.i));
		}

		r.tan0 = std::stod(str_v[7]);
		r.tan1 = std::stod(str_v[14]);
		r.dal0 = std::stod(str_v[8]);
		r.dal1 = std::stod(str_v[15]);
		r.dl0 = std::stod(str_v[9]);
		r.dl1 = std::stod(str_v[16]);
		r.vph0 = std::stod(str_v[10]);
		r.vph0 = std::stod(str_v[17]);
		r.vph00 = std::stod(str_v[11]);
		r.vph11 = std::stod(str_v[18]);
		res.push_back(r);
	}
	std::cout << ret.size() << " " << res.size() << std::endl;
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

void extraporate(std::string output, std::vector<mfile0::M_Chain>& chain, std::map<Key, UnitList>& ve, std::map<int, double>& zmap, std::map<UnitList, Output>& tmp) {
	std::ofstream ofs(output);

	int num = 0;
	int count = 0;
	Key k; int plz;
	double d;
	Output op;
	std::multimap<int, Prediction> map;
	std::set<int> set;
	std::multimap<int, Output> out;
	for (auto itr = chain.begin(); itr != chain.end(); itr++) {
		for (auto itr2 = itr->basetracks.begin(); itr2 != itr->basetracks.end(); itr2++) {
			op.org.ax = itr2->ax;
			op.org.ay = itr2->ay;
			op.org.x = itr2->x;
			op.org.y = itr2->y;
			op.org.z = itr2->z;
			op.org.PL = itr2->pos / 10;
			op.org.vph = itr2->ph % 10000;
			op.org.rawid = itr2->rawid;
			op.org.flg = -1;
			map.insert(std::make_pair(k.pl = itr2->pos / 10, op.org));

			k.pl = itr2->pos / 10;
			k.rid = itr2->rawid;

			if (ve.count(k) == 0)continue;
			op.u = ve.find(k)->second;
			if (k.pl % 2 == 1) {
				// \nu --> 54 76 98
				op.prd.PL = k.pl + 2;
				plz = zmap.find(k.pl + 2)->second;
				// 7-->9Ç÷ÇÃäOë}
			}
			else {
				op.prd.PL = k.pl - 2;
				plz = zmap.find(k.pl - 2)->second;
			}

			op.prd = op.org;
			op.prd.x = op.org.x + op.org.ax * (op.org.z - plz);
			op.prd.x = op.org.y + op.org.ay * (op.org.z - plz);

			out.insert(std::make_pair(op.prd.PL, op));
			set.insert(op.prd.PL);
		}
	}
	chain.clear();

	double pos_lat_all[2]; double pos_rad_all[2]; double ang_lat_all[2]; double ang_rad_all[2];
	//arrowane[0]+ [1]*rad
	pos_lat_all[0] = 500;
	pos_lat_all[1] = 0;
	ang_lat_all[0] = 0.015;
	ang_lat_all[1] = 0;

	double dl = 201, dal = 20;
	double x, y, angle;
	Prediction sel;
	int flg;
	for (auto itr = set.begin(); itr != set.end(); itr++) {
		auto p = map.equal_range(*itr);//btrk
		auto q = out.equal_range(*itr);//ve
		for (auto v_itr = q.first; v_itr != q.second; v_itr++) {//ve
			flg = 0;
			for (auto ir = p.first; ir != p.second; ir++) {//btrk
				angle = sqrt(v_itr->second.prd.ax * v_itr->second.prd.ax + v_itr->second.prd.ay * v_itr->second.prd.ay);
				if (fabs(v_itr->second.prd.x - ir->second.x) > 500)continue;
				if (fabs(v_itr->second.prd.y - ir->second.y) > 500)continue;
				if (sqrt((v_itr->second.prd.x - ir->second.x) * v_itr->second.prd.ay - (v_itr->second.prd.y - ir->second.y) * v_itr->second.prd.ax) / angle > pos_lat_all[0] + pos_lat_all[1] * angle)continue;
				if (sqrt((v_itr->second.prd.ax - ir->second.ax) * v_itr->second.prd.ay - (v_itr->second.prd.ay - ir->second.ay) * v_itr->second.prd.ax) / angle > ang_lat_all[0] + ang_lat_all[1] * angle)continue;

				//Ç»ÇÒÇ©1trackëIÇ‘->hitÇ…ìÀÇ¡çûÇﬁ
				if (sqrt((v_itr->second.prd.x - ir->second.x) * v_itr->second.prd.ay - (v_itr->second.prd.y - ir->second.y) * v_itr->second.prd.ax) / angle < dl) {
					flg = flg + 1;
				}				if (sqrt((v_itr->second.prd.ax - ir->second.ax) * v_itr->second.prd.ay - (v_itr->second.prd.ay - ir->second.ay) * v_itr->second.prd.ax) / angle < dl) {
					flg = flg + 2;
				}
				if (flg > 1)
				{
					dl = sqrt((v_itr->second.prd.x - ir->second.x) * v_itr->second.prd.ay - (v_itr->second.prd.y - ir->second.y) * v_itr->second.prd.ax) / angle;
					dal = sqrt((v_itr->second.prd.ax - ir->second.ax) * v_itr->second.prd.ay - (v_itr->second.prd.ay - ir->second.ay) * v_itr->second.prd.ax) / angle;
					v_itr->second.hit = ir->second;
				}
			}
			tmp.insert(std::make_pair(v_itr->second.u, v_itr->second));

		}
	}

	for (auto itr = out.begin(); itr != out.end(); itr++) {
		ofs << std::right << std::fixed << std::setfill(' ')
			<< std::setw(3) << std::setprecision(0) << itr->second.u.gid << " "
			<< std::setw(7) << std::setprecision(0) << itr->second.u.cid << " "
			<< std::setw(3) << std::setprecision(0) << itr->second.u.venum << " "
			// org
			<< std::setw(3) << std::setprecision(0) << itr->second.org.PL << " "
			<< std::setw(7) << std::setprecision(0) << itr->second.org.rawid << " "
			<< std::setw(3) << std::setprecision(0) << itr->second.org.vph << " "
			<< std::setw(7) << std::setprecision(4) << itr->second.org.ax << " "
			<< std::setw(7) << std::setprecision(4) << itr->second.org.ay << " "
			<< std::setw(8) << std::setprecision(1) << itr->second.org.x << " "
			<< std::setw(8) << std::setprecision(1) << itr->second.org.y << " "
			<< std::setw(2) << std::setprecision(1) << itr->second.org.z << " "
			//pred
			<< std::setw(3) << std::setprecision(0) << itr->second.prd.PL << " "
			<< std::setw(7) << std::setprecision(0) << itr->second.prd.rawid << " "
			<< std::setw(3) << std::setprecision(0) << itr->second.prd.vph << " "
			<< std::setw(7) << std::setprecision(4) << itr->second.prd.ax << " "
			<< std::setw(7) << std::setprecision(4) << itr->second.prd.ay << " "
			<< std::setw(8) << std::setprecision(1) << itr->second.prd.x << " "
			<< std::setw(8) << std::setprecision(1) << itr->second.prd.y << " "
			<< std::setw(2) << std::setprecision(1) << itr->second.prd.z << " "
			//hit
			<< std::setw(3) << std::setprecision(0) << itr->second.hit.PL << " "
			<< std::setw(7) << std::setprecision(0) << itr->second.hit.rawid << " "
			<< std::setw(3) << std::setprecision(0) << itr->second.hit.vph << " "
			<< std::setw(7) << std::setprecision(4) << itr->second.hit.ax << " "
			<< std::setw(7) << std::setprecision(4) << itr->second.hit.ay << " "
			<< std::setw(8) << std::setprecision(1) << itr->second.hit.x << " "
			<< std::setw(8) << std::setprecision(1) << itr->second.hit.y << " "
			<< std::setw(2) << std::setprecision(1) << itr->second.hit.z << " "
			<< std::endl;
	}


	out.clear();
	map.clear();
	set.clear();

}
void output_results(std::string output, std::map<UnitList, Output>& tmp, std::vector<Results>& mp){
	std::ofstream ofs(output);

	double ax, ay, x, y, dax, day, tan;
	for(auto itr=mp.begin();itr!=mp.end();itr++)	{	
		auto f = tmp.find(itr->i);
		if (f != tmp.end()) {
			ax = f->second.hit.ax;
			ay = f->second.hit.ay;
			x = f->second.hit.x;
			y = f->second.hit.y;
			dax = f->second.org.ax;
			day = f->second.org.ay;
			tan = sqrt(ax * ax + ay * ay);

			if (itr->i.k0.pl < 1) {
				itr->tan0 = sqrt(ax * ax + ay * ay);
				itr->dal0 = (day * ax - dax * ay) / tan;
				itr->dl0 = ((x- f->second.org.x) * ax - (y- f->second.org.y) * ay) / tan;
				if (tan < 0.0001) {
					itr->dal0 = dax;
					itr->dl0 = 0;
				}
				itr->i.k0.pl = f->second.hit.PL;
				itr->vph00 = f->second.hit.vph;
			}

			if (itr->i.k1.pl < 1) {
				itr->tan1 = sqrt(ax * ax + ay * ay);
				itr->dal1 = (day * ax - dax * ay) / tan;
				itr->dl1 = ((x - f->second.org.x) * ax - (y - f->second.org.y) * ay) / tan;
				if (tan < 0.0001) {
					itr->dal1 = dax;
					itr->dl1 = 0;
				}
				itr->i.k1.pl = f->second.hit.PL;
				itr->vph11 = f->second.hit.vph;
			}

		}
		ofs << std::right << std::fixed //<< std::setfill(' ')
		<< std::setw(2) << std::setprecision(0) << itr->i.ecc << " "
		<< std::setw(5) << std::setprecision(0) << itr->i.gid << " "
		<< std::setw(3) << std::setprecision(0) << itr->i.cid << " "
		<< std::setw(2) << std::setprecision(0) << itr->i.venum << " "
		<< std::setw(3) << std::setprecision(0) << itr->i.k0.pl << " "
		<< std::setw(10) << std::setprecision(0) << itr->i.k0.rid << " "
		<< std::setw(8) << std::setprecision(4) << itr->tan0 << " "
		<< std::setw(8) << std::setprecision(4) << itr->dal0 << " "
		<< std::setw(8) << std::setprecision(1) << itr->dl0 << " "
		<< std::setw(4) << std::setprecision(0) << itr->vph0 << " "
		<< std::setw(4) << std::setprecision(0) << itr->vph00 << " "
		<< std::setw(3) << std::setprecision(0) << itr->i.k1.pl << " "
		<< std::setw(10) << std::setprecision(0) << itr->i.k1.rid << " "
		<< std::setw(8) << std::setprecision(4) << itr->tan1 << " "
		<< std::setw(8) << std::setprecision(4) << itr->dal1 << " "
		<< std::setw(8) << std::setprecision(1) << itr->dl1 << " "
		<< std::setw(4) << std::setprecision(0) << itr->vph1 << " "
		<< std::setw(4) << std::setprecision(0) << itr->vph11 << std::endl;



	};

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
