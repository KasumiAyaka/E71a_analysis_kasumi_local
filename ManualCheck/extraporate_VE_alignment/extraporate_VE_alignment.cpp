#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

#include <filesystem>
#include <set>

#include <omp.h>

class align_param {
public:
	int id, signal, ix, iy;
	//視野中心
	double x, y, z;
	//parameter(9);
	double dx, dy, dz, x_rot, y_rot, z_rot, x_shrink, y_shrink, z_shrink, yx_shear, zx_shear, zy_shear;

};
class align_param2 {
public:
	align_param* corr_p[3];
	//3点の視野中心の重心(回転中心)
	double x, y, z;
	//parameter(9);
	double dx, dy, dz, x_rot, y_rot, z_rot, x_shrink, y_shrink, z_shrink, yx_shear, zx_shear, zy_shear;
public:
	//3つのparameterから計算
	void Calc_9param();
	void Calc_9param_shift();

};


class output_format_micro {
public:
	int pos, view, imager, zone, isg, ph, vph, px;
};
class output_format_base {
public:
	int pl, rawid;
	double ax, ay, x, y, z;
	output_format_micro m[2];
};

class output_format_link {
public:
	output_format_base b[2];
	double dax, day, dx, dy, dar, dal, dr, dl;
	void Calc_difference();

};

vxx::base_track_t pick_target_track(std::string file_in_ECC, int pl, int rawid);
std::vector<vxx::base_track_t> near_tracks(std::string file_in_ECC, int pl, vxx::base_track_t& base, double search_area, double search_angle);
double get_nominal_gap(std::string file_in_ECC, int pl0, int pl1);

std::vector<align_param> read_ali_param(std::string filename, bool output);
std::pair<vxx::base_track_t, align_param2>track_affineparam_correspondence(vxx::base_track_t& base, std::vector <align_param2>& param);
align_param2 search_param(std::vector<align_param*>& param, vxx::base_track_t& base, std::multimap<int, align_param2>& triangles);
double select_triangle_vale(align_param2 param, vxx::base_track_t& base);
void trans_base(vxx::base_track_t& base, align_param2* param);
std::vector<vxx::base_track_t> extra_track_all(std::vector<vxx::base_track_t>& base, align_param2& param);

std::vector<output_format_link> read_linket(std::string filename);
std::vector<output_format_link> match_baseid_extract(std::vector<vxx::base_track_t>& base, std::vector<output_format_link>& link, std::vector<std::pair<int, int>>& rawid_connect);

void output_file(std::string filename, vxx::base_track_t& ex_base, std::vector<vxx::base_track_t>& ali_extra, std::vector<vxx::base_track_t>& ali_connect, std::vector<std::pair<int, int>>& rawid_connect);
std::vector< align_param> search_param(vxx::base_track_t& base, std::vector<align_param>& param);
vxx::base_track_t extra_track(vxx::base_track_t base, align_param2& param);
std::vector<vxx::base_track_t> match_tracks(std::string file_in_ECC, int pl, std::vector<output_format_link>& link);
vxx::base_track_t extra_track(vxx::base_track_t base, align_param2& param);
std::vector<vxx::base_track_t> near_tracks_connect(std::string file_in_ECC, int pl, vxx::base_track_t& base, double search_area, double search_angle);
void output_file(std::string filename, vxx::base_track_t& base, std::vector<vxx::base_track_t>& ali);



int main(int argc, char** argv) {
	if (argc != 6) {
		fprintf(stderr, "usage prg extra-base-pl extra-base-rwaid to-pl ECC-path file-out\n");
		exit(1);
	}
	int pl = std::stoi(argv[1]);
	int rawid = std::stoi(argv[2]);
	int to_pl = std::stoi(argv[3]);
	std::string file_in_ECC = argv[4];
	std::string file_out_name = argv[5];

	if (pl > to_pl) {
		fprintf(stderr, "not supported PL%03d --> PL%03d\n", pl, to_pl);
		fprintf(stderr, "only PL%03d <-- PL%03d is OK\n", pl, to_pl);
		exit(1);
	}
	else if (pl == to_pl) {
		//外挿前のtrack
		vxx::base_track_t traget_base = pick_target_track(file_in_ECC, pl, rawid);
		//x(y)+-5000um以内 ax(ay)+-0.2未満の飛跡を抽出
		//距離でsort済み
		std::vector<vxx::base_track_t> ali_tracks = near_tracks_connect(file_in_ECC, pl, traget_base, 5000, 0.2);
		output_file(file_out_name, traget_base, ali_tracks);
		return 0;
	}
	//外挿前のtrack
	vxx::base_track_t traget_base = pick_target_track(file_in_ECC, pl, rawid);
	//x(y)+-5000um以内 ax(ay)+-0.2未満の飛跡を抽出
	//距離でsort済み
	std::vector<vxx::base_track_t> ali_tracks = near_tracks(file_in_ECC, pl, traget_base, 5000, 0.2);

	//補正マップの読み込み
	std::stringstream file_in_align;
	file_in_align << file_in_ECC << "\\Area0\\0\\align\\fine\\ali_"
		<< std::setw(3) << std::setfill('0') << pl << "_"
		<< std::setw(3) << std::setfill('0') << to_pl << "_interpolation.txt";
	std::vector<align_param> ali_param = read_ali_param(file_in_align.str(), false);

	printf("fin\n");
	std::vector<align_param> param_3 = search_param(traget_base, ali_param);
	align_param2 param;
	param.corr_p[0] = &(param_3[0]);
	param.corr_p[1] = &(param_3[1]);
	param.corr_p[2] = &(param_3[2]);
	//= search_param(traget_base, ali_param);
	//printf("fin\n");
	//printf("%d %d %d\n", param.x, param.y, param.z);
	//for (int i = 0; i < 3; i++) {
	//	printf("%3d %3d %8.1lf %8.1lf %8.1lf %8.4lf %8.4lf %8.4lf\n",
	//		param.corr_p[i]->ix,
	//		param.corr_p[i]->iy,
	//		param.corr_p[i]->x,
	//		param.corr_p[i]->y,
	//		param.corr_p[i]->z,
	//		param.corr_p[i]->dx,
	//		param.corr_p[i]->dy,
	//		param.corr_p[i]->dz);
	//}

	//外挿
	vxx::base_track_t traget_base_ex = extra_track(traget_base, param);
	std::vector<vxx::base_track_t> ali_tracks_ex = extra_track_all(ali_tracks, param);

	printf("fin\n");

	//linkletを使ってつながる飛跡を抽出
	std::stringstream file_in_link;
	file_in_link << file_in_ECC << "\\Area0\\0\\linklet_3d\\link_"
		<< std::setw(3) << std::setfill('0') << pl << "_"
		<< std::setw(3) << std::setfill('0') << to_pl << ".sel.bin";
	std::vector<output_format_link> linklet = read_linket(file_in_link.str());

	std::vector<std::pair<int, int>> rawid_connect;
	linklet = match_baseid_extract(ali_tracks, linklet, rawid_connect);
	std::vector<vxx::base_track_t> ali_tracks_connect = match_tracks(file_in_ECC, to_pl, linklet);
	printf("fin\n");

	//traget_base_ex
	//traget_base
	//ali_tracks
	//ali_tracks_ex
	//ali_tracks_connect
	//rawid_connect
	output_file(file_out_name, traget_base_ex, ali_tracks_ex, ali_tracks_connect, rawid_connect);

}
vxx::base_track_t pick_target_track(std::string file_in_ECC, int pl, int rawid) {
	std::stringstream file_in_base;
	file_in_base << file_in_ECC << "\\Area0\\PL" << std::setw(3) << std::setfill('0') << pl
		<< "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";

	std::vector<vxx::base_track_t> base;
	vxx::BvxxReader br;
	std::array<int, 2> index = { rawid, rawid + 1 };//1234<=rawid<=5678であるようなものだけを読む。
	base = br.ReadAll(file_in_base.str(), pl, 0, vxx::opt::index = index);
	if (base.size() != 1) {
		fprintf(stderr, "PL%03d rawid=%d not found\n", pl, rawid);
		exit(1);
	}
	return *base.begin();

}

std::vector<vxx::base_track_t> near_tracks(std::string file_in_ECC, int pl, vxx::base_track_t& base, double search_area, double search_angle) {

	std::stringstream file_in_base;
	file_in_base << file_in_ECC << "\\Area0\\PL" << std::setw(3) << std::setfill('0') << pl
		<< "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";

	std::vector<vxx::base_track_t> base_all;
	vxx::BvxxReader br;
	std::vector<vxx::CutArea> area;
	area.push_back(vxx::CutArea(base.x - search_area, base.x + search_area, base.y - search_area, base.y + search_area, base.ax - search_angle, base.ax + search_angle, base.ay - search_angle, base.ay + search_angle));
	//std::array<int, 2> index = { rawid, rawid + 1 };//1234<=rawid<=5678であるようなものだけを読む。
	base_all = br.ReadAll(file_in_base.str(), pl, 0, vxx::opt::a = area);

	std::multimap<double, vxx::base_track_t> base_dist;
	double dist;
	for (auto itr = base_all.begin(); itr != base_all.end(); itr++) {
		dist = sqrt(pow(itr->x - base.x, 2) + pow(itr->y - base.y, 2));
	}
	std::vector<vxx::base_track_t> ret;
	for (auto itr = base_all.begin(); itr != base_all.end(); itr++) {
		ret.push_back(*itr);
	}
	if (ret.size() < 2) {
		fprintf(stderr, "near track not found\n");
		exit(1);
	}
	return ret;
}

std::vector<vxx::base_track_t> near_tracks_connect(std::string file_in_ECC, int pl, vxx::base_track_t& base, double search_area, double search_angle) {

	std::stringstream file_in_base;
	file_in_base << file_in_ECC << "\\Area0\\PL" << std::setw(3) << std::setfill('0') << pl
		<< "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.connect.vxx";
	std::vector<vxx::base_track_t> base_all;
	vxx::BvxxReader br;
	std::vector<vxx::CutArea> area;
	area.push_back(vxx::CutArea(base.x - search_area, base.x + search_area, base.y - search_area, base.y + search_area, base.ax - search_angle, base.ax + search_angle, base.ay - search_angle, base.ay + search_angle));
	//std::array<int, 2> index = { rawid, rawid + 1 };//1234<=rawid<=5678であるようなものだけを読む。
	base_all = br.ReadAll(file_in_base.str(), pl, 0, vxx::opt::a = area);

	std::multimap<double, vxx::base_track_t> base_dist;
	double dist;
	for (auto itr = base_all.begin(); itr != base_all.end(); itr++) {
		dist = sqrt(pow(itr->x - base.x, 2) + pow(itr->y - base.y, 2));
	}
	std::vector<vxx::base_track_t> ret;
	for (auto itr = base_all.begin(); itr != base_all.end(); itr++) {
		ret.push_back(*itr);
	}
	if (ret.size() < 2) {
		fprintf(stderr, "near track not found\n");
		exit(1);
	}
	return ret;
}

std::vector<vxx::base_track_t> match_tracks(std::string file_in_ECC, int pl, std::vector<output_format_link>& link) {

	std::stringstream file_in_base;
	file_in_base << file_in_ECC << "\\Area0\\PL" << std::setw(3) << std::setfill('0') << pl
		<< "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";

	std::vector<vxx::base_track_t> base_all;
	vxx::BvxxReader br;
	base_all = br.ReadAll(file_in_base.str(), pl, 0);


	std::set<int> rawid;
	for (auto itr = link.begin(); itr != link.end(); itr++) {
		rawid.insert(itr->b[1].rawid);

	}
	std::vector<vxx::base_track_t> ret;
	for (auto itr = base_all.begin(); itr != base_all.end(); itr++) {
		if (rawid.count(itr->rawid) == 0)continue;
		ret.push_back(*itr);
	}
	return ret;


}

double get_nominal_gap(std::string file_in_ECC, int pl0, int pl1) {
	//gap nominal read
	std::stringstream structure_path;
	structure_path << file_in_ECC << "\\st\\st.dat";
	chamber1::Chamber chamber;
	chamber1::read_structure(structure_path.str(), chamber);
	std::map<int, double> z_map = chamber1::base_z_convert(chamber);
	if (z_map.count(pl0) + z_map.count(pl1) != 2) {
		fprintf(stderr, "nominal gap not found PL%03d - PL%03d\n", pl0, pl1);
		exit(1);
	}
	double z0 = z_map.at(pl0);
	double z1 = z_map.at(pl1);

	return z0 - z1;

}



std::vector<align_param> read_ali_param(std::string filename, bool output) {

	std::vector<align_param> ret;
	align_param param_tmp;
	std::ifstream ifs(filename);

	while (ifs >> param_tmp.id >> param_tmp.ix >> param_tmp.iy >> param_tmp.signal
		>> param_tmp.x >> param_tmp.y >> param_tmp.z
		>> param_tmp.x_rot >> param_tmp.y_rot >> param_tmp.z_rot
		>> param_tmp.x_shrink >> param_tmp.y_shrink >> param_tmp.z_shrink
		>> param_tmp.yx_shear >> param_tmp.zx_shear >> param_tmp.zy_shear
		>> param_tmp.dx >> param_tmp.dy >> param_tmp.dz) {
		ret.push_back(param_tmp);
		//printf("ix %d iy%d\n", param_tmp.ix, param_tmp.iy);

	}
	if (output == 1) {
		fprintf(stderr, "%s input finish\n", filename.c_str());
	}
	if (ret.size() == 0) {
		fprintf(stderr, "%s alignment miss!\n", filename.c_str());
		exit(1);
	}
	return ret;

}

std::vector <align_param2 >DelaunayDivide(std::vector <align_param >& corr) {

	//delaunay分割
	std::vector<double> x, y;
	for (auto itr = corr.begin(); itr != corr.end(); itr++) {
		x.push_back(itr->x);
		y.push_back(itr->y);
	}

	delaunay::DelaunayTriangulation DT(x, y); // (std::vector<double> x, std::vector<double> y, uint32_t seed_)
	DT.execute(); // (double min_delta = 1e-6, double max_delta = 1e-5, int max_miss_count = 30)
	std::vector<delaunay::Edge> edge = DT.get_edges();

	std::multimap<int, int> edge_map;

	for (auto itr = edge.begin(); itr != edge.end(); itr++) {
		edge_map.insert(std::make_pair(std::min(itr->first, itr->second), std::max(itr->first, itr->second)));

	}
	std::set<std::tuple<int, int, int>>triangle;
	std::set<int> vertex;
	for (auto itr = edge_map.begin(); itr != edge_map.end(); itr++) {
		//itr->firstの点=aを通る三角形の探索
		vertex.clear();
		auto range = edge_map.equal_range(itr->first);
		//aを通りitr->secondの点=bに行く。bのsetを作成
		for (auto res = range.first; res != range.second; res++) {
			vertex.insert(res->second);
		}
		//bを通る線分の探索
		for (auto itr2 = vertex.begin(); itr2 != vertex.end(); itr2++) {
			if (edge_map.count(*itr2) == 0)continue;
			auto range2 = edge_map.equal_range(*itr2);
			//bを通る線分の中からaから始まる線分を探す
			for (auto res = range2.first; res != range2.second; res++) {
				if (vertex.count(res->second) == 1) {
					triangle.insert(std::make_tuple(itr->first, *itr2, res->second));
				}
			}

		}
	}

	std::vector <align_param2 > ret;
	for (auto itr = triangle.begin(); itr != triangle.end(); itr++) {
		//printf("delaunay triangle %d %d %d\n", std::get<0>(*itr), std::get<1>(*itr), std::get<2>(*itr));
		align_param2 param;
		param.corr_p[0] = &(corr[std::get<0>(*itr)]);
		param.corr_p[1] = &(corr[std::get<1>(*itr)]);
		param.corr_p[2] = &(corr[std::get<2>(*itr)]);
		param.x = 0;
		param.y = 0;
		param.z = 0;
		param.z_shrink = 0;
		for (int i = 0; i < 3; i++) {
			param.x += param.corr_p[i]->x;
			param.y += param.corr_p[i]->y;
			param.z += param.corr_p[i]->z;
			param.z_shrink += param.corr_p[i]->z_shrink;
		}
		param.x = param.x / 3;
		param.y = param.y / 3;
		param.z = param.z / 3;
		param.z_shrink = param.z_shrink / 3;
		param.zx_shear = 0;
		param.zy_shear = 0;

		param.Calc_9param();

		ret.push_back(param);
	}

	return ret;

}
void align_param2::Calc_9param() {


	double bp[3][3], ap[3][3], cos_z, sin_z;
	for (int i = 0; i < 3; i++) {
		bp[i][0] = corr_p[i]->x;
		bp[i][1] = corr_p[i]->y;
		bp[i][2] = corr_p[i]->z;

		cos_z = cos(corr_p[i]->z_rot);
		sin_z = sin(corr_p[i]->z_rot);


		ap[i][0] = corr_p[i]->x_shrink * cos_z * (corr_p[i]->x) - corr_p[i]->y_shrink * sin_z * (corr_p[i]->y) + corr_p[i]->dx;
		ap[i][1] = corr_p[i]->x_shrink * sin_z * (corr_p[i]->x) + corr_p[i]->y_shrink * cos_z * (corr_p[i]->y) + corr_p[i]->dy;
		ap[i][2] = corr_p[i]->z + corr_p[i]->dz;
		//printf("bp%d %8.1lf %8.1lf %8.1lf\n",i, bp[i][0], bp[i][1], bp[i][2]);
		//printf("ap%d %8.1lf %8.1lf %8.1lf\n", i,ap[i][0], ap[i][1], ap[i][2]);
	}
	//apの位置ずれvectorを定義
	double dp[2][3];
	for (int i = 0; i < 3; i++) {
		dp[0][i] = ap[1][i] - ap[0][i];
		dp[1][i] = ap[2][i] - ap[0][i];
	}
	//printf("0-->1 x,y,z : %.1lf %.1lf %.1lf\n", dp[0][0], dp[0][1], dp[0][2]);
	//printf("0-->2 x,y,z : %.1lf %.1lf %.1lf\n", dp[1][0], dp[1][1], dp[1][2]);
	//法線vector
	double n_v[3];
	n_v[0] = (dp[0][1] * dp[1][2] - dp[0][2] * dp[1][1]);
	n_v[1] = (dp[0][2] * dp[1][0] - dp[0][0] * dp[1][2]);
	n_v[2] = (dp[0][0] * dp[1][1] - dp[0][1] * dp[1][0]);

	//std::cout << "normal vector" << std::endl;
	//for (int i = 0; i < 3; i++) {
	//	std::cout << std::setw(14) << std::fixed << std::setprecision(10) << n_v[i] << std::endl;
	//}

	x_rot = atan(n_v[1] / n_v[2]);
	n_v[1] = cos(x_rot) * n_v[1] - sin(x_rot) * n_v[2];
	n_v[2] = sin(x_rot) * n_v[1] + cos(x_rot) * n_v[2];
	//std::cout << "normal vector" << std::endl;
	//for (int i = 0; i < 3; i++) {
	//	std::cout << std::setw(14) << std::fixed << std::setprecision(10) << n_v[i] << std::endl;
	//}
	y_rot = atan(-1 * n_v[0] / n_v[2]);
	n_v[0] = cos(y_rot) * n_v[0] + sin(y_rot) * n_v[2];
	n_v[2] = -1 * sin(y_rot) * n_v[0] + cos(y_rot) * n_v[2];

	//std::cout << "normal vector" << std::endl;
	//for (int i = 0; i < 3; i++) {
	//	std::cout << std::setw(14) << std::fixed << std::setprecision(10) << n_v[i] << std::endl;
	//}

	//printf("x rot:%.6lf\n", x_rot);
	//printf("y rot:%.6lf\n", y_rot);


	matrix_3D::matrix_33 x_rot_mat(0, x_rot), y_rot_mat(1, y_rot);
	matrix_3D::vector_3D ap_v[3];
	for (int i = 0; i < 3; i++) {
		ap_v[i].x = ap[i][0];
		ap_v[i].y = ap[i][1];
		ap_v[i].z = ap[i][2];
	}
	for (int i = 0; i < 3; i++) {
		ap_v[i].matrix_multiplication(x_rot_mat);
		ap_v[i].matrix_multiplication(y_rot_mat);
	}
	//for (int i = 0; i < 3; i++) {
	//	printf("point %d\n", i);
	//	printf("\t %.2lf %.2lf %.2lf\n", bp[i][0], bp[i][1], bp[i][2]);
	//	printf("\t %.2lf %.2lf %.2lf\n", ap_v[i].x, ap_v[i].y, ap_v[i].z);
	//}
	dz = (ap_v[0].z - bp[0][2] + ap_v[1].z - bp[1][2] + ap_v[2].z - bp[2][2]) / 3;
	//printf("dz=%.2lf\n", dz);
	//3元方程式を解く
	double a[2][3][3] = { { {bp[0][0],bp[0][1],1},{bp[1][0],bp[1][1],1},{bp[2][0],bp[2][1],1} },  { {bp[0][0],bp[0][1],1},{bp[1][0],bp[1][1],1},{bp[2][0],bp[2][1],1} } };
	double b[2][3] = { { ap_v[0].x,ap_v[1].x,ap_v[2].x },{ ap_v[0].y,ap_v[1].y,ap_v[2].y } };
	double c[2][3] = { {1, 1, 1},{1,1,1} };
	//gauss(a[0], b[0], c[0]);
	//gauss(a[1], b[1], c[1]);
	GaussJorden(a[0], b[0], c[0]);
	GaussJorden(a[1], b[1], c[1]);
	z_rot = atan(c[1][0] / c[0][0]);
	x_shrink = c[0][0] / cos(z_rot);
	y_shrink = (c[0][0] * c[1][1] - c[0][1] * c[1][0]) / (c[0][0] * cos(z_rot) + c[1][0] * sin(z_rot));
	yx_shear = (c[0][1] * cos(z_rot) + c[1][1] * sin(z_rot)) / (c[0][0] * cos(z_rot) + c[1][0] * sin(z_rot));

	dx = c[0][2];
	dy = c[1][2];
	matrix_3D::vector_3D dr;
	dr.x = c[0][2];
	dr.y = c[1][2];
	dr.z = dz;

	x_rot = x_rot * -1;
	y_rot = y_rot * -1;
	matrix_3D::matrix_33 x_rot_mat_inv(0, x_rot), y_rot_mat_inv(1, y_rot);


	dr.matrix_multiplication(y_rot_mat_inv);
	dr.matrix_multiplication(x_rot_mat_inv);

	dx = dr.x;
	dy = dr.y;
	dz = dr.z;

	//printf("x rot: %.6lf\n",x_rot);
	//printf("y rot: %.6lf\n",y_rot);
	//printf("z rot: %.6lf\n",z_rot);
	//printf("x shrink: %.6lf\n", x_shrink);
	//printf("y shrink: %.6lf\n", y_shrink);
	//printf("z shrink: %.6lf\n", z_shrink);
	//printf("x shift: %.5lf\n", dx);
	//printf("y shift: %.5lf\n", dy);
	//printf("z shift: %.5lf\n", dz);
	//printf("yx shear: %.6lf\n", yx_shear);
	//printf("zx shear: %.6lf\n", zx_shear);
	//printf("zy shear: %.6lf\n", zy_shear);

	//std::vector< matrix_3D::vector_3D >point,point_after;
	//for (int i = 0; i < 3; i++) {
	//	matrix_3D::vector_3D p;
	//	p.x = corr_p[i]->x;
	//	p.y = corr_p[i]->y;
	//	p.z = corr_p[i]->z;
	//	point.push_back(p);
	//	p.x = corr_p[i]->x + corr_p[i]->dx;
	//	p.y = corr_p[i]->y + corr_p[i]->dy;
	//	p.z = corr_p[i]->z + corr_p[i]->dz;
	//	point_after.push_back(p);
	//}
	//trans_9para(point, *this);
	//for (auto p : point_after) {
	//	printf("x:%10.1lf y:%10.1lf z:%10.1lf\n", p.x, p.y, p.z);
	//}
}

void align_param2::Calc_9param_shift() {


	double bp[3][3], ap[3][3], cos_z, sin_z;
	for (int i = 0; i < 3; i++) {
		bp[i][0] = corr_p[i]->x;
		bp[i][1] = corr_p[i]->y;
		bp[i][2] = corr_p[i]->z;

		ap[i][0] = corr_p[i]->x + corr_p[i]->dx;
		ap[i][1] = corr_p[i]->y + corr_p[i]->dy;
		ap[i][2] = corr_p[i]->z + corr_p[i]->dz;

		//printf("bp%d %8.1lf %8.1lf %8.1lf\n",i, bp[i][0], bp[i][1], bp[i][2]);
		//printf("ap%d %8.1lf %8.1lf %8.1lf\n", i,ap[i][0], ap[i][1], ap[i][2]);
	}
	//apの位置ずれvectorを定義
	double dp[2][3];
	for (int i = 0; i < 3; i++) {
		dp[0][i] = ap[1][i] - ap[0][i];
		dp[1][i] = ap[2][i] - ap[0][i];
	}
	//printf("0-->1 x,y,z : %.1lf %.1lf %.1lf\n", dp[0][0], dp[0][1], dp[0][2]);
	//printf("0-->2 x,y,z : %.1lf %.1lf %.1lf\n", dp[1][0], dp[1][1], dp[1][2]);
	//法線vector
	double n_v[3];
	n_v[0] = (dp[0][1] * dp[1][2] - dp[0][2] * dp[1][1]);
	n_v[1] = (dp[0][2] * dp[1][0] - dp[0][0] * dp[1][2]);
	n_v[2] = (dp[0][0] * dp[1][1] - dp[0][1] * dp[1][0]);

	//std::cout << "normal vector" << std::endl;
	//for (int i = 0; i < 3; i++) {
	//	std::cout << std::setw(14) << std::fixed << std::setprecision(10) << n_v[i] << std::endl;
	//}

	x_rot = atan(n_v[1] / n_v[2]);
	n_v[1] = cos(x_rot) * n_v[1] - sin(x_rot) * n_v[2];
	n_v[2] = sin(x_rot) * n_v[1] + cos(x_rot) * n_v[2];
	//std::cout << "normal vector" << std::endl;
	//for (int i = 0; i < 3; i++) {
	//	std::cout << std::setw(14) << std::fixed << std::setprecision(10) << n_v[i] << std::endl;
	//}
	y_rot = atan(-1 * n_v[0] / n_v[2]);
	n_v[0] = cos(y_rot) * n_v[0] + sin(y_rot) * n_v[2];
	n_v[2] = -1 * sin(y_rot) * n_v[0] + cos(y_rot) * n_v[2];

	//std::cout << "normal vector" << std::endl;
	//for (int i = 0; i < 3; i++) {
	//	std::cout << std::setw(14) << std::fixed << std::setprecision(10) << n_v[i] << std::endl;
	//}

	//printf("x rot:%.6lf\n", x_rot);
	//printf("y rot:%.6lf\n", y_rot);


	matrix_3D::matrix_33 x_rot_mat(0, x_rot), y_rot_mat(1, y_rot);
	matrix_3D::vector_3D ap_v[3];
	for (int i = 0; i < 3; i++) {
		ap_v[i].x = ap[i][0];
		ap_v[i].y = ap[i][1];
		ap_v[i].z = ap[i][2];
	}
	for (int i = 0; i < 3; i++) {
		ap_v[i].matrix_multiplication(x_rot_mat);
		ap_v[i].matrix_multiplication(y_rot_mat);
	}
	//for (int i = 0; i < 3; i++) {
	//	printf("point %d\n", i);
	//	printf("\t %.2lf %.2lf %.2lf\n", bp[i][0], bp[i][1], bp[i][2]);
	//	printf("\t %.2lf %.2lf %.2lf\n", ap_v[i].x, ap_v[i].y, ap_v[i].z);
	//}
	dz = (ap_v[0].z - bp[0][2] + ap_v[1].z - bp[1][2] + ap_v[2].z - bp[2][2]) / 3;
	//printf("dz=%.2lf\n", dz);
	//3元方程式を解く
	double a[2][3][3] = { { {bp[0][0],bp[0][1],1},{bp[1][0],bp[1][1],1},{bp[2][0],bp[2][1],1} },  { {bp[0][0],bp[0][1],1},{bp[1][0],bp[1][1],1},{bp[2][0],bp[2][1],1} } };
	double b[2][3] = { { ap_v[0].x,ap_v[1].x,ap_v[2].x },{ ap_v[0].y,ap_v[1].y,ap_v[2].y } };
	double c[2][3] = { {1, 1, 1},{1,1,1} };
	//gauss(a[0], b[0], c[0]);
	//gauss(a[1], b[1], c[1]);
	GaussJorden(a[0], b[0], c[0]);
	GaussJorden(a[1], b[1], c[1]);
	z_rot = atan(c[1][0] / c[0][0]);
	x_shrink = c[0][0] / cos(z_rot);
	y_shrink = (c[0][0] * c[1][1] - c[0][1] * c[1][0]) / (c[0][0] * cos(z_rot) + c[1][0] * sin(z_rot));
	yx_shear = (c[0][1] * cos(z_rot) + c[1][1] * sin(z_rot)) / (c[0][0] * cos(z_rot) + c[1][0] * sin(z_rot));

	dx = c[0][2];
	dy = c[1][2];
	matrix_3D::vector_3D dr;
	dr.x = c[0][2];
	dr.y = c[1][2];
	dr.z = dz;

	x_rot = x_rot * -1;
	y_rot = y_rot * -1;
	matrix_3D::matrix_33 x_rot_mat_inv(0, x_rot), y_rot_mat_inv(1, y_rot);


	dr.matrix_multiplication(y_rot_mat_inv);
	dr.matrix_multiplication(x_rot_mat_inv);

	dx = dr.x;
	dy = dr.y;
	dz = dr.z;

	//printf("x rot: %.6lf\n",x_rot);
	//printf("y rot: %.6lf\n",y_rot);
	//printf("z rot: %.6lf\n",z_rot);
	//printf("x shrink: %.6lf\n", x_shrink);
	//printf("y shrink: %.6lf\n", y_shrink);
	//printf("z shrink: %.6lf\n", z_shrink);
	//printf("x shift: %.5lf\n", dx);
	//printf("y shift: %.5lf\n", dy);
	//printf("z shift: %.5lf\n", dz);
	//printf("yx shear: %.6lf\n", yx_shear);
	//printf("zx shear: %.6lf\n", zx_shear);
	//printf("zy shear: %.6lf\n", zy_shear);

	//std::vector< matrix_3D::vector_3D >point,point_after;
	//for (int i = 0; i < 3; i++) {
	//	matrix_3D::vector_3D p;
	//	p.x = corr_p[i]->x;
	//	p.y = corr_p[i]->y;
	//	p.z = corr_p[i]->z;
	//	point.push_back(p);
	//	p.x = corr_p[i]->x + corr_p[i]->dx;
	//	p.y = corr_p[i]->y + corr_p[i]->dy;
	//	p.z = corr_p[i]->z + corr_p[i]->dz;
	//	point_after.push_back(p);
	//}
	//trans_9para(point, *this);
	//for (auto p : point_after) {
	//	printf("x:%10.1lf y:%10.1lf z:%10.1lf\n", p.x, p.y, p.z);
	//}
}
//
//void GaussJorden(double in[3][3], double b[3], double c[3]) {
//
//
//	double a[3][4];
//	for (int i = 0; i < 3; i++) {
//		for (int j = 0; j < 4; j++) {
//			if (j < 3) {
//				a[i][j] = in[i][j];
//			}
//			else {
//				a[i][j] = b[i];
//			}
//		}
//	}
//	int N = 3;
//	double p, d;         // ピボット係数、ピボット行ｘ係数
//	double max, dummy;   // 最大絶対値、入れ替え時ダミー
//	int s;
//
//	//元の連立方程式をコンソール出力
//   //for (int i = 0; i < N; i++) {
//   //	for (int j = 0; j < N; j++)
//   //		printf("%+fx%d ", a[i][j], j + 1);
//   //	printf("= %+f\n", a[i][N]);
//   //}
//
//	for (int k = 0; k < N; k++) {
//		// 行入れ替え
//		max = 0; s = k;
//		for (int j = k; j < N; j++) {
//			if (fabs(a[j][k]) > max) {
//				max = fabs(a[j][k]);
//				s = j;
//			}
//		}
//		if (max == 0) {
//			printf("解けない！");
//			exit(1);
//		}
//		for (int j = 0; j <= N; j++) {
//			dummy = a[k][j];
//			a[k][j] = a[s][j];
//			a[s][j] = dummy;
//		}
//
//		// ピボット係数
//		p = a[k][k];
//
//		// ピボット行を p で除算
//		for (int j = k; j < N + 1; j++)
//			a[k][j] /= p;
//
//		// ピボット列の掃き出し
//		for (int i = 0; i < N; i++) {
//			if (i != k) {
//				d = a[i][k];
//				for (int j = k; j < N + 1; j++)
//					a[i][j] -= d * a[k][j];
//			}
//		}
//	}
//
//	// 結果出力
//	for (int k = 0; k < N; k++) {
//		c[k] = a[k][N];
//		//printf("x%d = %f\n", k + 1, a[k][N]);
//	}
//}
//
std::vector< align_param> search_param(vxx::base_track_t& base, std::vector<align_param>& param) {
	std::vector<align_param> param_shift;

	for (auto itr = param.begin(); itr != param.end(); itr++) {
		double x_tmp, y_tmp, z_tmp;
		x_tmp = itr->x;
		y_tmp = itr->y;
		z_tmp = itr->z;
		double  cos_z, sin_z;

		cos_z = cos(itr->z_rot);
		sin_z = sin(itr->z_rot);

		itr->x = itr->x_shrink * cos_z * (x_tmp)-itr->y_shrink * sin_z * (y_tmp)+itr->dx;
		itr->y = itr->x_shrink * sin_z * (x_tmp)+itr->y_shrink * cos_z * (y_tmp)+itr->dy;
		itr->z = z_tmp + itr->dz;

		itr->dx = x_tmp - itr->x;
		itr->dy = y_tmp - itr->y;
		itr->dz = z_tmp - itr->z;
		//param.x+param.dx --> 元のx

		param_shift.push_back(*itr);
	}


	std::vector <align_param2 >param2 = DelaunayDivide(param_shift);
	std::pair<vxx::base_track_t, align_param2>param_pair = track_affineparam_correspondence(base, param2);
	//printf("%d %d %d\n", param_pair.second.x, param_pair.second.y, param_pair.second.z);

	//for (int i = 0; i < 3; i++) {
	//	printf("%3d %3d %8.1lf %8.1lf %8.1lf %8.4lf %8.4lf %8.4lf\n",
	//		param_pair.second.corr_p[i]->ix,
	//		param_pair.second.corr_p[i]->iy,
	//		param_pair.second.corr_p[i]->x,
	//		param_pair.second.corr_p[i]->y,
	//		param_pair.second.corr_p[i]->z,
	//		param_pair.second.corr_p[i]->dx,
	//		param_pair.second.corr_p[i]->dy,
	//		param_pair.second.corr_p[i]->dz);
	//}
	align_param2 return_param;
	std::vector< align_param> ret;
	ret.push_back(*(param_pair.second.corr_p[0]));
	ret.push_back(*(param_pair.second.corr_p[1]));
	ret.push_back(*(param_pair.second.corr_p[2]));
	return_param = (param_pair.second);
	for (int i = 0; i < 3; i++) {
		return_param.corr_p[i]->dx = param_pair.second.corr_p[i]->dx;
		return_param.corr_p[i]->dy = param_pair.second.corr_p[i]->dy;
		return_param.corr_p[i]->dz = param_pair.second.corr_p[i]->dz;
		return_param.corr_p[i]->id = param_pair.second.corr_p[i]->id;
		return_param.corr_p[i]->ix = param_pair.second.corr_p[i]->ix;
		return_param.corr_p[i]->iy = param_pair.second.corr_p[i]->iy;
		return_param.corr_p[i]->signal = param_pair.second.corr_p[i]->signal;
		return_param.corr_p[i]->x = param_pair.second.corr_p[i]->x;
		return_param.corr_p[i]->y = param_pair.second.corr_p[i]->y;
		return_param.corr_p[i]->z = param_pair.second.corr_p[i]->z;
		return_param.corr_p[i]->x_rot = param_pair.second.corr_p[i]->x_rot;
		return_param.corr_p[i]->y_rot = param_pair.second.corr_p[i]->y_rot;
		return_param.corr_p[i]->z_rot = param_pair.second.corr_p[i]->z_rot;
		return_param.corr_p[i]->x_shrink = param_pair.second.corr_p[i]->x_shrink;
		return_param.corr_p[i]->y_shrink = param_pair.second.corr_p[i]->y_shrink;
		return_param.corr_p[i]->z_shrink = param_pair.second.corr_p[i]->z_shrink;
		return_param.corr_p[i]->yx_shear = param_pair.second.corr_p[i]->yx_shear;
		return_param.corr_p[i]->zx_shear = param_pair.second.corr_p[i]->zx_shear;
		return_param.corr_p[i]->zy_shear = param_pair.second.corr_p[i]->zy_shear;

	}
	//return_param.corr_p[0] = param_pair.second.corr_p[0];
	//return_param.corr_p[1] = param_pair.second.corr_p[1];
	//return_param.corr_p[2] = param_pair.second.corr_p[2];

	//for (int i = 0; i < 3; i++) {
	//	printf("%3d %3d %8.1lf %8.1lf %8.1lf %8.4lf %8.4lf %8.4lf\n",
	//		return_param.corr_p[i]->ix,
	//		return_param.corr_p[i]->iy,
	//		return_param.corr_p[i]->x,
	//		return_param.corr_p[i]->y,
	//		return_param.corr_p[i]->z,
	//		return_param.corr_p[i]->dx,
	//		return_param.corr_p[i]->dy,
	//		return_param.corr_p[i]->dz);
	//}

	return ret;
}

//basetrack-alignment mapの対応
std::pair<vxx::base_track_t, align_param2>track_affineparam_correspondence(vxx::base_track_t& base, std::vector <align_param2>& param) {

	//local alignの視野中心を取り出して、位置でhash
	//local alignの視野中心の作るdelaunay三角形をmapで対応

	std::map<int, align_param*> view_center;
	std::multimap<int, align_param2>triangles;
	double xmin = 999999, ymin = 999999, hash = 2000;
	for (auto itr = param.begin(); itr != param.end(); itr++) {
		for (int i = 0; i < 3; i++) {
			view_center.insert(std::make_pair(itr->corr_p[i]->id, (itr->corr_p[i])));
			triangles.insert(std::make_pair(itr->corr_p[i]->id, (*itr)));
			xmin = std::min(itr->corr_p[i]->x, xmin);
			ymin = std::min(itr->corr_p[i]->y, ymin);
		}
	}
	std::multimap<std::pair<int, int>, align_param*> view_center_hash;
	std::pair<int, int>id;
	for (auto itr = view_center.begin(); itr != view_center.end(); itr++) {
		id.first = int((itr->second->x - xmin) / hash);
		id.second = int((itr->second->y - ymin) / hash);
		view_center_hash.insert(std::make_pair(id, itr->second));
	}

	std::vector < std::pair<vxx::base_track_t*, align_param2>> ret;
	std::vector<align_param*> param_cand;
	int loop = 0, ix, iy, count = 0;
	ix = (base.x - xmin) / hash;
	iy = (base.x - ymin) / hash;
	loop = 1;
	while (true) {
		param_cand.clear();
		for (int iix = ix - loop; iix <= ix + loop; iix++) {
			for (int iiy = iy - loop; iiy <= iy + loop; iiy++) {
				id.first = iix;
				id.second = iiy;
				if (view_center_hash.count(id) != 0) {
					auto range = view_center_hash.equal_range(id);
					for (auto res = range.first; res != range.second; res++) {
						param_cand.push_back(res->second);
					}
				}
			}
		}
		if (param_cand.size() > 2)break;
		loop++;
	}
	align_param2 param2 = search_param(param_cand, base, triangles);
	return std::make_pair(base, param2);

}
align_param2 search_param(std::vector<align_param*>& param, vxx::base_track_t& base, std::multimap<int, align_param2>& triangles) {
	//三角形内部
	//最近接三角形
	double dist = 0;
	std::map<double, align_param* > dist_map;
	//align_paramを近い順にsort
	for (auto itr = param.begin(); itr != param.end(); itr++) {
		dist = ((*itr)->x - base.x) * ((*itr)->x - base.x) + ((*itr)->y - base.y) * ((*itr)->y - base.y);
		dist_map.insert(std::make_pair(dist, (*itr)));
	}

	double sign[3];
	bool flg = false;
	int id;

	align_param2 ret = triangles.begin()->second;
	for (auto itr = dist_map.begin(); itr != dist_map.end(); itr++) {
		if (itr != dist_map.begin())continue;


		//corrmapのID
		id = itr->second->id;
		if (triangles.count(id) == 0) {
			fprintf(stderr, "alignment triangle ID=%d not found\n", id);
			exit(1);
		}
		//idの属する三角形を探索
		auto range = triangles.equal_range(id);
		for (auto itr2 = range.first; itr2 != range.second; itr2++) {
			sign[0] = (itr2->second.corr_p[1]->x - itr2->second.corr_p[0]->x) * (base.y - itr2->second.corr_p[1]->y) - (itr2->second.corr_p[1]->y - itr2->second.corr_p[0]->y) * (base.x - itr2->second.corr_p[1]->x);
			sign[1] = (itr2->second.corr_p[2]->x - itr2->second.corr_p[1]->x) * (base.y - itr2->second.corr_p[2]->y) - (itr2->second.corr_p[2]->y - itr2->second.corr_p[1]->y) * (base.x - itr2->second.corr_p[2]->x);
			sign[2] = (itr2->second.corr_p[0]->x - itr2->second.corr_p[2]->x) * (base.y - itr2->second.corr_p[0]->y) - (itr2->second.corr_p[0]->y - itr2->second.corr_p[2]->y) * (base.x - itr2->second.corr_p[0]->x);
			//printf("point %.lf,%.1lf\n", base.x, base.y);
			//printf("triangle %.1lf %.1lf\n", itr2->second->corr_p[0]->x, itr2->second->corr_p[0]->y);
			//printf("triangle %.1lf %.1lf\n", itr2->second->corr_p[1]->x, itr2->second->corr_p[1]->y);
			//printf("triangle %.1lf %.1lf\n", itr2->second->corr_p[2]->x, itr2->second->corr_p[2]->y);
			//printf("sign %.1lf %1.lf %.1lf\n", sign[0], sign[1], sign[2]);
			//printf("  signbit %d %d %d\n", std::signbit(sign[0]), std::signbit(sign[1]), std::signbit(sign[2]));
			//printf("n signbit %d %d %d\n", !std::signbit(sign[0]), !std::signbit(sign[1]), !std::signbit(sign[2]));
			//printf("judge %d\n", (std::signbit(sign[0]) && std::signbit(sign[1]) && std::signbit(sign[2])) || (!std::signbit(sign[0]) && !std::signbit(sign[1]) && !std::signbit(sign[2])));
			//printf("\n");

			//符号が3つとも一致でtrue
			if ((std::signbit(sign[0]) && std::signbit(sign[1]) && std::signbit(sign[2])) || (!std::signbit(sign[0]) && !std::signbit(sign[1]) && !std::signbit(sign[2]))) {
				ret = itr2->second;
				flg = true;
				break;
			}
		}
		if (flg)break;
	}
	if (flg) {
		//printf("point in trianlge\n");
		return ret;
	}

	//distが最小になるcorrmapをとってくる
	dist = -1;
	for (auto itr = dist_map.begin(); itr != dist_map.end(); itr++) {
		//corrmapのID
		id = itr->second->id;
		if (triangles.count(id) == 0) {
			fprintf(stderr, "alignment triangle ID=%d not found\n", id);
			exit(1);
		}
		//idの属する三角形を探索
		auto range = triangles.equal_range(id);
		for (auto itr2 = range.first; itr2 != range.second; itr2++) {
			if (dist<0 || dist>select_triangle_vale(itr2->second, base)) {
				dist = select_triangle_vale(itr2->second, base);
				ret = itr2->second;
			}
		}
	}
	//printf("point not in trianlge\n");
	return ret;
}
double select_triangle_vale(align_param2 param, vxx::base_track_t& base) {
	double x, y;
	double dist = 0;
	x = (param.corr_p[0]->x + param.corr_p[1]->x + param.corr_p[2]->x) / 3;
	y = (param.corr_p[0]->y + param.corr_p[1]->y + param.corr_p[2]->y) / 3;
	dist = (base.x - x) * (base.x - x) + (base.y - y) * (base.y - y);
	return dist;
}

std::vector<vxx::base_track_t> extra_track_all(std::vector<vxx::base_track_t>& base, align_param2& param) {
	std::vector<vxx::base_track_t> ret;
	for (auto itr = base.begin(); itr != base.end(); itr++) {
		ret.push_back(extra_track(*itr, param));
	}
	return ret;
}
vxx::base_track_t extra_track(vxx::base_track_t base, align_param2& param) {
	double n[3];
	base.z = 0;
	n[0] = (param.corr_p[1]->y - param.corr_p[0]->y) * (param.corr_p[2]->z - param.corr_p[0]->z)
		- (param.corr_p[2]->y - param.corr_p[0]->y) * (param.corr_p[1]->z - param.corr_p[0]->z);
	n[1] = (param.corr_p[1]->z - param.corr_p[0]->z) * (param.corr_p[2]->x - param.corr_p[0]->x)
		- (param.corr_p[2]->z - param.corr_p[0]->z) * (param.corr_p[1]->x - param.corr_p[0]->x);
	n[2] = (param.corr_p[1]->x - param.corr_p[0]->x) * (param.corr_p[2]->y - param.corr_p[0]->y)
		- (param.corr_p[2]->x - param.corr_p[0]->x) * (param.corr_p[1]->y - param.corr_p[0]->y);

	//printf("normal %g %g %g\n", n[0], n[1], n[2]);
	double extra_dist;
	extra_dist = (n[0] * (param.corr_p[0]->x - base.x) + n[1] * (param.corr_p[0]->y - base.y) + n[2] * (param.corr_p[0]->z - base.z))
		/
		(n[0] * base.ax + n[1] * base.ay + n[2]);

	//printf("ex dist %g\n", extra_dist);
	vxx::base_track_t extra_base = base;
	extra_base.x = base.x + extra_dist * base.ax;
	extra_base.y = base.y + extra_dist * base.ay;
	extra_base.z = base.z + extra_dist;
	//printf("ex base %g %g %g\n", extra_base.x, extra_base.y, extra_base.z);

	//for (int i = 0; i < 3; i++) {
	//	printf("%g %g %g %g %g %g\n"
	//		, param.corr_p[i]->x
	//		, param.corr_p[i]->y
	//		, param.corr_p[i]->z
	//		, param.corr_p[i]->dx
	//		, param.corr_p[i]->dy
	//		, param.corr_p[i]->dz);
	//}
//逆変換で戻す
	param.Calc_9param_shift();

	trans_base(extra_base, &param);
	//printf("ex base inv %g %g %g\n", extra_base.x, extra_base.y, extra_base.z);

	return extra_base;
}
void trans_base(vxx::base_track_t& base, align_param2* param) {

	matrix_3D::matrix_33 x_rot_mat(0, param->x_rot), y_rot_mat(1, param->y_rot), z_rot_mat(2, param->z_rot), all_trans(0, 0), shear_mat(0, 0), shrink_mat(0, 0);

	shrink_mat.val[0][0] *= param->x_shrink;
	shrink_mat.val[1][1] *= param->y_shrink;
	//shrink_mat.val[2][2] *= param->z_shrink;
	shear_mat.val[0][1] = param->yx_shear;
	shear_mat.val[0][2] = param->zx_shear;
	shear_mat.val[1][2] = param->zy_shear;

	matrix_3D::vector_3D shift, center;
	center.x = param->x;
	center.y = param->y;
	center.z = param->z;
	shift.x = param->dx;
	shift.y = param->dy;
	shift.z = param->dz;

	all_trans.matrix_multiplication(shear_mat);
	all_trans.matrix_multiplication(shrink_mat);
	all_trans.matrix_multiplication(z_rot_mat);
	all_trans.matrix_multiplication(y_rot_mat);
	all_trans.matrix_multiplication(x_rot_mat);

	//all_trans.Print();
	matrix_3D::vector_3D base_p0, base_p1;
	double base_nominal_thick = 210;

	param->z_shrink = (param->corr_p[0]->z_shrink + param->corr_p[1]->z_shrink + param->corr_p[2]->z_shrink) / 3;
	base_p0.x = base.x;
	base_p0.y = base.y;
	base_p0.z = base.z;

	base_p1.x = base.x + base.ax * base_nominal_thick * param->z_shrink;
	base_p1.y = base.y + base.ay * base_nominal_thick * param->z_shrink;
	base_p1.z = base.z + base_nominal_thick * param->z_shrink;
	//base_p0 = matrix_3D::addition(base_p0, matrix_3D::const_multiple(center, -1));
	//base_p1 = matrix_3D::addition(base_p1, matrix_3D::const_multiple(center, -1));

	//変換の実行
	base_p0.matrix_multiplication(all_trans);
	base_p0 = matrix_3D::addition(base_p0, shift);
	base_p1.matrix_multiplication(all_trans);
	base_p1 = matrix_3D::addition(base_p1, shift);

	//原点をもとに戻す
	//base_p0 = matrix_3D::addition(base_p0, center);
	//base_p1 = matrix_3D::addition(base_p1, center);

	base.x = base_p0.x;
	base.y = base_p0.y;
	base.z = base_p0.z;

	//printf("ax:%.4lf --> %.4lf\n", (*itr)->ax, (base_p0.x - base_p1.x) / (base_p0.z - base_p1.z));
	//printf("ay:%.4lf --> %.4lf\n", (*itr)->ay, (base_p0.y - base_p1.y) / (base_p0.z - base_p1.z));

	base.ax = (base_p0.x - base_p1.x) / (base_p0.z - base_p1.z);
	base.ay = (base_p0.y - base_p1.y) / (base_p0.z - base_p1.z);


}

std::vector<output_format_link> read_linket(std::string filename) {
	std::ifstream ifs(filename, std::ios::binary);
	//filesize取得
	ifs.seekg(0, std::ios::end);
	int64_t eofpos = ifs.tellg();
	ifs.clear();
	ifs.seekg(0, std::ios::beg);
	int64_t begpos = ifs.tellg();
	int64_t nowpos = ifs.tellg();
	int64_t size2 = eofpos - begpos;
	int64_t GB = size2 / (1000 * 1000 * 1000);
	int64_t MB = (size2 - GB * 1000 * 1000 * 1000) / (1000 * 1000);
	int64_t KB = (size2 - GB * 1000 * 1000 * 1000 - MB * 1000 * 1000) / (1000);
	if (GB > 0) {
		std::cout << "FILE size :" << GB << "." << MB << " [GB]" << std::endl;
	}
	else {
		std::cout << "FILE size :" << MB << "." << KB << " [MB]" << std::endl;
	}
	std::vector<output_format_link> ret;
	int64_t count = 0;
	output_format_link l;
	while (ifs.read((char*)&l, sizeof(output_format_link))) {
		if (count % 10000 == 0) {
			nowpos = ifs.tellg();
			auto size1 = nowpos - begpos;
			std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
		}
		count++;
		ret.push_back(l);
	}
	return ret;

}

std::vector<output_format_link> match_baseid_extract(std::vector<vxx::base_track_t>& base, std::vector<output_format_link>& link, std::vector<std::pair<int, int>>& rawid_connect) {
	std::vector<output_format_link> ret;
	std::set<int> rawid;
	for (auto itr = base.begin(); itr != base.end(); itr++) {
		rawid.insert(itr->rawid);
	}
	for (auto itr = link.begin(); itr != link.end(); itr++) {
		if (rawid.count(itr->b[0].rawid) == 0)continue;
		rawid_connect.push_back(std::make_pair(itr->b[0].rawid, itr->b[1].rawid));
		ret.push_back(*itr);
	}
	return ret;
}

void output_file(std::string filename, vxx::base_track_t& ex_base, std::vector<vxx::base_track_t>& ali_extra, std::vector<vxx::base_track_t>& ali_connect, std::vector<std::pair<int, int>>& rawid_connect) {

	std::ofstream ofs(filename.c_str());


	ofs << std::right << std::fixed
		<< std::setw(3) << std::setprecision(0) << ex_base.pl << " "
		<< std::setw(10) << std::setprecision(0) << ex_base.rawid << " "
		<< std::setw(7) << std::setprecision(0) << ex_base.m[0].ph << " "
		<< std::setw(7) << std::setprecision(0) << ex_base.m[1].ph << " "
		<< std::setw(8) << std::setprecision(4) << ex_base.ax << " "
		<< std::setw(8) << std::setprecision(4) << ex_base.ay << " "
		<< std::setw(8) << std::setprecision(1) << ex_base.x << " "
		<< std::setw(8) << std::setprecision(1) << ex_base.y << std::endl;

	std::multimap<int, int>rawid_connect_map;
	for (auto itr = rawid_connect.begin(); itr != rawid_connect.end(); itr++) {
		rawid_connect_map.insert(*itr);
	}
	std::map<int, vxx::base_track_t> connect_base_map;
	for (auto itr = ali_connect.begin(); itr != ali_connect.end(); itr++) {
		connect_base_map.insert(std::make_pair(itr->rawid, *itr));
	}

	std::vector<std::pair<vxx::base_track_t, vxx::base_track_t >> ali;
	for (auto itr = ali_extra.begin(); itr != ali_extra.end(); itr++) {
		if (rawid_connect_map.count(itr->rawid) == 0)continue;
		auto range = rawid_connect_map.equal_range(itr->rawid);
		for (auto res = range.first; res != range.second; res++) {
			if (connect_base_map.count(res->second) == 0)continue;
			auto find_track = connect_base_map.find(res->second);
			ali.push_back(std::make_pair(*itr, find_track->second));
		}
	}
	int count = 0;
	for (auto itr = ali.begin(); itr != ali.end(); itr++) {
		if (count > 10)continue;
		count++;
		ofs << std::right << std::fixed << std::endl;
		ofs << std::right << std::fixed
			<< std::setw(3) << std::setprecision(0) << itr->first.pl << " "
			<< std::setw(10) << std::setprecision(0) << itr->first.rawid << " "
			<< std::setw(7) << std::setprecision(0) << itr->first.m[0].ph << " "
			<< std::setw(7) << std::setprecision(0) << itr->first.m[1].ph << " "
			<< std::setw(8) << std::setprecision(4) << itr->first.ax << " "
			<< std::setw(8) << std::setprecision(4) << itr->first.ay << " "
			<< std::setw(8) << std::setprecision(1) << itr->first.x << " "
			<< std::setw(8) << std::setprecision(1) << itr->first.y << std::endl;
		ofs << std::right << std::fixed
			<< std::setw(3) << std::setprecision(0) << itr->second.pl << " "
			<< std::setw(10) << std::setprecision(0) << itr->second.rawid << " "
			<< std::setw(7) << std::setprecision(0) << itr->second.m[0].ph << " "
			<< std::setw(7) << std::setprecision(0) << itr->second.m[1].ph << " "
			<< std::setw(8) << std::setprecision(4) << itr->second.ax << " "
			<< std::setw(8) << std::setprecision(4) << itr->second.ay << " "
			<< std::setw(8) << std::setprecision(1) << itr->second.x << " "
			<< std::setw(8) << std::setprecision(1) << itr->second.y << std::endl;
	}
}
void output_file(std::string filename, vxx::base_track_t& base, std::vector<vxx::base_track_t>& ali) {

	std::ofstream ofs(filename.c_str());


	ofs << std::right << std::fixed
		<< std::setw(3) << std::setprecision(0) << base.pl << " "
		<< std::setw(10) << std::setprecision(0) << base.rawid << " "
		<< std::setw(7) << std::setprecision(0) << base.m[0].ph << " "
		<< std::setw(7) << std::setprecision(0) << base.m[1].ph << " "
		<< std::setw(8) << std::setprecision(4) << base.ax << " "
		<< std::setw(8) << std::setprecision(4) << base.ay << " "
		<< std::setw(8) << std::setprecision(1) << base.x << " "
		<< std::setw(8) << std::setprecision(1) << base.y << std::endl << std::endl;

	int count = 0;
	for (auto itr = ali.begin(); itr != ali.end(); itr++) {
		if (count > 10)continue;
		count++;
		ofs << std::right << std::fixed
			<< std::setw(3) << std::setprecision(0) << itr->pl << " "
			<< std::setw(10) << std::setprecision(0) << itr->rawid << " "
			<< std::setw(7) << std::setprecision(0) << itr->m[0].ph << " "
			<< std::setw(7) << std::setprecision(0) << itr->m[1].ph << " "
			<< std::setw(8) << std::setprecision(4) << itr->ax << " "
			<< std::setw(8) << std::setprecision(4) << itr->ay << " "
			<< std::setw(8) << std::setprecision(1) << itr->x << " "
			<< std::setw(8) << std::setprecision(1) << itr->y << std::endl;
	}
}