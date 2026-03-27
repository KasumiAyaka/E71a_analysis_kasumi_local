// 2025/12/3
// kasumi

#pragma comment(lib,"FILE_structure.lib")
#include <cstdint> 
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>
#include <FILE_structure.hpp>

struct Key {
    int pl, rid;
};
bool operator<(const Key& lhs, const Key& rhs) {
    return std::tie(lhs.pl, lhs.rid) < std::tie(rhs.pl, rhs.rid);
}

struct UnitList {
    int ecc, gid, cid, venum;
    Key k0,k1;
};
bool operator<(const UnitList& lhs, const UnitList& rhs) {
    return std::tie(lhs.ecc, lhs.gid, lhs.cid, lhs.venum,lhs.k0, lhs.k1) < std::tie(rhs.ecc, rhs.gid, rhs.cid, rhs.venum, rhs.k0, rhs.k1);
}




class corrmap {
public:

    int id, pos[2], ix, iy;
    double x_area[2], y_area[2], position[6], angle[6], dz, signal, background, sn;
    double rms_x, rms_y, rms_ax, rms_ay, cx, cy;


};

class output_base {
public:
    int rawid;
    double ax, ay, x, y;
};
class output_format {
public:
    bool hit_flg;
    output_base ex, aft, bef;
};

struct Material {
    static const double water;
    static const double iron;
    static const double film;
    static const double pack;
    static const double angle_shrink;
};
struct EventTrackInfo {
    int gid, cid, vpl, pl, rid;
};

std::vector<vxx::base_track_t> read_base_id(std::string file_in_ECC, int pl, std::set<int>rawid);
std::vector<vxx::base_track_t> read_base_id(std::string file_in_ECC, int pl);
std::vector<corrmap>read_corrmap(std::string filename);
void calc_invers_corr_area(std::vector<corrmap>& corr);
std::vector<std::pair<vxx::base_track_t, corrmap>> correspond_corrmap_fwd(std::vector<vxx::base_track_t>& base, std::vector < corrmap>& corr);
std::vector<std::pair<vxx::base_track_t, corrmap>> correspond_corrmap_bwd(std::vector<vxx::base_track_t>& base, std::vector < corrmap>& corr);

const double Material::water = 2300.0;
const double Material::iron = 500.0;
const double Material::film = 350.0;
const double Material::pack = 109.0;
const double Material::angle_shrink = 0.989;

void set_basetrack(std::vector<Momentum_recon::Event_information>& momch, std::set<int>& eid, std::multimap<int, EventTrackInfo>& eve);
void extraporate_fwd_trk(int pl0, int pl1, vxx::base_track_t& base, corrmap& corr);
void extraporate_bwd_trk(int pl0, int pl1, vxx::base_track_t& base, corrmap& corr);
std::vector<UnitList> Read_VE_list(std::string input);
void base_matching(std::vector<vxx::base_track_t>& base0, std::vector<vxx::base_track_t>& bese_ref,
    double pos_lat_all[2], double pos_rad_all[2], double ang_lat_all[2], double ang_rad_all[2], std::vector<output_format>& out_base);
void output_basetracks(std::string filename, std::vector<output_format>& out);

int main(int argc, char** argv) {
    if (argc < 5) {
        fprintf(stderr, "usage:prg event.momch ECC_Area_path corrmap-path_ve output-path \n");
        fprintf(stderr, "usage:prg event.momch ECC_Area_path corrmap-path_ve output-path #ecc(except ecc5)\n");
        exit(1);
    }

    std::string input_txt, corr_path, corr_path_ve, output_path, ECC_Area_path, file_in_base[2], file_in_corr, file_out;

    input_txt = argv[1];
    ECC_Area_path = argv[2];
    corr_path_ve = argv[3];
    output_path = argv[4];
    int ecc = 5;
    if (argc == 6) {
        ecc = std::stoi(argv[5]);
    }

    corr_path = ECC_Area_path + "\\Area0\\0\\align";

    std::vector<UnitList> unit = Read_VE_list(input_txt);


    std::stringstream name, ss;
    std::vector<corrmap>corr;
    std::vector<std::pair<vxx::base_track_t, corrmap>>base_pair;
    std::vector<vxx::base_track_t>base_b[4];
    std::vector<vxx::base_track_t>base_f[4];
    std::vector<output_format>out;

    std::set<int>pllist;
    std::multimap<int, int>trklist;
    for (auto itr = unit.begin(); itr != unit.end(); itr++) {
        if (itr->k0.pl <= 0 && itr->k1.pl > 0) {
            pllist.insert(itr->k1.pl);
            trklist.insert(std::make_pair(itr->k1.pl, itr->k1.rid));
        }
        if (itr->k1.pl <= 0 && itr->k0.pl > 0) {
            pllist.insert(itr->k0.pl);
            trklist.insert(std::make_pair(itr->k0.pl, itr->k0.rid));
        }
    }

    double pos_lat_all[2];
    double pos_rad_all[2]; double ang_lat_all[2]; double ang_rad_all[2];
    
    //arrowane[0]+ [1]*rad
    pos_lat_all[0] = 30;
    pos_lat_all[1] = 0;
    pos_rad_all[0] = 40;
    pos_rad_all[1] = 30;
    ang_lat_all[0] = 0.05;
    ang_lat_all[1] = 0;
    ang_rad_all[0] = 0.1;
    ang_rad_all[1] = 0.1;
    //base1_trns --> 変換後
    //base[1] --> 変換前

    for (auto itr = pllist.begin(); itr != pllist.end(); itr++) {
        std::set<int>tkl;
        auto p = trklist.equal_range(*itr);
        for (auto q = p.first; q != p.second; q++) {
            tkl.insert(q->second);
        }

        if (*itr % 2 == 0) {
            // pl0-1 <-- pl0 pl1 への外挿:   32 5<--4 76 通常alignmentの
            base_b[0] = read_base_id(ECC_Area_path, *itr, tkl);
            ss << corr_path + "\\corrmap-align-" << std::setw(3) << std::setfill('0') << *itr << "-" << std::setw(3) << std::setfill('0') << *itr + 1 << ".lst";
            file_in_corr = ss.str();
            std::cout << file_in_corr << std::endl;
            corr = read_corrmap(file_in_corr);
            base_pair = correspond_corrmap_bwd(base_b[0], corr);
            for (auto itr0 = base_pair.begin(); itr0 != base_pair.end(); itr++) {
                extraporate_bwd_trk(*itr, *itr + 1, itr0->first, itr0->second);
                base_b[1].push_back(itr0->first);
            }
            ss.str("");
            corr.clear();
            base_pair.clear();

            //VE chk pl0-1 <-- pl0 pl1 への外挿:   32<--54 76 

            if (ecc == 5) {
                ss << corr_path_ve + "\\f_pl" << std::setw(3) << std::setfill('0') << *itr + 1 << "_" << std::setw(3) << std::setfill('0') << *itr - 2 << "\\corrmap-" << std::setw(3) << std::setfill('0') << *itr + 1 << "-" << std::setw(3) << std::setfill('0') << *itr - 2 << "-l.txt";
            }
            file_in_corr = ss.str();
            std::cout << file_in_corr << std::endl;
            corr = read_corrmap(file_in_corr);
            base_pair = correspond_corrmap_bwd(base_b[1], corr);
            for (auto itr0 = base_pair.begin(); itr0 != base_pair.end(); itr++) {
                extraporate_bwd_trk(*itr + 1, *itr - 2, itr0->first, itr0->second);
                base_b[2].push_back(itr0->first);
            }
            ss.str("");
            corr.clear();
            base_pair.clear();
            base_b[3] = read_base_id(ECC_Area_path, *itr - 2);

        }
        else {
            // pl0 pl1-->pl1+1 への外挿 : 32 5-->4 76 
            base_f[0] = read_base_id(ECC_Area_path, *itr, tkl);// pl0:downstream //5

            ss << corr_path + "\\corrmap-align-" << std::setw(3) << std::setfill('0') << *itr << "-" << std::setw(3) << std::setfill('0') << *itr - 1 << ".lst";
            file_in_corr = ss.str();
            std::cout << file_in_corr << std::endl;
            corr = read_corrmap(file_in_corr);
            base_pair = correspond_corrmap_fwd(base_f[0], corr);
            for (auto itr0 = base_pair.begin(); itr0 != base_pair.end(); itr0++) {
                extraporate_fwd_trk(*itr, *itr - 1, itr0->first, itr0->second); //5/4
                base_f[1].push_back(itr0->first);
            }
            ss.str("");
            corr.clear();
            base_pair.clear();

            // pl0 pl1-->pl1+1 への外挿 : 32 54--> 76 
            ss << corr_path_ve + "\\f_pl" << std::setw(3) << std::setfill('0') << *itr - 1 << "_" << std::setw(3) << std::setfill('0') << *itr + 3 << "\\corrmap-" << std::setw(3) << std::setfill('0') << *itr - 1 << "-" << std::setw(3) << std::setfill('0') << *itr+3 << "-l.txt";
            file_in_corr = ss.str();
            std::cout << file_in_corr << std::endl;
            corr = read_corrmap(file_in_corr);
            base_pair = correspond_corrmap_fwd(base_f[1], corr);
            for (auto itr0 = base_pair.begin(); itr0 != base_pair.end(); itr0++) {
                extraporate_fwd_trk(*itr, *itr + 2, itr0->first, itr0->second);//5,7
                base_f[2].push_back(itr0->first);
            }
            ss.str("");
            corr.clear();
            base_pair.clear();
            base_f[3] = read_base_id(ECC_Area_path, *itr + 2);
        }

        base_matching(base_f[2], base_f[3], pos_lat_all, pos_rad_all, ang_lat_all, ang_rad_all,out);
        base_matching(base_f[2], base_f[3], pos_lat_all, pos_rad_all, ang_lat_all, ang_rad_all, out);
    }
    output_basetracks(file_out, out);


}



std::vector<UnitList> Read_VE_list(std::string input) {
    std::vector<UnitList> ret;
    std::ifstream ifs(input);
    if (!ifs) {
        std::cerr << "File open error!" << std::endl;
        exit(0);
    }

    std::string str;					//1strein into
    std::vector<std::string> str_v;		//input 1 ward
    std::string buffer;

    UnitList u;
    while (std::getline(ifs, str)) {
        // gid cid pid nseg npl upl dpl
        str_v = StringSplit_with_tab(str);
        //std::cout << str << std::endl;
        u.ecc == std::stoi(str_v[1]);
            u.gid == std::stoi(str_v[2]);
        u.cid == std::stoi(str_v[3]);
        u.venum == std::stoi(str_v[4]);
        u.k0.pl == std::stoi(str_v[5]);
        u.k0.rid == std::stoi(str_v[6]);
        u.k1.pl == std::stoi(str_v[12]);
        u.k1.rid == std::stoi(str_v[13]);
        ret.push_back(u);
    }
    return ret;
}

std::vector<vxx::base_track_t> read_base_id(std::string file_in_ECC, int pl, std::set<int>rawid) {

    std::stringstream file_in_base;
    file_in_base << file_in_ECC << "\\Area0\\PL" << std::setw(3) << std::setfill('0') << pl
        << "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";

    std::vector<vxx::base_track_t> base_all;
    vxx::BvxxReader br;
    std::vector<vxx::CutArea> area;
    base_all = br.ReadAll(file_in_base.str(), pl, 0);
    std::vector<vxx::base_track_t> ret;
    for (auto itr = base_all.begin(); itr != base_all.end(); itr++) {
        if (rawid.count(itr->rawid) == 0)continue;
        itr->z = 0;
        ret.push_back(*itr);
    }
    return ret;
}
std::vector<vxx::base_track_t> read_base_id(std::string file_in_ECC, int pl) {

    std::stringstream file_in_base;
    file_in_base << file_in_ECC << "\\Area0\\PL" << std::setw(3) << std::setfill('0') << pl
        << "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";

    std::vector<vxx::base_track_t> base_all;
    vxx::BvxxReader br;
    std::vector<vxx::CutArea> area;
    base_all = br.ReadAll(file_in_base.str(), pl, 0);
    std::vector<vxx::base_track_t> ret;
    for (auto itr = base_all.begin(); itr != base_all.end(); itr++) {
        itr->z = 0;
        ret.push_back(*itr);
    }
    return ret;
}

std::vector<corrmap>read_corrmap(std::string filename) {
    std::vector<corrmap>ret;
    std::ifstream ifs(filename);

    corrmap corr;
    double buf_d[12];
    int count = 0;
    while (ifs >> corr.id >> corr.pos[0] >> corr.pos[1]
        >> corr.x_area[0] >> corr.x_area[1] >> corr.y_area[0] >> corr.y_area[1]
        >> corr.position[0] >> corr.position[1] >> corr.position[2] >> corr.position[3] >> corr.position[4] >> corr.position[5]
        >> corr.angle[0] >> corr.angle[1] >> corr.angle[2] >> corr.angle[3] >> corr.angle[4] >> corr.angle[5]
        >> corr.dz >> corr.signal >> corr.background >> corr.sn >> corr.rms_x >> corr.rms_y
        >> buf_d[0] >> buf_d[1]
        >> corr.rms_ax >> corr.rms_ay
        >> buf_d[2]
        >> corr.ix >> corr.iy
        >> buf_d[3] >> buf_d[4] >> buf_d[5] >> buf_d[6] >> buf_d[7] >> buf_d[8] >> buf_d[9] >> buf_d[10] >> buf_d[11]

        ) {
        if (count % 1000 == 0) {
            printf("\r read corrmap --> %d", count);

        }
        count++;
        ret.push_back(corr);
    }
    printf("\r read corrmap --> %d\n", count);
    return ret;
}
//ここまで読み込み

std::vector<std::pair<vxx::base_track_t, corrmap>> correspond_corrmap_fwd(std::vector<vxx::base_track_t>& base, std::vector <corrmap>& corr) {
    // VE alignmentで下流方向に外挿
    calc_invers_corr_area(corr);    // corrmapのareaはpos0系で書かれているため，pl1のtrkに適用するparamを選ぶための変換


    double x_min = 0, y_min = 0;
    for (int i = 0; i < corr.size(); i++) {
        if (i == 0) {
            x_min = corr[i].cx;
            y_min = corr[i].cy;
        }
        x_min = std::min(x_min, corr[i].cx);
        y_min = std::min(y_min, corr[i].cy);
    }

    double hash_size = 2000;
    std::multimap<std::pair<int, int>, corrmap> corr_hash;
    std::pair<int, int>id;
    for (int i = 0; i < corr.size(); i++) {
        id.first = (corr[i].cx - x_min) / hash_size;
        id.second = (corr[i].cy - y_min) / hash_size;
        corr_hash.insert(std::make_pair(id, corr[i]));
    }

    int ix, iy, loop_num;
    std::vector<std::pair<vxx::base_track_t, corrmap>>ret;
    std::vector<corrmap>corr_v;
    double dist;
    corrmap param;
    for (int i = 0; i < base.size(); i++) {
        corr_v.clear();
        loop_num = 0;

        ix = (base[i].x - x_min) / hash_size;
        iy = (base[i].y - y_min) / hash_size;

        while (corr_v.size() < 9) {
            // 適用するcorrmapの探索
            for (int iix = -1 * loop_num; iix <= loop_num; iix++) {
                for (int iiy = -1 * loop_num; iiy <= loop_num; iiy++) {
                    if (abs(iix) != loop_num && abs(iiy) != loop_num)continue;

                    id.first = ix + iix;
                    id.second = iy + iiy;
                    if (corr_hash.count(id) == 0)continue;
                    auto range = corr_hash.equal_range(id);
                    for (auto res = range.first; res != range.second; res++) {
                        corr_v.push_back(res->second);

                    }

                }

            }
            if (loop_num > 100)break;
            loop_num++;
        }
        // trackに最も近いcorrmapを選ぶ
        if (corr_v.size() == 0)continue;
        for (int j = 0; j < corr_v.size(); j++) {
            if (j == 0) {
                dist = pow(base[i].x - corr_v[j].cx, 2) + pow(base[i].y - corr_v[j].cy, 2);
                param = corr_v[j];
            }
            if (dist > pow(base[i].x - corr_v[j].cx, 2) + pow(base[i].y - corr_v[j].cy, 2)) {
                dist = pow(base[i].x - corr_v[j].cx, 2) + pow(base[i].y - corr_v[j].cy, 2);
                param = corr_v[j];
            }
        }
        ret.push_back(std::make_pair(base[i], param));
    }

    return ret;

}
std::vector<std::pair<vxx::base_track_t, corrmap>> correspond_corrmap_bwd(std::vector<vxx::base_track_t>& base, std::vector <corrmap>& corr) {
    // VE alignmentで上流方向に外挿=pl0 --> pl1の逆変換
    // 使うparamはpl0の位置でのでよい

    double x_min = 0, y_min = 0;
    for (int i = 0; i < corr.size(); i++) {
        if (i == 0) {
            x_min = corr[i].cx;
            y_min = corr[i].cy;
        }
        x_min = std::min(x_min, corr[i].cx);
        y_min = std::min(y_min, corr[i].cy);
    }

    double hash_size = 2000;
    std::multimap<std::pair<int, int>, corrmap> corr_hash;
    std::pair<int, int>id;
    for (int i = 0; i < corr.size(); i++) {
        id.first = (corr[i].cx - x_min) / hash_size;
        id.second = (corr[i].cy - y_min) / hash_size;
        corr_hash.insert(std::make_pair(id, corr[i]));
    }

    int ix, iy, loop_num;
    std::vector<std::pair<vxx::base_track_t, corrmap>>ret;
    std::vector<corrmap>corr_v;
    double dist;
    corrmap param;
    for (int i = 0; i < base.size(); i++) {
        corr_v.clear();
        loop_num = 0;

        ix = (base[i].x - x_min) / hash_size;
        iy = (base[i].y - y_min) / hash_size;

        while (corr_v.size() < 9) {
            // 適用するcorrmapの探索
            for (int iix = -1 * loop_num; iix <= loop_num; iix++) {
                for (int iiy = -1 * loop_num; iiy <= loop_num; iiy++) {
                    if (abs(iix) != loop_num && abs(iiy) != loop_num)continue;

                    id.first = ix + iix;
                    id.second = iy + iiy;
                    if (corr_hash.count(id) == 0)continue;
                    auto range = corr_hash.equal_range(id);
                    for (auto res = range.first; res != range.second; res++) {
                        corr_v.push_back(res->second);

                    }

                }

            }
            if (loop_num > 100)break;
            loop_num++;
        }
        // trackに最も近いcorrmapを選ぶ
        if (corr_v.size() == 0)continue;
        for (int j = 0; j < corr_v.size(); j++) {
            if (j == 0) {
                dist = pow(base[i].x - corr_v[j].cx, 2) + pow(base[i].y - corr_v[j].cy, 2);
                param = corr_v[j];
            }
            if (dist > pow(base[i].x - corr_v[j].cx, 2) + pow(base[i].y - corr_v[j].cy, 2)) {
                dist = pow(base[i].x - corr_v[j].cx, 2) + pow(base[i].y - corr_v[j].cy, 2);
                param = corr_v[j];
            }
        }
        ret.push_back(std::make_pair(base[i], param));
    }

    return ret;

}


void calc_invers_corr_area(std::vector<corrmap>& corr) {
    double tmp_x, tmp_y, denominator;
    for (int i = 0; i < corr.size(); i++) {
        // cernter of corrmap-area
        corr[i].cx = (corr[i].x_area[0] + corr[i].x_area[1]) / 2;
        corr[i].cy = (corr[i].y_area[0] + corr[i].y_area[1]) / 2;

        tmp_x = corr[i].cx;
        tmp_y = corr[i].cy;

        //視野中心の逆変換（座標の変換） pl71のpl70
        tmp_x = tmp_x - corr[i].position[4];
        tmp_y = tmp_y - corr[i].position[5];
        denominator = corr[i].position[0] * corr[i].position[3] - corr[i].position[1] * corr[i].position[2];

        corr[i].cx = (tmp_x * corr[i].position[3] - tmp_y * corr[i].position[1]) / denominator;
        corr[i].cy = (-1 * tmp_x * corr[i].position[2] + tmp_y * corr[i].position[0]) / denominator;
    }
    // pl1系の視野中心をpl0系に変換。pl1の座標でこの変換後後の視野中心と最も近いtrkがそのparamを適用するtrk
    //変換対象の track に最も近い中心座標を持つ区画を選び、補正値を適用する。
    //各区画の領域情報(xmin, xmax, ymin, ymax)は pos0 の座標系なので、これを pos1 の座標系に変換する。この変換にはその区画の afp を逆方向に適用する。
    //NETSCAN 内部では(xmin, ymin) (xmax, ymax) を pos1 の座標系に変換し、これらを対角線上の２頂点とする矩形領域を pos1 座標系での領域と近似している。
    //補正の適用方向は順方向、つまり、位置xyに対しては afp を、角度に対しては aft をそのまま適用し、z は z + dz とする。
}



// partner track of fwd direction
void extraporate_fwd_trk(int pl0, int pl1, vxx::base_track_t& base, corrmap& corr) {
    // 場合分けめんどいのでpl0-->vpl,pl1-->外挿先の#plにする
    // pl1:upstream
    // pl1-->pl0への外挿, pl1系-->pl0系への変換(PL071-->PL070, PL074-->pl077)
    double tan, ex_x, ex_y;
    Material m;
    double zgap, dz, rotation;

    //32 5-->4 76
    // norminalなzgap
    if (pl1 - pl0 == -1) {
        // iron
        zgap = m.iron;
    }
    else if (pl1 - pl0 == -3) {
        //iron + film + pack + pack
        //zgap = m.iron + m.film + m.pack + m.pack;
        zgap = m.film + m.pack + m.pack;
    }

    dz = zgap + corr.dz;
    //std::cout << " * zgap = " << zgap << " ; dz = " << dz << " [um]" << std::endl;
    // extraporation
    tan = sqrt(base.ax * base.ax + base.ay * base.ay);
    ex_x = base.x + base.ax * dz;
    ex_y = base.y + base.ay * dz;

    double tmp_x, tmp_y, denominator;


    // trf
    tmp_x = ex_x * corr.position[0] + ex_y * corr.position[1];
    tmp_y = ex_x * corr.position[2] + ex_y * corr.position[3];
    ex_x = tmp_x + corr.position[4];
    ex_y = tmp_y + corr.position[5];
    base.x = ex_x;
    base.y = ex_y;

    rotation = atan(corr.position[2] / corr.position[0]);
    base.ax = (base.ax * cos(rotation) - base.ay * sin(rotation)) * m.angle_shrink + corr.angle[4];
    base.ay = (base.ax * sin(rotation) + base.ay * cos(rotation)) * m.angle_shrink + corr.angle[5];

    base.pl = pl1;
    base.rawid = 0;
    base.m[0].ph = 0;
    base.m[1].ph = 0;

}

// partner track of bwd direction
void extraporate_bwd_trk(int pl0, int pl1, vxx::base_track_t& base, corrmap& corr) {
    // 場合分けめんどいのでpl0-->vpl,pl1-->外挿先の#plにする
    // pl1:upstream
    // pl0-->pl1への外挿, pl0系-->pl1系への変換(PL070->pl071, pl071->pl068)
    double tan, ex_x, ex_y;
    Material m;
    double zgap, dz, rotation;

    //32 5 < --4 76
    // norminalなzgap
    if (pl1 - pl0 == 1) {
        // iron
        zgap = m.iron;
    }
    else if (pl1 - pl0 == -3) {
        //iron + film + pack + pack
        zgap = m.film + m.pack + m.pack;
    }
    else 
    {
        std::cout << "order x" << std::endl;
    }

    dz = zgap + corr.dz;

    // extraporation
    tan = sqrt(base.ax * base.ax + base.ay * base.ay);
    ex_x = base.x - base.ax * dz;
    ex_y = base.y - base.ay * dz;

    double tmp_x, tmp_y, denominator;


    // inverse trf
    tmp_x = ex_x - corr.position[4];
    tmp_y = ex_y - corr.position[5];
    denominator = corr.position[0] * corr.position[3] - corr.position[1] * corr.position[2];

    ex_x = (tmp_x * corr.position[3] - tmp_y * corr.position[1]) / denominator;
    ex_y = (-1 * tmp_x * corr.position[2] + tmp_y * corr.position[0]) / denominator;
    base.x = ex_x;
    base.y = ex_y;

    rotation = atan(corr.position[2] / corr.position[0]);
    base.ax = (base.ax * cos(rotation) - base.ay * sin(rotation)) * m.angle_shrink + corr.angle[4];
    base.ay = (base.ax * sin(rotation) + base.ay * cos(rotation)) * m.angle_shrink + corr.angle[5];

    base.pl = pl1;
    base.rawid = 0;
    base.m[0].ph = 0;
    base.m[1].ph = 0;

}

void set_basetrack(std::vector<Momentum_recon::Event_information>& momch, std::set<int>& eid, std::multimap<int, EventTrackInfo>& eve) {
    EventTrackInfo btk = { 0 };
    int d_pl;
    for (auto ev : momch) {
        eid.insert(ev.groupid);
        btk.gid = ev.groupid;
        btk.vpl = ev.vertex_pl;
        for (auto c : ev.chains) {
            btk.cid = c.chainid;
            d_pl = 133;
            if (c.base.rbegin()->pl - c.base.begin()->pl + 1 < 15)continue;
            for (auto b : c.base) {
                if (abs(b.pl - btk.vpl) < d_pl) {
                    btk.pl = b.pl;
                    btk.rid = b.rawid;
                    d_pl = abs(btk.pl - btk.vpl);
                }
            }
            eve.insert(std::make_pair(ev.groupid, btk));
        }
    }
}

void base_matching(std::vector<vxx::base_track_t>& base0, std::vector<vxx::base_track_t>& bese_ref,
    double pos_lat_all[2], double pos_rad_all[2], double ang_lat_all[2], double ang_rad_all[2], std::vector<output_format>& out_base) {

    std::cout << "check:" << pos_lat_all[0] << " " << pos_lat_all[1] << std::endl;
    double x_min = 0, y_min = 0;
    for (int i = 0; i < bese_ref.size(); i++) {
        if (i == 0) {
            x_min = bese_ref[i].x;
            y_min = bese_ref[i].y;
        }
        x_min = std::min(bese_ref[i].x, x_min);
        y_min = std::min(bese_ref[i].y, y_min);
    }

    double hash_size = 2000;
    std::multimap<std::pair<int, int>, vxx::base_track_t>base1_hash;
    std::pair<int, int>id;
    double dz_mean = 0;
    for (int i = 0; i < bese_ref.size(); i++) {
        id.first = (bese_ref[i].x - x_min) / hash_size;
        id.second = (bese_ref[i].y - y_min) / hash_size;
        base1_hash.insert(std::make_pair(id, bese_ref[i]));
        dz_mean += bese_ref[i].z;
    }
    dz_mean /= bese_ref.size();
    printf("dz mean=%g\n", dz_mean);


    double angle, ex_x, ey_y;
    double ix, iy;
    double diff_pos_rad, diff_pos_lat, diff_ang_rad, diff_ang_lat;
    std::vector<vxx::base_track_t>base_buf;
    bool flg = false;
    int count = 0;
    for (int i = 0; i < base0.size(); i++) {
        if (count % 10000 == 0) {
            fprintf(stderr, "\r base matchin %lld/%lld(%4.1lf%%)", count, base0.size(), count * 100. / base0.size());
        }
        count++;

        base_buf.clear();
        output_format out;
        out.ex.rawid = base0[i].rawid;
        out.ex.ax = base0[i].ax;
        out.ex.ay = base0[i].ay;
        out.ex.x = base0[i].x;
        out.ex.y = base0[i].y;
        out.hit_flg = false;
        out.bef.rawid = -1;
        out.bef.ax = 0;
        out.bef.ay = 0;
        out.bef.x = 0;
        out.bef.y = 0;



        //angle = sqrt(base0[i].ax * base0[i].ax + base0[i].ay * base0[i].ay);
        //ex_x = base0[i].x + base0[i].ax * dz_mean;
        //ey_y = base0[i].y + base0[i].ay * dz_mean;
        //ix = (ex_x - x_min) / hash_size;
        //iy = (ey_y - y_min) / hash_size;
        //out.ex.x = ex_x;
        //out.ex.y = ey_y;


        for (int iix = -1; iix <= 1; iix++) {
            for (int iiy = -1; iiy <= 1; iiy++) {
                id.first = ix + iix;
                id.second = iy + iiy;
                if (base1_hash.count(id) == 0)continue;
                auto range = base1_hash.equal_range(id);
                for (auto res = range.first; res != range.second; res++) {
                    base_buf.push_back(res->second);

                }
            }

        }

        for (int j = 0; j < base_buf.size(); j++) {
            //ex_x = base0[i].x + base0[i].ax * base_buf[j].z;
            //ey_y = base0[i].y + base0[i].ay * base_buf[j].z;

            if (angle > 0.01) {
                diff_pos_rad = ((ex_x - base_buf[j].x) * base0[i].ax + (ey_y - base_buf[j].y) * base0[i].ay) / angle;
                diff_pos_lat = ((ex_x - base_buf[j].x) * base0[i].ay - (ey_y - base_buf[j].y) * base0[i].ax) / angle;
                diff_ang_rad = ((base0[i].ax - base_buf[j].ax) * base0[i].ax + (base0[i].ay - base_buf[j].ay) * base0[i].ay) / angle;
                diff_ang_lat = ((base0[i].ax - base_buf[j].ax) * base0[i].ay - (base0[i].ay - base_buf[j].ay) * base0[i].ax) / angle;

            }
            else {
                diff_pos_rad = (ey_y - base_buf[j].y);
                diff_pos_lat = (ex_x - base_buf[j].x);
                diff_ang_rad = base0[i].ay - base_buf[j].ay;
                diff_ang_lat = base0[i].ax - base_buf[j].ax;
            }
            if (fabs(diff_pos_rad) > pos_rad_all[0] + pos_rad_all[1] * angle)continue;
            if (fabs(diff_pos_lat) > pos_lat_all[0] + pos_lat_all[1] * angle)continue;
            if (fabs(diff_ang_rad) > ang_rad_all[0] + ang_rad_all[1] * angle)continue;
            if (fabs(diff_ang_lat) > ang_lat_all[0] + ang_lat_all[1] * angle)continue;
            out.ex.x = ex_x;
            out.ex.y = ey_y;

            out.hit_flg = true;
            out.aft.rawid = base_buf[j].rawid;
            out.aft.ax = base_buf[j].ax;
            out.aft.ay = base_buf[j].ay;
            out.aft.x = base_buf[j].x;
            out.aft.y = base_buf[j].y;
            break;
        }
        out_base.push_back(out);
    }
    fprintf(stderr, "\r base matching %lld/%lld(%4.1lf%%)\n", count, base0.size(), count * 100. / base0.size());


}
void output_basetracks(std::string filename, std::vector<output_format>& out) {
    std::ofstream ofs(filename);
    int count = 0;
    for (auto itr = out.begin(); itr != out.end(); itr++) {
        if (count % 10000 == 0) {
            fprintf(stderr, "\r write result %d/%d(%4.1lf%%)", count, out.size(), count * 100. / out.size());
        }
        count++;
        ofs << std::right << std::fixed
            << std::setw(10) << std::setprecision(0) << itr->ex.rawid << " "
            << std::setw(8) << std::setprecision(4) << itr->ex.ax << " "
            << std::setw(8) << std::setprecision(4) << itr->ex.ay << " "
            << std::setw(10) << std::setprecision(1) << itr->ex.x << " "
            << std::setw(10) << std::setprecision(1) << itr->ex.y << " "
            << std::setw(2) << std::setprecision(1) << itr->hit_flg << " "
            << std::setw(10) << std::setprecision(0) << itr->aft.rawid << " "
            << std::setw(8) << std::setprecision(4) << itr->aft.ax << " "
            << std::setw(8) << std::setprecision(4) << itr->aft.ay << " "
            << std::setw(10) << std::setprecision(1) << itr->aft.x << " "
            << std::setw(10) << std::setprecision(1) << itr->aft.y << " "
            << std::setw(10) << std::setprecision(0) << itr->bef.rawid << " "
            << std::setw(8) << std::setprecision(4) << itr->bef.ax << " "
            << std::setw(8) << std::setprecision(4) << itr->bef.ay << " "
            << std::setw(10) << std::setprecision(1) << itr->bef.x << " "
            << std::setw(10) << std::setprecision(1) << itr->bef.y << std::endl;
    }
    fprintf(stderr, "\r write result %d/%d(%4.1lf%%)\n", count, out.size(), count * 100. / out.size());

}

