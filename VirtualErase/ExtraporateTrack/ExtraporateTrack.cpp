// extrapolationと同じ
//I:\NINJA\E71a\work\kasumi\work\test\x64\Release\exp.exe
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

class micro_track_subset_t {
public:
    double ax, ay;
    double z;
    int ph;
    int pos, col, row, zone, isg;
    int64_t rawid;
};

class base_track_t {
public:
    double ax, ay;
    double x, y, z;
    int pl;
    int isg, zone;
    int dmy;    // In ROOT, you will have to add this member because CINT does not handle 8byte alignment. 
    int64_t rawid;
    micro_track_subset_t m[2];
};

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

std::vector<base_track_t> read_base(std::string filename);
std::vector<vxx::base_track_t> read_base_id(std::string file_in_ECC, int pl, std::set<int>rawid);
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


int main(int argc, char** argv) {
    if (argc < 5) {
        fprintf(stderr, "usage:prg event.momch ECC_Area_path corrmap-path_ve output-path \n");
        fprintf(stderr, "usage:prg event.momch ECC_Area_path corrmap-path_ve output-path #ecc(except ecc5)\n");
        fprintf(stderr, "I : \NINJA\E71a\ECC5 I : \NINJA\E71a\work\kasumi\ECC5\VirtualErase\n");
        fprintf(stderr, " T : \NINJA\E71a\ECC2 K : \NINJA\E71a\work\zkasumi\VirtualErase\ECC2\01_Alingment\\corrmap\n");
        exit(1);
    }

    std::string in_momch, corr_path, corr_path_ve, output_path, ECC_Area_path, file_in_base[2], file_in_corr, file_out;

    in_momch = argv[1];
    ECC_Area_path = argv[2];
    corr_path_ve = argv[3];
    output_path = argv[4];
    int ecc = 5;
    if (argc == 6) {
        ecc = std::stoi(argv[5]);
    }

    corr_path = ECC_Area_path + "\\Area0\\0\\align";

    std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);
    std::multimap<int, EventTrackInfo> eve;
    std::set<int>eid;
    set_basetrack(momch, eid, eve);

    int vpl;
    for (auto itr = eid.begin(); itr != eid.end(); itr++) {
        std::set<int>f;
        std::set<int>b;
        std::map<int, int>chainid;
        auto itr1 = eve.equal_range(*itr);
        for (auto itr2 = itr1.first; itr2 != itr1.second; itr2++) {
            vpl = itr2->second.vpl;
            chainid.insert(std::make_pair(itr2->second.rid, itr2->second.cid));
            if (itr2->second.pl > itr2->second.vpl) {
                // pl - vpl: bwd-->外挿後に順変換
                b.insert(itr2->second.rid);
            }
            else {
                // fwd-->外挿後に逆変換
                f.insert(itr2->second.rid);
            }
        }
        std::cout << " \t EVENT " << *itr << std::endl;
        std::cout << " * the number of fwd trk " << f.size() << std::endl;
        std::cout << " * the number of bwd trk " << b.size() << std::endl;
        std::stringstream name, ss;
        std::vector<corrmap>corr;
        std::vector<std::pair<vxx::base_track_t, corrmap>>base_pair;
        std::vector<vxx::base_track_t>base_b[4];
        std::vector<vxx::base_track_t>base_f[4];
        if (b.size() != 0) {
            base_b[0] = read_base_id(ECC_Area_path, vpl + 1, b);// pl0:downstream

            ss << corr_path + "\\corrmap-align-" << std::setw(3) << std::setfill('0') << vpl + 1 << "-" << std::setw(3) << std::setfill('0') << vpl + 2 << ".lst";
            file_in_corr = ss.str();
            std::cout << file_in_corr << std::endl;
            corr = read_corrmap(file_in_corr);
            base_pair = correspond_corrmap_bwd(base_b[0], corr);
            for (auto itr = base_pair.begin(); itr != base_pair.end(); itr++) {
                extraporate_bwd_trk(vpl + 1, vpl + 2, itr->first, itr->second);
                base_b[1].push_back(itr->first);
            }
            ss.str("");
            corr.clear();
            base_pair.clear();

            // vpl - 1
            if(ecc==5){
            ss << corr_path_ve + "\\f_pl" << std::setw(3) << std::setfill('0') << vpl + 2 << "_" << std::setw(3) << std::setfill('0') << vpl - 1 << "\\corrmap-" << std::setw(3) << std::setfill('0') << vpl + 2 << "-" << std::setw(3) << std::setfill('0') << vpl - 1 << "-l.txt";
            }
            else {
                ss << corr_path_ve << "\\corrmap-" << std::setw(3) << std::setfill('0') << vpl + 2 << "-" << std::setw(3) << std::setfill('0') << vpl - 1 << "-l.txt";
            }
            file_in_corr = ss.str();
            std::cout << file_in_corr << std::endl;
            corr = read_corrmap(file_in_corr);
            base_pair = correspond_corrmap_bwd(base_b[1], corr);
            for (auto itr = base_pair.begin(); itr != base_pair.end(); itr++) {
                extraporate_bwd_trk(vpl + 1, vpl - 1, itr->first, itr->second);
                base_b[2].push_back(itr->first);
            }
            ss.str("");
            corr.clear();
            base_pair.clear();

            // vpl 
            ss << corr_path + "\\corrmap-align-" << std::setw(3) << std::setfill('0') << vpl - 1 << "-" << std::setw(3) << std::setfill('0') << vpl << ".lst";
            file_in_corr = ss.str();
            std::cout << file_in_corr << std::endl;
            corr = read_corrmap(file_in_corr);
            base_pair = correspond_corrmap_bwd(base_b[2], corr);
            for (auto itr = base_pair.begin(); itr != base_pair.end(); itr++) {
                extraporate_bwd_trk(vpl + 1, vpl, itr->first, itr->second);
                base_b[3].push_back(itr->first);
            }
            ss.str("");
            corr.clear();
            base_pair.clear();


        }
        if (f.size() != 0) {
            base_f[0] = read_base_id(ECC_Area_path, vpl, f);// pl0:downstream

            // vpl - 1
            ss << corr_path + "\\corrmap-align-" << std::setw(3) << std::setfill('0') << vpl - 1 << "-" << std::setw(3) << std::setfill('0') << vpl << ".lst";
            file_in_corr = ss.str();
            std::cout << file_in_corr << std::endl;
            corr = read_corrmap(file_in_corr);
            base_pair = correspond_corrmap_fwd(base_f[0], corr);
            for (auto itr = base_pair.begin(); itr != base_pair.end(); itr++) {
                extraporate_fwd_trk(vpl, vpl - 1, itr->first, itr->second);
                base_f[1].push_back(itr->first);
            }
            ss.str("");
            corr.clear();
            base_pair.clear();

            // vpl + 2
            if (ecc == 5) {
                ss << corr_path_ve + "\\f_pl" << std::setw(3) << std::setfill('0') << vpl + 2 << "_" << std::setw(3) << std::setfill('0') << vpl - 1 << "\\corrmap-" << std::setw(3) << std::setfill('0') << vpl + 2 << "-" << std::setw(3) << std::setfill('0') << vpl - 1 << "-l.txt";
            }
            else {
                ss << corr_path_ve << "\\corrmap-" << std::setw(3) << std::setfill('0') << vpl + 2 << "-" << std::setw(3) << std::setfill('0') << vpl - 1 << "-l.txt";
            }
            file_in_corr = ss.str();
            std::cout << file_in_corr << std::endl;
            corr = read_corrmap(file_in_corr);
            base_pair = correspond_corrmap_fwd(base_f[1], corr);
            for (auto itr = base_pair.begin(); itr != base_pair.end(); itr++) {
                extraporate_fwd_trk(vpl, vpl + 2, itr->first, itr->second);
                base_f[2].push_back(itr->first);
            }
            ss.str("");
            corr.clear();
            base_pair.clear();

            // vpl + 1
            ss << corr_path + "\\corrmap-align-" << std::setw(3) << std::setfill('0') << vpl + 1 << "-" << std::setw(3) << std::setfill('0') << vpl + 2 << ".lst";
            file_in_corr = ss.str();
            std::cout << file_in_corr << std::endl;
            corr = read_corrmap(file_in_corr);
            base_pair = correspond_corrmap_fwd(base_f[2], corr);
            for (auto itr = base_pair.begin(); itr != base_pair.end(); itr++) {
                extraporate_fwd_trk(vpl, vpl + 1, itr->first, itr->second);
                base_f[3].push_back(itr->first);
            }
            ss.str("");
            corr.clear();
            base_pair.clear();


        }
        // set folder name
        name << output_path << "\\eve" << std::setw(5) << std::setfill('0') << *itr << ".txt";
        std::string output = name.str();
        std::ofstream ofs(output);
        ofs << "bwd " << std::endl;
        for (int itr = 0; itr != base_b[0].size(); itr++) {
           // ofs << std::next(eve.equal_range(base_b[0][itr].rawid).first,itr)->second.gid << " " << std::next(eve.equal_range(base_b[0][itr].rawid).first, itr)->second.vpl << " " << std::next(eve.equal_range(base_b[0][itr].rawid).first, itr)->second.cid << std::endl;
            ofs<<chainid.find(base_b[0][itr].rawid)->second << std::endl;
            ofs << std::right << std::fixed
                << std::setw(3) << std::setprecision(0) << base_b[0][itr].pl << " "
                << std::setw(10) << std::setprecision(0) << base_b[0][itr].rawid << " "
                << std::setw(7) << std::setprecision(0) << base_b[0][itr].m[0].ph << " "
                << std::setw(7) << std::setprecision(0) << base_b[0][itr].m[1].ph << " "
                << std::setw(8) << std::setprecision(4) << base_b[0][itr].ax << " "
                << std::setw(8) << std::setprecision(4) << base_b[0][itr].ay << " "
                << std::setw(8) << std::setprecision(1) << base_b[0][itr].x << " "
                << std::setw(8) << std::setprecision(1) << base_b[0][itr].y << " "
                << std::setw(8) << std::setprecision(1) << base_b[0][itr].z << std::endl;
            ofs << std::right << std::fixed
                << std::setw(3) << std::setprecision(0) << base_b[1][itr].pl << " "
                << std::setw(10) << std::setprecision(0) << base_b[1][itr].rawid << " "
                << std::setw(7) << std::setprecision(0) << base_b[1][itr].m[0].ph << " "
                << std::setw(7) << std::setprecision(0) << base_b[1][itr].m[1].ph << " "
                << std::setw(8) << std::setprecision(4) << base_b[1][itr].ax << " "
                << std::setw(8) << std::setprecision(4) << base_b[1][itr].ay << " "
                << std::setw(8) << std::setprecision(1) << base_b[1][itr].x << " "
                << std::setw(8) << std::setprecision(1) << base_b[1][itr].y << " "
                << std::setw(8) << std::setprecision(1) << base_b[1][itr].z << std::endl;
            ofs << std::right << std::fixed
                << std::setw(3) << std::setprecision(0) << base_b[2][itr].pl << " "
                << std::setw(10) << std::setprecision(0) << base_b[2][itr].rawid << " "
                << std::setw(7) << std::setprecision(0) << base_b[2][itr].m[0].ph << " "
                << std::setw(7) << std::setprecision(0) << base_b[2][itr].m[1].ph << " "
                << std::setw(8) << std::setprecision(4) << base_b[2][itr].ax << " "
                << std::setw(8) << std::setprecision(4) << base_b[2][itr].ay << " "
                << std::setw(8) << std::setprecision(1) << base_b[2][itr].x << " "
                << std::setw(8) << std::setprecision(1) << base_b[2][itr].y << " "
                << std::setw(8) << std::setprecision(1) << base_b[2][itr].z << std::endl;
            ofs << std::right << std::fixed
                << std::setw(3) << std::setprecision(0) << base_b[3][itr].pl << " "
                << std::setw(10) << std::setprecision(0) << base_b[3][itr].rawid << " "
                << std::setw(7) << std::setprecision(0) << base_b[3][itr].m[0].ph << " "
                << std::setw(7) << std::setprecision(0) << base_b[3][itr].m[1].ph << " "
                << std::setw(8) << std::setprecision(4) << base_b[3][itr].ax << " "
                << std::setw(8) << std::setprecision(4) << base_b[3][itr].ay << " "
                << std::setw(8) << std::setprecision(1) << base_b[3][itr].x << " "
                << std::setw(8) << std::setprecision(1) << base_b[3][itr].y << " "
                << std::setw(8) << std::setprecision(1) << base_b[3][itr].z << std::endl;
        }
        ofs << "fwd" << std::endl;
        for (int itr = 0; itr != base_f[0].size(); itr++) {
            ofs << chainid.find(base_f[0][itr].rawid)->second << std::endl;
            ofs << std::right << std::fixed
                << std::setw(3) << std::setprecision(0) << base_f[0][itr].pl << " "
                << std::setw(10) << std::setprecision(0) << base_f[0][itr].rawid << " "
                << std::setw(7) << std::setprecision(0) << base_f[0][itr].m[0].ph << " "
                << std::setw(7) << std::setprecision(0) << base_f[0][itr].m[1].ph << " "
                << std::setw(8) << std::setprecision(4) << base_f[0][itr].ax << " "
                << std::setw(8) << std::setprecision(4) << base_f[0][itr].ay << " "
                << std::setw(8) << std::setprecision(1) << base_f[0][itr].x << " "
                << std::setw(8) << std::setprecision(1) << base_f[0][itr].y << " "
                << std::setw(8) << std::setprecision(1) << base_f[0][itr].z << std::endl;
            ofs << std::right << std::fixed
                << std::setw(3) << std::setprecision(0) << base_f[1][itr].pl << " "
                << std::setw(10) << std::setprecision(0) << base_f[1][itr].rawid << " "
                << std::setw(7) << std::setprecision(0) << base_f[1][itr].m[0].ph << " "
                << std::setw(7) << std::setprecision(0) << base_f[1][itr].m[1].ph << " "
                << std::setw(8) << std::setprecision(4) << base_f[1][itr].ax << " "
                << std::setw(8) << std::setprecision(4) << base_f[1][itr].ay << " "
                << std::setw(8) << std::setprecision(1) << base_f[1][itr].x << " "
                << std::setw(8) << std::setprecision(1) << base_f[1][itr].y << " "
                << std::setw(8) << std::setprecision(1) << base_f[1][itr].z << std::endl;
            ofs << std::right << std::fixed
                << std::setw(3) << std::setprecision(0) << base_f[2][itr].pl << " "
                << std::setw(10) << std::setprecision(0) << base_f[2][itr].rawid << " "
                << std::setw(7) << std::setprecision(0) << base_f[2][itr].m[0].ph << " "
                << std::setw(7) << std::setprecision(0) << base_f[2][itr].m[1].ph << " "
                << std::setw(8) << std::setprecision(4) << base_f[2][itr].ax << " "
                << std::setw(8) << std::setprecision(4) << base_f[2][itr].ay << " "
                << std::setw(8) << std::setprecision(1) << base_f[2][itr].x << " "
                << std::setw(8) << std::setprecision(1) << base_f[2][itr].y << " "
                << std::setw(8) << std::setprecision(1) << base_f[2][itr].z << std::endl;
            ofs << std::right << std::fixed
                << std::setw(3) << std::setprecision(0) << base_f[3][itr].pl << " "
                << std::setw(10) << std::setprecision(0) << base_f[3][itr].rawid << " "
                << std::setw(7) << std::setprecision(0) << base_f[3][itr].m[0].ph << " "
                << std::setw(7) << std::setprecision(0) << base_f[3][itr].m[1].ph << " "
                << std::setw(8) << std::setprecision(4) << base_f[3][itr].ax << " "
                << std::setw(8) << std::setprecision(4) << base_f[3][itr].ay << " "
                << std::setw(8) << std::setprecision(1) << base_f[3][itr].x << " "
                << std::setw(8) << std::setprecision(1) << base_f[3][itr].y << " "
                << std::setw(8) << std::setprecision(1) << base_f[3][itr].z << std::endl;
        }
        f.clear();
        b.clear();
        for (int i = 0; i < 4; i++) {
            base_b[i].clear();
            base_f[i].clear();
        }
    }


}

std::vector<base_track_t> read_base(std::string filename) {
    std::vector<base_track_t> ret;

    //binary modeで開く
    std::ifstream ifs(filename, std::ios::binary);
    int count = 0;
    base_track_t base;
    while (ifs.read((char*)&base, sizeof(base))) {
        if (count % 10000 == 0) {
            printf("\r read base --> %d", count);
        }
        count++;
        ret.push_back(base);

    }
    printf("\r read base --> %d fin\n", count);
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
    // pl1-->pl0への外挿, pl1系-->pl0系への変換
    double tan, ex_x, ex_y;
    Material m;
    double zgap, dz, rotation;

    // norminalなzgap
    if (pl1 - pl0 == -1) {
        // iron
        zgap = m.iron;
    }
    else if (pl1 - pl0 == 2) {
        //iron + film + pack + pack
        //zgap = m.iron + m.film + m.pack + m.pack;
        zgap = m.film + m.pack + m.pack;
    }
    else if (pl1 - pl0 == 1) {
        // iron + file + pack + pack + film + iron
        zgap = m.film + m.iron;
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
    // pl0-->pl1への外挿, pl0系-->pl1系への変換
    double tan, ex_x, ex_y;
    Material m;
    double zgap, dz, rotation;

    // norminalなzgap
    if (pl1 - pl0 == 1) {
        // iron
        zgap = m.iron;
    }
    else if (pl1 - pl0 == -2) {
        //iron + film + pack + pack
        zgap = m.film + m.pack + m.pack;
    }
    else if (pl1 - pl0 == -1) {
        // iron + file + pack + pack + film + iron
        zgap =m.film + m.iron;
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