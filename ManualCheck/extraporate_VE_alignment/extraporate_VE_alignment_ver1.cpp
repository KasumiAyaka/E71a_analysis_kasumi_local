#include <cstdint> 
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>

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

std::vector<base_track_t> read_base(std::string filename);
std::vector<corrmap>read_corrmap(std::string filename);
void calc_invers_corr_area(std::vector<corrmap>& corr);
std::vector<std::pair<base_track_t, corrmap>> correspond_corrmap(std::vector<base_track_t>& base, std::vector < corrmap>& corr);
std::vector<base_track_t>base_trans(std::vector<std::pair<base_track_t, corrmap>>& base_pair, double dz_nominal);
std::vector<output_format>  base_matching(std::vector <base_track_t>& base0, std::vector<base_track_t>& base1,
    double pos_lat_all[2], double pos_rad_all[2], double ang_lat_all[2], double ang_rad_all[2]);
void pick_up_base_inf(std::vector<output_format>& out, std::vector <base_track_t>& base);
void output_basetracks(std::string filename, std::vector<output_format>& out);


int main(int argc, char** argv) {
    if (argc != 6) {
        fprintf(stderr, "usage:prg file-in-base-bin0 file-in-bin1 file-in-corr file-out　dz\n");
        //exp.exe   ..\..\..\b070.bbt   ..\..\..\b071.bbt   ..\..\..\corrmap-070-071-l1.txt
        //argv0       1                    2                  3
        exit(1);
    }

    std::string file_in_base[2], file_in_corr, file_out;
    file_in_base[0] = argv[1];
    file_in_base[1] = argv[2];
    file_in_corr = argv[3];
    file_out = argv[4];
    double dz_nominal = -(350 + 500); 
    dz_nominal = std::stod(argv[5]);

    std::vector<base_track_t>base[2];
    base[0] = read_base(file_in_base[0]);
    base[1] = read_base(file_in_base[1]);

    std::vector<corrmap>corr;
    corr = read_corrmap(file_in_corr);

    std::vector<std::pair<base_track_t, corrmap>>base_pair = correspond_corrmap(base[1], corr);

    
    std::vector<base_track_t> base1_trns = base_trans(base_pair, dz_nominal);

    double pos_lat_all[2];
    double pos_rad_all[2]; 
    double ang_lat_all[2]; 
    double ang_rad_all[2];
    //arrowane[0]+ [1]*rad
    //3つを　固定→1つを変化x４回

    pos_lat_all[0] = 40;
    pos_lat_all[1] = 0;
    pos_rad_all[0] = 40;
    pos_rad_all[1] = 40;
    ang_lat_all[0] = 0.05;
    ang_lat_all[1] = 0.05;
    ang_rad_all[0] = 0.1;
    ang_rad_all[1] = 0.1;
    /*
    pos_lat_all[0] = 30;
    pos_lat_all[1] = 0;
    pos_rad_all[0] = 40;
    pos_rad_all[1] = 30;
    ang_lat_all[0] = 0.05;
    ang_lat_all[1] = 0;
    ang_rad_all[0] = 0.1;
    ang_rad_all[1] = 0.1;
      ang_rad_all[1] = 0.1;
    */


    //base1_trns --> 変換後
    //base[1] --> 変換前
    std::vector<output_format> out_base = base_matching(base[0], base1_trns, pos_lat_all, pos_rad_all, ang_lat_all, ang_rad_all);
    pick_up_base_inf(out_base, base[1]);
    output_basetracks(file_out, out_base);
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

std::vector<std::pair<base_track_t, corrmap>> correspond_corrmap(std::vector<base_track_t>& base, std::vector < corrmap>& corr) {
    calc_invers_corr_area(corr);

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
    std::vector<std::pair<base_track_t, corrmap>>ret;
    std::vector<corrmap>corr_v;
    double dist;
    corrmap param;
    for (int i = 0; i < base.size(); i++) {
        corr_v.clear();
        loop_num = 0;

        ix = (base[i].x - x_min) / hash_size;
        iy = (base[i].y - y_min) / hash_size;

        while (corr_v.size() < 9) {
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
        corr[i].cx = (corr[i].x_area[0] + corr[i].x_area[1]) / 2;
        corr[i].cy = (corr[i].y_area[0] + corr[i].y_area[1]) / 2;

        tmp_x = corr[i].cx;
        tmp_y = corr[i].cy;


        tmp_x = tmp_x - corr[i].position[4];
        tmp_y = tmp_y - corr[i].position[5];
        denominator = corr[i].position[0] * corr[i].position[3] - corr[i].position[1] * corr[i].position[2];

        corr[i].cx = (tmp_x * corr[i].position[3] - tmp_y * corr[i].position[1]) / denominator;
        corr[i].cy = (-1 * tmp_x * corr[i].position[2] + tmp_y * corr[i].position[0]) / denominator;
        //視野中心の逆変換（座標の変換） pl71のpl70
    }
}

std::vector<base_track_t>base_trans(std::vector<std::pair<base_track_t, corrmap>>& base_pair, double dz_nominal) {
    std::vector<base_track_t> ret;
    corrmap param;
    double rotation;
    double angle_shrink = 0.995;   //0.989;//大きい
    double dz_tune = -3; // 平均からのずれ
    for (int i = 0; i < base_pair.size(); i++) {
        base_track_t base;
        param = base_pair[i].second;
        base = base_pair[i].first;

        base.x = base_pair[i].first.x * param.position[0] + base_pair[i].first.y * param.position[1] + param.position[4];
        base.y = base_pair[i].first.x * param.position[2] + base_pair[i].first.y * param.position[3] + param.position[5];


        rotation = atan(param.position[2] / param.position[0]);
        base.ax = (base_pair[i].first.ax * cos(rotation) - base_pair[i].first.ay * sin(rotation)) * angle_shrink + param.angle[4];
        base.ay = (base_pair[i].first.ax * sin(rotation) + base_pair[i].first.ay * cos(rotation)) * angle_shrink + param.angle[5];

        base.z = dz_nominal + param.dz + dz_tune;
        ret.push_back(base);
    }
    return ret;
}


std::vector<output_format>  base_matching(std::vector <base_track_t>& base0, std::vector<base_track_t>& base1,
    double pos_lat_all[2], double pos_rad_all[2], double ang_lat_all[2], double ang_rad_all[2]) {

    double x_min = 0, y_min = 0;
    for (int i = 0; i < base1.size(); i++) {
        if (i == 0) {
            x_min = base1[i].x;
            y_min = base1[i].y;
        }
        x_min = std::min(base1[i].x, x_min);
        y_min = std::min(base1[i].y, y_min);
    }

    double hash_size = 2000;
    std::multimap<std::pair<int, int>, base_track_t>base1_hash;
    std::pair<int, int>id;
    double dz_mean = 0;
    for (int i = 0; i < base1.size(); i++) {
        id.first = (base1[i].x - x_min) / hash_size;
        id.second = (base1[i].y - y_min) / hash_size;
        base1_hash.insert(std::make_pair(id, base1[i]));
        dz_mean += base1[i].z;
    }
    dz_mean /= base1.size();
    printf("dz mean=%g\n", dz_mean);

    std::vector<output_format> out_base;

    double angle, ex_x, ey_y;
    double ix, iy;
    double diff_pos_rad, diff_pos_lat, diff_ang_rad, diff_ang_lat;
    std::vector<base_track_t>base_buf;
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
        out.aft.rawid = -1;
        out.aft.ax = 0;
        out.aft.ay = 0;
        out.aft.x = 0;
        out.aft.y = 0;



        angle = sqrt(base0[i].ax * base0[i].ax + base0[i].ay * base0[i].ay);
        ex_x = base0[i].x + base0[i].ax * dz_mean;
        ey_y = base0[i].y + base0[i].ay * dz_mean;
        ix = (ex_x - x_min) / hash_size;
        iy = (ey_y - y_min) / hash_size;
        out.ex.x = ex_x;
        out.ex.y = ey_y;


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
            ex_x = base0[i].x + base0[i].ax * base_buf[j].z;
            ey_y = base0[i].y + base0[i].ay * base_buf[j].z;

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
    fprintf(stderr, "\r base matchin %lld/%lld(%4.1lf%%)\n", count, base0.size(), count * 100. / base0.size());

    return out_base;

}

void pick_up_base_inf(std::vector<output_format>& out, std::vector <base_track_t>& base) {
    std::map<int, base_track_t> base_map;
    for (auto itr = base.begin(); itr != base.end(); itr++) {
        base_map.insert(std::make_pair(itr->rawid, *itr));
    }
    for (auto itr = out.begin(); itr != out.end(); itr++) {
        if (!itr->hit_flg)continue;
        if (base_map.count(itr->aft.rawid) == 0)continue;
        auto base_before = base_map.at(itr->aft.rawid);
        itr->bef.rawid = base_before.rawid;
        itr->bef.ax = base_before.ax;
        itr->bef.ay = base_before.ay;
        itr->bef.x = base_before.x;
        itr->bef.y = base_before.y;
    }
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
            // 外挿先
            << std::setw(10) << std::setprecision(0) << itr->ex.rawid << " "
            << std::setw(8) << std::setprecision(4) << itr->ex.ax << " "
            << std::setw(8) << std::setprecision(4) << itr->ex.ay << " "
            << std::setw(10) << std::setprecision(1) << itr->ex.x << " "
            << std::setw(10) << std::setprecision(1) << itr->ex.y << " "
            << std::setw(2) << std::setprecision(1) << itr->hit_flg << " "
            // 座標変換後
            << std::setw(10) << std::setprecision(0) << itr->aft.rawid << " "
            << std::setw(8) << std::setprecision(4) << itr->aft.ax << " "
            << std::setw(8) << std::setprecision(4) << itr->aft.ay << " "
            << std::setw(10) << std::setprecision(1) << itr->aft.x << " "
            << std::setw(10) << std::setprecision(1) << itr->aft.y << " "
            // before 座標変換
            << std::setw(10) << std::setprecision(0) << itr->bef.rawid << " "
            << std::setw(8) << std::setprecision(4) << itr->bef.ax << " "
            << std::setw(8) << std::setprecision(4) << itr->bef.ay << " "
            << std::setw(10) << std::setprecision(1) << itr->bef.x << " "
            << std::setw(10) << std::setprecision(1) << itr->bef.y << std::endl;
    }
    fprintf(stderr, "\r write result %d/%d(%4.1lf%%)\n", count, out.size(), count * 100. / out.size());

}