

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <random>
#include <cstdlib>
#include <string>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <map>

struct Config {
    int L = 64;                    
    int threshold = 4;             
    long long total_grains = 500000; 
    long long transient = 50000;   
    int snapshot_interval = 50000; 
    unsigned seed = 42;            
    std::string output_dir = ".";  
};

Config parse_args(int argc, char* argv[]) {
    Config cfg;
    for (int i = 1; i < argc; i += 2) {
        std::string key = argv[i];
        if (i + 1 >= argc) break;
        std::string val = argv[i + 1];
        if (key == "--L") cfg.L = std::stoi(val);
        else if (key == "--threshold") cfg.threshold = std::stoi(val);
        else if (key == "--grains") cfg.total_grains = std::stoll(val);
        else if (key == "--transient") cfg.transient = std::stoll(val);
        else if (key == "--seed") cfg.seed = std::stoul(val);
        else if (key == "--outdir") cfg.output_dir = val;
    }
    return cfg;
}


class Sandpile {
public:
    int L, threshold;
    std::vector<std::vector<int>> grid;
    std::mt19937 rng;

    
    static constexpr int dx[] = {0, 0, 1, -1};
    static constexpr int dy[] = {1, -1, 0, 0};

    Sandpile(int L, int threshold, unsigned seed)
        : L(L), threshold(threshold), grid(L, std::vector<int>(L, 0)), rng(seed) {}

    struct AvalancheResult {
        long long size;      
        int duration;         
        int unique_sites;     
        int max_radius;      
    };

    AvalancheResult add_grain() {
        std::uniform_int_distribution<int> dist(0, L - 1);
        int x0 = dist(rng), y0 = dist(rng);
        grid[x0][y0]++;
        return relax(x0, y0);
    }

    AvalancheResult relax(int x0, int y0) {
        AvalancheResult res{0, 0, 0, 0};
        if (grid[x0][y0] < threshold) return res;

        std::vector<std::vector<bool>> toppled(L, std::vector<bool>(L, false));
        std::vector<std::pair<int,int>> current;

        for (int i = 0; i < L; i++)
            for (int j = 0; j < L; j++)
                if (grid[i][j] >= threshold)
                    current.push_back({i, j});

        while (!current.empty()) {
            res.duration++;
            std::vector<std::pair<int,int>> next;

            for (size_t ci = 0; ci < current.size(); ci++) {
                int x = current[ci].first;
                int y = current[ci].second;
                if (grid[x][y] < threshold) continue;
                // Topple
                grid[x][y] -= threshold;
                res.size++;
                toppled[x][y] = true;

                int dist_from_origin = std::abs(x - x0) + std::abs(y - y0);
                res.max_radius = std::max(res.max_radius, dist_from_origin);

                for (int d = 0; d < 4; d++) {
                    int nx = x + dx[d], ny = y + dy[d];
                    if (nx >= 0 && nx < L && ny >= 0 && ny < L) {
                        grid[nx][ny]++;
                        if (grid[nx][ny] >= threshold)
                            next.push_back({nx, ny});
                    }
                    
                }
            }
            current = next;
        }

        for (int i = 0; i < L; i++)
            for (int j = 0; j < L; j++)
                if (toppled[i][j]) res.unique_sites++;

        return res;
    }


    double avg_height() const {
        double sum = 0;
        for (auto& row : grid) for (int v : row) sum += v;
        return sum / (L * L);
    }
};

constexpr int Sandpile::dx[];
constexpr int Sandpile::dy[];


int main(int argc, char* argv[]) {
    Config cfg = parse_args(argc, argv);
    std::cout << "=== BTW Sandpile: Earthquake Fault Dynamics ===" << std::endl;
    std::cout << "Lattice: " << cfg.L << "x" << cfg.L
              << "  Threshold: " << cfg.threshold
              << "  Grains: " << cfg.total_grains << std::endl;

    Sandpile sp(cfg.L, cfg.threshold, cfg.seed);

    
    std::ofstream f_sizes(cfg.output_dir + "/avalanche_sizes.csv");
    std::ofstream f_dur(cfg.output_dir + "/avalanche_durations.csv");
    std::ofstream f_conn(cfg.output_dir + "/connectivity_stats.csv");
    std::ofstream f_snap(cfg.output_dir + "/grid_snapshots.csv");

    f_sizes << "grain,size\n";
    f_dur   << "grain,duration\n";
    f_conn  << "grain,size,duration,unique_sites,max_radius\n";

    long long avalanche_count = 0;

    for (long long t = 1; t <= cfg.total_grains; t++) {
        auto res = sp.add_grain();

        
        if (t > cfg.transient && res.size > 0) {
            f_sizes << t << "," << res.size << "\n";
            f_dur   << t << "," << res.duration << "\n";
            f_conn  << t << "," << res.size << "," << res.duration
                    << "," << res.unique_sites << "," << res.max_radius << "\n";
            avalanche_count++;
        }

       
        if (t % cfg.snapshot_interval == 0) {
            f_snap << "# grain=" << t << " avg_height=" << sp.avg_height() << "\n";
            for (int i = 0; i < cfg.L; i++) {
                for (int j = 0; j < cfg.L; j++) {
                    if (j > 0) f_snap << ",";
                    f_snap << sp.grid[i][j];
                }
                f_snap << "\n";
            }
        }

        if (t % 100000 == 0)
            std::cout << "  Progress: " << t << "/" << cfg.total_grains
                      << "  avg_height=" << sp.avg_height() << std::endl;
    }

    f_sizes.close(); f_dur.close(); f_conn.close(); f_snap.close();

    std::cout << "Done. Avalanches recorded: " << avalanche_count << std::endl;
    return 0;
}
