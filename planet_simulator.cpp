#define _USE_MINGW_ANSI_STDIO 1
#include <bits/stdc++.h>
#include <time.h>
using namespace std;

const long double eps = 1e-14;

// physics constants
const long double G = 6.67e-11;

char* file_name = (char *) "orbit";
FILE *file;
bool log_ = true;

bool slow_ = false;

struct Dot{
    long double x, y;
    Dot(): x(0), y(0){}
    Dot(long double __x, long double __y): x(__x), y(__y){}
    inline Dot friend operator+(Dot __x, Dot __y) {return Dot(__x.x + __y.x, __x.y + __y.y);}
    inline Dot friend operator-(Dot __x, Dot __y) {return Dot(__x.x - __y.x, __x.y - __y.y);}
    inline Dot friend operator*(long double __k, Dot __x) {return Dot(__x.x * __k, __x.y * __k);}
    inline Dot friend operator*(Dot __x, long double __k) {return Dot(__x.x * __k, __x.y * __k);}
    inline Dot friend operator/(Dot __x, long double __k) {return Dot(__x.x / __k, __x.y / __k);}
};

Dot abs(Dot __d) {return {abs(__d.x), abs(__d.y)};}

template<typename _Tp = int>
inline _Tp pow2(_Tp __x) {return __x*__x;}

long double geo_dis(Dot __a, Dot __b) {return sqrt(pow2(__a.x - __b.x) + pow2(__a.y - __b.y));}

struct Planet{
    int id;
    Dot p, v;
    long double mass;
};

const int N = 32;
Planet planets[N];

long long timestep = 4*60*60ll, log_gap = 4*24*60*60;

void init_args(int argc, const char * argv[]){
    for(int i = 1; i < argc; i++){
        if(i < argc - 1){
            if(!strcmp(argv[i], "-o")){
                file_name = (char *) argv[++i];
            }
            else if(!strcmp(argv[i], "-step")){
                sscanf(argv[++i], "%lld", &timestep);
            }
            else if(!strcmp(argv[i], "-log_gap")){
                sscanf(argv[++i], "%lld", &log_gap);
            }
        }
        if(!strcmp(argv[i], "-nolog")){
            log_ = false;
        }
        else if(!strcmp(argv[i], "-slow")){
            slow_ = true;
        }
    }
}

void adapt_p_verlet(vector<Dot>& accelerations) {
    for (int i = 0; i < N; i++) {
        planets[i].p = planets[i].p + planets[i].v * timestep + 0.5 * accelerations[i] * pow2(timestep);
    }
}

vector<Dot> compute_accelerations() {
    vector<Dot> accelerations(N, Dot(0, 0));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j) {
                Dot direction = planets[j].p - planets[i].p;
                long double distance = geo_dis(planets[i].p, planets[j].p);
                accelerations[i] = accelerations[i] + G * planets[j].mass * direction / pow2(distance) / distance;
            }
        }
    }
    return accelerations;
}

void adapt_v_verlet(const vector<Dot>& old_accelerations, const vector<Dot>& new_accelerations) {
    for (int i = 0; i < N; i++) {
        planets[i].v = planets[i].v + 0.5 * (old_accelerations[i] + new_accelerations[i]) * timestep;
    }
}

void iterate(long long T) {
    fprintf(file, "%lld\n", T / log_gap);
    long long time = 0, gaps = -1;
    vector<Dot> accelerations = compute_accelerations();
    while ((time += timestep) < T) {
        adapt_p_verlet(accelerations);
        vector<Dot> new_accelerations = compute_accelerations();
        adapt_v_verlet(accelerations, new_accelerations);
        accelerations = new_accelerations;

        if (time / log_gap != gaps) {
            gaps = time / log_gap;
            if (slow_) this_thread::sleep_for((chrono::microseconds(1)));
            if (log_) printf("%lld : %lld\n", time / log_gap, time);
            for (Planet plt : planets) {
                fprintf(file, "%.Lf\n%.Lf\n", plt.p.x / 1e9, plt.p.y / 1e9);
            }
        }
    }
}

int main(int argc, const char * argv[]){
    init_args(argc, argv);
    file = fopen(file_name, "w");
    if(file == nullptr){
        return -1;
    }

    planets[0].id = 1;
    planets[0].mass = 2e30;
    planets[0].p = Dot(0, 0);
    planets[0].v = Dot(0, 0);

    planets[1].id = 2;
    planets[1].mass = 6e24;
    planets[1].p = Dot(1.5e11, 0);
    planets[1].v = Dot(0, 2.978e4);
    
    // for(int i = 2; i < N; i++){
    //     planets[i].id = i + 1;
    //     planets[i].mass = 6e2;
    //     planets[i].p = Dot(1.5e11 + 1e10 * cos(2 * M_PI / (N - 1) * (i - 2)), 1e10 * sin(2 * M_PI / (N - 1) * (i - 2)));
    //     planets[i].v = Dot(-3e4 * sin(2 * M_PI / (N - 1) * (i - 2)), 3e4 * cos(2 * M_PI / (N - 1) * (i - 2)));
    // }

    for(int i = 2; i < N; i++){
        planets[i].id = i + 1;
        planets[i].mass = 6e2;
        planets[i].p = Dot(1.5e11 + 1e9 * cos(2 * M_PI / (N - 1) * (i - 2)), 1e9 * sin(2 * M_PI / (N - 1) * (i - 2)));
        planets[i].v = Dot(0*-3e3 * sin(2 * M_PI / (N - 1) * (i - 2)), 3e3 * cos(2 * M_PI / (N - 1) * (i - 2)) + 3e4);
    }

    fprintf(file, "%d\n", N);
    iterate(1ll*20*365*24*60*60); // 365,246,060
    return 0;
}
