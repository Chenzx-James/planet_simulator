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
    inline long double friend operator*(Dot __x, Dot __y) {return __x.x * __y.x + __x.y * __y.y;}
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

const int N = 4;
Planet planets[N];

const long long time_step = 10*60ll, log_gap = 4*24*60*60;

void init_args(int argc, const char * argv[]){
    for(int i = 1; i < argc; i++){
        if(i < argc - 1){
            if(!strcmp(argv[i], "-o")){
                file_name = (char *) argv[++i];
            }
            else if(!strcmp(argv[i], "-step")){
                sscanf(argv[++i], "%lld", time_step);
            }
            else if(!strcmp(argv[i], "-log_gap")){
                sscanf(argv[++i], "%lld", log_gap);
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

void adapt_p(){
    for(Planet *plt = planets; plt - planets < N; plt++){
        plt->p = plt->p + Dot(plt->v.x, plt->v.y) * time_step;
    }
}

void adapt_v(){
    for(Planet *plt = planets; plt - planets < N; plt++){
        for(Planet *x = planets; x - planets < N; x++){
            if(x != plt){
                plt->v = plt->v + G*x->mass*(x->p-plt->p)/(pow2(plt->p.x - x->p.x) + pow2(plt->p.y - x->p.y))/geo_dis(plt->p, x->p)* time_step;
            }
        }
//        while(sqrtl(plt->v * plt->v) > 2e5) plt->v = plt->v / 2;
    }
}

void iterate(long long T){
    fprintf(file, "%lld\n", T / log_gap);
    long long time = 0, gaps = -1;
    while((time += time_step) < T){
        adapt_p();
        adapt_v();
        if(time / log_gap != gaps){
            gaps = time / log_gap;
            if(slow_) this_thread::sleep_for((chrono::microseconds(1)));
        	if(log_) printf("%lld : %lld\n", time / log_gap, time);
            fprintf(file, "%lld\n", time);
            for(Planet plt : planets){
                fprintf(file, "%Lf\n%Lf\n", plt.p.x, plt.p.y);
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
    planets[0].v = Dot(-1e3, 9e2);
    planets[1].id = 2;
    planets[1].mass = 6e24;
    planets[1].p = Dot(1.5e11, 0);
    planets[1].v = Dot(0, 2.978e4);
    planets[2].id = 3;
    planets[2].mass = 3e28;
    planets[2].p = Dot(0, 1.0e11);
    planets[2].v = Dot(3e4, 0);
    planets[3].id = 4;
    planets[3].mass = 9e28;
    planets[3].p = Dot(1.8e11, -6576466);
    planets[3].v = Dot(0, -3.5e4);
    iterate(30*365*24*60*60); // 365,246,060
    system("pause");
    return 0;
}
