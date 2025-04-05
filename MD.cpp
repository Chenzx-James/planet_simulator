#include <bits/stdc++.h>
#define SimpleOutput 1
#define Iterations 1000
#define LogStep 1

FILE *b;

// 基本常量
const long double kB = 1.38064852e-23; // 玻尔兹曼常数
const long double NA = 6.02214076e23;  // 阿伏伽德罗常数
const long double e = 1.6e-19; // 元素电荷(1.6021766208（98）×10-19 C)
const long double epsilon = 1.65e-21; // 势阱深度 (J)
const long double sigma = 3.4e-10;    // 碰撞直径 (m)
const long double k_e = 8.9875517923e9; // 库仑常数 (N·m²/C²)

// 原子/分子类
class Atom {
public:
    std::string element;
    long double mass;    // 质量 (kg)
    long double charge;  // 电荷
    std::vector<long double> position; // 位置 (x,y,z)
    std::vector<long double> velocity; // 速度 (vx,vy,vz)
    std::vector<long double> force;    // 受力 (fx,fy,fz)
    
    Atom(std::string el, long double m, long double q) 
        : element(el), mass(m), charge(q), position(3,0), velocity(3,0), force(3,0) {}
    
    Atom() 
        : element(0), mass(0), charge(0), position(3,0), velocity(3,0), force(3,0) {}
};

// 分子系统类
class MolecularSystem {
private:
    std::vector<Atom> atoms;
    long double temperature; // 温度 (K)
    long double timeStep;    // 时间步长 (s)
    
public:
    MolecularSystem(long double temp, long double dt) : temperature(temp), timeStep(dt) {}
    
    void addAtom(const Atom& atom) {
        atoms.push_back(atom);
    }

    size_t getAtomCount() const {
        return atoms.size();
    }
    
    // 计算力场 (简化的Lennard-Jones势能)
    void calculateForces() {
        // 重置所有力为零
        for (auto& atom : atoms) {
            std::fill(atom.force.begin(), atom.force.end(), 0.0);
        }
        
        // 计算原子间作用力
        for (size_t i = 0; i < atoms.size(); ++i) {
            for (size_t j = i+1; j < atoms.size(); ++j) {
                long double dx = atoms[j].position[0] - atoms[i].position[0];
                long double dy = atoms[j].position[1] - atoms[i].position[1];
                long double dz = atoms[j].position[2] - atoms[i].position[2];
                
                long double r = sqrt(dx*dx + dy*dy + dz*dz);
                
                if (r > 0) {
                    // Lennard-Jones力
                    long double ljForce = 24 * epsilon * (2*pow(sigma/r, 13) - pow(sigma/r, 7)) / r;
                    
                    // 静电力 (库仑力)
                    long double electrostaticForce = -k_e * atoms[i].charge * atoms[j].charge / (r * r);
                    
                    // 总力
                    // no lj
                    long double f = electrostaticForce;
                    f = fmin(f, 1e-7); // 限制最大力
                    // printf("Force of %s and %s: %e %e\n", atoms[i].element.c_str(), atoms[j].element.c_str(), ljForce, electrostaticForce);
                    
                    // 分配力到各分量
                    atoms[i].force[0] += f * dx/r;
                    atoms[i].force[1] += f * dy/r;
                    atoms[i].force[2] += f * dz/r;
                    
                    atoms[j].force[0] -= f * dx/r;
                    atoms[j].force[1] -= f * dy/r;
                    atoms[j].force[2] -= f * dz/r;
                }
            }
        }
    }
    
    // 更新位置和速度 (Verlet算法)
    void integrate() {
        for (auto& atom : atoms) {
            // 更新位置: r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2
            for (int i = 0; i < 3; ++i) {
                atom.position[i] += atom.velocity[i] * timeStep + 
                                   0.5 * (atom.force[i]/atom.mass) * timeStep * timeStep;
            }
            
            // 保存当前力用于下一步
            std::vector<long double> oldForce = atom.force;
            
            // 计算新力
            calculateForces();
            
            // 更新速度: v(t+dt) = v(t) + 0.5*(a(t)+a(t+dt))*dt
            for (int i = 0; i < 3; ++i) {
                atom.velocity[i] += 0.5 * (oldForce[i] + atom.force[i]) / atom.mass * timeStep;
            }
        }
    }
    
    // 温度控制 (简单的速度缩放)
    void thermostat() {
        long double totalKE = 0.0;
        for (const auto& atom : atoms) {
            long double v2 = atom.velocity[0]*atom.velocity[0] + 
                         atom.velocity[1]*atom.velocity[1] + 
                         atom.velocity[2]*atom.velocity[2];
            totalKE += 0.5 * atom.mass * v2;
        }
        
        long double currentTemp = 2.0 * totalKE / (3.0 * kB * atoms.size());
        long double scaleFactor = sqrt(temperature / currentTemp);
        
        for (auto& atom : atoms) {
            for (int i = 0; i < 3; ++i) {
                atom.velocity[i] *= scaleFactor;
            }
        }
    }
    
    // 运行模拟
    void run(int steps) {
        // 初始化随机速度
        // std::random_device rd;
        // std::mt19937 gen(rd());
        // std::normal_distribution<long double> dist(0.0, sqrt(kB * temperature));
        
        // for (auto& atom : atoms) {
        //     for (int i = 0; i < 3; ++i) {
        //         atom.velocity[i] = dist(gen) / sqrt(atom.mass);
        //     }
        // }
        
        // 主循环
        for (int step = 0; step < steps; ++step) {
            integrate();
            thermostat();
            // 输出当前状态到文件x.txt
            if(step%LogStep == 0){
                if (SimpleOutput){
                    for (const auto& atom : atoms) {
                        printf("%.0Lf\n%.0Lf\n", atom.position[0]*1e11*2,
                                            atom.position[1]*1e11*2);
                    }
                    // 在b.txt中输出atom:force
                    // for (const auto& atom : atoms) {
                    //     fprintf(b, "Atom %s: Force (%LE)\n", atom.element.c_str(), sqrt(atom.force[0]*atom.force[0] + atom.force[1]*atom.force[1] + atom.force[2]*atom.force[2]));
                    // }
                    // fprintf(b, "\n");
                }
                else{
                    printf("Step %d:\n", step);
                    for (const auto& atom : atoms) {
                        printf("Atom %s: Position (%Lf, %Lf, %Lf)\n", atom.element.c_str(), 
                            atom.position[0]*1e10, atom.position[1]*1e10, atom.position[2]*1e10);
                        printf("Velocity (%Lf, %Lf, %Lf)\n",
                            atom.velocity[0]*1e10, atom.velocity[1]*1e10, atom.velocity[2]*1e10);
                        printf("Force (%Lf, %Lf, %Lf)\n",
                            atom.force[0]*1e10, atom.force[1]*1e10, atom.force[2]*1e10);
                    }
                    printf("\n");
                }
            }
            // print force
            // 可以在这里添加反应检测和处理的代码
        }
    }
};

int main() {
    freopen("x.txt", "w", stdout);
    // 打开b.txt
    b = fopen("b.txt", "w");
    // 创建系统: 温度300K，时间步长1fs
    MolecularSystem system(300.0, 1e-15);
    
    // 随机10个Oxygen和20个Hydrogen原子
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<long double> dist(-10e-10, 10e-10);
    std::uniform_real_distribution<long double> dist1(-1, 1);
    for (int i = 0; i < 10; ++i) {
        Atom o("O", 2.656e-26, -2*e); // 氧原子质量
        o.position[0] = dist(gen);
        o.position[1] = dist(gen);
        o.position[2] = dist(gen);
        system.addAtom(o);
        Atom h1("H", 1.673e-27, e); // 氢原子质量
        h1.position[0] = dist(gen);
        h1.position[1] = dist(gen);
        h1.position[2] = dist(gen);
        system.addAtom(h1);
        Atom h2("H", 1.673e-27, e); // 氢原子质量
        h2.position[0] = dist(gen);
        h2.position[1] = dist(gen);
        h2.position[2] = dist(gen);
        system.addAtom(h2);
    }

    // 添加一些原子 (示例: 氧和氢)
    Atom o1("O", 2.656e-26, -2*e); // 氧原子质量
    Atom h1("H", 1.673e-27, e); // 氢原子质量
    Atom h2("H", 1.673e-27, e);

    // 添加CO2分子
    Atom c1("C", 1.993e-26, 4*e); // 碳原子质量
    Atom o2("O", 2.656e-26, -2*e); // 氧原子质量
    Atom o3("O", 2.656e-26, -2*e); // 氧原子质量
    
    // 设置初始位置 (水分子构型)
    o1.position = {0.0, 0.0, 0.0};
    h1.position = {0.957e-10, 0.0, 0.0};
    h2.position = {-0.239e-10, 0.927e-10, 0.0};


    c1.position = {0.0, 0.0, 0.0};
    o2.position = {-1.16e-10, 0.0, 0.0};
    o3.position = {1.16e-10, 0, 0.0};
    
    system.addAtom(o1);
    system.addAtom(h1);
    system.addAtom(h2);
    
    // system.addAtom(c1);
    // system.addAtom(o2);
    // system.addAtom(o3);

    // 运行1000步模拟
    printf("%d\n", int(system.getAtomCount()));
    printf("%d\n", Iterations / LogStep);
    system.run(Iterations);
    
    return 0;
}