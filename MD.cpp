#include <bits/stdc++.h>
#define SimpleOutput 1
#define Iterations 10000
#define LogStep 1

// 基本常量
const double kB = 1.38064852e-23; // 玻尔兹曼常数
const double NA = 6.02214076e23;  // 阿伏伽德罗常数
const double e = 1.6e-19; // 元素电荷(1.6021766208（98）×10-19 C)

// 原子/分子类
class Atom {
public:
    std::string element;
    double mass;    // 质量 (kg)
    double charge;  // 电荷
    std::vector<double> position; // 位置 (x,y,z)
    std::vector<double> velocity; // 速度 (vx,vy,vz)
    std::vector<double> force;    // 受力 (fx,fy,fz)
    
    Atom(std::string el, double m, double q) 
        : element(el), mass(m), charge(q), position(3,0), velocity(3,0), force(3,0) {}
};

// 分子系统类
class MolecularSystem {
private:
    std::vector<Atom> atoms;
    double temperature; // 温度 (K)
    double timeStep;    // 时间步长 (s)
    
public:
    MolecularSystem(double temp, double dt) : temperature(temp), timeStep(dt) {}
    
    void addAtom(const Atom& atom) {
        atoms.push_back(atom);
    }

    size_t getAtomCount() const {
        return atoms.size();
    }
    
    // 计算力场 (简化的Lennard-Jones势能)
    void calculateForces() {
        const double epsilon = 1.65e-21; // 势阱深度 (J)
        const double sigma = 3.4e-10;    // 碰撞直径 (m)
        const double k_e = 8.9875517923e9; // 库仑常数 (N·m²/C²)
        
        // 重置所有力为零
        for (auto& atom : atoms) {
            std::fill(atom.force.begin(), atom.force.end(), 0.0);
        }
        
        // 计算原子间作用力
        for (size_t i = 0; i < atoms.size(); ++i) {
            for (size_t j = i+1; j < atoms.size(); ++j) {
                double dx = atoms[j].position[0] - atoms[i].position[0];
                double dy = atoms[j].position[1] - atoms[i].position[1];
                double dz = atoms[j].position[2] - atoms[i].position[2];
                
                double r = sqrt(dx*dx + dy*dy + dz*dz);
                
                if (r > 0) {
                    // Lennard-Jones力
                    double ljForce = 24 * epsilon * (2*pow(sigma/r, 13) - pow(sigma/r, 7)) / r;
                    
                    // 静电力 (库仑力)
                    double electrostaticForce = -k_e * atoms[i].charge * atoms[j].charge / (r * r);
                    
                    // 总力
                    double f = ljForce + electrostaticForce;
                    f = fmin(f, 1e-15); // 限制最大力
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
            std::vector<double> oldForce = atom.force;
            
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
        double totalKE = 0.0;
        for (const auto& atom : atoms) {
            double v2 = atom.velocity[0]*atom.velocity[0] + 
                         atom.velocity[1]*atom.velocity[1] + 
                         atom.velocity[2]*atom.velocity[2];
            totalKE += 0.5 * atom.mass * v2;
        }
        
        double currentTemp = 2.0 * totalKE / (3.0 * kB * atoms.size());
        double scaleFactor = sqrt(temperature / currentTemp);
        
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
        // std::normal_distribution<double> dist(0.0, sqrt(kB * temperature));
        
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
                        printf("%lf\n%lf\n", atom.position[0]*1e11*2,
                                            atom.position[1]*1e11*2);
                    }
                }
                else{
                    printf("Step %d:\n", step);
                    for (const auto& atom : atoms) {
                        printf("Atom %s: Position (%lf, %lf, %lf)\n", atom.element.c_str(), 
                            atom.position[0]*1e10, atom.position[1]*1e10, atom.position[2]*1e10);
                        printf("Velocity (%lf, %lf, %lf)\n",
                            atom.velocity[0]*1e10, atom.velocity[1]*1e10, atom.velocity[2]*1e10);
                        printf("Force (%lf, %lf, %lf)\n",
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
    // 创建系统: 温度300K，时间步长1fs
    MolecularSystem system(300.0, 1e-15);
    
    // 添加一些原子 (示例: 氧和氢)
    Atom o1("O", 2.656e-26, -2*e); // 氧原子质量
    Atom h1("H", 1.673e-27, e); // 氢原子质量
    Atom h2("H", 1.673e-27, e);

    // 添加CO2分子
    Atom c1("C", 1.993e-26, 0); // 碳原子质量
    Atom o2("O", 2.656e-26, -2*e); // 氧原子质量
    Atom o3("O", 2.656e-26, -2*e); // 氧原子质量
    
    // 设置初始位置 (水分子构型)
    o1.position = {0.0, 0.0, 0.0};
    h1.position = {0.957e-10, 0.0, 0.0};
    h2.position = {-0.239e-10, 0.927e-10, 0.0};

    c1.position = {0.0, 0.0, 0.0};
    o2.position = {-1.16e-10, 0.0, 0.0};
    o3.position = {1.16e-10, 0, 0.0};
    
    // system.addAtom(o1);
    // system.addAtom(h1);
    // system.addAtom(h2);
    
    system.addAtom(c1);
    system.addAtom(o2);
    system.addAtom(o3);

    // 运行1000步模拟
    printf("%d\n", int(system.getAtomCount()));
    printf("%d\n", Iterations / LogStep);
    system.run(Iterations);
    
    return 0;
}