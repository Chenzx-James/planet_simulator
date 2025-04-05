#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#define LogStep 1
#define Iterations 1000

FILE *file;

// 三维向量
struct Vec3 {
    double x, y, z;
    Vec3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
    Vec3 operator+(const Vec3& other) const { return Vec3(x + other.x, y + other.y, z + other.z); }
    Vec3 operator-(const Vec3& other) const { return Vec3(x - other.x, y - other.y, z - other.z); }
    Vec3 operator*(double s) const { return Vec3(x * s, y * s, z * s); }
    Vec3& operator+=(const Vec3& other) { x += other.x; y += other.y; z += other.z; return *this; }
    Vec3& operator-=(const Vec3& other) { x -= other.x; y -= other.y; z -= other.z; return *this; }
    double dot(const Vec3& other) const { return x * other.x + y * other.y + z * other.z; }
    Vec3 cross(const Vec3& other) const { return Vec3(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x); }
    double norm() const { return std::sqrt(x*x + y*y + z*z); }
    Vec3 normalized() const { double n = norm(); return Vec3(x/n, y/n, z/n); }
};

// 四元数
struct Quaternion {
    double w, x, y, z;
    Quaternion(double w = 1, double x = 0, double y = 0, double z = 0) : w(w), x(x), y(y), z(z) {}
    Quaternion conjugate() const { return Quaternion(w, -x, -y, -z); }
    Quaternion operator*(const Quaternion& other) const {
        return Quaternion(
            w*other.w - x*other.x - y*other.y - z*other.z,
            w*other.x + x*other.w + y*other.z - z*other.y,
            w*other.y - x*other.z + y*other.w + z*other.x,
            w*other.z + x*other.y - y*other.x + z*other.w
        );
    }
    Vec3 rotate(const Vec3& v) const {
        Quaternion qv(0, v.x, v.y, v.z);
        Quaternion res = (*this) * qv * this->conjugate();
        return Vec3(res.x, res.y, res.z);
    }
    void integrate(const Vec3& omega, double dt) {
        double angle = omega.norm() * dt;
        if (angle > 1e-10) {
            Vec3 axis = omega * (1.0 / angle);
            double s = std::sin(angle / 2);
            Quaternion dq(std::cos(angle / 2), axis.x * s, axis.y * s, axis.z * s);
            *this = dq * *this;
        } else {
            Quaternion dq(1, omega.x * dt / 2, omega.y * dt / 2, omega.z * dt / 2);
            *this = dq * *this;
        }
        normalize();
    }
    void normalize() {
        double len = std::sqrt(w*w + x*x + y*y + z*z);
        w /= len; x /= len; y /= len; z /= len;
    }
};

// 水分子结构
struct Water {
    Vec3 position;      // 质心位置（O原子）
    Vec3 velocity;
    Quaternion orient;  // 取向
    Vec3 angularVel;
    Vec3 inertia;       // 主转动惯量

    static constexpr double mass = 18.0;                // 总质量 (O:16, H2:2)
    static constexpr double chargeO = -0.82;            // O电荷 (e)
    static constexpr double chargeH = 0.41;             // H电荷 (e)
    static constexpr double sigma = 3.166;               // LJ参数 (Å)
    static constexpr double epsilon = 0.1553;            // kcal/mol
    static constexpr double bondLen = 1.0;               // O-H键长 (Å)
    static constexpr double angle = 109.47 * M_PI/180;   // 键角 (rad)

    Vec3 h1Rel, h2Rel;  // H相对质心的位置

    Water() {
        // 初始化H的相对位置
        h1Rel = Vec3(bondLen * std::sin(angle/2), 0, bondLen * std::cos(angle/2));
        h2Rel = Vec3(bondLen * std::sin(angle/2), 0, -bondLen * std::cos(angle/2));
        // 计算转动惯量
        double mH = 1.0;
        inertia.x = mH * (h1Rel.y*h1Rel.y + h1Rel.z*h1Rel.z) + mH * (h2Rel.y*h2Rel.y + h2Rel.z*h2Rel.z);
        inertia.y = mH * (h1Rel.x*h1Rel.x + h1Rel.z*h1Rel.z) + mH * (h2Rel.x*h2Rel.x + h2Rel.z*h2Rel.z);
        inertia.z = mH * (h1Rel.x*h1Rel.x + h1Rel.y*h1Rel.y) + mH * (h2Rel.x*h2Rel.x + h2Rel.y*h2Rel.y);
    }

    Vec3 h1Pos() const { return position + orient.rotate(h1Rel); }
    Vec3 h2Pos() const { return position + orient.rotate(h2Rel); }
};

class Simulation {
    std::vector<Water> waters;
    double boxSize;
    double cutoff;
    double temp;
    double timeStep;

    // 物理常数
    const double kCoulomb = 332.06371;  // kcal·Å/(mol·e²)
    const double kB = 0.0019872041;     // kcal/(mol·K)

public:
    Simulation(int nWaters, double box, double cut, double t, double dt)
        : boxSize(box), cutoff(cut), temp(t), timeStep(dt) {
        // 初始化水分子位置
        int perSide = static_cast<int>(std::ceil(std::cbrt(nWaters)));
        double spacing = box / perSide;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> velDist(0.0, std::sqrt(kB * temp / Water::mass));

        for (int i = 0; i < nWaters; ++i) {
            Water w;
            w.position = Vec3(
                ((i % perSide) * spacing + boxSize) / 3,
                (((i / perSide) % perSide) * spacing + boxSize) / 3,
                ((i / (perSide * perSide)) * spacing + boxSize) / 3
            );
            // 初始速度
            w.velocity = Vec3(velDist(gen)*1, velDist(gen)*1, velDist(gen)*1);
            waters.push_back(w);
        }
    }

    int analyzeHBonds(const std::vector<Water>& waters, double maxDist = 3.0, double minAngle = 150.0) {
        int hbonds = 0;
        for (size_t i = 0; i < waters.size(); ++i) {
            const Water& wi = waters[i];
            Vec3 oi = wi.position;
            Vec3 h1i = wi.h1Pos(), h2i = wi.h2Pos();
            for (size_t j = 0; j < waters.size(); ++j) {
                if (i == j) continue;
                const Water& wj = waters[j];
                Vec3 oj = wj.position;
                // O-O距离
                double ooDist = (oi - oj).norm();
                if (ooDist > maxDist) continue;
                // 检查H是否指向O
                for (const Vec3& hi : {h1i, h2i}) {
                    Vec3 oh = hi - oi;
                    Vec3 oo = oj - oi;
                    double angle = std::acos(oh.dot(oo) / (oh.norm() * oo.norm())) * 180 / M_PI;
                    if (angle >= minAngle) {
                        hbonds++;
                    }
                }
            }
        }
        return hbonds;
        printf("Hydrogen bonds: %d\n", hbonds);
    }

    void applyPBC(Vec3& pos) {
        pos.x = std::fmod(pos.x + boxSize, boxSize);
        pos.y = std::fmod(pos.y + boxSize, boxSize);
        pos.z = std::fmod(pos.z + boxSize, boxSize);
    }

    void computeForces(std::vector<Vec3>& forcesO, std::vector<Vec3>& forcesH1, std::vector<Vec3>& forcesH2) {
        std::fill(forcesO.begin(), forcesO.end(), Vec3());
        std::fill(forcesH1.begin(), forcesH1.end(), Vec3());
        std::fill(forcesH2.begin(), forcesH2.end(), Vec3());

        for (size_t i = 0; i < waters.size(); ++i) {
            const Water& wi = waters[i];
            Vec3 oi = wi.position;
            Vec3 h1i = wi.h1Pos();
            Vec3 h2i = wi.h2Pos();

            for (size_t j = i + 1; j < waters.size(); ++j) {
                const Water& wj = waters[j];
                Vec3 oj = wj.position;
                Vec3 h1j = wj.h1Pos();
                Vec3 h2j = wj.h2Pos();

                // LJ (O-O)
                Vec3 deltaOO = oi - oj;
                double rOO = deltaOO.norm();
                if (rOO < cutoff) {
                    double r6 = std::pow(Water::sigma / rOO, 6);
                    double ljForce = 4 * Water::epsilon * (12 * r6*r6 / rOO - 6 * r6 / rOO) / rOO;
                    forcesO[i] += deltaOO * (ljForce / rOO);
                    forcesO[j] -= deltaOO * (ljForce / rOO);
                }

                // Coulomb interactions
                auto addCoulomb = [&](const Vec3& r1, double q1, const Vec3& r2, double q2, Vec3& f1, Vec3& f2) {
                    Vec3 dr = r1 - r2;
                    double r = dr.norm();
                    if (r < cutoff) {
                        double force = kCoulomb * q1 * q2 / (r * r);
                        Vec3 f = dr * (force / r);
                        f1 += f;
                        f2 -= f;
                    }
                };

                // Oi与所有其他原子
                addCoulomb(oi, Water::chargeO, oj, Water::chargeO, forcesO[i], forcesO[j]);
                addCoulomb(oi, Water::chargeO, h1j, Water::chargeH, forcesO[i], forcesH1[j]);
                addCoulomb(oi, Water::chargeO, h2j, Water::chargeH, forcesO[i], forcesH2[j]);
                addCoulomb(h1i, Water::chargeH, oj, Water::chargeO, forcesH1[i], forcesO[j]);
                addCoulomb(h1i, Water::chargeH, h1j, Water::chargeH, forcesH1[i], forcesH1[j]);
                addCoulomb(h1i, Water::chargeH, h2j, Water::chargeH, forcesH1[i], forcesH2[j]);
                addCoulomb(h2i, Water::chargeH, oj, Water::chargeO, forcesH2[i], forcesO[j]);
                addCoulomb(h2i, Water::chargeH, h1j, Water::chargeH, forcesH2[i], forcesH1[j]);
                addCoulomb(h2i, Water::chargeH, h2j, Water::chargeH, forcesH2[i], forcesH2[j]);
            }
        }
    }

    void integrate() {
        std::vector<Vec3> forcesO(waters.size()), forcesH1(waters.size()), forcesH2(waters.size());
        computeForces(forcesO, forcesH1, forcesH2);

        for (size_t i = 0; i < waters.size(); ++i) {
            Water& w = waters[i];
            // 平动
            Vec3 totalForce = forcesO[i] + forcesH1[i] + forcesH2[i];
            Vec3 accel = totalForce * (1.0 / Water::mass);
            w.velocity += accel * (timeStep / 2);
            w.position += w.velocity * timeStep;
            applyPBC(w.position);

            // 计算力矩
            Vec3 rH1 = w.orient.rotate(w.h1Rel);
            Vec3 rH2 = w.orient.rotate(w.h2Rel);
            Vec3 torque = rH1.cross(forcesH1[i]) + rH2.cross(forcesH2[i]);

            // 转换到本体坐标系
            Quaternion inv = w.orient.conjugate();
            Vec3 torqueBody = inv.rotate(torque);
            Vec3 omegaBody = inv.rotate(w.angularVel);

            // 欧拉方程计算角加速度
            Vec3 alphaBody;
            alphaBody.x = (torqueBody.x - (omegaBody.y * omegaBody.z * (w.inertia.z - w.inertia.y))) / w.inertia.x;
            alphaBody.y = (torqueBody.y - (omegaBody.z * omegaBody.x * (w.inertia.x - w.inertia.z))) / w.inertia.y;
            alphaBody.z = (torqueBody.z - (omegaBody.x * omegaBody.y * (w.inertia.y - w.inertia.x))) / w.inertia.z;

            Vec3 alpha = w.orient.rotate(alphaBody);

            // 更新角速度
            w.angularVel += alpha * (timeStep / 2);
            w.orient.integrate(w.angularVel, timeStep);
            w.angularVel += alpha * (timeStep / 2);
        }
    }

    void integrateStep1() {
        std::vector<Vec3> forcesO_prev(waters.size()), forcesH1_prev(waters.size()), forcesH2_prev(waters.size());
        // 第一次力计算
        computeForces(forcesO_prev, forcesH1_prev, forcesH2_prev);

        // 更新平动速度半步
        for (size_t i = 0; i < waters.size(); ++i) {
            Water& w = waters[i];
            Vec3 totalForce = forcesO_prev[i] + forcesH1_prev[i] + forcesH2_prev[i];
            Vec3 accel = totalForce * (1.0 / Water::mass);
            w.velocity += accel * (timeStep / 2); // 平动速度半步更新
        }

        // 更新旋转速度半步
        for (size_t i = 0; i < waters.size(); ++i) {
            Water& w = waters[i];
            // 计算力矩
            Vec3 rH1 = w.orient.rotate(w.h1Rel);
            Vec3 rH2 = w.orient.rotate(w.h2Rel);
            Vec3 torque = rH1.cross(forcesH1_prev[i]) + rH2.cross(forcesH2_prev[i]);
            // 转换到本体坐标系
            Quaternion inv = w.orient.conjugate();
            Vec3 torqueBody = inv.rotate(torque);
            Vec3 omegaBody = inv.rotate(w.angularVel);
            // 计算角加速度
            Vec3 alphaBody(
                (torqueBody.x - (omegaBody.y * omegaBody.z * (w.inertia.z - w.inertia.y))) / w.inertia.x,
                (torqueBody.y - (omegaBody.z * omegaBody.x * (w.inertia.x - w.inertia.z))) / w.inertia.y,
                (torqueBody.z - (omegaBody.x * omegaBody.y * (w.inertia.y - w.inertia.x))) / w.inertia.z
            );
            Vec3 alpha = w.orient.rotate(alphaBody);
            // 角速度半步更新
            w.angularVel += alpha * (timeStep / 2);
        }
    }

    void integrateStep2() {
        // 更新位置和方向
        for (size_t i = 0; i < waters.size(); ++i) {
            Water& w = waters[i];
            w.position += w.velocity * timeStep; // 平动位置更新
            applyPBC(w.position);
            w.orient.integrate(w.angularVel, timeStep); // 方向更新
        }

        // 第二次力计算
        std::vector<Vec3> forcesO_new(waters.size()), forcesH1_new(waters.size()), forcesH2_new(waters.size());
        computeForces(forcesO_new, forcesH1_new, forcesH2_new);

        // 更新平动速度后半步
        for (size_t i = 0; i < waters.size(); ++i) {
            Water& w = waters[i];
            Vec3 totalForce = forcesO_new[i] + forcesH1_new[i] + forcesH2_new[i];
            Vec3 accel = totalForce * (1.0 / Water::mass);
            w.velocity += accel * (timeStep / 2); // 平动速度后半步更新
        }
    }

    void applyThermostat() {
        double totalKE = 0;
        for (auto& w : waters) {
            totalKE += 0.5 * Water::mass * w.velocity.dot(w.velocity);
            Vec3 omegaBody = w.orient.conjugate().rotate(w.angularVel);
            totalKE += 0.5 * (w.inertia.x * omegaBody.x*omegaBody.x + w.inertia.y * omegaBody.y*omegaBody.y + w.inertia.z * omegaBody.z*omegaBody.z);
        }
        double currentTemp = totalKE / (3 * waters.size() * kB);
        double scale = std::sqrt(temp / currentTemp);
        for (auto& w : waters) {
            w.velocity = w.velocity * scale;
            w.angularVel = w.angularVel * scale;
        }
    }

    void run() {
        fprintf(file, "%lu\n", waters.size()*3);
        fprintf(file, "%.2lf\n", boxSize);
        fprintf(file, "%d\n", Iterations/LogStep-1);
        for (int step = 0; step < Iterations; ++step) {
            // integrate();
            integrateStep1();  // 第一步：更新速度半步
            integrateStep2();  // 第二步：更新位置+方向，重新计算力，完成速度更新
            if (step % 1 == 0) applyThermostat();

            // 输出轨迹
            if (step % LogStep == 0) {
                for (auto& w : waters) {
                    double scale = 1.0;
                    fprintf(file, "%.2lf\n%.2lf\n%.2lf\n", (w.position.x-boxSize/2)*scale, (w.position.y-boxSize/2)*scale, (w.position.z-boxSize/2)*scale);
                    Vec3 h1 = w.h1Pos(), h2 = w.h2Pos();
                    fprintf(file, "%.2lf\n%.2lf\n%.2lf\n", (h1.x-boxSize/2)*scale, (h1.y-boxSize/2)*scale, (h1.z-boxSize/2)*scale);
                    fprintf(file, "%.2lf\n%.2lf\n%.2lf\n", (h2.x-boxSize/2)*scale, (h2.y-boxSize/2)*scale, (h2.z-boxSize/2)*scale);
                }
                fprintf(file, "%d\n", analyzeHBonds(waters));
            }
        }
    }
};

int main() {
    // freopen("t", "w", stdout);
    file = fopen("t", "w");
    // 300K=26°C
    // 400K=127°C
    Simulation sim(10, 5.0, 10.0, 300.0, 0.001); // 100水分子，30Å盒子，10Å截断，300K，0.001ps步长
    sim.run();
    return 0;
}
