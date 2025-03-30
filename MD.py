import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from ase.md import VelocityVerlet

# 创建H2 + O反应体系
system = Atoms("H2O", positions=[[0, 0, 0], [0.7, 0, 0], [0, 0.5, 1]])
system.calc = EMT()  # 使用经验力场
dyn = VelocityVerlet(system, timestep=1e-6)  # 时间步长1fs

# 运行模拟
for step in range(10000):
    dyn.run(10)
    print(f"Step {step}: Energy = {system.get_potential_energy():.3f} eV")
