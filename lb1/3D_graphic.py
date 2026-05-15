import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Загружаем данные из файла
def load_data(filename):
    """Загрузка данных из файла"""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.split()
            if len(parts) >= 4:
                t = float(parts[0])
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                data.append([t, x, y, z])
    
    return np.array(data)

# Загружаем данные
print("Введите название файла:")
filename = input()
data = load_data(filename)
t = data[:, 0]
x = data[:, 1]
y = data[:, 2]
z = data[:, 3]

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Рисуем траекторию линией
ax.plot(x, y, z, 'b-', linewidth=0.8, alpha=0.7, label='Траектория')

# Отмечаем начало и конец
ax.scatter(x[0], y[0], z[0], c='green', s=50, marker='o', label='Начальная точка')
ax.scatter(x[-1], y[-1], z[-1], c='red', s=50, marker='o', label='Конечная точка')

ax.plot([0, 0], [0, 0], [0, z.max()], 'k-', linewidth=2)  # ось Z
ax.plot([0, 0], [0, y.max()], [0, 0], 'k-', linewidth=2)  # ось Y
ax.plot([0, x.max()], [0, 0], [0, 0], 'k-', linewidth=2)  # ось X

ax.set_xlabel('X', fontsize=15)
ax.set_ylabel('Y', fontsize=15)
ax.set_zlabel('Z', fontsize=15)
ax.set_title('3D график Аттрактора Рёслера', fontsize=14)
ax.legend(fontsize=12)

# Настраиваем угол обзора
ax.view_init(elev=30, azim=60)

plt.tight_layout()
plt.savefig('rossler_3d.png', dpi=150, bbox_inches='tight')
plt.show()