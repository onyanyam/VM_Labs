import numpy as np
import matplotlib.pyplot as plt

def load_data(filename):
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.split()
            t = float(parts[0])
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])
            data.append([t, x, y, z])
    
    return np.array(data)

# загружаем данные
print("Введите название файла:")
filename = input().strip()
data = load_data(filename)
t = data[:, 0]
x = data[:, 1]
y = data[:, 2]
z = data[:, 3]

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

# График x(t)
ax1.plot(t, x, 'b-', linewidth=1, alpha=0.8)
ax1.set_ylabel('x(t)', fontsize=12, color='b')
ax1.tick_params(axis='y', labelcolor='b')
ax1.grid(True, alpha=0.3)
ax1.set_title('Аттрактор Рёсслера', fontsize=14)

# График y(t)
ax2.plot(t, y, 'r-', linewidth=1, alpha=0.8)
ax2.set_ylabel('y(t)', fontsize=12, color='r')
ax2.tick_params(axis='y', labelcolor='r')
ax2.grid(True, alpha=0.3)

# График z(t)
ax3.plot(t, z, 'g-', linewidth=1, alpha=0.8)
ax3.set_ylabel('z(t)', fontsize=12, color='g')
ax3.tick_params(axis='y', labelcolor='g')
ax3.grid(True, alpha=0.3)
ax3.set_xlabel('Время, t', fontsize=12)

plt.tight_layout()
plt.show()

# три отдельных окна (альтернатива)
print("\nОткрыть каждый график отдельно? (y/n):")
choice = input().strip().lower()

if choice == 'y':
    # График x(t)
    plt.figure(figsize=(10, 5))
    plt.plot(t, x, 'b-', linewidth=1, alpha=0.8)
    plt.xlabel('Время, t', fontsize=12)
    plt.ylabel('x(t)', fontsize=12)
    plt.title('Аттрактор Рёсслера - x(t)', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    
    # График y(t)
    plt.figure(figsize=(10, 5))
    plt.plot(t, y, 'r-', linewidth=1, alpha=0.8)
    plt.xlabel('Время, t', fontsize=12)
    plt.ylabel('y(t)', fontsize=12)
    plt.title('Аттрактор Рёсслера - y(t)', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    
    # График z(t)
    plt.figure(figsize=(10, 5))
    plt.plot(t, z, 'g-', linewidth=1, alpha=0.8)
    plt.xlabel('Время, t', fontsize=12)
    plt.ylabel('z(t)', fontsize=12)
    plt.title('Аттрактор Рёсслера - z(t)', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

# Фазовые портреты
print("\nПостроить фазовые портреты (плоскости)? (y/n):")
choice2 = input().strip().lower()

if choice2 == 'y':
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # XY проекция
    axes[0].plot(x, y, 'b-', linewidth=0.8, alpha=0.7)
    axes[0].set_xlabel('x', fontsize=12)
    axes[0].set_ylabel('y', fontsize=12)
    axes[0].set_title('Плоскость XY', fontsize=12)
    axes[0].grid(True, alpha=0.3)
    axes[0].set_aspect('equal')
    
    # XZ проекция
    axes[1].plot(x, z, 'r-', linewidth=0.8, alpha=0.7)
    axes[1].set_xlabel('x', fontsize=12)
    axes[1].set_ylabel('z', fontsize=12)
    axes[1].set_title('Плоскость XZ', fontsize=12)
    axes[1].grid(True, alpha=0.3)
    axes[1].set_aspect('equal')
    
    # YZ проекция
    axes[2].plot(y, z, 'g-', linewidth=0.8, alpha=0.7)
    axes[2].set_xlabel('y', fontsize=12)
    axes[2].set_ylabel('z', fontsize=12)
    axes[2].set_title('Плоскость YZ', fontsize=12)
    axes[2].grid(True, alpha=0.3)
    axes[2].set_aspect('equal')
    
    plt.suptitle('Фазовые портреты аттрактора Рёсслера', fontsize=14)
    plt.tight_layout()
    plt.show()