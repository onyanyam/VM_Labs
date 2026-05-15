import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection

def load_data(filename):
    try:
        data = np.loadtxt(filename, comments='#')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        
        t = data[:, 0]
        x = data[:, 1]
        y = data[:, 2]
        z = data[:, 3]
        
        return t, x, y, z
    except Exception as e:
        print(f"Ошибка при загрузке файла: {e}")
        return None, None, None, None

def create_3d_animation_gradient_trail(filename, skip=1, interval=10, linewidth=1.5, alpha=0.8):
    # Загрузка данных
    t, x, y, z = load_data(filename)
    
    if t is None:
        return
    
    # Прореживание данных для ускорения
    x = x[::skip]
    y = y[::skip]
    z = z[::skip]
    t = t[::skip]
    
    print(f"Загружено {len(x)} точек")
    print(f"Диапазон x: [{x.min():.2f}, {x.max():.2f}]")
    print(f"Диапазон y: [{y.min():.2f}, {y.max():.2f}]")
    print(f"Диапазон z: [{z.min():.2f}, {z.max():.2f}]")
    
    # Создание фигуры
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Настройка внешнего вида
    ax.set_xlabel('X', fontsize=12, labelpad=10)
    ax.set_ylabel('Y', fontsize=12, labelpad=10)
    ax.set_zlabel('Z', fontsize=12, labelpad=10)
    ax.set_title('Rossler Attractor - Gradient Trail', fontsize=14, pad=20)
    
    # Текущая точка
    current_point, = ax.plot([], [], [], 'ro', markersize=8, alpha=0.9,
                             markeredgecolor='black', markeredgewidth=1)
    
    # Коллекция для линий с градиентом
    trail_collection = None
    
    # Текст с информацией
    info_text = ax.text2D(0.02, 0.95, '', transform=ax.transAxes, 
                         fontsize=10, bbox=dict(boxstyle="round,pad=0.3", 
                                               facecolor='white', alpha=0.8))
    
    # Прогресс
    progress_text = ax.text2D(0.02, 0.02, '', transform=ax.transAxes, 
                             fontsize=10, bbox=dict(boxstyle="round,pad=0.3",
                                                   facecolor='white', alpha=0.8))
    
    # Установка границ с отступом
    margin = 0.1
    ax.set_xlim(x.min() - margin*(x.max()-x.min()), x.max() + margin*(x.max()-x.min()))
    ax.set_ylim(y.min() - margin*(y.max()-y.min()), y.max() + margin*(y.max()-y.min()))
    ax.set_zlim(z.min() - margin*(z.max()-z.min()), z.max() + margin*(z.max()-z.min()))
    
    # Начальный угол обзора
    ax.view_init(elev=25, azim=-60)
    
    def init():
        current_point.set_data([], [])
        current_point.set_3d_properties([])
        info_text.set_text('')
        progress_text.set_text('')
        return current_point, info_text, progress_text
    
    def update(frame):
        nonlocal trail_collection
        
        # Удаляем старую коллекцию
        if trail_collection is not None:
            trail_collection.remove()
        
        # Создаем сегменты для Line3DCollection
        points = np.column_stack((x[:frame+1], y[:frame+1], z[:frame+1]))
        
        if len(points) > 1:
            # Создаем сегменты (точка i и i+1)
            segments = np.array([points[i:i+2] for i in range(len(points)-1)])
            
            # Цвета для сегментов (от синего к красному через весь спектр)
            colors = plt.cm.plasma(np.linspace(0, 1, len(segments)))
            
            # Создаем коллекцию линий
            trail_collection = Line3DCollection(segments, colors=colors, 
                                               linewidths=linewidth, alpha=alpha)
            ax.add_collection3d(trail_collection)
        
        # Текущая точка
        current_point.set_data([x[frame]], [y[frame]])
        current_point.set_3d_properties([z[frame]])
        
        # Обновление информации
        progress = (frame / len(x)) * 100
        info_text.set_text(f'Time: {t[frame]:.2f}\n'
                          f'Position: ({x[frame]:.3f}, {y[frame]:.3f}, {z[frame]:.3f})')
        progress_text.set_text(f'Progress: {progress:.1f}% ({frame+1}/{len(x)})')
        
        # Обновление заголовка
        ax.set_title(f'Аттрактор Рёслера - Время: {t[frame]:.2f}', fontsize=14, pad=20)
        
        return current_point, info_text, progress_text
    
    # Создание анимации
    anim = FuncAnimation(fig, update, frames=len(x), init_func=init,
                        interval=interval, blit=False, repeat=True)
    
    plt.tight_layout()
    plt.show()
    
    return anim

def main():
    
    # Выбор файла
    filename = input("Название файла:").strip()
    if not filename:
        filename = "result.txt"
    
    # Параметры анимации
    print("\nНастройки анимации:")
    skip = int(input("Шаг прореживания (1 - все точки, 10 - каждую 10-ю): ") or "1")
    interval = int(input("Интервал между кадрами (мс): ") or "10")
    linewidth = 1.5
    alpha = 0.8
    
    print("\nЗапуск анимации...")
    create_3d_animation_gradient_trail(filename, skip, interval, linewidth, alpha)

if __name__ == "__main__":
    main()