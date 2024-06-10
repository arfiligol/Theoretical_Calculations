import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# 設置參數
length = 10  # 波導長度
c = 1  # 波速
frequency = 1  # 頻率
wavelength = c / frequency  # 波長
k = 2 * np.pi / wavelength  # 波數
omega = 2 * np.pi * frequency  # 角頻率
x = np.linspace(0, length, 10000)  # 空間點
t = np.linspace(0, 40, 10000)  # 時間點，增加到 4 以顯示更多周期

# 計算波形
def incident_wave(x, t):
    return np.sin(k * x - omega * t)

def reflected_wave(x, t):
    return np.sin(k * (length - x) - omega * t)

# 創建動畫
fig, ax = plt.subplots()
line1, = ax.plot(x, np.zeros_like(x), label='Incident Wave', color='blue')
line2, = ax.plot(x, np.zeros_like(x), label='Reflected Wave', color='red')
line3, = ax.plot(x, np.zeros_like(x), label='Standing Wave', color='green')
ax.set_ylim(-2, 2)
ax.set_xlim(0, length)
ax.set_xlabel('Position')
ax.set_ylabel('Amplitude')
ax.set_title('Wave Propagation and Standing Wave Formation')
ax.legend()

def update(frame):
    t_current = t[frame]
    y1 = incident_wave(x, t_current)
    y2 = reflected_wave(x, t_current) * (x <= t_current * c)
    y3 = y1 + y2
    line1.set_ydata(y1 * (x <= t_current * c))
    line2.set_ydata(y2)
    line3.set_ydata(y3)
    return line1, line2, line3

ani = animation.FuncAnimation(fig, update, frames=len(t), blit=True)

# 使用 ffmpeg 保存動畫
Writer = animation.writers['ffmpeg']
writer = Writer(fps=120, metadata=dict(artist='Me'), bitrate=1800)
ani.save('wave_propagation_and_reflection.mp4', writer=writer)
