import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
import math
import matplotlib.pyplot as plt
from tkinter import Label
from PIL import Image, ImageTk
from tkinter import Canvas,NW
from io import BytesIO
# CL=[1.0672, 1.0672, 1.0672, 1.0672, 1.0672,1.033, 1.0988, 1.0799, 1.0799, 1.0799,1.0799, 1.0799, 1.0799, 1.0799, 1.0799,1.1303, 1.1303, 1.1303, 1.1303, 1.1303,1.1303]
# a = [5.5,5.5,5.5,5.5,5,5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5]
# Cd = [0.02601, 0.02601, 0.02601, 0.02601, 0.02601,0.02552, 0.02574, 0.01357, 0.01357, 0.01357,0.01357, 0.01357, 0.01357, 0.01357, 0.01357,0.01357, 0.01125, 0.01125, 0.01125, 0.01125,0.01125, 0.01125]
# 


def calculate_air_density(altitude):
    rho_0 = 1.225  # 海平面上的标准空气密度，单位：kg/m^3
    L = 0.0065     # 大气温度随海拔变化的标准温度梯度，单位：K/m
    T_0 = 288.15   # 海平面上的标准温度，单位：K
    g = 9.80665    # 重力加速度，单位：m/s^2
    R = 287.058     # 气体常数，单位：J/(kg·K)
    
    rho = rho_0 * ((T_0 - (L * altitude)) / T_0) ** (g / (R * L))
    return rho

def calculate(rated_power, rated_wind_speed, blade_number, wind_utilization_coefficient, transmission_efficiency, tip_speed_ratio, cutoff_number, air_density, a, CL, Cd):
    rotor_diameter = (8 * rated_power / (math.pi * air_density * rated_wind_speed ** 3 * wind_utilization_coefficient * transmission_efficiency)) ** (1/2)
    hub_radius = rotor_diameter * 0.05
    blade_length = rotor_diameter / 2
    firstfile = blade_length * 0.15
    per_section = (blade_length - firstfile) / (cutoff_number - 1)
    n = (30 * rated_wind_speed * tip_speed_ratio) / (math.pi * rotor_diameter / 2)
    angle_speed = (2 * math.pi * n) / 60

    AI = [0] * cutoff_number
    BI = [0] * cutoff_number
    p_SUM = 0
    CI = [0] * cutoff_number
    B = [0] * cutoff_number
    pcr_ri = [0] * cutoff_number
    Speedper = [0] * cutoff_number
    Fai = [0] * cutoff_number
    DT = [0] * cutoff_number

    for i in range(cutoff_number):
        per_ri = firstfile + i * per_section
        pcr_ri[i] = per_ri
        speedper = tip_speed_ratio * per_ri / blade_length
        Speedper[i] = speedper
        fai = (1/3) * math.atan(speedper) + (math.pi / 3)
        Ki = (((speedper * speedper) + 1)) ** (1/2) * math.cos(fai)
        Ai = (1 - Ki) / 2
        Hi = math.sqrt(1 + ((1 - Ki) / (speedper * speedper)))
        Bi = (Hi - 1) / 2
        AI[i] = Ai
        BI[i] = Bi
        ifei = math.atan((1 - Ai) / (speedper * (1 + Bi))) * (180 / math.pi)
        Fai[i] = ifei

        Re = (2 * rated_wind_speed * per_ri * 0.9) / 0.00001698
        Ci = (per_ri * 8 * math.pi * Ai * math.sin(ifei * math.pi / 180) * math.sin(ifei * math.pi / 180)) / (blade_number * CL[i] * (1 - Ai) * math.cos(ifei * math.pi / 180) * math.cos(ifei * math.pi / 180))
        CI[i] = Ci
        B[i] = ifei - a[i]

        vx = rated_wind_speed * (1 - Ai)
        vy = angle_speed * per_ri * (1 + Bi)
        W = math.sqrt((vx ** 2) + (vy ** 2))
        dt = (CL[i] * math.sin(ifei * math.pi / 180)) - (Cd[i] * math.cos(ifei * math.pi / 180))
        dx = (CL[i] * math.sin(ifei * math.pi / 180)) + (Cd[i] * math.cos(ifei * math.pi / 180))
        DT[i] = dt

        if i == 0:
            dm = 0.5 * blade_number * air_density * Ci * W * W * dt * ((per_ri * per_ri) / 2)
        else:
            dm = 0.5 * blade_number * air_density * Ci * W * W * dt * (((per_ri * per_ri) - (pcr_ri[i - 1] * pcr_ri[i - 1])) / 2)
        Pi = dm * angle_speed
        p_SUM = p_SUM + Pi

    P_power = p_SUM * 0.9
    ideal_power = 0.5 * air_density * rated_wind_speed ** 3 * math.pi * (rotor_diameter / 2) ** 2
    Cp = P_power / ideal_power
    return pcr_ri, Fai, AI, BI, CI, DT, B, Cp,rotor_diameter

def verify(V, rated_power, rated_wind_speed, blade_number, wind_utilization_coefficient, transmission_efficiency, tip_speed_ratio, cutoff_number, air_density, a, CL, Cd):
    pcr_ri, Fai, AI, BI, CI, DT, B, Cp ,rotor_diameter= calculate(
        rated_power, rated_wind_speed, blade_number, wind_utilization_coefficient,
        transmission_efficiency, tip_speed_ratio, cutoff_number, air_density, a, CL, Cd
    )

    Angle_Speed = [0] * len(V)
    dM = [0] * len(V)
    and_Speed = [[0 for _ in range(cutoff_number)] for _ in range(len(V))]
    DM = [[0 for _ in range(cutoff_number)] for _ in range(len(V))]
    nomalangle_speed = [0] * len(V)

    for i in range(len(V)):
        for j in range(cutoff_number):
            angle_speed = ((1 - AI[j]) * V[i]) / (math.tan(Fai[j] * (math.pi / 180)) * (1 + BI[j]) * pcr_ri[j])
            Angle_Speed[i] += angle_speed
            if j == cutoff_number - 1:
                Angle_Speed[i] /= cutoff_number
                nomalangle_speed[i] = Angle_Speed[i]

    n_ram = [0] * len(V)
    ratio = [0] * len(V)
    p = [0] * len(V)
    cP = [0] * len(V)
    shoot = [0] * len(V)
    Dm = [0] * len(V)

    for i in range(len(V)):
        if i > rated_wind_speed:
            Angle_Speed[i] = Angle_Speed[rated_wind_speed + 1]

    for i in range(len(V)):
        n_ram[i] = (60 * Angle_Speed[i]) / (2 * math.pi)
        ratio[i] = Angle_Speed[i] * pcr_ri[-1] / V[i]

    for i in range(len(V)):
        for j in range(cutoff_number):
            and_Speed[i][j] = math.sqrt((V[i] * (1 - AI[j])) ** 2 + (Angle_Speed[i] * pcr_ri[j] * (1 + BI[j])) ** 2)
            if j == 0:
                DM[i][j] = 0.5 * blade_number * air_density * CI[j] * and_Speed[i][j] ** 2 * DT[j] * ((pcr_ri[j] ** 2) / 2)
            else:
                DM[i][j] = 0.5 * blade_number * air_density * CI[j] * and_Speed[i][j] ** 2 * DT[j] * ((pcr_ri[j] ** 2 - pcr_ri[j - 1] ** 2) / 2)
            dM[i] += DM[i][j]
        if i > rated_wind_speed:
            dM[i] = dM[rated_wind_speed + 1]

        p[i] = dM[i] * Angle_Speed[i]*1.1
        Dm[i] = p[i] / Angle_Speed[i]
        cP[i] = p[i] / (0.5 * air_density * V[i] ** 3 * math.pi * (pcr_ri[-1]) ** 2) * transmission_efficiency
        shoot[i] = (2 * p[i]) / nomalangle_speed[i]

    data1 = {'x': pcr_ri, 'y': B}
    data2 = {'x': pcr_ri, 'y': CI}
    data3 = {'x': V, 'y': cP}
    data4 = {'x': V, 'y': n_ram}
    data5 = {'x': V, 'y': shoot}
    data6 = {'x': V, 'y': p}
    data7 = {'x': V, 'y': ratio}
    data8 = {'x': V, 'y': Dm}

    return data1, data2, data3, data4, data5, data6, data7, data8

# 创建Tkinter应用程序窗口
class WindTurbineApp:
    def __init__(self, root):
        self.root = root
        self.root.title("风力计模拟性能曲线 by Alexyin")
        self.current_offset =0
        self.create_widgets()


    def create_widgets(self):
        self.frame = ttk.Frame(self.root)
        self.frame.pack(padx=10, pady=1, fill=tk.BOTH, expand=True)

        # 输入参数
        self.create_input_fields()

        # 计算按钮
        self.calc_button = ttk.Button(self.frame, text="计算", command=self.calculate_and_plot)
        self.calc_button.grid(row=13, columnspan=2, pady=10)

       
        # 绘图区域
        self.fig, self.axs = plt.subplots(4, 2, figsize=(9, 9))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canvas.get_tk_widget().grid(row=0, column=2, rowspan=15, pady=20)
        label_text ="风力发电机组各个性能曲线"
        label =Label(self.frame,text=label_text)
        label.grid(row=0,column=2,pady=30)

        
        self.radius_value = tk.StringVar()
        self.radius_value.set("风轮直径显示区域")  # 初始值设为0.0，你可以根据需要修改
        self.radius_display = ttk.Label(self.frame, textvariable=self.radius_value)
        self.radius_display.grid(row=0, column=2, sticky=tk.W, pady=60)
        # 显示署名
        signature_text = "Made by Alexyin"
        signature_label = Label(self.frame, text=signature_text)
        signature_label.grid(row=14, column=2, pady=10, sticky="se", padx=10)


    def create_input_fields(self):
        fields = [
            ("额定功率 (W):", "rated_power"),
            ("额定风速 (m/s):", "rated_wind_speed"),
            ("叶片数:", "blade_number"),
            ("风能利用系数:", "wind_utilization_coefficient"),
            ("传动损失:", "transmission_efficiency"),
            ("叶尖速比:", "tip_speed_ratio"),
            ("海拔 (m):", "altitude"),
            ("截面数:", "cutoff_number"),
            ("切入风速(m/s):", "start_wind"),
            ("切出风速 (m/s):", "end_wind"),
            ("安装角:", "a"),
            ("升力系数:", "CL"),
            ("阻力系数:", "Cd"),
        ]

        self.entries = {}
        for i, (label_text, var_name) in enumerate(fields):
            label = ttk.Label(self.frame, text=label_text)
            label.grid(row=i, column=0, sticky=tk.W, pady=2)

            entry = ttk.Entry(self.frame)
            entry.grid(row=i, column=1, pady=2)
            self.entries[var_name] = entry

    def calculate_and_plot(self):
        try:
            rated_power = float(self.entries["rated_power"].get())
            rated_wind_speed = int(self.entries["rated_wind_speed"].get())
            blade_number = int(self.entries["blade_number"].get())
            wind_utilization_coefficient = float(self.entries["wind_utilization_coefficient"].get())
            transmission_efficiency = float(self.entries["transmission_efficiency"].get())
            tip_speed_ratio = float(self.entries["tip_speed_ratio"].get())
            altitude = float(self.entries["altitude"].get())
            cutoff_number = int(self.entries["cutoff_number"].get())
            start_wind = int(self.entries["start_wind"].get())
            end_wind = int(self.entries["end_wind"].get())

            V = [start_wind + (0.5 * i) for i in range((end_wind - start_wind) * 2 + 1)]
            
            a = [float(val.strip()) for val in self.entries["a"].get().split(',')]
            CL = [float(val.strip()) for val in self.entries["CL"].get().split(',')]
            Cd = [float(val.strip()) for val in self.entries["Cd"].get().split(',')]
            

            air_density = calculate_air_density(altitude)
            rotor_diameter = (8 * rated_power / (math.pi * air_density * rated_wind_speed ** 3 * wind_utilization_coefficient * transmission_efficiency)) ** (1/2)

            self.radius_value.set("叶轮直径为: {:.2f}m".format(rotor_diameter))
            data1, data2, data3, data4, data5, data6, data7, data8 = verify(
                V, rated_power, rated_wind_speed, blade_number, wind_utilization_coefficient,
            transmission_efficiency, tip_speed_ratio, cutoff_number, air_density, a, CL, Cd,)

         # 清空之前的绘图
            for ax in self.axs.flatten():
                ax.clear()
    
         # 绘制新的数据
            self.axs[0, 0].plot(data1['x'], data1['y'])
            self.axs[0, 0].set_xlabel('R(m)')
            self.axs[0, 0].set_ylabel('Install Angle(°)')

            self.axs[0, 1].plot(data2['x'], data2['y'])
            self.axs[0, 1].set_xlabel('R(m)')
            self.axs[0, 1].set_ylabel('Chord Length(m)')

            self.axs[1, 0].plot(data3['x'], data3['y'])
            self.axs[1, 0].set_xlabel('V(m/s)')
            self.axs[1, 0].set_ylabel('CP')

            self.axs[1, 1].plot(data4['x'], data4['y'])
            self.axs[1, 1].set_xlabel('V(m/s)')
            self.axs[1, 1].set_ylabel('Rev Speed(r/min)')

            self.axs[2, 0].plot(data5['x'], data5['y'])
            self.axs[2, 0].set_xlabel('V(m/s)')
            self.axs[2, 0].set_ylabel('Thrust(N)')

            self.axs[2, 1].plot(data6['x'], data6['y'])
            self.axs[2, 1].set_xlabel('V(m/s)')
            self.axs[2, 1].set_ylabel('P(w)')

            self.axs[3, 0].plot(data7['x'], data7['y'])
            self.axs[3, 0].set_xlabel('V(m/s)')
            self.axs[3, 0].set_ylabel('Tip Speed')

            self.axs[3, 1].plot(data8['x'], data8['y'])
            self.axs[3, 1].set_xlabel('V(m/s)')
            self.axs[3, 1].set_ylabel('Triv Torque(N)')

            self.canvas.draw()
    
        except ValueError as e:
            print("Invalid input:", e)


# 启动Tkinter应用程序
if __name__ == "__main__":
    root = tk.Tk()
    app = WindTurbineApp(root)
    root.mainloop()
