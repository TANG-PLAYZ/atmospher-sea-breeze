import math as m
import numpy as np
import matplotlib.pyplot as plt

# 比熱
Co = 3890   #海洋比熱(J/kg K)
Cl = 2500   #陸地比熱(J/kg K)

# 質量
Mo = 3e8*27.75**2   #海洋質量 (kg)
Ml = 3e8*27.75**2   #陸地質量 (kg)

# 空氣密度 1.293x(實際壓力/101325)x(273.15/實際絕對溫度)
RHOo = 1.2  #近海洋空氣密度 (kg/m3)
RHOl = 1.2  #近陸地空氣密度 (kg/m3)

# 時間
t = 600                     #模擬間隔 (s)
times = int(86400*30/t)     #模擬次數
tick = np.arange(1,times,1) #計次

# 輻射量
S = 1366        #短波輻射 (W/m2)
theta = 0       #輻射入射角度 (rad)
Sr = 0          #地表接收輻射 (W/m2)
a = 0.3         #反射率
e = 5.67e-8     #史帝芬-波茲曼常數 (W/m2 K4)

# 畫圖用陣列
Toa=np.arange(1.,times,1.)
Tla=np.arange(1.,times,1.)
To=np.arange(1.,times,1.)
Tl=np.arange(1.,times,1.)
Top0=np.arange(1.,times,1.)
Top1=np.arange(1.,times,1.)
Ta2=np.arange(1.,times,1.)
dT=np.zeros(times-1)
Gh=np.zeros(times-1)

# 初始溫度(K)
T02=280     #海洋第二層大氣
T12=280     #陸地第二層大氣
T00=290.    #海洋第一層大氣
T01=290.    #陸地第一層大氣
T10=310.    #海洋地表
T11=310.    #陸地地表

# 地表熱容量變化
dQo=0   #海
dQl=0   #陸

# 對流公式:q=hA(Ta-Ts)
ho=0.6  #水散熱係數
hl=2    #岩石散熱係數
h=10    #自然對流常數

# 地表面積
Ao=27750**2 #海
Al=27750**2 #陸

for i in range(times):
    #計算地表接收輻射
    Sr = S * m.sin(2*m.pi*theta/86400) #計算入射角度
    if (Sr<0):
        Sr, S0 = 0, 0
    else:
        S0 = 1366
    
    #第二層大氣
    T02=round(((T00**4)/2)**0.25,4) #海洋
    T12=round(((T01**4)/2)**0.25,4) #陸地

    dQo = round(Sr*(1-a)*Ao-Ao*e*T10**4+Ao*e*T00**4-ho*(T10-T00)*Ao,4)*t        #計算水域熱量變化(以初始溫度時為零)
    T00=round(((T10**4)/2+ho*(T10-T00)/e/2-h*(T00-T01)/e/2+(T02**4)/2)**0.25,4) #以熱量平恆計算水域上方大氣溫度
    
    T10+= round(dQo/Mo/Co,4)    #以熱量變化計算水域溫度變化
    
    dQl = round(Sr*(1-a)*Al-Al*e*T11**4+Al*e*T01**4-hl*(T11-T01)*Al,4)*t        #計算陸地熱量變化(以初始溫度時為零)
    
    T01=round(((T11**4)/2+hl*(T11-T01)/e/2-h*(T01-T00)/e/2+(T12**4)/2)**0.25,4) #以熱量平恆計算陸地上方大氣溫度
    T11+= round(dQl/Ml/Cl,4)    #以熱量變化計算陸地溫度變化

    theta+=1*t
    #將結果記錄在陣列中
    Toa[i-1]=T00
    Tla[i-1]=T01
    To[i-1]=T10
    Tl[i-1]=T11
    Ta2[i-1]=T12
    Top0[i-1]=T02
    dT[i-1]=T00-T01
    #計算氣壓梯度力(加速度) Gh=-(R*dT)/(M*dn) *M:莫耳質量; dT:溫度差; dn:兩地距離
    Gh[i-1]=-((8.314462618*dT[i-1])/(0.0289647*27750))  #假設兩地距離27.75km


ws = np.zeros(times-1)
dt = t
Gh[:3310] = 0
ws = np.cumsum(Gh) * dt    #計算風速(m/s) Gh=dv/dt

#plt.plot(tick,Toa,'b.')
#plt.plot(tick,Tla,'g.')
#plt.plot(tick,dT,'c.')
#plt.plot(tick,Gh,'y.')
plt.plot(tick,ws,'k.')
#plt.title('Horizontal Pressure Gradient Force',fontsize=15)
plt.title('Simulation of Land-sea Winds',fontsize=15)
plt.xlabel('Time(10min)',fontsize=12)
plt.ylabel('Velocity(m/s)',fontsize=12)
plt.legend(['Wind speed'])
#plt.xlim([6000,times])
plt.show()
