import numpy as np
from numpy import *
from math import *
import sys, os,  codecs
import sympy
from scipy import sparse

def main():
    #path = input('input1.txt')
    #try:
    file = open('input1.txt', 'r')# encoding='utf-8'
    #except:
        #print('文件不存在，请检查，程序结束运行。\n')
    data = file.readlines()
    print(data)
    time = 0
    line = data[time].strip('\r\n')#去除换行符
    print(line)
    temps = line.split(' ', line.count(' '))#根据空格分割
    print(temps)
    gannum = int(temps[0])#计算杆数目、节点数目、位移数目、受力数目
    diannum = int(temps[1])
    weiyinum = int(temps[2])
    linum = int(temps[3])
    #print(gannum)

    #建立杆点矩阵G
    G = np.zeros((gannum, 4))
    for i in range(diannum+1,diannum+gannum+1):
        line = data[i].strip('\r\n')
        temps1 = line.split(' ', line.count(' '))
        print(temps1)
        for j in range(0, 4):
            G[i-diannum-1][j]=int(temps1[j])
    print("G")
    print(G)

    #建立点位移矩阵B
    B = np.zeros((diannum, 5))
    for i in range(1,diannum+1):
        line = data[i].strip('\r\n')
        temps2 = line.split(' ', line.count(' '))
        for j in range(0,5):
            B[i-1][j]=int(temps2[j])
    print("B")
    print(B)

    #杆长、位移编码lamuda建立、杆转角、EA，EL
    pai=math.pi
    lamuda = np.zeros((gannum, 6))
    l = np.zeros((gannum))
    a = np.zeros((gannum))
    EA = np.zeros((gannum))
    EI = np.zeros((gannum))
    for i in range(1,gannum+1):
        print("B[int(G[i - 1][1]) - 1][2]")
        print(int(B[int(G[i - 1][0]) - 1][2]))
        print(B[int(G[i - 1][1]) - 1][3])
        print(B[int(G[i - 1][1]) - 1][4])
        lamuda[i - 1] = [int(B[int(G[i - 1][0]) - 1][2]), int(B[int(G[i - 1][0]) - 1][3]), int(B[int(G[i - 1][0]) - 1][4]),
                         int(B[int(G[i - 1][1]) - 1][2]), int(B[int(G[i - 1][1]) - 1][3]), int(B[int(G[i - 1][1]) - 1][4])]
        print("lamuda")
        print(lamuda)

        l[i-1]=sqrt((B[int(G[i-1][0])-1][0]-B[int(G[i-1][1])-1][0])**2+
                    (B[int(G[i-1][0])-1][1]-B[int(G[i-1][1])-1][1])**2)

        #print("B[int(G[i - 1][0]) - 1][0] - B[int(G[i - 1][1]) - 1][0]")
        #print((B[int(G[i - 1][0]) - 1][1] - B[int(G[i - 1][1]) - 1][1])/(B[int(G[i - 1][0]) - 1][0] - B[int(G[i - 1][1]) - 1][0]))
        #print(atan2(0))
        if  int(B[int(G[i - 1][0]) - 1][0] - B[int(G[i - 1][1]) - 1][0])!=0:
            a[i - 1] = math.atan2((B[int(G[i - 1][0]) - 1][1] - B[int(G[i - 1][1]) - 1][1]),
                             (B[int(G[i - 1][1]) - 1][0] - B[int(G[i - 1][0]) - 1][0]))*180/pai
        if (B[int(G[i - 1][0]) - 1][0] - B[int(G[i - 1][1]) - 1][0])==0 and (
                B[int(G[i - 1][0]) - 1][1] - B[int(G[i - 1][1]) - 1][1]) < 0:
            a[i-1] = 90
        if (B[int(G[i - 1][0]) - 1][0] - B[int(G[i - 1][1]) - 1][0]) == 0 and (
                B[int(G[i - 1][0]) - 1][1] - B[int(G[i - 1][1]) - 1][1]) > 0:
            a[i - 1] = -90
        #print("a[i - 1]")
        #print(a[i-1])
        EA[i-1]=G[i-1][2]
        EI[i - 1] = G[i - 1][3]
    print("l")
    print(l)
    print("a")
    print(a)
    print("lamuda")
    print(lamuda)
    print("EA")
    print(EA)
    print("EI")
    print(EI)

    #建立单元元素矩阵
    ke = np.zeros((gannum, 6, 6))
    for i in range(0, gannum):
        ke[i]=k66(EA[i],EI[i],l[i],a[i])#调用下方单元矩阵建立函数k66
    print("ke")
    print(ke)

    #组装总刚度矩阵
    K=np.zeros((weiyinum,weiyinum))
    for i in range (0, weiyinum):
        for j in range (0, weiyinum):
            for n in range (0, gannum):
                for m in range (0,6):
                    for s in range (0,6):
                        if i == lamuda[n][m]-1 and j==lamuda[n][s]-1:
                            K[i][j]=K[i][j]+ke[n][m][s]
    print("K")
    print(K)

    #建立外荷载
    #输入1 1 1 30：1/2为节点力/固端力，1/2为节点/杆编号，1(q)(Fx)/2(F)/3为第几种力/位移编码，30为大小KN/KN/m(向下为正）
    Fe = np.zeros((gannum,6))
    Pe = np.zeros((gannum, 6))
    PE = np.zeros((weiyinum))
    for i in range (0,linum):
        line = data[i+1+gannum+diannum].strip('\r\n')
        temps3 = line.split(' ', line.count(' '))
        if int(temps3[0])==2:
            if int(temps3[2])==1:
                Fe[int(temps3[1])-1]=[0,int(temps3[3])*l[int(temps3[1])-1]/2, int(temps3[3])*l[int(temps3[1])-1]*l[int(temps3[1])-1]/12,
                                    0,int(temps3[3])*l[int(temps3[1])-1]/2, -int(temps3[3])*l[int(temps3[1])-1]*l[int(temps3[1])-1]/12]
            if int(temps3[2])==2:
                Fe[int(temps3[1])-1]=[0,int(temps3[3])/2,int(temps3[3])*l[int(temps3[1])-1]/8,
                                    0,int(temps3[3])/2,-int(temps3[3])*l[int(temps3[1])-1]/8]
                print("int(temps3[1]")
                print(int(temps3[1]))
    for i in range(0,gannum):
        for j in range(0,6):
            PE[int(lamuda[i][j])-1]=PE[int(lamuda[i][j])-1]-Fe[i][j]#建立总荷载PE
            #print(PE[int(lamuda[i][j])-1])
    Pe=-Fe
    Pj=np.zeros((gannum, 6))
    PJ=np.zeros((weiyinum))
    for i in range (0,linum):
        line = data[i+1+gannum+diannum].strip('\r\n')
        temps4 = line.split(' ', line.count(' '))
        if int(temps4[0]) == 1:
            for i in range (0,gannum):
                for j in range (0,6):
                    if lamuda[i][j]==int(temps4[2]):
                        Pj[i][j]=int(temps4[3])
                        PJ[int(temps4[2])-1]=int(temps4[3])#建立总荷载PJ
    print("Fe")
    print(Fe)
    print("Pj")
    print(Pj)
    print("PJ")
    print(PJ)
    print("PE")
    print(PE)
    P=PJ+PE
    DERTA=np.zeros((weiyinum))
    DERTA=np.dot(P, np.linalg.inv(K))#计算总位移DERTA
    print("总位移DERTA")
    print(DERTA)
    derta=np.zeros((gannum, 6))
    for i in range (0,gannum):
        for j in range (0,6):
            if int(lamuda[i][j])!=0:
                derta[i][j]=DERTA[int(lamuda[i][j])-1]#计算节点位移derta
    print("节点位移derta")
    print(derta)

    #计算节点力
    F=np.zeros((gannum,6))
    for i in range(0,gannum) :
        F[i]=np.dot(derta[i],ke[i])+Fe[i]#总体坐标系下的方程
        print(i+1)
        print("杆的杆端力")
        print(F[i])
    #print(F)
    print('计算完成\n')


# 建立单元刚度矩阵
def k66(EA, EI, l, a):
    pai=math.pi
    T=np.mat([[cos(a*pai/180), sin(a*pai/180), 0, 0, 0, 0],
        [-sin(a*pai/180), cos(a*pai/180), 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, cos(a*pai/180), sin(a*pai/180), 0],
        [0, 0, 0, -sin(a*pai/180), cos(a*pai/180), 0],
        [0, 0, 0, 0, 0, 1]])
    #print("T")
    #print(T)
    k6=np.mat([[EA/l, 0, 0, -EA/l, 0, 0],
        [0, (12*EI)/(l**3), (6*EI)/(l**2), 0, -(12*EI)/(l**3), (6*EI)/(l**2)],
        [0, (6*EI)/(l**2), (4*EI)/l, 0, -(6*EI)/(l**2), (2*EI)/l],
        [-EA / l, 0, 0, EA / l, 0, 0],
        [0, -(12 * EI) / (l **3), -(6 * EI) / (l ** 2), 0, (12 * EI) / (l ** 3), -(6 * EI) / (l ** 2)],
        [0, (6 * EI) / (l ** 2), (2 * EI) / l, 0, -(6 * EI) / (l ** 2), (4 * EI) / l]])
    ke1=T.T @ k6
    ke=ke1 @ T#根据角度进行变换
    return (ke)



if __name__ == '__main__':
    main()